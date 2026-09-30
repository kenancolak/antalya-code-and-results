// Antalya endpoint land-cover audit and retained-core screening
//
// The canonical 2017 and 2025 endpoint maps are the analysis inputs used to
// reconstruct retained cores and endpoint-threshold sensitivity. An optional
// diagnostic branch rebuilds the Random Forest classifications from the
// currently available Sentinel-2 catalogue and private training assets.
// Historical validation matrices are archived separately in
// data/derived/classification_validation_confusion_matrices.csv because a
// current-catalogue rebuild is not guaranteed to reproduce the historical
// map/validation pairing exactly.
//
// The temporally consistent SCL-based p90/p95/p97 screen is implemented in
// 03_scene_meteorology_annual_scl.js.

// 0. Analysis settings
var RUN_UID = '20260422_unified_canonical_v02';

var RF_TREES = 500;
var RF_BAG_FRAC = 0.632;
var RF_SEED = 42;
var CLOUD_PCT_MAX = 40;

var VALID_FRAC17_MIN      = 0.99;
var VALID_FRAC25_MIN      = 0.99;
var COMMON_FRAC_FULL_MIN  = 0.98;
var COMMON_FRAC_UNION_MIN = 0.98;

var CHANGE_THR_DEFAULT = 0.10;
var PURITY_THR_DEFAULT = 0.85;

var changeThrList = ee.List([0.05, 0.10, 0.15, 0.20]);
var purityThrList = ee.List([0.70, 0.75, 0.80, 0.85, 0.90, 0.95]);
var MIN_UU = 10;
var MIN_FF = 10;

// Routine verification uses the canonical endpoint images. Set the diagnostic
// flag to true only to audit a current-catalogue Random Forest reconstruction.
var RUN_CURRENT_CATALOG_RECONSTRUCTION_DIAGNOSTIC = false;
var PRINT_THRESHOLD_SWEEP = false;
var CREATE_RECONSTRUCTION_EXPORT_TASKS = false;

var USE_INLAND_WATER_FILTER = false;
// var inlandWaterFC = ee.FeatureCollection('projects/tikkenan/assets/inlandwater');

var SHOW_NDVI = false;
var SHOW_LULC = false;
var SHOW_CORES = false;
var EXPORT_FOLDER = 'GEE_Exports';
var RUN_TAG = 'sinir_denemefish5_train3aCanonical_' + RUN_UID;

// 1. Input assets and map context
var boundaryFC     = ee.FeatureCollection('projects/tikkenan/assets/sinir');
var gridFC      = ee.FeatureCollection('projects/tikkenan/assets/denemefish5');
var trainingFC    = ee.FeatureCollection('projects/tikkenan/assets/train3a_balanced_canonical_v01');
var validation2017FC   = ee.FeatureCollection('projects/tikkenan/assets/valpoints200_17');
var validation2025FC   = ee.FeatureCollection('projects/tikkenan/assets/valpoints200_25');
var canonicalLulc2017Raw = ee.Image('projects/tikkenan/assets/final_lulc2017raw');
var canonicalLulc2025Raw = ee.Image('projects/tikkenan/assets/final_lulc2025raw');

var aoiBounds = boundaryFC.geometry().bounds();

Map.centerObject(boundaryFC, 10);
Map.addLayer(
  boundaryFC.style({color: 'red', fillColor: '00000000', width: 2}),
  {},
  'Boundary',
  false
);
Map.addLayer(
  gridFC.style({color: 'yellow', fillColor: '00000000', width: 1}),
  {},
  'Grid',
  false
);

// 2. Standardise validation labels
var val2017 = validation2017FC.map(function(f) {
  return f.set('lc', ee.Number(f.get('LC')).toInt());
}).filterBounds(boundaryFC);

var val2025 = validation2025FC.map(function(f) {
  return f.set('lc', ee.Number(f.get('LC')).toInt());
}).filterBounds(boundaryFC);

print('Boundary feature count:', boundaryFC.size());
print('Grid feature count:', gridFC.size());
print('Diagnostic training asset size:', trainingFC.size());
print('Diagnostic training asset histogram:', trainingFC.aggregate_histogram('lc'));
print('Validation 2017 points:', val2017.size());
print('Validation 2025 points:', val2025.size());
print('Validation 2017 lc histogram:', val2017.aggregate_histogram('lc'));
print('Validation 2025 lc histogram:', val2025.aggregate_histogram('lc'));

// 3. Sentinel-2 preprocessing
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');

function maskS2clouds(image) {
  var qa = image.select('QA60');
  var cloudBitMask  = 1 << 10;
  var cirrusBitMask = 1 << 11;

  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
    .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask)
    .divide(10000)
    .select(['B2','B3','B4','B5','B6','B7','B8','B11','B12'])
    .copyProperties(image, ['system:time_start']);
}

function addIndices(image) {
  var ndvi = image.normalizedDifference(['B8','B4']).rename('NDVI');
  var ndbi = image.normalizedDifference(['B11','B8']).rename('NDBI');
  var ndwi = image.normalizedDifference(['B3','B8']).rename('NDWI');

  var tcb = image.expression(
    '0.25*B2 + 0.3*B3 + 0.33*B4 + 0.1*B8 + 0.6*B11 + 0.6*B12', {
      B2: image.select('B2'),
      B3: image.select('B3'),
      B4: image.select('B4'),
      B8: image.select('B8'),
      B11: image.select('B11'),
      B12: image.select('B12')
    }
  ).rename('TCB');

  var bsi = image.expression(
    '((B11 + B4) - (B8 + B2)) / ((B11 + B4) + (B8 + B2))', {
      B11: image.select('B11'),
      B4: image.select('B4'),
      B8: image.select('B8'),
      B2: image.select('B2')
    }
  ).rename('BSI');

  return image.addBands([ndvi, ndbi, ndwi, tcb, bsi]);
}

var bands = [
  'B2','B3','B4','B5','B6','B7','B8','B11','B12',
  'NDVI','NDBI','NDWI','TCB','BSI'
];

function makeSummerComposite(year) {
  year = ee.Number(year);
  var start = ee.Date.fromYMD(year, 6, 1);
  var end   = start.advance(3, 'month');

  var comp = s2
    .filterBounds(boundaryFC)
    .filterDate(start, end)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_PCT_MAX))
    .map(maskS2clouds)
    .map(addIndices)
    .median()
    .clipToCollection(boundaryFC)
    .set('year', year);

  return comp.select(bands);
}

// 4. Endpoint summer composites
var comp2017 = makeSummerComposite(2017);
var comp2025 = makeSummerComposite(2025);

if (SHOW_NDVI) {
  Map.addLayer(
    comp2017.select('NDVI'),
    {min: -0.2, max: 0.8, palette: ['blue', 'white', 'green']},
    'NDVI 2017',
    false
  );
  Map.addLayer(
    comp2025.select('NDVI'),
    {min: -0.2, max: 0.8, palette: ['blue', 'white', 'green']},
    'NDVI 2025',
    false
  );
}

// 5. Random Forest classification model
var rf = ee.Classifier.smileRandomForest({
  numberOfTrees: RF_TREES,
  bagFraction: RF_BAG_FRAC,
  seed: RF_SEED
}).train({
  features: trainingFC,
  classProperty: 'lc',
  inputProperties: bands
});

// 6. Endpoint land-cover maps
var lulc2017_raw = comp2017.classify(rf).rename('LULC').unmask(0).toInt16().set('year', 2017);
var lulc2025_raw = comp2025.classify(rf).rename('LULC').unmask(0).toInt16().set('year', 2025);

var kernel = ee.Kernel.square({radius: 1});
var lulc2017_sm = lulc2017_raw.focal_mode({kernel: kernel, iterations: 1}).toInt16();
var lulc2025_sm = lulc2025_raw.focal_mode({kernel: kernel, iterations: 1}).toInt16();

var vis = {
  min: 1,
  max: 5,
  palette: ['d73027','1a9850','fee08b','74add1','a6611a']
};

if (SHOW_LULC) {
  Map.addLayer(lulc2017_raw.updateMask(lulc2017_raw.neq(0)), vis, 'LULC 2017 RAW', false);
  Map.addLayer(lulc2025_raw.updateMask(lulc2025_raw.neq(0)), vis, 'LULC 2025 RAW', false);
}

// 7. Independent accuracy assessment
function accuracyReport(lulcImg, valPoints, label) {
  var val = lulcImg.sampleRegions({
    collection: valPoints,
    properties: ['lc'],
    scale: 10,
    geometries: false,
    tileScale: 4
  })
  .filter(ee.Filter.neq('LULC', 0))
  .filter(ee.Filter.notNull(['lc', 'LULC']));

  print('Validation sample size (' + label + '):', val.size());

  var cm = val.errorMatrix('lc', 'LULC');
  print('Confusion Matrix (' + label + '):', cm);
  print('OA (' + label + '):', cm.accuracy());
  print('Kappa (' + label + '):', cm.kappa());

  return val;
}

if (RUN_CURRENT_CATALOG_RECONSTRUCTION_DIAGNOSTIC) {
  print(
    'Diagnostic only: current-catalogue validation may differ from the archived historical matrices.'
  );
  accuracyReport(lulc2017_raw, val2017, '2017 CURRENT-CATALOGUE RAW');
  accuracyReport(lulc2025_raw, val2025, '2025 CURRENT-CATALOGUE RAW');
}

// 8. Grid-level endpoint stability and purity
var lulc17 = canonicalLulc2017Raw.rename('LULC17').toInt16();
var lulc25 = canonicalLulc2025Raw.rename('LULC25').toInt16();

var valid17 = lulc17.gt(0).rename('valid17').unmask(0);
var valid25 = lulc25.gt(0).rename('valid25').unmask(0);

var commonValid = valid17.and(valid25).rename('commonValid').unmask(0);
var unionValid  = valid17.or(valid25).rename('unionValid').unmask(0);

var u17 = lulc17.eq(1).rename('u17').unmask(0);
var f17 = lulc17.eq(2).rename('f17').unmask(0);
var a17 = lulc17.eq(3).rename('a17').unmask(0);
var g17 = lulc17.eq(4).rename('g17').unmask(0);
var b17 = lulc17.eq(5).rename('b17').unmask(0);

var u25 = lulc25.eq(1).rename('u25').unmask(0);
var f25 = lulc25.eq(2).rename('f25').unmask(0);
var a25 = lulc25.eq(3).rename('a25').unmask(0);
var g25 = lulc25.eq(4).rename('g25').unmask(0);
var b25 = lulc25.eq(5).rename('b25').unmask(0);

var chg = lulc17.neq(lulc25).and(commonValid).rename('chg').unmask(0);

var statsImg = ee.Image.cat([
  valid17, valid25, commonValid, unionValid,
  u17, f17, a17, g17, b17,
  u25, f25, a25, g25, b25,
  chg
]).toFloat();

var gridStats0 = statsImg.reduceRegions({
  collection: gridFC.filterBounds(boundaryFC),
  reducer: ee.Reducer.mean(),
  scale: 10,
  crs: 'EPSG:32636',
  tileScale: 4
});

var gridStats1;
if (USE_INLAND_WATER_FILTER) {
  gridStats1 = gridStats0.map(function(f) {
    var intersects = inlandWaterFC.filterBounds(f.geometry()).size().gt(0);
    return f.set('isInlandWater', ee.Number(intersects));
  });
} else {
  gridStats1 = gridStats0.map(function(f) {
    return f.set('isInlandWater', 0);
  });
}

function shortLabel(code) {
  code = ee.Number(code);
  return ee.String(
    ee.Algorithms.If(code.eq(1), 'U',
      ee.Algorithms.If(code.eq(2), 'F',
        ee.Algorithms.If(code.eq(3), 'A',
          ee.Algorithms.If(code.eq(4), 'G', 'B')
        )
      )
    )
  );
}

function addCoreBaseFields(f) {
  var vals17 = ee.List([
    ee.Number(f.get('u17')),
    ee.Number(f.get('f17')),
    ee.Number(f.get('a17')),
    ee.Number(f.get('g17')),
    ee.Number(f.get('b17'))
  ]);

  var vals25 = ee.List([
    ee.Number(f.get('u25')),
    ee.Number(f.get('f25')),
    ee.Number(f.get('a25')),
    ee.Number(f.get('g25')),
    ee.Number(f.get('b25'))
  ]);

  var purity17 = ee.Number(vals17.reduce(ee.Reducer.max()));
  var purity25 = ee.Number(vals25.reduce(ee.Reducer.max()));

  var dom17 = ee.Number(vals17.indexOf(purity17)).add(1);
  var dom25 = ee.Number(vals25.indexOf(purity25)).add(1);
  var domSame = dom17.eq(dom25);
  var changeRatio = ee.Number(f.get('chg'));

  var validFrac17 = ee.Number(f.get('valid17'));
  var validFrac25 = ee.Number(f.get('valid25'));
  var commonFracFull = ee.Number(f.get('commonValid'));
  var unionFrac = ee.Number(f.get('unionValid'));

  var commonFracUnion = ee.Number(
    ee.Algorithms.If(
      unionFrac.gt(0),
      commonFracFull.divide(unionFrac),
      0
    )
  );

  var enoughSupport = validFrac17.gte(VALID_FRAC17_MIN)
    .and(validFrac25.gte(VALID_FRAC25_MIN))
    .and(commonFracFull.gte(COMMON_FRAC_FULL_MIN))
    .and(commonFracUnion.gte(COMMON_FRAC_UNION_MIN));

  var notInlandWater = ee.Number(f.get('isInlandWater')).eq(0);
  var eligible = enoughSupport.and(notInlandWater);

  return f.set({
    'dom17': dom17,
    'dom25': dom25,
    'dom17_lbl': shortLabel(dom17),
    'dom25_lbl': shortLabel(dom25),
    'purity17': purity17,
    'purity25': purity25,
    'changeRatio': changeRatio,
    'domSame': ee.Number(domSame),
    'validFrac17': validFrac17,
    'validFrac25': validFrac25,
    'unionFrac': unionFrac,
    'commonFrac_full': commonFracFull,
    'commonFrac_union': commonFracUnion,
    'enoughSupport': ee.Number(enoughSupport),
    'notInlandWater': ee.Number(notInlandWater),
    'eligible': ee.Number(eligible)
  });
}

var gridStats = gridStats1.map(addCoreBaseFields);

function addDefaultCoreFields(f) {
  var dom17 = ee.Number(f.get('dom17'));

  var keep = ee.Number(f.get('eligible')).eq(1)
    .and(ee.Number(f.get('domSame')).eq(1))
    .and(ee.Number(f.get('changeRatio')).lte(CHANGE_THR_DEFAULT))
    .and(ee.Number(f.get('purity17')).gte(PURITY_THR_DEFAULT))
    .and(ee.Number(f.get('purity25')).gte(PURITY_THR_DEFAULT));

  var coreLabel = ee.String(
    ee.Algorithms.If(
      keep,
      shortLabel(dom17).cat(shortLabel(dom17)),
      'NA'
    )
  );

  var coreCode = ee.Number(
    ee.Algorithms.If(
      keep,
      dom17.multiply(11),
      0
    )
  );

  return f.set({
    'run_uid': RUN_UID,
    'keep_default': ee.Number(keep),
    'core_default': coreLabel,
    'coreCode_default': coreCode
  });
}

var gridCoreDefault = gridStats.map(addDefaultCoreFields);

var uuFC = gridCoreDefault.filter(ee.Filter.eq('core_default', 'UU'));
var ffFC = gridCoreDefault.filter(ee.Filter.eq('core_default', 'FF'));
var aaFC = gridCoreDefault.filter(ee.Filter.eq('core_default', 'AA'));
var ggFC = gridCoreDefault.filter(ee.Filter.eq('core_default', 'GG'));
var bbFC = gridCoreDefault.filter(ee.Filter.eq('core_default', 'BB'));

print('Default retained core counts (canonical endpoint assets):');
print('UU:', uuFC.size());
print('FF:', ffFC.size());
print('AA:', aaFC.size());
print('GG:', ggFC.size());
print('BB:', bbFC.size());
print('Eligible grid count:', gridStats.filter(ee.Filter.eq('eligible', 1)).size());

if (SHOW_CORES) {
  Map.addLayer(
    uuFC.style({color:'red', fillColor:'ff000055', width:1}),
    {},
    'UU cores',
    true
  );
}

// 9. Sensitivity of endpoint criteria
var gridLite = gridStats.select([
  'eligible', 'domSame', 'changeRatio', 'purity17', 'purity25', 'dom17'
]);

var endpointSweepList = changeThrList.map(function(chThr) {
  chThr = ee.Number(chThr);

  return purityThrList.map(function(pThr) {
    pThr = ee.Number(pThr);

    var screened = gridLite
      .filter(ee.Filter.eq('eligible', 1))
      .filter(ee.Filter.eq('domSame', 1))
      .filter(ee.Filter.lte('changeRatio', chThr))
      .filter(ee.Filter.gte('purity17', pThr))
      .filter(ee.Filter.gte('purity25', pThr));

    var nUU = screened.filter(ee.Filter.eq('dom17', 1)).size();
    var nFF = screened.filter(ee.Filter.eq('dom17', 2)).size();
    var nAA = screened.filter(ee.Filter.eq('dom17', 3)).size();
    var nGG = screened.filter(ee.Filter.eq('dom17', 4)).size();
    var nBB = screened.filter(ee.Filter.eq('dom17', 5)).size();

    var feasible = ee.Number(nUU).gte(MIN_UU).and(ee.Number(nFF).gte(MIN_FF));

    return ee.Feature(null, {
      'run_uid': RUN_UID,
      'changeThr': chThr,
      'purityThr': pThr,
      'nRetained': screened.size(),
      'nUU': nUU,
      'nFF': nFF,
      'nAA': nAA,
      'nGG': nGG,
      'nBB': nBB,
      'feasible': ee.Number(feasible)
    });
  });
});

var endpointSweepFC = ee.FeatureCollection(ee.List(endpointSweepList).flatten());

// 10. Endpoint checks and exports
print(
  'Endpoint retained-core counts:',
  gridCoreDefault.filter(ee.Filter.neq('core_default', 'NA'))
    .aggregate_histogram('core_default')
);
if (PRINT_THRESHOLD_SWEEP) {
  print('Endpoint threshold sweep:', endpointSweepFC);
}

if (CREATE_RECONSTRUCTION_EXPORT_TASKS) {
  Export.image.toDrive({
    image: lulc2017_raw,
    description: 'Antalya_LULC_2017_CURRENT_CATALOGUE_RAW_' + RUN_TAG,
    folder: EXPORT_FOLDER,
    fileNamePrefix: 'Antalya_LULC_2017_CURRENT_CATALOGUE_RAW_' + RUN_TAG,
    region: aoiBounds,
    scale: 10,
    crs: 'EPSG:32636',
    maxPixels: 1e13
  });

  Export.image.toDrive({
    image: lulc2025_raw,
    description: 'Antalya_LULC_2025_CURRENT_CATALOGUE_RAW_' + RUN_TAG,
    folder: EXPORT_FOLDER,
    fileNamePrefix: 'Antalya_LULC_2025_CURRENT_CATALOGUE_RAW_' + RUN_TAG,
    region: aoiBounds,
    scale: 10,
    crs: 'EPSG:32636',
    maxPixels: 1e13
  });

  Export.image.toDrive({
    image: lulc2017_sm,
    description: 'Antalya_LULC_2017_CURRENT_CATALOGUE_SMOOTH_' + RUN_TAG,
    folder: EXPORT_FOLDER,
    fileNamePrefix: 'Antalya_LULC_2017_CURRENT_CATALOGUE_SMOOTH_' + RUN_TAG,
    region: aoiBounds,
    scale: 10,
    crs: 'EPSG:32636',
    maxPixels: 1e13
  });

  Export.image.toDrive({
    image: lulc2025_sm,
    description: 'Antalya_LULC_2025_CURRENT_CATALOGUE_SMOOTH_' + RUN_TAG,
    folder: EXPORT_FOLDER,
    fileNamePrefix: 'Antalya_LULC_2025_CURRENT_CATALOGUE_SMOOTH_' + RUN_TAG,
    region: aoiBounds,
    scale: 10,
    crs: 'EPSG:32636',
    maxPixels: 1e13
  });
}

var gridCoreDefaultExport = gridCoreDefault.map(function(f) {
  return ee.Feature(null, f.toDictionary(f.propertyNames()));
});

Export.table.toDrive({
  collection: gridCoreDefaultExport,
  description: 'EndpointCore_GridStats_Default_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'EndpointCore_GridStats_Default_' + RUN_TAG,
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: endpointSweepFC,
  description: 'EndpointCore_ThresholdSweep_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'EndpointCore_ThresholdSweep_' + RUN_TAG,
  fileFormat: 'CSV'
});
