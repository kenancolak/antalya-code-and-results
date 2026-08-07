// REPOSITORY STATUS
// Archived copy of the analysis workflow used in the study.
// Provided to document analytical logic, parameters, and export structure.
// This is not a turnkey reproducibility script: required study-created GEE
// assets remain private and exact rerunning from another account is not claimed.
// The analytical operations below have not been revised for this archive.
//
// Antalya land-cover classification and stability screening
//
// This script produces the 2017 and 2025 endpoint land-cover maps,
// identifies endpoint retained cores, and evaluates annual spectral stability.
// Run this script before the thermal analysis in 02_thermal_dlst.js.

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

var YEARS = ee.List.sequence(2017, 2025);
var PERCENTILES = ee.List([90, 95, 97]);
var DEFAULT_ANCHOR_PCT = 95;
var UNSTABLE_RATIO_DEFAULT = 0.10;
var unstableRatioThrList = ee.List([0.05, 0.10, 0.15, 0.20]);

var USE_BSI_DIAGNOSTIC = false;
var USE_INLAND_WATER_FILTER = false;
// var inlandWaterFC = ee.FeatureCollection('projects/tikkenan/assets/inlandwater');

var SHOW_NDVI = false;
var SHOW_LULC = false;
var SHOW_CORES = false;
var SHOW_ANNUAL = false;
var SHOW_MAXJUMP = false;

var EXPORT_FOLDER = 'GEE_Exports';
var RUN_TAG = 'sinir_denemefish5_train3aCanonical_' + RUN_UID;

// 1. Input assets and map context
var boundaryFC     = ee.FeatureCollection('projects/tikkenan/assets/sinir');
var gridFC      = ee.FeatureCollection('projects/tikkenan/assets/denemefish5');
var trainingFC    = ee.FeatureCollection('projects/tikkenan/assets/train3a_balanced_canonical_v01');
var validation2017FC   = ee.FeatureCollection('projects/tikkenan/assets/valpoints200_17');
var validation2025FC   = ee.FeatureCollection('projects/tikkenan/assets/valpoints200_25');

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
print('Canonical trainingFC training size:', trainingFC.size());
print('Canonical trainingFC training histogram:', trainingFC.aggregate_histogram('lc'));
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

// 7. Accuracy assessment using reserved validation data
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

var valRes2017 = accuracyReport(lulc2017_raw, val2017, '2017 RAW');
var valRes2025 = accuracyReport(lulc2025_raw, val2025, '2025 RAW');

// 8. Grid-level endpoint stability and purity
var lulc17 = lulc2017_raw.rename('LULC17').toInt16();
var lulc25 = lulc2025_raw.rename('LULC25').toInt16();

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

print('Default retained core counts (canonical asset chain):');
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

// 10. Annual Sentinel-2 composites
function makeAnnualIndexComposite(year) {
  year = ee.Number(year);
  var start = ee.Date.fromYMD(year, 6, 1);
  var end   = start.advance(3, 'month');

  var img = s2
    .filterBounds(boundaryFC)
    .filterDate(start, end)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_PCT_MAX))
    .map(maskS2clouds)
    .map(addIndices)
    .median()
    .clipToCollection(boundaryFC)
    .set('year', year);

  var keepBands = ee.List(['NDVI', 'NDBI', 'TCB']);
  keepBands = ee.List(
    ee.Algorithms.If(USE_BSI_DIAGNOSTIC, keepBands.add('BSI'), keepBands)
  );

  return img.select(keepBands);
}

var annualIC = ee.ImageCollection.fromImages(
  YEARS.map(function(y) {
    return makeAnnualIndexComposite(y);
  })
);

print('Annual summer composite collection:', annualIC);

// 11. Consecutive-year spectral differences
var DIFF_YEARS = ee.List.sequence(2017, 2024);

var annualDiffIC = ee.ImageCollection.fromImages(
  DIFF_YEARS.map(function(y) {
    y = ee.Number(y);

    var img1 = ee.Image(annualIC.filter(ee.Filter.eq('year', y)).first());
    var img2 = ee.Image(annualIC.filter(ee.Filter.eq('year', y.add(1))).first());

    var diff = img2.subtract(img1).abs();
    var out = diff.select(['NDVI', 'NDBI', 'TCB'], ['dNDVI', 'dNDBI', 'dTCB']);

    out = ee.Image(
      ee.Algorithms.If(
        USE_BSI_DIAGNOSTIC,
        out.addBands(diff.select('BSI').rename('dBSI')),
        out
      )
    );

    return out.set({
      'startYear': y,
      'endYear': y.add(1)
    });
  })
);

print('Annual diff collection:', annualDiffIC);

// 12. Maximum annual spectral change
var maxJump = annualDiffIC.max();

var maxJumpMain = maxJump.select(
  ['dNDVI', 'dNDBI', 'dTCB'],
  ['maxdNDVI', 'maxdNDBI', 'maxdTCB']
);

maxJumpMain = ee.Image(
  ee.Algorithms.If(
    USE_BSI_DIAGNOSTIC,
    maxJumpMain.addBands(maxJump.select('dBSI').rename('maxdBSI')),
    maxJumpMain
  )
);

if (SHOW_MAXJUMP) {
  Map.addLayer(
    maxJumpMain.select('maxdNDVI'),
    {min: 0, max: 0.5, palette: ['white','orange','red']},
    'max |ΔNDVI|',
    false
  );
}

// 13. Calibration of annual-change thresholds
var anchorFC = gridCoreDefault.filter(ee.Filter.eq('keep_default', 1));
print('Anchor grid count:', anchorFC.size());

var anchorSample = maxJumpMain.clipToCollection(anchorFC).sample({
  region: anchorFC.geometry(),
  scale: 10,
  numPixels: 50000,
  seed: 42,
  geometries: false,
  tileScale: 4
});

print('Anchor sample size:', anchorSample.size());

var ndviDict = anchorSample.reduceColumns({
  reducer: ee.Reducer.percentile(PERCENTILES),
  selectors: ['maxdNDVI']
});

var ndbiDict = anchorSample.reduceColumns({
  reducer: ee.Reducer.percentile(PERCENTILES),
  selectors: ['maxdNDBI']
});

var tcbDict = anchorSample.reduceColumns({
  reducer: ee.Reducer.percentile(PERCENTILES),
  selectors: ['maxdTCB']
});

var thresholdFC = ee.FeatureCollection(
  PERCENTILES.map(function(p) {
    p = ee.Number(p);
    var key = ee.String('p').cat(p.format('%.0f'));

    return ee.Feature(null, {
      'run_uid': RUN_UID,
      'anchor_pct': p,
      'thr_ndvi': ee.Number(ndviDict.get(key)),
      'thr_ndbi': ee.Number(ndbiDict.get(key)),
      'thr_tcb': ee.Number(tcbDict.get(key))
    });
  })
);

print('Annual threshold calibration table:', thresholdFC);

// 14. Default annual-change thresholds
var defaultThrFeat = ee.Feature(
  thresholdFC.filter(ee.Filter.eq('anchor_pct', DEFAULT_ANCHOR_PCT)).first()
);

var THR_NDVI = ee.Number(defaultThrFeat.get('thr_ndvi'));
var THR_NDBI = ee.Number(defaultThrFeat.get('thr_ndbi'));
var THR_TCB  = ee.Number(defaultThrFeat.get('thr_tcb'));

print('Default annual-change thresholds:');
print('anchor_pct:', DEFAULT_ANCHOR_PCT);
print('thr_ndvi:', THR_NDVI);
print('thr_ndbi:', THR_NDBI);
print('thr_tcb :', THR_TCB);

// 15. Annual stability classification
var unstableNDVI = maxJumpMain.select('maxdNDVI').gt(THR_NDVI).rename('unstNDVI').unmask(0);
var unstableNDBI = maxJumpMain.select('maxdNDBI').gt(THR_NDBI).rename('unstNDBI').unmask(0);
var unstableTCB  = maxJumpMain.select('maxdTCB').gt(THR_TCB).rename('unstTCB').unmask(0);
var unstableAny  = unstableNDVI.or(unstableNDBI).or(unstableTCB).rename('unstAny').unmask(0);

var annualMetricsImg = ee.Image.cat([
  maxJumpMain,
  unstableNDVI.toFloat(),
  unstableNDBI.toFloat(),
  unstableTCB.toFloat(),
  unstableAny.toFloat()
]);

var annualGrid0 = annualMetricsImg.reduceRegions({
  collection: gridCoreDefault,
  reducer: ee.Reducer.mean(),
  scale: 10,
  crs: 'EPSG:32636',
  tileScale: 4
});

function addAnnualFlags(f) {
  var unstableRatio = ee.Number(f.get('unstAny'));

  var urbanDominant = ee.Number(f.get('dom17')).eq(1)
    .and(ee.Number(f.get('dom25')).eq(1));

  var annualKeep = ee.Number(f.get('eligible')).eq(1)
    .and(urbanDominant)
    .and(unstableRatio.lte(UNSTABLE_RATIO_DEFAULT));

  var confirmedStableUrban = annualKeep
    .and(ee.String(f.get('core_default')).compareTo('UU').eq(0));

  return f.set({
    'run_uid': RUN_UID,
    'default_anchor_pct': DEFAULT_ANCHOR_PCT,
    'default_thr_ndvi': THR_NDVI,
    'default_thr_ndbi': THR_NDBI,
    'default_thr_tcb': THR_TCB,
    'default_unstableRatio_thr': UNSTABLE_RATIO_DEFAULT,
    'urbanDominant_annualBase': ee.Number(urbanDominant),
    'unstableRatio_default': unstableRatio,
    'annualKeep_default': ee.Number(annualKeep),
    'confirmedStableUrban_default': ee.Number(confirmedStableUrban)
  });
}

var annualGridDefault = annualGrid0.map(addAnnualFlags);

print(
  'Annual stable urban count (default):',
  annualGridDefault.filter(ee.Filter.eq('annualKeep_default', 1)).size()
);

print(
  'Confirmed stable urban count (default intersection with UU):',
  annualGridDefault.filter(ee.Filter.eq('confirmedStableUrban_default', 1)).size()
);

// 15a. Endpoint and annually stable urban populations

var endpointUU_all = annualGridDefault
  .filter(ee.Filter.eq('core_default', 'UU'));

var spectralConfirmedUU = annualGridDefault
  .filter(ee.Filter.eq('confirmedStableUrban_default', 1));

var spectralRejectedUU = annualGridDefault
  .filter(ee.Filter.eq('core_default', 'UU'))
  .filter(ee.Filter.eq('confirmedStableUrban_default', 0));

print('Diagnostic | Endpoint UU total:', endpointUU_all.size());
print('Diagnostic | Spectrally confirmed UU:', spectralConfirmedUU.size());
print('Diagnostic | Spectrally rejected UU:', spectralRejectedUU.size());

print(
  'Spectrally rejected UU grids:',
  spectralRejectedUU.select([
    'core_default',
    'keep_default',
    'changeRatio',
    'purity17',
    'purity25',
    'unstableRatio_default',
    'annualKeep_default',
    'confirmedStableUrban_default'
  ])
);

// 15b. Spatial diagnostic for rejected urban cores

// Purpose:
// - Show unstable pixels across the full retained anchor population.
// - Show only unstable pixels, without a white raster background.
// - Draw rejected UU grids with red outline.
// - Draw approved anchor grids with black outline.

// All retained endpoint cores: expected 83 grids
var anchorAll = annualGridDefault
  .filter(ee.Filter.eq('keep_default', 1));

// Rejected UU grids:
// endpoint UU, but not confirmed by annual spectral stability.
// Expected: 4 grids.
var rejectedUU = anchorAll
  .filter(ee.Filter.eq('core_default', 'UU'))
  .filter(ee.Filter.eq('confirmedStableUrban_default', 0));

// Approved anchor grids:
// all 83 anchor grids except the rejected UU grids.
// Because ee.Filter.not() does not exist in GEE JS API,
// we use the logical equivalent:
// approved = non-UU anchors OR confirmed UU anchors.
var approvedAnchor = anchorAll.filter(
  ee.Filter.or(
    ee.Filter.neq('core_default', 'UU'),
    ee.Filter.eq('confirmedStableUrban_default', 1)
  )
);

// Unstable pixels only within the full 83 anchor grid population.
// selfMask() hides all zero-valued stable pixels.
// Therefore, no white/black raster background is drawn.
var unstablePixels_anchorAll = unstableAny
  .selfMask()
  .clipToCollection(anchorAll);

// Count checks
print('Diagnostic | Anchor grids total:', anchorAll.size());
print('Diagnostic | Approved anchor grids:', approvedAnchor.size());
print('Diagnostic | Rejected UU grids:', rejectedUU.size());

// 1) Unstable pixels within all 83 anchors: red pixels only
Map.addLayer(
  unstablePixels_anchorAll,
  {
    min: 1,
    max: 1,
    palette: ['ff0000']
  },
  'Unstable pixels within all 83 anchor grids',
  true
);

// 2) Approved anchor grids: black outline, no fill
Map.addLayer(
  approvedAnchor.style({
    color: '000000',
    fillColor: '00000000',
    width: 2
  }),
  {},
  'Approved anchor grids: black outline',
  true
);

// 3) Rejected UU grids: red outline, no fill
Map.addLayer(
  rejectedUU.style({
    color: 'ff0000',
    fillColor: '00000000',
    width: 3
  }),
  {},
  'Rejected UU grids: red outline',
  true
);

// Compact count panel for visual verification.

var uuPanel = ui.Panel({
  style: {
    position: 'top-left',
    padding: '8px',
    width: '260px'
  }
});

uuPanel.add(ui.Label({
  value: 'UU stability check',
  style: {
    fontWeight: 'bold',
    fontSize: '14px'
  }
}));

endpointUU_all.size().evaluate(function(n) {
  uuPanel.add(ui.Label('Endpoint UU cores: ' + n));
});

spectralConfirmedUU.size().evaluate(function(n) {
  uuPanel.add(ui.Label('Spectrally confirmed UU: ' + n));
});

spectralRejectedUU.size().evaluate(function(n) {
  uuPanel.add(ui.Label('Rejected UU grids: ' + n));
});

Map.add(uuPanel);

// 16. Sensitivity of annual-stability criteria
var annualSweepList = PERCENTILES.map(function(p) {
  p = ee.Number(p);
  var feat = ee.Feature(thresholdFC.filter(ee.Filter.eq('anchor_pct', p)).first());

  var thrNDVI = ee.Number(feat.get('thr_ndvi'));
  var thrNDBI = ee.Number(feat.get('thr_ndbi'));
  var thrTCB  = ee.Number(feat.get('thr_tcb'));

  var uNDVI = maxJumpMain.select('maxdNDVI').gt(thrNDVI).unmask(0);
  var uNDBI = maxJumpMain.select('maxdNDBI').gt(thrNDBI).unmask(0);
  var uTCB  = maxJumpMain.select('maxdTCB').gt(thrTCB).unmask(0);

  var uAny = uNDVI.or(uNDBI).or(uTCB).rename('unstAnyP').unmask(0).toFloat();

  var annualP0 = uAny.reduceRegions({
    collection: gridCoreDefault,
    reducer: ee.Reducer.mean(),
    scale: 10,
    crs: 'EPSG:32636',
    tileScale: 4
  });

  var annualP = annualP0.map(function(f) {
    var val = ee.Number(
      ee.Algorithms.If(
        f.propertyNames().contains('unstAnyP'),
        f.get('unstAnyP'),
        ee.Algorithms.If(
          f.propertyNames().contains('mean'),
          f.get('mean'),
          -9999
        )
      )
    );
    return f.set('unstableRatioP', val);
  });

  return unstableRatioThrList.map(function(ur) {
    ur = ee.Number(ur);

    var stableUrban = annualP
      .filter(ee.Filter.eq('eligible', 1))
      .filter(ee.Filter.eq('dom17', 1))
      .filter(ee.Filter.eq('dom25', 1))
      .filter(ee.Filter.lte('unstableRatioP', ur));

    var intersection = stableUrban.filter(ee.Filter.eq('core_default', 'UU'));

    return ee.Feature(null, {
      'run_uid': RUN_UID,
      'anchor_pct': p,
      'thr_ndvi': thrNDVI,
      'thr_ndbi': thrNDBI,
      'thr_tcb': thrTCB,
      'unstableRatio_thr': ur,
      'is_default_combo': ee.Number(p.eq(DEFAULT_ANCHOR_PCT).and(ur.eq(UNSTABLE_RATIO_DEFAULT))),
      'nAnnualStableUrban': stableUrban.size(),
      'nIntersection_UU': intersection.size()
    });
  });
});

var annualSweepTables = ee.FeatureCollection(ee.List(annualSweepList).flatten());

// 17. Run summary
var finalSummary = ee.FeatureCollection([
  ee.Feature(null, {
    'run_uid': RUN_UID,
    'UU_default': uuFC.size(),
    'FF_default': ffFC.size(),
    'AA_default': aaFC.size(),
    'GG_default': ggFC.size(),
    'BB_default': bbFC.size(),
    'anchor_count': anchorFC.size(),
    'keep_default_count': gridCoreDefault.filter(ee.Filter.eq('keep_default', 1)).size(),
    'default_anchor_pct': DEFAULT_ANCHOR_PCT,
    'default_thr_ndvi': THR_NDVI,
    'default_thr_ndbi': THR_NDBI,
    'default_thr_tcb': THR_TCB,
    'default_unstableRatio_thr': UNSTABLE_RATIO_DEFAULT,
    'annual_stable_default': annualGridDefault.filter(ee.Filter.eq('annualKeep_default', 1)).size(),
    'confirmed_stable_default': annualGridDefault.filter(ee.Filter.eq('confirmedStableUrban_default', 1)).size()
  })
]);

print('Run summary:', finalSummary);

// 18. Confirmed urban-core export
var annualGridWithCoords = annualGridDefault.map(function(f) {
  var c_ll = f.geometry().centroid(1);
  var xy_ll = ee.List(c_ll.coordinates());

  var c_utm = c_ll.transform('EPSG:32636', 1);
  var xy_utm = ee.List(c_utm.coordinates());

  return f.set({
    'GridID': ee.String(f.get('system:index')),
    'lon': ee.Number(xy_ll.get(0)),
    'lat': ee.Number(xy_ll.get(1)),
    'x_utm': ee.Number(xy_utm.get(0)),
    'y_utm': ee.Number(xy_utm.get(1))
  });
});

var confirmedStableUrbanFC = annualGridWithCoords
  .filter(ee.Filter.eq('confirmedStableUrban_default', 1))
  .map(function(f) {
    return ee.Feature(f.geometry(), {
      run_uid:  f.get('run_uid'),
      GridID:   f.get('GridID'),
      core_def: f.get('core_default'),
      keep_def: f.get('keep_default'),
      dom17:    f.get('dom17_lbl'),
      dom25:    f.get('dom25_lbl'),
      pur17:    f.get('purity17'),
      pur25:    f.get('purity25'),
      chg_rt:   f.get('changeRatio'),
      anc_pct:  f.get('default_anchor_pct'),
      th_ndvi:  f.get('default_thr_ndvi'),
      th_ndbi:  f.get('default_thr_ndbi'),
      th_tcb:   f.get('default_thr_tcb'),
      ur_thr:   f.get('default_unstableRatio_thr'),
      urb_dom:  f.get('urbanDominant_annualBase'),
      ur_def:   f.get('unstableRatio_default'),
      ann_keep: f.get('annualKeep_default'),
      conf_stb: f.get('confirmedStableUrban_default'),
      lon:      f.get('lon'),
      lat:      f.get('lat'),
      x_utm:    f.get('x_utm'),
      y_utm:    f.get('y_utm')
    });
  });

// 19. Export tables and rasters
Export.image.toDrive({
  image: unstableAny.toByte().clipToCollection(boundaryFC),
  description: 'Annual_UnstablePixels_Default_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'Annual_UnstablePixels_Default_' + RUN_TAG,
  region: aoiBounds,
  scale: 10,
  crs: 'EPSG:32636',
  maxPixels: 1e13
});

Export.image.toDrive({
  image: maxJumpMain.clipToCollection(boundaryFC),
  description: 'Annual_MaxJumpStack_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'Annual_MaxJumpStack_' + RUN_TAG,
  region: aoiBounds,
  scale: 10,
  crs: 'EPSG:32636',
  maxPixels: 1e13
});

Export.table.toDrive({
  collection: confirmedStableUrbanFC,
  description: 'ConfirmedStableUrban_Default_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'ConfirmedStableUrban_Default_' + RUN_TAG,
  fileFormat: 'SHP'
});

Export.image.toDrive({
  image: lulc2017_raw,
  description: 'Antalya_LULC_2017_RAW_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'Antalya_LULC_2017_RAW_' + RUN_TAG,
  region: aoiBounds,
  scale: 10,
  crs: 'EPSG:32636',
  maxPixels: 1e13
});

Export.image.toDrive({
  image: lulc2025_raw,
  description: 'Antalya_LULC_2025_RAW_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'Antalya_LULC_2025_RAW_' + RUN_TAG,
  region: aoiBounds,
  scale: 10,
  crs: 'EPSG:32636',
  maxPixels: 1e13
});

Export.image.toDrive({
  image: lulc2017_sm,
  description: 'Antalya_LULC_2017_SMOOTH_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'Antalya_LULC_2017_SMOOTH_' + RUN_TAG,
  region: aoiBounds,
  scale: 10,
  crs: 'EPSG:32636',
  maxPixels: 1e13
});

Export.image.toDrive({
  image: lulc2025_sm,
  description: 'Antalya_LULC_2025_SMOOTH_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'Antalya_LULC_2025_SMOOTH_' + RUN_TAG,
  region: aoiBounds,
  scale: 10,
  crs: 'EPSG:32636',
  maxPixels: 1e13
});

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

Export.table.toDrive({
  collection: thresholdFC,
  description: 'AnnualSpectral_ThresholdCalibration_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'AnnualSpectral_ThresholdCalibration_' + RUN_TAG,
  fileFormat: 'CSV'
});

var annualGridDefaultExport = annualGridDefault.map(function(f) {
  return ee.Feature(null, f.toDictionary(f.propertyNames()));
});

Export.table.toDrive({
  collection: annualGridDefaultExport,
  description: 'AnnualSpectral_GridStats_Default_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'AnnualSpectral_GridStats_Default_' + RUN_TAG,
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: annualSweepTables,
  description: 'AnnualSpectral_ThresholdSweep_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'AnnualSpectral_ThresholdSweep_' + RUN_TAG,
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: finalSummary,
  description: 'FINAL_Summary_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'FINAL_Summary_' + RUN_TAG,
  fileFormat: 'CSV'
});
