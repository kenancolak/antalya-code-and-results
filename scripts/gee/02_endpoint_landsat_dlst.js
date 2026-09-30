// Antalya endpoint land-surface-temperature analysis
//
// This script reconstructs the retained cores from the canonical endpoint
// land-cover maps audited in 01_endpoint_landcover_retained_cores.js, calculates 2017–2025 LST change
// on common valid Landsat pixels, and exports the primary and sensitivity tables.
//
// Class codes: 1 urban, 2 forest, 3 agriculture, 4 greenhouse, 5 bare surface.
// Projected analysis CRS: EPSG:32636.

// 0. Analysis settings
var RUN_UID = '20260426_tc_v03';

var boundaryFC = ee.FeatureCollection('projects/tikkenan/assets/sinir');
var gridFC  = ee.FeatureCollection('projects/tikkenan/assets/denemefish5');

// Update these paths if the endpoint maps are stored under different asset names.
var lulc2017_raw = ee.Image('projects/tikkenan/assets/final_lulc2017raw');
var lulc2025_raw = ee.Image('projects/tikkenan/assets/final_lulc2025raw');

var USE_INLAND_WATER_FILTER = false;
// var inlandWaterFC = ee.FeatureCollection('projects/tikkenan/assets/inlandwater');

var VALID_FRAC17_MIN      = 0.99;
var VALID_FRAC25_MIN      = 0.99;
var COMMON_FRAC_FULL_MIN  = 0.98;
var COMMON_FRAC_UNION_MIN = 0.98;

var CHANGE_THR_DEFAULT = 0.10;
var PURITY_THR_DEFAULT = 0.85;

// Threshold sensitivity lists for row-level thermal export
var changeThrList = ee.List([0.05, 0.10, 0.15, 0.20]);
var purityThrList = ee.List([0.70, 0.75, 0.80, 0.85, 0.90, 0.95]);

var MIN_UU = 10;
var MIN_FF = 10;

var USE_WRS_FILTER = true;
var WRS_PATH = 178;
var WRS_ROW  = 34;

var EXPORT_RETAINED_ONLY = true;
var LANDSAT_COMMON_FRAC_MIN = 0.00;

var SHOW_DFLT_LST = false;
var SHOW_DFLT_DLST = false;
var SHOW_DEFAULT_COMMON = false;
// Printing the 24-threshold x four-support row-level collection can exceed the
// Code Editor memory limit. Keep false for routine verification; the export
// task and the released sensitivity table remain available.
var PRINT_THRESHOLD_THERMAL_ROWS = false;

var EXPORT_FOLDER = 'GEE_Exports';
var RUN_TAG = 'sinir_denemefish5_train3aCanonical_' + RUN_UID;

// 1. Map context and input checks
Map.centerObject(boundaryFC, 10);
Map.addLayer(
  boundaryFC.style({color:'red', fillColor:'00000000', width:2}),
  {},
  'Boundary',
  false
);
Map.addLayer(
  gridFC.style({color:'yellow', fillColor:'00000000', width:1}),
  {},
  'Grid',
  false
);

print('Boundary feature count:', boundaryFC.size());
print('Grid feature count:', gridFC.size());

// 2. Reconstruct retained cores from endpoint land-cover maps
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

var gridCoreDefault0 = gridStats.map(addDefaultCoreFields);

var gridCoreDefault = gridCoreDefault0.map(function(f) {
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

print('Default retained core counts:');
print('UU:', gridCoreDefault.filter(ee.Filter.eq('core_default', 'UU')).size());
print('FF:', gridCoreDefault.filter(ee.Filter.eq('core_default', 'FF')).size());
print('AA:', gridCoreDefault.filter(ee.Filter.eq('core_default', 'AA')).size());
print('GG:', gridCoreDefault.filter(ee.Filter.eq('core_default', 'GG')).size());
print('BB:', gridCoreDefault.filter(ee.Filter.eq('core_default', 'BB')).size());

// 3. Landsat LST preparation
var L8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');
var L9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2');

function maybeFilterWRS(ic) {
  if (USE_WRS_FILTER) {
    return ic.filter(ee.Filter.eq('WRS_PATH', WRS_PATH))
             .filter(ee.Filter.eq('WRS_ROW', WRS_ROW));
  }
  return ic;
}

function prepLST(img) {
  var qa = img.select('QA_PIXEL');

  var dilatedBit = 1 << 1;
  var cirrusBit  = 1 << 2;
  var cloudBit   = 1 << 3;
  var shadowBit  = 1 << 4;
  var snowBit    = 1 << 5;

  var mask = qa.bitwiseAnd(dilatedBit).eq(0)
    .and(qa.bitwiseAnd(cirrusBit).eq(0))
    .and(qa.bitwiseAnd(cloudBit).eq(0))
    .and(qa.bitwiseAnd(shadowBit).eq(0))
    .and(qa.bitwiseAnd(snowBit).eq(0));

  var lstC = img.select('ST_B10')
    .multiply(0.00341802)
    .add(149.0)
    .subtract(273.15)
    .rename('LST');

  return lstC.updateMask(mask)
    .copyProperties(img, ['system:time_start', 'SPACECRAFT_ID']);
}

function getPreparedCollection(sensor, start, end) {
  var ic;
  if (sensor === 'L8') {
    ic = maybeFilterWRS(L8);
  } else if (sensor === 'L9') {
    ic = maybeFilterWRS(L9);
  } else {
    ic = maybeFilterWRS(L8).merge(maybeFilterWRS(L9));
  }

  return ic
    .filterBounds(boundaryFC)
    .filterDate(start, end)
    .map(prepLST);
}

function imageOnDate(dateStr, sensor) {
  var start = ee.Date(dateStr);
  var end   = start.advance(1, 'day');
  var ic = getPreparedCollection(sensor, start, end);
  return ee.Image(ic.first()).set('scene_date', dateStr);
}

function meanFromDateList(dateList, sensor) {
  var imgs = ee.ImageCollection.fromImages(
    ee.List(dateList).map(function(d) {
      return imageOnDate(ee.String(d), sensor);
    })
  );

  return imgs.mean().rename('LST');
}

function scenarioEndpointMeans(dateList2017, sensor2017, dateList2025, sensor2025) {
  var mean2017 = meanFromDateList(dateList2017, sensor2017);
  var mean2025 = meanFromDateList(dateList2025, sensor2025);

  var commonMask = mean2017.mask().and(mean2025.mask());

  var lst2017_common = mean2017.updateMask(commonMask).rename('LST2017');
  var lst2025_common = mean2025.updateMask(commonMask).rename('LST2025');
  var dlst_common = lst2025_common.subtract(lst2017_common).rename('dLST');
  var commonBand = commonMask.rename('commonValidLST').unmask(0).toFloat();

  return ee.Image.cat([lst2017_common, lst2025_common, dlst_common, commonBand]);
}

// 4. Thermal-support scenarios
var DEFAULT_2017 = ee.List(['2017-06-25', '2017-07-11', '2017-08-12']);
var DEFAULT_2025 = ee.List(['2025-06-23', '2025-07-17', '2025-08-10']);

var S1_2017 = ee.List(['2017-06-25', '2017-07-11', '2017-08-12', '2017-08-28']);
var S1_2025 = ee.List(['2025-06-23', '2025-07-17', '2025-08-10', '2025-08-26']);

var S2_2017 = ee.List(['2017-06-25', '2017-07-11', '2017-08-12', '2017-08-28']);
var S2_2025 = ee.List([
  '2025-06-07', '2025-06-15', '2025-06-23',
  '2025-07-01', '2025-07-17', '2025-07-25',
  '2025-08-10', '2025-08-26'
]);

var SENSCHK_2017 = ee.List(['2017-06-25', '2017-07-11', '2017-08-12', '2017-08-28']);
var SENSCHK_2025 = ee.List(['2025-06-15', '2025-07-01', '2025-07-17']);

// 5. Aggregate thermal images to the fixed grid
// The same reducer can return retained cores or the full analytical grid.
// retainedOnly = true  -> retained cores used in the primary comparison
// retainedOnly = false -> full grid used in threshold sensitivity checks

function scenarioToGridTable(img, scenarioName, retainedOnly) {
  var reduced = img.reduceRegions({
    collection: gridCoreDefault,
    reducer: ee.Reducer.mean(),
    scale: 30,
    crs: 'EPSG:32636',
    tileScale: 4
  });

  reduced = reduced.map(function(f) {
    var commonFrac = ee.Number(
      ee.Algorithms.If(
        f.propertyNames().contains('commonValidLST'),
        ee.Algorithms.If(f.get('commonValidLST'), f.get('commonValidLST'), 0),
        ee.Algorithms.If(f.get('mean'), f.get('mean'), 0)
      )
    );

    return f.set({
      'run_uid': RUN_UID,
      'scenario': scenarioName,
      'lst2017': f.get('LST2017'),
      'lst2025': f.get('LST2025'),
      'dlst': f.get('dLST'),
      'commonFrac_landsat': commonFrac
    });
  });

  var filtered = ee.FeatureCollection(
    ee.Algorithms.If(
      retainedOnly,
      reduced.filter(ee.Filter.eq('keep_default', 1)),
      reduced
    )
  );

  filtered = filtered.filter(
    ee.Filter.gte('commonFrac_landsat', LANDSAT_COMMON_FRAC_MIN)
  );

  filtered = filtered.map(function(f) {
    return ee.Feature(f.geometry(), {
      run_uid: f.get('run_uid'),
      GridID: f.get('GridID'),
      scenario: f.get('scenario'),

      // Default-core information
      core_default: f.get('core_default'),
      keep_default: f.get('keep_default'),
      coreCode_default: f.get('coreCode_default'),

      // Endpoint LULC / stability information
      eligible: f.get('eligible'),
      domSame: f.get('domSame'),
      dom17: f.get('dom17'),
      dom25: f.get('dom25'),
      dom17_lbl: f.get('dom17_lbl'),
      dom25_lbl: f.get('dom25_lbl'),
      purity17: f.get('purity17'),
      purity25: f.get('purity25'),
      changeRatio: f.get('changeRatio'),
      validFrac17: f.get('validFrac17'),
      validFrac25: f.get('validFrac25'),
      commonFrac_full: f.get('commonFrac_full'),
      commonFrac_union: f.get('commonFrac_union'),

      // Class fractions, useful for diagnostics
      u17: f.get('u17'),
      f17: f.get('f17'),
      a17: f.get('a17'),
      g17: f.get('g17'),
      b17: f.get('b17'),
      u25: f.get('u25'),
      f25: f.get('f25'),
      a25: f.get('a25'),
      g25: f.get('g25'),
      b25: f.get('b25'),

      // Thermal values
      lst2017: f.get('lst2017'),
      lst2025: f.get('lst2025'),
      dlst: f.get('dlst'),
      commonFrac_landsat: f.get('commonFrac_landsat'),

      // Coordinates
      lon: f.get('lon'),
      lat: f.get('lat'),
      x_utm: f.get('x_utm'),
      y_utm: f.get('y_utm')
    });
  });

  return filtered;
}

// 6. Primary thermal support
var defaultImg = scenarioEndpointMeans(DEFAULT_2017, 'L8', DEFAULT_2025, 'L8L9');

// Retained-core table used by the primary analysis.
var defaultTable = scenarioToGridTable(defaultImg, 'DEFAULT_3SCENE', true);

// Full analytical grid used for threshold sensitivity.
var defaultTable_all = scenarioToGridTable(defaultImg, 'DEFAULT_3SCENE', false);

print('Default retained-core scenario row count:', defaultTable.size());
print('Default full-grid scenario row count:', defaultTable_all.size());

// ArcGIS-compatible retained-core exports
// These exports use the primary thermal table defined above.

// Retained-core polygons
var retainedPolys0 = gridCoreDefault.filter(
  ee.Filter.inList('core_default', ['UU', 'FF', 'AA', 'BB'])
);

// Short field names compatible with the shapefile format
var retainedPolys1 = retainedPolys0.map(function(f) {
  var gid = ee.String(
    ee.Algorithms.If(
      f.propertyNames().contains('GridID'),
      f.get('GridID'),
      f.get('system:index')
    )
  );

  var core = ee.String(f.get('core_default'));

  return ee.Feature(f.geometry(), {
    runid: RUN_UID,
    GridID: gid,
    core: core,
    lbl: core.cat('_').cat(gid),
    keepdf: f.get('keep_default'),
    pur17: f.get('purity17'),
    pur25: f.get('purity25'),
    chgrt: f.get('changeRatio'),
    lon: f.get('lon'),
    lat: f.get('lat'),
    xutm: f.get('x_utm'),
    yutm: f.get('y_utm')
  });
});

// Join the primary thermal values by GridID
var joinFilter = ee.Filter.equals({
  leftField: 'GridID',
  rightField: 'GridID'
});

var joined = ee.Join.saveFirst('mtch').apply({
  primary: retainedPolys1,
  secondary: defaultTable,
  condition: joinFilter
});

var retainedPolys = ee.FeatureCollection(joined).map(function(f) {
  var m = ee.Feature(f.get('mtch'));

  return ee.Feature(f.geometry(), {
    runid: f.get('runid'),
    GridID: f.get('GridID'),
    core: f.get('core'),
    lbl: f.get('lbl'),
    keepdf: f.get('keepdf'),
    pur17: f.get('pur17'),
    pur25: f.get('pur25'),
    chgrt: f.get('chgrt'),
    lst17: m.get('lst2017'),
    lst25: m.get('lst2025'),
    dlst:  m.get('dlst'),
    cfrac: m.get('commonFrac_landsat'),
    lon: f.get('lon'),
    lat: f.get('lat'),
    xutm: f.get('xutm'),
    yutm: f.get('yutm')
  });
});

// Centroid layer for point-based spatial checks
var retainedPts = retainedPolys.map(function(f) {
  return ee.Feature(f.geometry().centroid(1), {
    runid: f.get('runid'),
    GridID: f.get('GridID'),
    core: f.get('core'),
    lbl: f.get('lbl'),
    keepdf: f.get('keepdf'),
    pur17: f.get('pur17'),
    pur25: f.get('pur25'),
    chgrt: f.get('chgrt'),
    lst17: f.get('lst17'),
    lst25: f.get('lst25'),
    dlst:  f.get('dlst'),
    cfrac: f.get('cfrac'),
    lon: f.get('lon'),
    lat: f.get('lat'),
    xutm: f.get('xutm'),
    yutm: f.get('yutm')
  });
});

// Map layers used to verify class membership
var uuArc = retainedPolys.filter(ee.Filter.eq('core', 'UU'));
var ffArc = retainedPolys.filter(ee.Filter.eq('core', 'FF'));
var aaArc = retainedPolys.filter(ee.Filter.eq('core', 'AA'));
var bbArc = retainedPolys.filter(ee.Filter.eq('core', 'BB'));

Map.addLayer(
  uuArc.style({color: 'red', fillColor: 'ff000055', width: 2}),
  {},
  'ArcGIS UU',
  true
);

Map.addLayer(
  ffArc.style({color: '008000', fillColor: '00800055', width: 2}),
  {},
  'ArcGIS FF',
  true
);

Map.addLayer(
  aaArc.style({color: 'ffa500', fillColor: 'ffa50055', width: 2}),
  {},
  'ArcGIS AA',
  true
);

Map.addLayer(
  bbArc.style({color: '8b4513', fillColor: '8b451355', width: 2}),
  {},
  'ArcGIS BB',
  true
);

Map.addLayer(
  retainedPts.style({
    color: '000000',
    pointSize: 3,
    pointShape: 'circle'
  }),
  {},
  'ArcGIS retained centroids',
  false
);

// Console checks
print('Retained polygons for ArcGIS:', retainedPolys.limit(10));
print('Retained centroids for ArcGIS:', retainedPts.limit(10));
print('UU GridIDs:', uuArc.aggregate_array('GridID'));
print('FF GridIDs:', ffArc.aggregate_array('GridID'));
print('AA GridIDs:', aaArc.aggregate_array('GridID'));
print('BB GridIDs:', bbArc.aggregate_array('GridID'));

// Polygon shapefile export
Export.table.toDrive({
  collection: retainedPolys,
  description: 'ArcGIS_RetainedCores_Polygons_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'ArcGIS_RetainedCores_Polygons_' + RUN_TAG,
  fileFormat: 'SHP'
});

// Centroid shapefile export
Export.table.toDrive({
  collection: retainedPts,
  description: 'ArcGIS_RetainedCores_Centroids_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'ArcGIS_RetainedCores_Centroids_' + RUN_TAG,
  fileFormat: 'SHP'
});

print('Default scenario row count:', defaultTable.size());

if (SHOW_DFLT_LST) {
  Map.addLayer(
    defaultImg.select('LST2017'),
    {min:20, max:45, palette:['blue','cyan','yellow','orange','red']},
    'Default LST 2017',
    false
  );
  Map.addLayer(
    defaultImg.select('LST2025'),
    {min:20, max:45, palette:['blue','cyan','yellow','orange','red']},
    'Default LST 2025',
    false
  );
}

if (SHOW_DFLT_DLST) {
  Map.addLayer(
    defaultImg.select('dLST'),
    {min:-1, max:3, palette:['blue','white','red']},
    'Default dLST',
    false
  );
}

if (SHOW_DEFAULT_COMMON) {
  Map.addLayer(
    defaultImg.select('commonValidLST'),
    {min:0, max:1, palette:['black','white']},
    'Default common valid LST',
    false
  );
}

// 7. Alternative thermal supports
var s1Img = scenarioEndpointMeans(S1_2017, 'L8', S1_2025, 'L8L9');
var s2Img = scenarioEndpointMeans(S2_2017, 'L8', S2_2025, 'L8L9');
var sensChkImg = scenarioEndpointMeans(SENSCHK_2017, 'L8', SENSCHK_2025, 'L8');

// Scenario tables restricted to retained cores.
var s1Table = scenarioToGridTable(s1Img, 'S1_4SCENE', true);
var s2Table = scenarioToGridTable(s2Img, 'S2_ALLCLEAR_MEAN', true);
var sensChkTable = scenarioToGridTable(sensChkImg, 'SENSOR_2025_L8_ONLY_ALLCLEAR', true);

var allScenariosTable = defaultTable
  .merge(s1Table)
  .merge(s2Table)
  .merge(sensChkTable);

print('All retained-core scenarios row count:', allScenariosTable.size());

// Full-grid scenario tables used in threshold sensitivity checks.
var s1Table_all = scenarioToGridTable(s1Img, 'S1_4SCENE', false);
var s2Table_all = scenarioToGridTable(s2Img, 'S2_ALLCLEAR_MEAN', false);
var sensChkTable_all = scenarioToGridTable(sensChkImg, 'SENSOR_2025_L8_ONLY_ALLCLEAR', false);

var allScenariosTable_all = defaultTable_all
  .merge(s1Table_all)
  .merge(s2Table_all)
  .merge(sensChkTable_all);
  
// 7a. Threshold-specific grid membership and thermal values
// This table records membership under each threshold combination.
// Each row retains the corresponding endpoint and thermal values,
// allowing the sensitivity results to be reconstructed outside GEE.

function makeEndpointThresholdThermalRows(scenTable) {
  var comboList = changeThrList.map(function(chThr) {
    chThr = ee.Number(chThr);

    return purityThrList.map(function(pThr) {
      pThr = ee.Number(pThr);

      var screened = scenTable
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

      var feasible = ee.Number(nUU).gte(MIN_UU)
        .and(ee.Number(nFF).gte(MIN_FF));

      return screened.map(function(f) {
        var dom = ee.Number(f.get('dom17'));
        var coreThr = shortLabel(dom).cat(shortLabel(dom));

        return ee.Feature(null, {
          run_uid: RUN_UID,
          scenario: f.get('scenario'),
          GridID: f.get('GridID'),

          // Threshold combination
          changeThr: chThr,
          purityThr: pThr,
          feasible: ee.Number(feasible),

          // Counts repeated for transparent diagnostics
          nUU_combo: nUU,
          nFF_combo: nFF,
          nAA_combo: nAA,
          nGG_combo: nGG,
          nBB_combo: nBB,

          // Membership under this threshold
          keep_thr: 1,
          core_thr: coreThr,
          dom17: f.get('dom17'),
          dom25: f.get('dom25'),
          dom17_lbl: f.get('dom17_lbl'),
          dom25_lbl: f.get('dom25_lbl'),

          // Endpoint metrics
          changeRatio: f.get('changeRatio'),
          purity17: f.get('purity17'),
          purity25: f.get('purity25'),
          validFrac17: f.get('validFrac17'),
          validFrac25: f.get('validFrac25'),
          commonFrac_full: f.get('commonFrac_full'),
          commonFrac_union: f.get('commonFrac_union'),

          // Default-core status, for comparison only
          core_default: f.get('core_default'),
          keep_default: f.get('keep_default'),

          // Class fractions
          u17: f.get('u17'),
          f17: f.get('f17'),
          a17: f.get('a17'),
          g17: f.get('g17'),
          b17: f.get('b17'),
          u25: f.get('u25'),
          f25: f.get('f25'),
          a25: f.get('a25'),
          g25: f.get('g25'),
          b25: f.get('b25'),

          // Thermal values
          lst2017: f.get('lst2017'),
          lst2025: f.get('lst2025'),
          dlst: f.get('dlst'),
          commonFrac_landsat: f.get('commonFrac_landsat'),

          // Coordinates
          lon: f.get('lon'),
          lat: f.get('lat'),
          x_utm: f.get('x_utm'),
          y_utm: f.get('y_utm')
        });
      });
    });
  });

  return ee.FeatureCollection(ee.List(comboList).flatten()).flatten();
}

var endpointThresholdThermalRows_default =
  makeEndpointThresholdThermalRows(defaultTable_all);

var endpointThresholdThermalRows_s1 =
  makeEndpointThresholdThermalRows(s1Table_all);

var endpointThresholdThermalRows_s2 =
  makeEndpointThresholdThermalRows(s2Table_all);

var endpointThresholdThermalRows_sens =
  makeEndpointThresholdThermalRows(sensChkTable_all);

var endpointThresholdThermalRows =
  endpointThresholdThermalRows_default
    .merge(endpointThresholdThermalRows_s1)
    .merge(endpointThresholdThermalRows_s2)
    .merge(endpointThresholdThermalRows_sens);

if (PRINT_THRESHOLD_THERMAL_ROWS) {
  print(
    'Endpoint threshold x thermal row-level table:',
    endpointThresholdThermalRows.limit(20)
  );

  print(
    'Endpoint threshold x thermal row count:',
    endpointThresholdThermalRows.size()
  );

  // The default threshold combination should reproduce the primary retained-core set.
  print(
    'Default threshold rows, DEFAULT_3SCENE only:',
    endpointThresholdThermalRows
      .filter(ee.Filter.eq('scenario', 'DEFAULT_3SCENE'))
      .filter(ee.Filter.eq('changeThr', CHANGE_THR_DEFAULT))
      .filter(ee.Filter.eq('purityThr', PURITY_THR_DEFAULT))
      .limit(20)
  );
}

print('All full-grid scenarios row count:', allScenariosTable_all.size());

// 8. Run summary
var thermalSummary = ee.FeatureCollection([
  ee.Feature(null, {
    'run_uid': RUN_UID,
    'UU_default': gridCoreDefault.filter(ee.Filter.eq('core_default', 'UU')).size(),
    'FF_default': gridCoreDefault.filter(ee.Filter.eq('core_default', 'FF')).size(),
    'AA_default': gridCoreDefault.filter(ee.Filter.eq('core_default', 'AA')).size(),
    'GG_default': gridCoreDefault.filter(ee.Filter.eq('core_default', 'GG')).size(),
    'BB_default': gridCoreDefault.filter(ee.Filter.eq('core_default', 'BB')).size(),
    'default_rows': defaultTable.size(),
    'all_scenarios_rows': allScenariosTable.size(),
    'landsat_common_frac_min': LANDSAT_COMMON_FRAC_MIN
  })
]);

print('Run summary:', thermalSummary);

// 9. Export tables and rasters
Export.table.toDrive({
  collection: defaultTable,
  description: 'Grid_dLST_Default_byCore_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'Grid_dLST_Default_byCore_' + RUN_TAG,
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: allScenariosTable,
  description: 'Grid_dLST_AllScenarios_byCore_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'Grid_dLST_AllScenarios_byCore_' + RUN_TAG,
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: thermalSummary,
  description: 'FINAL_ThermalSummary_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'FINAL_ThermalSummary_' + RUN_TAG,
  fileFormat: 'CSV'
});

Export.image.toDrive({
  image: defaultImg.select('dLST').clipToCollection(boundaryFC),
  description: 'Default_dLST_Raster_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'Default_dLST_Raster_' + RUN_TAG,
  region: boundaryFC.geometry().bounds(),
  scale: 30,
  crs: 'EPSG:32636',
  maxPixels: 1e13
});

Export.image.toDrive({
  image: defaultImg.select('commonValidLST').clipToCollection(boundaryFC),
  description: 'Default_CommonValidLST_Raster_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'Default_CommonValidLST_Raster_' + RUN_TAG,
  region: boundaryFC.geometry().bounds(),
  scale: 30,
  crs: 'EPSG:32636',
  maxPixels: 1e13
});

// Primary thermal support for the full grid
Export.table.toDrive({
  collection: defaultTable_all,
  description: 'Grid_dLST_Default_AllGrids_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'Grid_dLST_Default_AllGrids_' + RUN_TAG,
  fileFormat: 'CSV'
});

// Alternative thermal supports for the full grid
Export.table.toDrive({
  collection: allScenariosTable_all,
  description: 'Grid_dLST_AllScenarios_AllGrids_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'Grid_dLST_AllScenarios_AllGrids_' + RUN_TAG,
  fileFormat: 'CSV'
});

// Threshold-specific row-level thermal table
Export.table.toDrive({
  collection: endpointThresholdThermalRows,
  description: 'EndpointThreshold_dLST_RowLevel_AllScenarios_' + RUN_TAG,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'EndpointThreshold_dLST_RowLevel_AllScenarios_' + RUN_TAG,
  fileFormat: 'CSV'
});
