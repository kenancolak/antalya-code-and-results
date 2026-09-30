// Antalya revision: memory-safe yearly Landsat core tables and ERA5-Land tables.
//
// Why this version exists
// -----------------------
// The full 2017-2025 graph caused "Too many concurrent aggregations" and
// "User memory limit exceeded". This script preserves the data content but
// exports core statistics and meteorology in separate YEARLY tasks.
//
// Output
// ------
// 1) Nine CoreScene CSV files, one for each year from 2017 through 2025.
// 2) Nine SceneMeteo CSV files, one for each year from 2017 through 2025.
//
// Do not print or preview the output FeatureCollections. Start the batch export
// tasks from the Tasks tab. The tables will be joined later by scene_id in R.

// -----------------------------------------------------------------------------
// 0. Settings and source assets
// -----------------------------------------------------------------------------
var boundaryFC = ee.FeatureCollection('projects/tikkenan/assets/sinir');
var gridFC = ee.FeatureCollection('projects/tikkenan/assets/denemefish5');
var lulc2017Raw = ee.Image('projects/tikkenan/assets/final_lulc2017raw');
var lulc2025Raw = ee.Image('projects/tikkenan/assets/final_lulc2025raw');

var RUN_UID = '20260919_revision_yearsplit_v03';
var EXPORT_FOLDER = 'GEE_Exports';
// First run only the 2017 test. After both 2017 tasks finish successfully,
// replace this line with:
// var EXPORT_YEARS = [2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var EXPORT_YEARS = [2017];
var WRS_PATH = 178;
var WRS_ROW = 34;

var VALID_FRAC17_MIN = 0.99;
var VALID_FRAC25_MIN = 0.99;
var COMMON_FRAC_FULL_MIN = 0.98;
var COMMON_FRAC_UNION_MIN = 0.98;
var CHANGE_THR_DEFAULT = 0.10;
var PURITY_THR_DEFAULT = 0.85;
var CORE_ROW_VALID_MIN = 0.90;

var S0_2017 = ee.List(['2017-06-25', '2017-07-11', '2017-08-12']);
var S0_2025 = ee.List(['2025-06-23', '2025-07-17', '2025-08-10']);
var S1_2017 = ee.List(['2017-06-25', '2017-07-11', '2017-08-12', '2017-08-28']);
var S1_2025 = ee.List(['2025-06-23', '2025-07-17', '2025-08-10', '2025-08-26']);
var S2_2017 = ee.List(['2017-06-25', '2017-07-11', '2017-08-12', '2017-08-28']);
var S2_2025 = ee.List([
  '2025-06-07', '2025-06-15', '2025-06-23', '2025-07-01',
  '2025-07-17', '2025-07-25', '2025-08-10', '2025-08-26'
]);
var S3_2017 = ee.List(['2017-06-25', '2017-07-11', '2017-08-12', '2017-08-28']);
var S3_2025 = ee.List(['2025-06-15', '2025-07-01', '2025-07-17']);

var S0_ALL = S0_2017.cat(S0_2025);
var S1_ALL = S1_2017.cat(S1_2025);
var S2_ALL = S2_2017.cat(S2_2025);
var S3_ALL = S3_2017.cat(S3_2025);

function flag(condition) {
  return ee.Number(ee.Algorithms.If(condition, 1, 0));
}

function shortLabel(code) {
  code = ee.Number(code);
  return ee.String(ee.Algorithms.If(code.eq(1), 'U',
    ee.Algorithms.If(code.eq(2), 'F',
      ee.Algorithms.If(code.eq(3), 'A',
        ee.Algorithms.If(code.eq(4), 'G', 'B')))));
}

// -----------------------------------------------------------------------------
// 1. Reconstruct the canonical 47 UU and 17 FF cores
// -----------------------------------------------------------------------------
var lulc17 = lulc2017Raw.rename('LULC17').toInt16();
var lulc25 = lulc2025Raw.rename('LULC25').toInt16();

var valid17 = lulc17.gt(0).rename('valid17').unmask(0);
var valid25 = lulc25.gt(0).rename('valid25').unmask(0);
var commonValid = valid17.and(valid25).rename('commonValid').unmask(0);
var unionValid = valid17.or(valid25).rename('unionValid').unmask(0);

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

var changed = lulc17.neq(lulc25).and(commonValid).rename('changed').unmask(0);

var coreStatsImage = ee.Image.cat([
  valid17, valid25, commonValid, unionValid,
  u17, f17, a17, g17, b17,
  u25, f25, a25, g25, b25,
  changed
]).toFloat();

var gridStats0 = coreStatsImage.reduceRegions({
  collection: gridFC.filterBounds(boundaryFC),
  reducer: ee.Reducer.mean(),
  scale: 10,
  crs: 'EPSG:32636',
  tileScale: 8
});

function addCoreFields(f) {
  var vals17 = ee.List([
    ee.Number(f.get('u17')), ee.Number(f.get('f17')),
    ee.Number(f.get('a17')), ee.Number(f.get('g17')),
    ee.Number(f.get('b17'))
  ]);
  var vals25 = ee.List([
    ee.Number(f.get('u25')), ee.Number(f.get('f25')),
    ee.Number(f.get('a25')), ee.Number(f.get('g25')),
    ee.Number(f.get('b25'))
  ]);

  var purity17 = ee.Number(vals17.reduce(ee.Reducer.max()));
  var purity25 = ee.Number(vals25.reduce(ee.Reducer.max()));
  var dom17 = ee.Number(vals17.indexOf(purity17)).add(1);
  var dom25 = ee.Number(vals25.indexOf(purity25)).add(1);
  var commonFracFull = ee.Number(f.get('commonValid'));
  var unionFrac = ee.Number(f.get('unionValid'));
  var commonFracUnion = ee.Number(ee.Algorithms.If(
    unionFrac.gt(0), commonFracFull.divide(unionFrac), 0
  ));

  var eligible = ee.Number(f.get('valid17')).gte(VALID_FRAC17_MIN)
    .and(ee.Number(f.get('valid25')).gte(VALID_FRAC25_MIN))
    .and(commonFracFull.gte(COMMON_FRAC_FULL_MIN))
    .and(commonFracUnion.gte(COMMON_FRAC_UNION_MIN));

  var keep = eligible
    .and(dom17.eq(dom25))
    .and(ee.Number(f.get('changed')).lte(CHANGE_THR_DEFAULT))
    .and(purity17.gte(PURITY_THR_DEFAULT))
    .and(purity25.gte(PURITY_THR_DEFAULT));

  var label = ee.String(ee.Algorithms.If(
    keep, shortLabel(dom17).cat(shortLabel(dom17)), 'NA'
  ));

  var bare17 = ee.Number(f.get('b17'));
  var bare25 = ee.Number(f.get('b25'));
  var centroidUtm = f.geometry().centroid(1).transform('EPSG:32636', 1);
  var xy = ee.List(centroidUtm.coordinates());
  var xUtm = ee.Number(xy.get(0));
  var yUtm = ee.Number(xy.get(1));

  function blockId(sizeM) {
    sizeM = ee.Number(sizeM);
    return xUtm.divide(sizeM).floor().format('%d')
      .cat('_').cat(yUtm.divide(sizeM).floor().format('%d'));
  }

  return f.set({
    GridID: ee.String(f.get('system:index')),
    dom17: dom17,
    dom25: dom25,
    purity17: purity17,
    purity25: purity25,
    changeRatio: f.get('changed'),
    urbanFrac17: f.get('u17'),
    urbanFrac25: f.get('u25'),
    forestFrac17: f.get('f17'),
    forestFrac25: f.get('f25'),
    agricultureFrac17: f.get('a17'),
    agricultureFrac25: f.get('a25'),
    greenhouseFrac17: f.get('g17'),
    greenhouseFrac25: f.get('g25'),
    bareFrac17: bare17,
    bareFrac25: bare25,
    bareFracMax: bare17.max(bare25),
    bareFracIncrease: bare25.subtract(bare17),
    keep_default: flag(keep),
    core_default: label,
    x_utm: xUtm,
    y_utm: yUtm,
    block_600: blockId(600),
    block_1000: blockId(1000),
    block_1500: blockId(1500),
    block_2000: blockId(2000)
  });
}

var gridCoreAll = gridStats0.map(addCoreFields);
var analysisCores = gridCoreAll.filter(
  ee.Filter.inList('core_default', ['UU', 'FF'])
);
var uuCores = analysisCores.filter(ee.Filter.eq('core_default', 'UU'));
var ffCores = analysisCores.filter(ee.Filter.eq('core_default', 'FF'));

print('UU count; must be 47:', uuCores.size());
print('FF count; must be 17:', ffCores.size());

// -----------------------------------------------------------------------------
// 2. Landsat scene preparation
// -----------------------------------------------------------------------------
var L8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');
var L9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2');

function prepLandsat(img) {
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 1).eq(0)
    .and(qa.bitwiseAnd(1 << 2).eq(0))
    .and(qa.bitwiseAnd(1 << 3).eq(0))
    .and(qa.bitwiseAnd(1 << 4).eq(0))
    .and(qa.bitwiseAnd(1 << 5).eq(0));

  var lst = img.select('ST_B10')
    .multiply(0.00341802).add(149.0).subtract(273.15)
    .rename('LST_C').updateMask(clear);

  var alignedMask = lst.mask();
  var stQa = img.select('ST_QA').multiply(0.01)
    .rename('ST_QA_K').updateMask(alignedMask);
  var stAtran = img.select('ST_ATRAN').multiply(0.0001)
    .rename('ST_ATRAN').updateMask(alignedMask);
  var stEmis = img.select('ST_EMIS').multiply(0.0001)
    .rename('ST_EMIS').updateMask(alignedMask);
  var stCdist = img.select('ST_CDIST').multiply(0.01)
    .rename('ST_CDIST_km').updateMask(alignedMask);

  var opticalUnsaturated = img.select('QA_RADSAT').bitwiseAnd(
    (1 << 3) | (1 << 4) | (1 << 5)
  ).eq(0);
  var nir = img.select('SR_B5').multiply(0.0000275).add(-0.2);
  var red = img.select('SR_B4').multiply(0.0000275).add(-0.2);
  var swir1 = img.select('SR_B6').multiply(0.0000275).add(-0.2);
  var ndmiDen = nir.add(swir1);
  var ndviDen = nir.add(red);
  var opticalMask = clear.and(opticalUnsaturated).and(alignedMask);

  var ndmi = nir.subtract(swir1).divide(ndmiDen).rename('NDMI')
    .updateMask(ndmiDen.abs().gt(0.0001))
    .updateMask(opticalMask);
  var ndvi = nir.subtract(red).divide(ndviDen).rename('NDVI')
    .updateMask(ndviDen.abs().gt(0.0001))
    .updateMask(opticalMask);

  var landsatValid = alignedMask.rename('landsat_valid').unmask(0).toFloat();
  var opticalValid = opticalMask.rename('optical_valid').unmask(0).toFloat();

  return ee.Image.cat([
    lst, ndmi, ndvi, stQa, stAtran, stEmis, stCdist,
    landsatValid, opticalValid
  ]).copyProperties(img, img.propertyNames());
}

var scenes = L8.merge(L9)
  .filterBounds(boundaryFC)
  .filter(ee.Filter.eq('WRS_PATH', WRS_PATH))
  .filter(ee.Filter.eq('WRS_ROW', WRS_ROW))
  .filter(ee.Filter.eq('PROCESSING_LEVEL', 'L2SP'))
  .filterDate('2017-06-01', '2025-09-01')
  .filter(ee.Filter.calendarRange(6, 8, 'month'))
  .map(prepLandsat)
  .sort('system:time_start');

print('Total JJA scenes; must be 75:', scenes.size());

// -----------------------------------------------------------------------------
// 3. ERA5-Land exact-time preparation
// -----------------------------------------------------------------------------
var ERA = ee.ImageCollection('ECMWF/ERA5_LAND/HOURLY');
var HOUR_MS = 60 * 60 * 1000;
var INSTANT_BANDS = [
  'temperature_2m',
  'dewpoint_temperature_2m',
  'u_component_of_wind_10m',
  'v_component_of_wind_10m',
  'volumetric_soil_water_layer_1',
  'volumetric_soil_water_layer_2',
  'volumetric_soil_water_layer_3',
  'volumetric_soil_water_layer_4'
];

function floorHour(date) {
  date = ee.Date(date);
  return ee.Date(ee.Number(date.millis())
    .divide(HOUR_MS).floor().multiply(HOUR_MS));
}

function eraImageAtHour(date) {
  date = ee.Date(date);
  return ee.Image(ERA.filterDate(date, date.advance(1, 'minute')).first());
}

function interpolateInstantaneous(date) {
  date = ee.Date(date);
  var loDate = floorHour(date);
  var hiDate = loDate.advance(1, 'hour');
  var weight = ee.Number(date.millis()).subtract(loDate.millis()).divide(HOUR_MS);
  var lo = eraImageAtHour(loDate).select(INSTANT_BANDS);
  var hi = eraImageAtHour(hiDate).select(INSTANT_BANDS);
  return lo.add(hi.subtract(lo).multiply(weight));
}

function completedHourlySum(date, hours, bandName) {
  var endHour = floorHour(date);
  var firstEndHour = endHour.advance(
    ee.Number(hours).subtract(1).multiply(-1), 'hour'
  );
  return ERA.filterDate(firstEndHour, endHour.advance(1, 'minute'))
    .select(bandName).sum();
}

function saturationVapourPressure(tempC) {
  return tempC.expression(
    '0.6108 * exp((17.27 * T) / (T + 237.3))', {T: tempC}
  );
}

function weatherImageAt(date) {
  date = ee.Date(date);
  var instant = interpolateInstantaneous(date);
  var t2m = instant.select('temperature_2m').subtract(273.15).rename('t2m_C');
  var td2m = instant.select('dewpoint_temperature_2m')
    .subtract(273.15).rename('td2m_C');
  var rh = saturationVapourPressure(td2m)
    .divide(saturationVapourPressure(t2m)).multiply(100)
    .clamp(0, 100).rename('rh2m_pct');
  var vpd = saturationVapourPressure(t2m)
    .subtract(saturationVapourPressure(td2m)).max(0).rename('vpd_kPa');
  var u10 = instant.select('u_component_of_wind_10m').rename('u10_mps');
  var v10 = instant.select('v_component_of_wind_10m').rename('v10_mps');
  var wind = u10.pow(2).add(v10.pow(2)).sqrt().rename('wind10_mps');
  var sm1 = instant.select('volumetric_soil_water_layer_1').rename('sm_0_7');
  var sm2 = instant.select('volumetric_soil_water_layer_2').rename('sm_7_28');
  var sm3 = instant.select('volumetric_soil_water_layer_3').rename('sm_28_100');
  var sm4 = instant.select('volumetric_soil_water_layer_4').rename('sm_100_289');
  var smRoot100 = sm1.multiply(0.07).add(sm2.multiply(0.21))
    .add(sm3.multiply(0.72)).rename('sm_root_0_100');
  var ssrdPrev3h = completedHourlySum(
    date, 3, 'surface_solar_radiation_downwards_hourly'
  ).divide(3 * 60 * 60).rename('ssrd_prev3h_Wm2');
  var p7 = completedHourlySum(date, 7 * 24, 'total_precipitation_hourly')
    .multiply(1000).rename('precip_prev7d_mm');
  var p30 = completedHourlySum(date, 30 * 24, 'total_precipitation_hourly')
    .multiply(1000).rename('precip_prev30d_mm');

  var localDate = date.advance(3, 'hour');
  var localYear = ee.Number.parse(localDate.format('YYYY', 'UTC'));
  var localMonth = ee.Number.parse(localDate.format('MM', 'UTC'));
  var localDay = ee.Number.parse(localDate.format('dd', 'UTC'));
  var localMidnightUTC = ee.Date.fromYMD(localYear, localMonth, localDay)
    .advance(-3, 'hour');
  var localDayMax = ERA.filterDate(
    localMidnightUTC, localMidnightUTC.advance(1, 'day')
  ).select('temperature_2m').max().subtract(273.15)
    .rename('t2m_localday_max_C');

  return ee.Image.cat([
    t2m, td2m, rh, vpd, u10, v10, wind,
    sm1, sm2, sm3, sm4, smRoot100,
    ssrdPrev3h, p7, p30, localDayMax
  ]).toFloat();
}

// -----------------------------------------------------------------------------
// 4. Normalized scene-level ERA5 rows
// -----------------------------------------------------------------------------
function weatherRowForScene(img) {
  img = ee.Image(img);
  var date = ee.Date(img.get('system:time_start'));
  var dateString = date.format('YYYY-MM-dd', 'UTC');
  var weather = weatherImageAt(date);
  var d = ee.Dictionary(weather.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: boundaryFC.geometry(),
    scale: 11132,
    maxPixels: 1e8,
    tileScale: 4
  }));

  return ee.Feature(null, {
    run_uid: RUN_UID,
    scene_id: img.get('LANDSAT_PRODUCT_ID'),
    landsat_scene_id: img.get('LANDSAT_SCENE_ID'),
    sensor: img.get('SPACECRAFT_ID'),
    date_utc: dateString,
    acquisition_utc: date.format("YYYY-MM-dd'T'HH:mm:ss.SSS", 'UTC'),
    acquisition_trt: date.advance(3, 'hour')
      .format("YYYY-MM-dd'T'HH:mm:ss.SSS", 'UTC'),
    year: ee.Number.parse(date.format('YYYY', 'UTC')),
    doy: ee.Number.parse(date.format('D', 'UTC')),
    month: ee.Number.parse(date.format('MM', 'UTC')),
    scene_cloud_pct: img.get('CLOUD_COVER'),
    land_cloud_pct: img.get('CLOUD_COVER_LAND'),
    sun_elevation_deg: img.get('SUN_ELEVATION'),
    s0_member: flag(S0_ALL.contains(dateString)),
    s1_member: flag(S1_ALL.contains(dateString)),
    s2_current_member: flag(S2_ALL.contains(dateString)),
    s3_member: flag(S3_ALL.contains(dateString)),
    t2m_C_aoi: d.get('t2m_C'),
    td2m_C_aoi: d.get('td2m_C'),
    rh2m_pct_aoi: d.get('rh2m_pct'),
    vpd_kPa_aoi: d.get('vpd_kPa'),
    u10_mps_aoi: d.get('u10_mps'),
    v10_mps_aoi: d.get('v10_mps'),
    wind10_mps_aoi: d.get('wind10_mps'),
    sm_0_7_aoi: d.get('sm_0_7'),
    sm_7_28_aoi: d.get('sm_7_28'),
    sm_28_100_aoi: d.get('sm_28_100'),
    sm_100_289_aoi: d.get('sm_100_289'),
    sm_root_0_100_aoi: d.get('sm_root_0_100'),
    ssrd_prev3h_Wm2_aoi: d.get('ssrd_prev3h_Wm2'),
    precip_prev7d_mm_aoi: d.get('precip_prev7d_mm'),
    precip_prev30d_mm_aoi: d.get('precip_prev30d_mm'),
    t2m_localday_max_C_aoi: d.get('t2m_localday_max_C')
  });
}

function weatherRowsForYear(year) {
  var yearScenes = scenes.filter(ee.Filter.calendarRange(year, year, 'year'));
  var list = yearScenes.toList(yearScenes.size());
  return ee.FeatureCollection(list.map(function(obj) {
    return weatherRowForScene(ee.Image(obj));
  }));
}

// -----------------------------------------------------------------------------
// 5. Normalized core-by-scene Landsat rows
// -----------------------------------------------------------------------------
function coreRowsForScene(img) {
  img = ee.Image(img);
  var date = ee.Date(img.get('system:time_start'));
  var dateString = date.format('YYYY-MM-dd', 'UTC');

  var rows = img.select([
    'LST_C', 'NDMI', 'NDVI',
    'ST_QA_K', 'ST_ATRAN', 'ST_EMIS', 'ST_CDIST_km',
    'landsat_valid', 'optical_valid'
  ]).reduceRegions({
    collection: analysisCores,
    reducer: ee.Reducer.mean(),
    scale: 30,
    crs: 'EPSG:32636',
    tileScale: 8
  });

  return rows.map(function(f) {
    return ee.Feature(null, {
      run_uid: RUN_UID,
      scene_id: img.get('LANDSAT_PRODUCT_ID'),
      sensor: img.get('SPACECRAFT_ID'),
      date_utc: dateString,
      acquisition_utc: date.format("YYYY-MM-dd'T'HH:mm:ss.SSS", 'UTC'),
      year: ee.Number.parse(date.format('YYYY', 'UTC')),
      doy: ee.Number.parse(date.format('D', 'UTC')),
      month: ee.Number.parse(date.format('MM', 'UTC')),
      GridID: f.get('GridID'),
      core_default: f.get('core_default'),
      LST_C: f.get('LST_C'),
      NDMI: f.get('NDMI'),
      NDVI: f.get('NDVI'),
      ST_QA_K: f.get('ST_QA_K'),
      ST_ATRAN: f.get('ST_ATRAN'),
      ST_EMIS: f.get('ST_EMIS'),
      ST_CDIST_km: f.get('ST_CDIST_km'),
      core_scene_valid_frac: f.get('landsat_valid'),
      core_scene_optical_valid_frac: f.get('optical_valid'),
      core_row_valid_min_for_R: CORE_ROW_VALID_MIN,
      purity17: f.get('purity17'),
      purity25: f.get('purity25'),
      changeRatio: f.get('changeRatio'),
      urbanFrac17: f.get('urbanFrac17'),
      urbanFrac25: f.get('urbanFrac25'),
      forestFrac17: f.get('forestFrac17'),
      forestFrac25: f.get('forestFrac25'),
      agricultureFrac17: f.get('agricultureFrac17'),
      agricultureFrac25: f.get('agricultureFrac25'),
      greenhouseFrac17: f.get('greenhouseFrac17'),
      greenhouseFrac25: f.get('greenhouseFrac25'),
      bareFrac17: f.get('bareFrac17'),
      bareFrac25: f.get('bareFrac25'),
      bareFracMax: f.get('bareFracMax'),
      bareFracIncrease: f.get('bareFracIncrease'),
      x_utm: f.get('x_utm'),
      y_utm: f.get('y_utm'),
      block_600: f.get('block_600'),
      block_1000: f.get('block_1000'),
      block_1500: f.get('block_1500'),
      block_2000: f.get('block_2000'),
      s0_member: flag(S0_ALL.contains(dateString)),
      s1_member: flag(S1_ALL.contains(dateString)),
      s2_current_member: flag(S2_ALL.contains(dateString)),
      s3_member: flag(S3_ALL.contains(dateString))
    });
  });
}

function coreRowsForYear(year) {
  var yearScenes = scenes.filter(ee.Filter.calendarRange(year, year, 'year'));
  var list = yearScenes.toList(yearScenes.size());
  var nested = list.map(function(obj) {
    return coreRowsForScene(ee.Image(obj)).toList(analysisCores.size());
  });
  return ee.FeatureCollection(ee.List(nested).flatten());
}

// -----------------------------------------------------------------------------
// 6. Create independent yearly export tasks
// -----------------------------------------------------------------------------
EXPORT_YEARS.forEach(function(year) {
  Export.table.toDrive({
    collection: weatherRowsForYear(year),
    description: 'Antalya_SceneMeteo_' + year + '_' + RUN_UID,
    folder: EXPORT_FOLDER,
    fileNamePrefix: 'Antalya_SceneMeteo_' + year + '_' + RUN_UID,
    fileFormat: 'CSV'
  });

  Export.table.toDrive({
    collection: coreRowsForYear(year),
    description: 'Antalya_CoreScene_' + year + '_' + RUN_UID,
    folder: EXPORT_FOLDER,
    fileNamePrefix: 'Antalya_CoreScene_' + year + '_' + RUN_UID,
    fileFormat: 'CSV'
  });
});

print('Created two tasks for each year listed in EXPORT_YEARS.');
print('Do not print the task collections. Start exports from the Tasks tab.');
