// Antalya revision audit: Landsat scene inventory, exact-time ERA5-Land
// meteorology, Landsat forest-moisture/quality diagnostics, and a Sentinel-2
// annual stability/vegetation-support recheck (2017-2025 JJA).
//
// Purpose
//   1. Reconstruct the same 47 UU and 17 FF endpoint cores as the canonical
//      thermal script (changeRatio <= 0.10; purity >= 0.85).
//   2. Audit every Landsat 8/9 L2SP scene on WRS-2 path/row 178/034 in JJA.
//   3. Export exact acquisition-time meteorology and per-core LST/NDMI/NDVI.
//   4. Recompute, but do not silently replace, the historical AOI-clear rule.
//   5. Recheck annual spectral stability and vegetation support using SCL.
//
// Important interpretation
//   ERA5-Land is regional background forcing (~11 km), not 300 m microclimate.
//   Landsat NDMI is a canopy/vegetation-water proxy, not soil moisture.
//   The script does not select scenes using LST values.
//   The Sentinel-2 recheck deliberately does not use QA60 because the official
//   Earth Engine catalogue documents a QA60 masking gap from 2022-01-25 through
//   2024-02-28. SCL is used consistently for every year instead.

// Expected assets (same paths as the canonical Antalya scripts)
var boundaryFC = ee.FeatureCollection('projects/tikkenan/assets/sinir');
var gridFC = ee.FeatureCollection('projects/tikkenan/assets/denemefish5');
var lulc2017Raw = ee.Image('projects/tikkenan/assets/final_lulc2017raw');
var lulc2025Raw = ee.Image('projects/tikkenan/assets/final_lulc2025raw');

var RUN_UID = '20260920_revision_exact_v03';
var EXPORT_FOLDER = 'GEE_Exports';
var WRS_PATH = 178;
var WRS_ROW = 34;
var HISTORICAL_AOI_CLEAR_MIN = 0.95;
var MODEL_CORE_CLEAR_MIN = 0.95;
var CORE_ROW_VALID_MIN = 0.90; // applied later in R; exported here for audit

var VALID_FRAC17_MIN = 0.99;
var VALID_FRAC25_MIN = 0.99;
var COMMON_FRAC_FULL_MIN = 0.98;
var COMMON_FRAC_UNION_MIN = 0.98;
var CHANGE_THR_DEFAULT = 0.10;
var PURITY_THR_DEFAULT = 0.85;

// Sentinel-2 annual recheck. The same JJA window and scene-level cloud ceiling
// are used in every year. NDVI >= 0.20 retains the manuscript's vegetation-
// support definition; the revision replaces the single-year IVS by the median
// of annual support fractions, calculated later in R.
var S2_CLOUD_PCT_MAX = 40;
var VEGETATION_NDVI_THR = 0.20;
var YEARS = ee.List.sequence(2017, 2025);
var DIFF_YEARS = ee.List.sequence(2017, 2024);
var ANNUAL_PERCENTILES = ee.List([90, 95, 97]);
var ANNUAL_DEFAULT_PCT = 95;
var ANNUAL_UNSTABLE_RATIO_THR = 0.10;
var ANNUAL_FULL_SUPPORT_MIN = 0.90;
// Full scene/core previews trigger many simultaneous reducers in the Code
// Editor. Leave false for routine verification; exports are unaffected.
var PRINT_HEAVY_PREVIEWS = false;

// Current canonical scenario definitions. These are audit labels, not filters.
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

// -----------------------------------------------------------------------------
// 1. Reconstruct the endpoint core population exactly as in 02_endpoint_landsat_dlst.js
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
  tileScale: 4
});

function shortLabel(code) {
  code = ee.Number(code);
  return ee.String(ee.Algorithms.If(code.eq(1), 'U',
    ee.Algorithms.If(code.eq(2), 'F',
      ee.Algorithms.If(code.eq(3), 'A',
        ee.Algorithms.If(code.eq(4), 'G', 'B')))));
}

function flag(condition) {
  return ee.Number(ee.Algorithms.If(condition, 1, 0));
}

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

  // Reviewer-facing contamination diagnostics. These are determined only by
  // the two endpoint maps, before any thermal outcome is examined.
  var bare17 = ee.Number(f.get('b17'));
  var bare25 = ee.Number(f.get('b25'));
  var bareMax = bare17.max(bare25);
  var bareIncrease = bare25.subtract(bare17);

  // Fixed spatial identifiers for inference. The same UTM origin and block
  // membership are used for every date; neither depends on LST.
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
    bareFracMax: bareMax,
    bareFracIncrease: bareIncrease,
    commonFrac_full: commonFracFull,
    commonFrac_union: commonFracUnion,
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
var analysisCores = gridCoreAll.filter(ee.Filter.inList('core_default', ['UU', 'FF']));
var uuCores = analysisCores.filter(ee.Filter.eq('core_default', 'UU'));
var ffCores = analysisCores.filter(ee.Filter.eq('core_default', 'FF'));
var coreGeometry = analysisCores.geometry();
var ffCentroids = ffCores.map(function(f) {
  return ee.Feature(f.geometry().centroid(1), {GridID: f.get('GridID')});
});

print('Reconstructed UU count (must be 47):', uuCores.size());
print('Reconstructed FF count (must be 17):', ffCores.size());

// -----------------------------------------------------------------------------
// 2. Landsat preparation: canonical LST mask + aligned optical diagnostics
// -----------------------------------------------------------------------------
var L8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');
var L9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2');

function prepLandsat(img) {
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 1).eq(0) // dilated cloud
    .and(qa.bitwiseAnd(1 << 2).eq(0))    // cirrus
    .and(qa.bitwiseAnd(1 << 3).eq(0))    // cloud
    .and(qa.bitwiseAnd(1 << 4).eq(0))    // cloud shadow
    .and(qa.bitwiseAnd(1 << 5).eq(0));   // snow/ice

  var lst = img.select('ST_B10')
    .multiply(0.00341802).add(149.0).subtract(273.15)
    .rename('LST_C')
    .updateMask(clear);

  // USGS Collection 2 surface-temperature diagnostic layers. These do not
  // alter LST; they allow us to test whether endpoint/group differences could
  // be associated with retrieval uncertainty, atmospheric transmittance,
  // emissivity, or proximity to cloud.
  var stQa = img.select('ST_QA').multiply(0.01)
    .rename('ST_QA_K').updateMask(clear);
  var stAtran = img.select('ST_ATRAN').multiply(0.0001)
    .rename('ST_ATRAN').updateMask(clear);
  var stEmis = img.select('ST_EMIS').multiply(0.0001)
    .rename('ST_EMIS').updateMask(clear);
  var stCdist = img.select('ST_CDIST').multiply(0.01)
    .rename('ST_CDIST_km').updateMask(clear);

  // Collection 2 SR scaling must precede ratio calculation because of offset.
  // Saturated optical pixels are excluded from NDMI/NDVI, while the canonical
  // thermal mask itself remains unchanged.
  var opticalUnsaturated = img.select('QA_RADSAT').bitwiseAnd(
    (1 << 3) | (1 << 4) | (1 << 5)
  ).eq(0);
  var nir = img.select('SR_B5').multiply(0.0000275).add(-0.2);
  var red = img.select('SR_B4').multiply(0.0000275).add(-0.2);
  var swir1 = img.select('SR_B6').multiply(0.0000275).add(-0.2);
  var ndmiDen = nir.add(swir1);
  var ndviDen = nir.add(red);
  var ndmi = nir.subtract(swir1).divide(ndmiDen).rename('NDMI')
    .updateMask(ndmiDen.abs().gt(0.0001));
  var ndvi = nir.subtract(red).divide(ndviDen).rename('NDVI')
    .updateMask(ndviDen.abs().gt(0.0001));

  // Optical support is a quality-controlled subset of the thermal/QA support.
  var alignedMask = lst.mask();
  var opticalMask = clear.and(opticalUnsaturated).and(alignedMask);
  ndmi = ndmi.updateMask(opticalMask).updateMask(ndmi.gte(-1).and(ndmi.lte(1)));
  ndvi = ndvi.updateMask(opticalMask).updateMask(ndvi.gte(-1).and(ndvi.lte(1)));
  var landsatValid = alignedMask.rename('landsat_valid').unmask(0).toFloat();
  var opticalValid = opticalMask.rename('optical_valid').unmask(0).toFloat();

  stQa = stQa.updateMask(alignedMask);
  stAtran = stAtran.updateMask(alignedMask);
  stEmis = stEmis.updateMask(alignedMask);
  stCdist = stCdist.updateMask(alignedMask);

  return ee.Image.cat([
    lst, ndmi, ndvi, stQa, stAtran, stEmis, stCdist,
    landsatValid, opticalValid
  ])
    .copyProperties(img, img.propertyNames());
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

print('All L8/L9 JJA candidate scenes, 2017-2025:', scenes.size());

// -----------------------------------------------------------------------------
// 3. Exact-time ERA5-Land helpers
// -----------------------------------------------------------------------------
var ERA = ee.ImageCollection('ECMWF/ERA5_LAND/HOURLY');
var HOUR_MS = 60 * 60 * 1000;
var INSTANT_BANDS = [
  'temperature_2m', 'dewpoint_temperature_2m',
  'u_component_of_wind_10m', 'v_component_of_wind_10m',
  'volumetric_soil_water_layer_1',
  'volumetric_soil_water_layer_2',
  'volumetric_soil_water_layer_3',
  'volumetric_soil_water_layer_4'
];

function floorHour(date) {
  date = ee.Date(date);
  return ee.Date(ee.Number(date.millis()).divide(HOUR_MS).floor().multiply(HOUR_MS));
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

// Sum the last N complete ERA5-Land hourly increments ending at floor(t).
// This is strictly antecedent: no post-acquisition precipitation is included.
function completedHourlySum(date, hours, bandName) {
  var endHour = floorHour(date);
  var firstEndHour = endHour.advance(ee.Number(hours).subtract(1).multiply(-1), 'hour');
  return ERA.filterDate(firstEndHour, endHour.advance(1, 'minute'))
    .select(bandName).sum();
}

function saturationVapourPressure(tempC) {
  // Use tempC as the expression's primary image so the ERA5-Land projection is
  // retained (rather than inheriting ee.Image.constant's default projection).
  return tempC.expression(
    '0.6108 * exp((17.27 * T) / (T + 237.3))', {T: tempC}
  );
}

function weatherImageAt(date) {
  date = ee.Date(date);
  var instant = interpolateInstantaneous(date);
  var t2m = instant.select('temperature_2m').subtract(273.15).rename('t2m_C');
  var td2m = instant.select('dewpoint_temperature_2m').subtract(273.15).rename('td2m_C');
  var rh = saturationVapourPressure(td2m)
    .divide(saturationVapourPressure(t2m)).multiply(100)
    .clamp(0, 100).rename('rh2m_pct');
  var vpd = saturationVapourPressure(t2m)
    .subtract(saturationVapourPressure(td2m)).max(0).rename('vpd_kPa');
  var u10 = instant.select('u_component_of_wind_10m').rename('u10_mps');
  var v10 = instant.select('v_component_of_wind_10m').rename('v10_mps');
  var wind = u10.pow(2)
    .add(v10.pow(2))
    .sqrt().rename('wind10_mps');

  var sm1 = instant.select('volumetric_soil_water_layer_1').rename('sm_0_7');
  var sm2 = instant.select('volumetric_soil_water_layer_2').rename('sm_7_28');
  var sm3 = instant.select('volumetric_soil_water_layer_3').rename('sm_28_100');
  var sm4 = instant.select('volumetric_soil_water_layer_4').rename('sm_100_289');
  var smRoot100 = sm1.multiply(0.07).add(sm2.multiply(0.21))
    .add(sm3.multiply(0.72)).rename('sm_root_0_100');

  // Earth Engine supplies an hourly disaggregation of the ERA5-Land
  // accumulated radiation field. Use the three fully completed antecedent
  // hours; this avoids both the daily reset and post-acquisition information.
  var ssrdPrev3h = completedHourlySum(
    date, 3, 'surface_solar_radiation_downwards_hourly'
  ).divide(3 * 60 * 60).rename('ssrd_prev3h_Wm2');

  var p7 = completedHourlySum(date, 7 * 24, 'total_precipitation_hourly')
    .multiply(1000).rename('precip_prev7d_mm');
  var p30 = completedHourlySum(date, 30 * 24, 'total_precipitation_hourly')
    .multiply(1000).rename('precip_prev30d_mm');

  // Local-day maximum (TRT = UTC+3 throughout the study period).
  var localDate = date.advance(3, 'hour');
  var localYear = ee.Number.parse(localDate.format('YYYY', 'UTC'));
  var localMonth = ee.Number.parse(localDate.format('MM', 'UTC'));
  var localDay = ee.Number.parse(localDate.format('dd', 'UTC'));
  var localMidnightUTC = ee.Date.fromYMD(localYear, localMonth, localDay)
    .advance(-3, 'hour');
  var localDayMax = ERA.filterDate(localMidnightUTC, localMidnightUTC.advance(1, 'day'))
    .select('temperature_2m').max().subtract(273.15).rename('t2m_localday_max_C');

  return ee.Image.cat([
    t2m, td2m, rh, vpd, u10, v10, wind, sm1, sm2, sm3, sm4, smRoot100,
    ssrdPrev3h, p7, p30, localDayMax
  ]).toFloat();
}

function regionalMean(image, geometry) {
  return image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: geometry,
    scale: 11132,
    maxPixels: 1e8,
    tileScale: 2
  });
}

// A 300 m FF polygon is smaller than an ERA5-Land pixel and can contain no
// coarse-grid pixel centre. Sample the containing ERA pixel at each FF centroid
// and then take an equal-core mean; do not pretend this is 300 m soil moisture.
function equalPointMean(image, points, bandNames) {
  bandNames = ee.List(bandNames);
  var sampled = image.select(bandNames).sampleRegions({
    collection: points,
    properties: ['GridID'],
    scale: 11132,
    geometries: false,
    tileScale: 2
  });
  var values = bandNames.map(function(b) {
    return sampled.aggregate_mean(ee.String(b));
  });
  return ee.Dictionary.fromLists(bandNames, values);
}

function landsatGroupMean(image, features, bandNames) {
  return ee.Dictionary(image.select(bandNames).reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: features.geometry(),
    scale: 30,
    crs: 'EPSG:32636',
    maxPixels: 1e8,
    tileScale: 4
  }));
}

// -----------------------------------------------------------------------------
// 4. One audit row per scene
// -----------------------------------------------------------------------------
function sceneAuditFeature(img) {
  img = ee.Image(img);
  var date = ee.Date(img.get('system:time_start'));
  var dateString = date.format('YYYY-MM-dd', 'UTC');
  var year = ee.Number.parse(date.format('YYYY', 'UTC'));
  var doy = ee.Number.parse(date.format('D', 'UTC'));

  var aoiClear = ee.Number(img.select('landsat_valid').reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: boundaryFC.geometry(),
    scale: 30,
    crs: 'EPSG:32636',
    maxPixels: 1e9,
    tileScale: 4
  }).get('landsat_valid'));

  var coreClear = ee.Number(img.select('landsat_valid').reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: coreGeometry,
    scale: 30,
    crs: 'EPSG:32636',
    maxPixels: 1e8,
    tileScale: 4
  }).get('landsat_valid'));

  var weather = weatherImageAt(date);
  var aoiWeather = ee.Dictionary(regionalMean(weather, boundaryFC.geometry()));
  var forestWeather = equalPointMean(
    weather, ffCentroids,
    ['sm_0_7', 'sm_7_28', 'sm_28_100', 'sm_100_289', 'sm_root_0_100']
  );
  var diagnosticBands = [
    'LST_C', 'NDMI', 'NDVI', 'ST_QA_K', 'ST_ATRAN', 'ST_EMIS',
    'ST_CDIST_km'
  ];
  var uuLandsat = landsatGroupMean(img, uuCores, diagnosticBands);
  var ffLandsat = landsatGroupMean(img, ffCores, diagnosticBands);

  var currentS2 = flag(S2_ALL.contains(dateString));
  var recomputedS2 = flag(aoiClear.gte(HISTORICAL_AOI_CLEAR_MIN));

  return ee.Feature(null, {
    run_uid: RUN_UID,
    scene_id: img.get('LANDSAT_PRODUCT_ID'),
    landsat_scene_id: img.get('LANDSAT_SCENE_ID'),
    system_index: img.get('system:index'),
    sensor: img.get('SPACECRAFT_ID'),
    date_utc: dateString,
    acquisition_utc: date.format("YYYY-MM-dd'T'HH:mm:ss.SSS", 'UTC'),
    acquisition_trt: date.advance(3, 'hour').format("YYYY-MM-dd'T'HH:mm:ss.SSS", 'UTC'),
    year: year,
    doy: doy,
    month: ee.Number.parse(date.format('MM', 'UTC')),
    wrs_path: img.get('WRS_PATH'),
    wrs_row: img.get('WRS_ROW'),
    scene_cloud_pct: img.get('CLOUD_COVER'),
    land_cloud_pct: img.get('CLOUD_COVER_LAND'),
    sun_elevation_deg: img.get('SUN_ELEVATION'),
    aoi_clear_frac: aoiClear,
    core_clear_frac: coreClear,
    historical_aoi_clear_min: HISTORICAL_AOI_CLEAR_MIN,
    model_core_clear_min: MODEL_CORE_CLEAR_MIN,
    s0_member: flag(S0_ALL.contains(dateString)),
    s1_member: flag(S1_ALL.contains(dateString)),
    s2_current_member: currentS2,
    s2_recomputed_aoi95: recomputedS2,
    s2_membership_mismatch: flag(currentS2.neq(recomputedS2)),
    s3_member: flag(S3_ALL.contains(dateString)),
    model_scene_keep_core95: flag(coreClear.gte(MODEL_CORE_CLEAR_MIN)),
    t2m_C_aoi: aoiWeather.get('t2m_C'),
    td2m_C_aoi: aoiWeather.get('td2m_C'),
    rh2m_pct_aoi: aoiWeather.get('rh2m_pct'),
    vpd_kPa_aoi: aoiWeather.get('vpd_kPa'),
    u10_mps_aoi: aoiWeather.get('u10_mps'),
    v10_mps_aoi: aoiWeather.get('v10_mps'),
    wind10_mps_aoi: aoiWeather.get('wind10_mps'),
    sm_0_7_aoi: aoiWeather.get('sm_0_7'),
    sm_7_28_aoi: aoiWeather.get('sm_7_28'),
    sm_28_100_aoi: aoiWeather.get('sm_28_100'),
    sm_100_289_aoi: aoiWeather.get('sm_100_289'),
    sm_root_0_100_aoi: aoiWeather.get('sm_root_0_100'),
    ssrd_prev3h_Wm2_aoi: aoiWeather.get('ssrd_prev3h_Wm2'),
    precip_prev7d_mm_aoi: aoiWeather.get('precip_prev7d_mm'),
    precip_prev30d_mm_aoi: aoiWeather.get('precip_prev30d_mm'),
    t2m_localday_max_C_aoi: aoiWeather.get('t2m_localday_max_C'),
    sm_0_7_ff_background: forestWeather.get('sm_0_7'),
    sm_7_28_ff_background: forestWeather.get('sm_7_28'),
    sm_28_100_ff_background: forestWeather.get('sm_28_100'),
    sm_100_289_ff_background: forestWeather.get('sm_100_289'),
    sm_root_0_100_ff_background: forestWeather.get('sm_root_0_100'),
    lst_C_uu: uuLandsat.get('LST_C'),
    lst_C_ff: ffLandsat.get('LST_C'),
    ndmi_uu: uuLandsat.get('NDMI'),
    ndmi_ff: ffLandsat.get('NDMI'),
    ndvi_uu: uuLandsat.get('NDVI'),
    ndvi_ff: ffLandsat.get('NDVI'),
    st_qa_K_uu: uuLandsat.get('ST_QA_K'),
    st_qa_K_ff: ffLandsat.get('ST_QA_K'),
    st_atran_uu: uuLandsat.get('ST_ATRAN'),
    st_atran_ff: ffLandsat.get('ST_ATRAN'),
    st_emis_uu: uuLandsat.get('ST_EMIS'),
    st_emis_ff: ffLandsat.get('ST_EMIS'),
    st_cdist_km_uu: uuLandsat.get('ST_CDIST_km'),
    st_cdist_km_ff: ffLandsat.get('ST_CDIST_km')
  });
}

var sceneList = scenes.toList(scenes.size());
var sceneAudit = ee.FeatureCollection(sceneList.map(function(obj) {
  return sceneAuditFeature(ee.Image(obj));
}));

var endpointAudit = sceneAudit.filter(ee.Filter.inList('year', [2017, 2025]));
if (PRINT_HEAVY_PREVIEWS) {
  print('Endpoint scene audit:', endpointAudit.sort('acquisition_utc'));
  print(
    'S2 hard-coded versus recomputed AOI-clear mismatches (must be inspected):',
    endpointAudit.filter(ee.Filter.eq('s2_membership_mismatch', 1)).size()
  );
}

// -----------------------------------------------------------------------------
// 5. One row per retained UU/FF core and scene
// -----------------------------------------------------------------------------
function coreRowsForScene(img) {
  img = ee.Image(img);
  var date = ee.Date(img.get('system:time_start'));
  var dateString = date.format('YYYY-MM-dd', 'UTC');
  var rows = img.select([
    'LST_C', 'NDMI', 'NDVI', 'ST_QA_K', 'ST_ATRAN', 'ST_EMIS',
    'ST_CDIST_km', 'landsat_valid', 'optical_valid'
  ]).reduceRegions({
    collection: analysisCores,
    reducer: ee.Reducer.mean(),
    scale: 30,
    crs: 'EPSG:32636',
    tileScale: 4
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

var nestedCoreRows = sceneList.map(function(obj) {
  return coreRowsForScene(ee.Image(obj)).toList(analysisCores.size());
});
var coreSceneRows = ee.FeatureCollection(ee.List(nestedCoreRows).flatten());

if (PRINT_HEAVY_PREVIEWS) {
  print('Core-scene row count before R validity filtering:', coreSceneRows.size());
  print('Core-scene preview:', coreSceneRows.limit(10));
}

// -----------------------------------------------------------------------------
// 6. Sentinel-2 annual vegetation support and spectral-stability recheck
// -----------------------------------------------------------------------------
var S2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
var S2_REFLECTANCE_BANDS = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12'];

function maskS2WithSCL(image) {
  var scl = image.select('SCL');

  // Conservatively retain only dark-area, vegetation, bare/urban and water
  // classes. Exclude no-data, saturated pixels, cloud shadow, SCL class 7
  // (low-probability cloud/unclassified), medium/high cloud, cirrus and snow.
  var clear = scl.eq(2)
    .or(scl.eq(4))
    .or(scl.eq(5))
    .or(scl.eq(6));

  return image.select(S2_REFLECTANCE_BANDS)
    .multiply(0.0001)
    .updateMask(clear)
    .copyProperties(image, ['system:time_start', 'CLOUDY_PIXEL_PERCENTAGE']);
}

function addS2AnnualIndices(image) {
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var ndbi = image.normalizedDifference(['B11', 'B8']).rename('NDBI');
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

  return image.addBands([ndvi, ndbi, tcb]);
}

function annualS2Source(year) {
  year = ee.Number(year);
  var start = ee.Date.fromYMD(year, 6, 1);
  var end = start.advance(3, 'month');
  return S2
    .filterBounds(boundaryFC)
    .filterDate(start, end)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', S2_CLOUD_PCT_MAX))
    .map(maskS2WithSCL)
    .map(addS2AnnualIndices);
}

function makeAnnualS2Image(year) {
  year = ee.Number(year);
  var source = annualS2Source(year);
  var med = source.median().clipToCollection(boundaryFC);
  var ndvi = med.select('NDVI');
  var vegSupport = ndvi.gte(VEGETATION_NDVI_THR)
    .rename('veg_support').updateMask(ndvi.mask());
  var valid = ndvi.mask().rename('s2_valid').unmask(0).toFloat();
  var nObs = source.select('B4').count().rename('s2_nobs').toFloat();

  return ee.Image.cat([
    med.select(['NDVI', 'NDBI', 'TCB']),
    vegSupport.toFloat(), valid, nObs
  ]).set({
    year: year,
    source_scene_count: source.size(),
    cloud_mask_method: 'SCL_2_4_5_6',
    vegetation_ndvi_threshold: VEGETATION_NDVI_THR
  });
}

var annualS2 = ee.ImageCollection.fromImages(YEARS.map(makeAnnualS2Image));
print('Annual S2 collection (SCL recheck):', annualS2);

function vegetationRowsForYear(year) {
  year = ee.Number(year);
  var image = ee.Image(annualS2.filter(ee.Filter.eq('year', year)).first());
  var rows = image.select([
    'NDVI', 'veg_support', 's2_valid', 's2_nobs'
  ]).reduceRegions({
    collection: analysisCores,
    reducer: ee.Reducer.mean(),
    scale: 10,
    crs: 'EPSG:32636',
    tileScale: 4
  });

  return rows.map(function(f) {
    var support = ee.Number(f.get('veg_support'));
    return ee.Feature(null, {
      run_uid: RUN_UID,
      GridID: f.get('GridID'),
      core_default: f.get('core_default'),
      year: year,
      annual_ndvi_mean: f.get('NDVI'),
      vegetation_support_frac: support,
      IVS_year: ee.Number(1).subtract(support),
      s2_valid_frac: f.get('s2_valid'),
      s2_valid_observations_mean: f.get('s2_nobs'),
      source_scene_count: image.get('source_scene_count'),
      cloud_mask_method: image.get('cloud_mask_method'),
      vegetation_ndvi_threshold: VEGETATION_NDVI_THR,
      x_utm: f.get('x_utm'),
      y_utm: f.get('y_utm')
    });
  });
}

var nestedVegetationRows = YEARS.map(function(year) {
  return vegetationRowsForYear(year).toList(analysisCores.size());
});
var annualVegetationRows = ee.FeatureCollection(
  ee.List(nestedVegetationRows).flatten()
);
print('Annual core vegetation-support rows:', annualVegetationRows.size());

var annualDiffs = ee.ImageCollection.fromImages(DIFF_YEARS.map(function(year) {
  year = ee.Number(year);
  var current = ee.Image(annualS2.filter(ee.Filter.eq('year', year)).first());
  var next = ee.Image(annualS2.filter(ee.Filter.eq('year', year.add(1))).first());
  return next.select(['NDVI', 'NDBI', 'TCB'])
    .subtract(current.select(['NDVI', 'NDBI', 'TCB']))
    .abs()
    .rename(['dNDVI', 'dNDBI', 'dTCB'])
    .set({start_year: year, end_year: year.add(1)});
}));

var annualPairCount = annualDiffs.select('dNDVI').count()
  .rename('annual_pair_count').toFloat();
var annualFullSupportMask = annualPairCount.eq(DIFF_YEARS.size());
var annualFullSupport = annualFullSupportMask
  .rename('annual_full_support').unmask(0).toFloat();
var annualMaxJump = annualDiffs.max().rename([
  'maxdNDVI', 'maxdNDBI', 'maxdTCB'
]).updateMask(annualFullSupportMask);
var anchorCores = gridCoreAll.filter(ee.Filter.eq('keep_default', 1));
// Use every fully supported 10 m anchor pixel. Do not draw a random sample:
// the exact revision calibration used 74,431 pixels and direct empirical
// percentiles, so sampling can move borderline cells between configurations.
var anchorGeometry = anchorCores.geometry();
var exactThresholds = annualMaxJump.reduceRegion({
  // maxRaw exceeds the 74,431-pixel calibration population, forcing direct
  // percentile calculation rather than Earth Engine's histogram approximation.
  reducer: ee.Reducer.percentile(
    [90, 95, 97], ['p90', 'p95', 'p97'], null, null, 100000
  ),
  geometry: anchorGeometry,
  scale: 10,
  crs: 'EPSG:32636',
  maxPixels: 10000000,
  tileScale: 4
});

var anchorPixelCount = annualMaxJump.select('maxdNDVI').reduceRegion({
  reducer: ee.Reducer.count(),
  geometry: anchorGeometry,
  scale: 10,
  crs: 'EPSG:32636',
  maxPixels: 10000000,
  tileScale: 4
}).get('maxdNDVI');

var annualThresholdTable = ee.FeatureCollection(ANNUAL_PERCENTILES.map(function(p) {
  p = ee.Number(p);
  var key = ee.String('p').cat(p.format('%.0f'));
  return ee.Feature(null, {
    run_uid: RUN_UID,
    cloud_mask_method: 'SCL_2_4_5_6',
    anchor_percentile: p,
    anchor_valid_pixel_count: anchorPixelCount,
    threshold_method: 'all_valid_anchor_pixels_direct_percentile',
    threshold_maxdNDVI: exactThresholds.get(
      ee.String('maxdNDVI_').cat(key)
    ),
    threshold_maxdNDBI: exactThresholds.get(
      ee.String('maxdNDBI_').cat(key)
    ),
    threshold_maxdTCB: exactThresholds.get(
      ee.String('maxdTCB_').cat(key)
    )
  });
}));

function unstableImageForPercentile(p, bandName) {
  var row = ee.Feature(
    annualThresholdTable.filter(ee.Filter.eq('anchor_percentile', p)).first()
  );
  var unstable = annualMaxJump.select('maxdNDVI')
    .gt(ee.Number(row.get('threshold_maxdNDVI')))
    .or(annualMaxJump.select('maxdNDBI')
      .gt(ee.Number(row.get('threshold_maxdNDBI'))))
    .or(annualMaxJump.select('maxdTCB')
      .gt(ee.Number(row.get('threshold_maxdTCB'))));
  return unstable.rename(bandName)
    .updateMask(annualFullSupportMask).toFloat();
}

var unstableP90 = unstableImageForPercentile(90, 'unstableRatio_p90');
var unstableP95 = unstableImageForPercentile(95, 'unstableRatio_p95');
var unstableP97 = unstableImageForPercentile(97, 'unstableRatio_p97');

var annualStabilityStack = ee.Image.cat([
  annualMaxJump, unstableP90, unstableP95, unstableP97,
  annualPairCount, annualFullSupport
]);

var annualCoreStats0 = annualStabilityStack.reduceRegions({
  collection: analysisCores,
  reducer: ee.Reducer.mean(),
  scale: 10,
  crs: 'EPSG:32636',
  tileScale: 4
});

var annualCoreStats = annualCoreStats0.map(function(f) {
  var isUrban = ee.String(f.get('core_default')).compareTo('UU').eq(0);
  var enoughAnnualSupport = ee.Number(f.get('annual_full_support'))
    .gte(ANNUAL_FULL_SUPPORT_MIN);
  var stableP90 = isUrban.and(enoughAnnualSupport).and(
    ee.Number(f.get('unstableRatio_p90')).lte(ANNUAL_UNSTABLE_RATIO_THR)
  );
  var stableP95 = isUrban.and(enoughAnnualSupport).and(
    ee.Number(f.get('unstableRatio_p95')).lte(ANNUAL_UNSTABLE_RATIO_THR)
  );
  var stableP97 = isUrban.and(enoughAnnualSupport).and(
    ee.Number(f.get('unstableRatio_p97')).lte(ANNUAL_UNSTABLE_RATIO_THR)
  );

  return ee.Feature(null, f.toDictionary(f.propertyNames())).set({
    run_uid: RUN_UID,
    cloud_mask_method: 'SCL_2_4_5_6',
    annual_unstable_ratio_threshold: ANNUAL_UNSTABLE_RATIO_THR,
    annual_full_support_min: ANNUAL_FULL_SUPPORT_MIN,
    annual_stable_urban_p90: flag(stableP90),
    annual_stable_urban_p95: flag(stableP95),
    annual_stable_urban_p97: flag(stableP97)
  });
});

print('Annual SCL thresholds:', annualThresholdTable);
print(
  'SCL p95-confirmed UU count:',
  annualCoreStats.filter(ee.Filter.eq('annual_stable_urban_p95', 1)).size()
);
print(
  'SCL p90-confirmed UU count (strict sensitivity):',
  annualCoreStats.filter(ee.Filter.eq('annual_stable_urban_p90', 1)).size()
);
print(
  'SCL p97-confirmed UU count:',
  annualCoreStats.filter(ee.Filter.eq('annual_stable_urban_p97', 1)).size()
);

// -----------------------------------------------------------------------------
// 7. Exports
// -----------------------------------------------------------------------------
Export.table.toDrive({
  collection: sceneAudit.sort('acquisition_utc'),
  description: 'Antalya_JJA_Scene_Meteo_Audit_' + RUN_UID,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'Antalya_JJA_Scene_Meteo_Audit_' + RUN_UID,
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: coreSceneRows,
  description: 'Antalya_JJA_CoreScene_LST_NDMI_QA_' + RUN_UID,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'Antalya_JJA_CoreScene_LST_NDMI_QA_' + RUN_UID,
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: annualVegetationRows,
  description: 'Antalya_Annual_Core_VegetationSupport_' + RUN_UID,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'Antalya_Annual_Core_VegetationSupport_' + RUN_UID,
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: annualThresholdTable,
  description: 'Antalya_Annual_Stability_SCL_Thresholds_' + RUN_UID,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'Antalya_Annual_Stability_SCL_Thresholds_' + RUN_UID,
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: annualCoreStats,
  description: 'Antalya_Annual_Stability_SCL_CoreStats_' + RUN_UID,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'Antalya_Annual_Stability_SCL_CoreStats_' + RUN_UID,
  fileFormat: 'CSV'
});
