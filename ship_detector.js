// ------------------------------------------------------------
// Bellingcat-style ship detection (ROI-based) + SDT-style visualization
// - Water forced to dark blue
// - Ships forced to bright white
// - Green vector detections filtered by minimum length
// ------------------------------------------------------------

// 0) INPUTS (EDIT THESE)
var TARGET_DATE = ee.Date('2025-07-10');   // target date (UTC)
var PICK_WITHIN_DAYS = 10;                 // widen if your area has sparse S1 coverage
var PASS = null;                          // 'ASCENDING' | 'DESCENDING' | null

// Water mask quality is critical for SDT-like output.
// Start high (80–95) for ports/coasts; lower only if deep offshore water.
var WATER_OCC_MIN = 0;                   // 0–100 (JRC occurrence)

// SDT-like simple intensity threshold on VV (in dB; COPERNICUS/S1_GRD is dB in GEE)
var THRESHOLD_DB = 0;                     // ships are typically very bright; 0 is a common start but you can also go below 0
var MIN_LENGTH_M = 50;                    // tune - min ship length (meters)

// Mask cleanup / smoothing
var COAST_ERODE_PX = 2;                   // erode water mask to avoid shoreline false positives
var MORPH_RADIUS_PX = 2;                  // smooth detections (0 off, 1–2 mild)
var MIN_PIXELS = 5;                       // remove small speckle blobs

// 1) GEOMETRY (ROI is your drawn polygon)
var area = ee.FeatureCollection(geometry).geometry().simplify(10).buffer(1);

Map.setOptions('HYBRID');
Map.centerObject(area, 13);

var outlineArea = ee.Image().paint(ee.FeatureCollection([ee.Feature(area)]), 1, 2);
Map.addLayer(outlineArea, {}, 'GEOMETRY');

// 2) WATER MASK (JRC Global Surface Water occurrence)   [oai_citation:2‡Google for Developers](https://developers.google.com/earth-engine/datasets/catalog/JRC_GSW1_4_GlobalSurfaceWater)
var gsw = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence');
var waterRaw = gsw.gte(WATER_OCC_MIN).clip(area);

// Erode away from coastlines to remove mixed land/water pixels (reduces false positives)
var water = waterRaw.focal_min(COAST_ERODE_PX);

Map.addLayer(waterRaw.selfMask(), {palette: ['0000ff']}, 'Water mask raw (debug)', false);
Map.addLayer(water.selfMask(), {palette: ['00aaff']}, 'Water mask eroded (debug)', false);

// 3) SENTINEL-1 COLLECTION (S1_GRD in GEE is already in dB)  [oai_citation:3‡bellingcat](https://www.bellingcat.com/resources/2023/05/11/peering-beyond-the-clouds-a-guide-to-bellingcats-ship-detection-tool/)
var start = TARGET_DATE.advance(-PICK_WITHIN_DAYS, 'day');
var end   = TARGET_DATE.advance(PICK_WITHIN_DAYS + 1, 'day');

var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(area)
  .filterDate(start, end)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.eq('productType', 'GRD'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'));

if (PASS !== null) s1 = s1.filter(ee.Filter.eq('orbitProperties_pass', PASS));

print('Candidates in ± window:', s1.size());
print('Candidate times:', s1.aggregate_array('system:time_start'));

// 4) PICK CLOSEST ACQUISITION (avoid empty-day problem)
var s1Closest = s1.map(function(im) {
  var d = ee.Number(im.date().difference(TARGET_DATE, 'day')).abs();
  return im.set('diffDays', d);
}).sort('diffDays');

// Guard against empty collection using evaluate()
s1Closest.size().evaluate(function(n) {
  if (n === 0) {
    print('No Sentinel-1 acquisitions found within ±' + PICK_WITHIN_DAYS +
          ' days of ' + TARGET_DATE.format('YYYY-MM-dd').getInfo() +
          '. Increase PICK_WITHIN_DAYS, remove PASS, or verify ROI.');
    return;
  }

  var img = ee.Image(s1Closest.first()).clip(area);
  var acqDate = img.date();
  print('Chosen acquisition date:', acqDate);

  // 5) VV masked to water
  var vv = img.select('VV').updateMask(water);

  // 6) SHIP MASK (binary) = VV > threshold
  var shipMask = vv.gt(THRESHOLD_DB).selfMask();

  // Optional mild morphology to tidy edges before vectorizing
  if (MORPH_RADIUS_PX > 0) {
    shipMask = shipMask.focal_max(MORPH_RADIUS_PX).focal_min(MORPH_RADIUS_PX); // closing
    shipMask = shipMask.focal_min(MORPH_RADIUS_PX).focal_max(MORPH_RADIUS_PX); // opening
  }

  // Remove small speckle blobs
  var ccount = shipMask.connectedPixelCount(100, true);
  shipMask = shipMask.updateMask(ccount.gte(MIN_PIXELS));

  // ------------------------------------------------------------
  // A) SDT-STYLE VISUALIZATION LAYER (THIS IS THE KEY CHANGE)
  // Water is forced to dark blue, ships forced to white.
  // ------------------------------------------------------------

  // Paint water as a constant dark colour (independent of VV stretch)
  var waterBlueRGB = ee.Image.constant(1)
    .updateMask(water)
    .visualize({palette: ['081d58']}); // dark blue

  // Paint ships as white on top
  var shipsWhiteRGB = ee.Image.constant(1)
    .updateMask(shipMask)
    .visualize({palette: ['ffffff']}); // white

  // Blend into a single RGB layer
  var sdtView = waterBlueRGB.blend(shipsWhiteRGB);

  Map.addLayer(sdtView, {}, 'SDT-style view (blue water, white ships)', true);

  // Optional: show VV for context (off by default; can look “grey” depending on sea state)
  // Keep it OFF unless you need context.
  Map.addLayer(vv.resample('bicubic'), {min: -30, max: -5}, 'VV (dB) context', false);

  // ------------------------------------------------------------
  // B) VECTOR DETECTIONS (green points), filtered by length proxy
  // ------------------------------------------------------------

  // Vectorize polygons (needed for size estimation)
  var polys = shipMask.reduceToVectors({
    geometry: area,
    scale: 10,
    geometryType: 'polygon',
    eightConnected: true,
    bestEffort: true,
    maxPixels: 1e8
  });

  // Length proxy = bounding-box diagonal (meters)
  function bboxDiagonalMeters(geom) {
    var b = geom.bounds(1);
    var coords = ee.List(b.coordinates().get(0));
    var p0 = ee.List(coords.get(0));
    var p2 = ee.List(coords.get(2));
    return ee.Geometry.Point(p0).distance(ee.Geometry.Point(p2));
  }

  // Convert polygons to centroid points with attributes
  var ships = polys.map(function(f) {
    var g = f.geometry();
    var len = bboxDiagonalMeters(g);
    return ee.Feature(g.centroid(1)).set({
      acq_date: acqDate.format('YYYY-MM-dd'),
      length_m: len
    });
  });

  var shipsFiltered = ships.filter(ee.Filter.gte('length_m', MIN_LENGTH_M));

  print('Ship candidates (raw):', ships.size());
  print('Ship candidates (>= ' + MIN_LENGTH_M + ' m):', shipsFiltered.size());
  print('Candidates sample:', shipsFiltered.limit(20));

  Map.addLayer(shipsFiltered, {color: '39ff14'}, 'Ship detections (points)', true);
});