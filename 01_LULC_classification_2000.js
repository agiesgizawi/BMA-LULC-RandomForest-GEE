// Original GEE script: https://code.earthengine.google.com/c339e392c0012ff4bcd2b2b51c24492a
// Note: link requires Google authentication; full code provided below for archival purposes

// =========================================================
// LULC Classification 2000 - Bandung Metropolitan Area
// Object-Based Image Analysis (OBIA) with Random Forest
// Landsat 7 Surface Reflectance, Collection 2
// Part of: BMA LULC change and green space conversion susceptibility study
// =========================================================

// ------------------------
// Define study area style
// ------------------------
var study_area_vis = {
  color: "ff0000",
  fillColor: "00000000",
  width: 2
};
Map.addLayer(study_area.style(study_area_vis));
Map.centerObject(study_area, 10);

// ------------------------
// Define year list
// ------------------------
var startDate = '2002-01-01';
var endDate = '2002-12-31';

// ------------------------
// Load Landsat 7 collection
// ------------------------
var l7_sr = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'); // Surface Reflectance
var l7_toa = ee.ImageCollection('LANDSAT/LE07/C02/T1_TOA'); // Top-of-Atmosphere

// ------------------------
// Cloud masking and scaling functions for SR
// ------------------------
// Apply scaling factors to Landsat 7 SR data
function applyScaleFactors(image) {
  var opticalBands = image.select('SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7').multiply(0.0000275).add(-0.2);
  // The QA_PIXEL band is not scaled
  return image.addBands(opticalBands, null, true);
}

// Mask clouds and cloud shadows using QA_PIXEL (Collection 2)
function maskL7srClouds(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 5);
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
               .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask);
}

// ------------------------
// Pan-sharpening function (for TOA data)
// ------------------------
function panSharpen(image, rgbBands, panBand) {
  var rgb = image.select(['B3','B2','B1']);
  var pan = image.select('B8');
  var hsv = rgb.rgbToHsv();
  var intensity = hsv.select('value');
  var blended = intensity.multiply(0.5).add(pan.multiply(0.5)); // 50% blend
  var panSharpened = ee.Image.cat(
    hsv.select('hue'),
    hsv.select('saturation'),
    blended
  ).hsvToRgb();
  return panSharpened;
}

// ------------------------
// Filter Landsat 7 collections
// ------------------------
var filter = ee.Filter.and(
  ee.Filter.bounds(study_area),
  ee.Filter.date(startDate, endDate),
  ee.Filter.lt('CLOUD_COVER', 50)
);

var filteredL7_sr = l7_sr.filter(filter);
var filteredL7_toa = l7_toa.filter(filter);

// ------------------------
// Create median composite from SR data for analysis
// ------------------------
var l7compositeMasked = filteredL7_sr
    .map(applyScaleFactors)
    .map(maskL7srClouds)
    .median()
    .clip(study_area);

// ------------------------
// Select and pan-sharpen the clearest TOA scene
// ------------------------
var bestImage_toa = filteredL7_toa.sort('CLOUD_COVER').first();
var sharpened = panSharpen(bestImage_toa, ['B3', 'B2', 'B1'], 'B8').clip(study_area);

// ------------------------
// Function to calculate cloud cover % per image
// ------------------------
var calcCloudCover = function(image) {
  var ones = ee.Image.constant(1).clip(study_area).rename('constant');

  var totalPixels = ones.updateMask(image.select('QA_PIXEL'))
                        .reduceRegion({
                          reducer: ee.Reducer.sum(),
                          geometry: study_area,
                          scale: 30,
                          maxPixels: 1e9
                        });

  var maskedImage = maskL7srClouds(image);

  var clearPixels = ones.updateMask(maskedImage.select('QA_PIXEL'))
                        .reduceRegion({
                          reducer: ee.Reducer.sum(),
                          geometry: study_area,
                          scale: 30,
                          maxPixels: 1e9
                        });

  var total = ee.Number(totalPixels.get('constant'));
  var clear = ee.Number(clearPixels.get('constant'));
  var cloudPct = total.subtract(clear).divide(total).multiply(100);

  return image.set('cloud_percentage', cloudPct);
};

// Apply to all images in the collection
var collectionWithCloud = l7_sr.filterBounds(study_area)
                            .filterDate(startDate, endDate)
                            .filterMetadata('CLOUD_COVER', 'less_than', 50)
                            .map(calcCloudCover);

print('Cloud cover % per image:', collectionWithCloud.aggregate_array('cloud_percentage'));

var avgCloudCover = collectionWithCloud.aggregate_mean('cloud_percentage');
print('Average cloud cover % across all images:', avgCloudCover);

// Convert image collection to a FeatureCollection with cloud percentage, for export
var cloudStats = collectionWithCloud.map(function(image) {
  var date = ee.Date(image.get('system:time_start')).format('YYYY-MM-dd');
  var cloudPct = image.get('cloud_percentage');
  return ee.Feature(null, {
    'date': date,
    'cloud_percentage': cloudPct
  });
});

print('Cloud stats table:', cloudStats);

Export.table.toDrive({
  collection: cloudStats,
  description: 'Landsat_Cloud_Percentage',
  fileFormat: 'CSV'
});

// ------------------------
// Band map (uses the SR composite for analysis)
// ------------------------
var bandMap = {
  NIR: l7compositeMasked.select('SR_B4'),
  SWIR: l7compositeMasked.select('SR_B5'),
  RED: l7compositeMasked.select('SR_B3'),
  GREEN: l7compositeMasked.select('SR_B2'),
  BLUE: l7compositeMasked.select('SR_B1')
};

// ------------------------
// Spectral indices: NDVI, NDBI, MNDWI, BSI
// ------------------------
var ndvi = l7compositeMasked.expression('(NIR - RED) / (NIR + RED)', bandMap).rename('NDVI');
var ndbi = l7compositeMasked.expression('(SWIR - NIR) / (SWIR + NIR)', bandMap).rename('NDBI');
var mndwi = l7compositeMasked.expression('(GREEN - SWIR) / (GREEN + SWIR)', bandMap).rename('MNDWI');
var bsi = l7compositeMasked.expression('((SWIR + RED) - (NIR + BLUE)) / ((SWIR + RED) + (NIR + BLUE))', bandMap).rename('BSI');

// ------------------------
// Visualizations
// ------------------------
var L7_SR_viz = {bands:['SR_B3', 'SR_B2', 'SR_B1'], min:[0], max:[0.3]};

Map.addLayer(l7compositeMasked, L7_SR_viz, 'Landsat 7 SR Composite');
Map.addLayer(sharpened, {min: 0, max: 0.3}, 'Pan-sharpened RGB');
Map.addLayer(ndvi, {min: -0.2, max: 1, palette: ['red', 'yellow', 'green']}, 'NDVI 2000');
Map.addLayer(ndbi, {min: -1, max: 1, palette: ['blue', 'white', 'red']}, 'NDBI 2000');
Map.addLayer(mndwi, {min: -1, max: 1, palette: ['red', 'white', 'blue']}, 'MNDWI 2000');
Map.addLayer(bsi, {min: -1, max: 1, palette: ['white', 'brown', 'black']}, 'BSI 2000');

// ------------------------
// OBIA Step: SNIC segmentation (uses the SR composite)
// ------------------------
var imageForSegmentation = l7compositeMasked
                           .select(['B3', 'B2', 'B1'])
                           .addBands(ndvi);

var seeds = ee.Algorithms.Image.Segmentation.seedGrid(2, 'square');

var segment = ee.Algorithms.Image.Segmentation.SNIC({
  image: l7compositeMasked,
  connectivity: 4,
  compactness: 0,
  seeds: seeds
});

Map.addLayer(segment.select('clusters').randomVisualizer(), {}, 'Clusters', false);
Map.addLayer(segment, {bands: ['SR_B3_mean', 'SR_B2_mean', 'SR_B1_mean'], min: 0, max: 0.3}, 'OBIA');

// Image for classification
var imageObject = segment.select(['SR_B.*_mean']);
var bandsName = imageObject.bandNames();

// ------------------------
// Merge training samples
// ------------------------
var trainingSamples = builtup.merge(water)
                             .merge(bare)
                             .merge(agricultural_land)
                             .merge(forest);

// ------------------------
// Sample training data from the object-based image
// ------------------------
var obiaSamples = imageObject.sampleRegions({
  collection: trainingSamples,
  properties: ['landcover'],
  scale: 30,
  tileScale: 2,
  geometries: true
}).randomColumn('random', 100);

var trainingSet = obiaSamples.filter(ee.Filter.lt('random', 0.8));
var testingSet = obiaSamples.filter(ee.Filter.gte('random', 0.8));

// ------------------------
// Train Random Forest classifier (100 trees)
// ------------------------
var obiaClassifier = ee.Classifier.smileRandomForest(100).train({
  features: trainingSet,
  classProperty: 'landcover',
  inputProperties: bandsName
});

// ------------------------
// Accuracy assessment
// ------------------------
var validated = testingSet.classify(obiaClassifier);
var cm = validated.errorMatrix('landcover', 'classification');
print('OBIA Confusion Matrix:', cm);
print('OBIA Overall Accuracy:', cm.accuracy());
print('OBIA Kappa Coefficient:', cm.kappa());

// ------------------------
// Classify the object-based image
// ------------------------
var obiaClassified2000 = imageObject.classify(obiaClassifier);

// ------------------------
// Visualization
// ------------------------
var classNames = [
  'Built-up', 'Water body', 'Bare land',
  'Agricultural land', 'Forest'
];

var classPalette = [
  'd7191c', // 0 - Built-up
  '2b83ba', // 1 - Water body
  'fdae61', // 2 - Bare land
  'ffffbf', // 3 - Agricultural land
  '008000'  // 4 - Forest
];

Map.addLayer(obiaClassified2000, {
  min: 0,
  max: 4,
  palette: classPalette
}, 'OBIA LULC 2000');

// ------------------------
// Export the final 2000 classification map
// ------------------------
Export.image.toAsset({
  image: obiaClassified2000,
  description: 'LULC_Map_2000',
  assetId: 'users/agiegizawi/LULC_2000',
  scale: 30,
  region: study_area,
  maxPixels: 1e13
});

// ------------------------
// Export Landsat 7 composite image to Google Drive
// ------------------------
Export.image.toDrive({
  image: l7compositeMasked,
  description: 'Landsat7_SR_Composite_2000',
  folder: 'GEE_Exports',
  scale: 30,
  region: study_area,
  fileFormat: 'GeoTIFF',
  maxPixels: 1e13
});

// ------------------------
// Legend UI
// ------------------------
var legend = ui.Panel({ style: { position: 'bottom-left', padding: '8px' } });
legend.add(ui.Label('Land Cover Classification', { fontWeight: 'bold' }));

for (var i = 0; i < classNames.length; i++) {
  var colorBox = ui.Label('', {
    backgroundColor: '#' + classPalette[i],
    padding: '8px',
    margin: '0 0 4px 0',
    width: '20px'
  });
  var description = ui.Label(classNames[i], { margin: '0 0 4px 6px' });
  legend.add(ui.Panel([colorBox, description], ui.Panel.Layout.Flow('horizontal')));
}
Map.add(legend);
