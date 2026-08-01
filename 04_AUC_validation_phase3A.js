// Original GEE script: https://code.earthengine.google.com/3ebd3b59ae4f19f5627e90979921d43c
// Note: link requires Google authentication; full code provided below for archival purposes

// =========================================================
// Phase 3A: Predictive Model Validation (Backcasting, 2000-2012)
// ROC Curve + AUC Calculation for Green Space to Urban Conversion
// Part of: BMA LULC change and green space conversion susceptibility study
// =========================================================

// =================================================================================
// Phase 3A: Predictive Model Validation Workflow (with ROC + AUC)
// =================================================================================

// --- Step 1: Define Inputs and Parameters ---

// 1. LULC Maps (Start and End of Validation Period)
var lulc_2000_asset   = 'projects/ee-agiegizawi/assets/LULC_2000';
var lulc_2012_asset   = 'projects/ee-agiegizawi/assets/LULC_2012';

// 2. Drivers (Must be the 2000 Separated Version)
// Note: this asset intentionally lives under a separate project path (no "ee-" prefix)
var drivers_2000_asset = 'projects/agiegizawi/assets/drivers_2000_separated_v1';

// 3. Study Area
var studyArea          = ee.FeatureCollection('projects/ee-agiegizawi/assets/ADM_BMA_48S');

// Load the images
var lulc_2000  = ee.Image(lulc_2000_asset);
var lulc_2012  = ee.Image(lulc_2012_asset);
var drivers_2000 = ee.Image(drivers_2000_asset);

// VERIFY BANDS: Print to console to ensure you see 'builtup_density' AND 'builtup_proxy'
print('Drivers 2000 Bands:', drivers_2000.bandNames());

// --- Step 2: Define Class Scheme ---
// LULC scheme in assets:
// 0 = Urban/Built-up
// 1 = Water Body
// 2 = Bare Land
// 3 = Agricultural Land
// 4 = Forest

// Redefine Green Space = Agriculture (3) + Forest (4)
var greenSpace2000 = lulc_2000.eq(3).or(lulc_2000.eq(4));
var greenSpace2012 = lulc_2012.eq(3).or(lulc_2012.eq(4));

// Urban/Built-up
var urban2000 = lulc_2000.eq(0);
var urban2012 = lulc_2012.eq(0);

// --- Step 3: Create Ground Truth Change Map ---
// Change: Green Space (2000) → Urban (2012)
var actualChange = greenSpace2000.and(urban2012);

// Stable Green Space: Green in 2000 and still green in 2012
var stableGreenSpace = greenSpace2000.and(greenSpace2012);

// Combine: change = 1, stable = 0
var groundTruth = actualChange.multiply(1).add(stableGreenSpace.multiply(0));

// Mask only original green space areas (focus analysis)
var analysisMask = greenSpace2000;
groundTruth = groundTruth.updateMask(analysisMask).rename('class');

// --- Step 4: Collect Training Data (WITH 70/30 SPLIT) ---
// Merge driver variables with ground truth
var trainingImage = drivers_2000.addBands(groundTruth);

// 1. Stratified sampling for balanced classes
var allSamples = trainingImage.stratifiedSample({
  numPoints: 0, // overridden by classPoints
  classBand: 'class',
  region: studyArea.geometry(),
  scale: 30,
  geometries: true,
  classValues: [0, 1],       // 0 = stable, 1 = change
  classPoints: [1000, 1000]  // equal points per class
}).randomColumn('random');

// 2. Split into Training (70%) and Testing (30%)
var split = 0.7;
var trainSet = allSamples.filter(ee.Filter.lt('random', split));
var testSet  = allSamples.filter(ee.Filter.gte('random', split));

print('Total Samples:', allSamples.size());
print('Training Set (70%):', trainSet.size());
print('Testing Set (30%):', testSet.size());

// --- Step 5: Train Random Forest Classifier in PROBABILITY mode ---
// Note: .bandNames() automatically picks up your new separate layers
var classifier = ee.Classifier.smileRandomForest({numberOfTrees: 100, seed:42})
  .setOutputMode('PROBABILITY').train({
  features: trainSet,
  classProperty: 'class',
  inputProperties: drivers_2000.bandNames()
});

// --- Step 6: Predict the Past Change (Validation) ---
var prediction = drivers_2000.classify(classifier, 'probability');

// --- Step 6.1: Reclassify Probability into 4 Susceptibility Classes ---
// Class 1 (Low):        0.00 - 0.25
// Class 2 (Moderate):   0.25 - 0.50
// Class 3 (High):       0.50 - 0.75
// Class 4 (Very High):  0.75 - 1.00

var susceptibility4 = ee.Image(0).byte()
  .where(prediction.lte(0.25), 1)
  .where(prediction.gt(0.25).and(prediction.lte(0.50)), 2)
  .where(prediction.gt(0.50).and(prediction.lte(0.75)), 3)
  .where(prediction.gt(0.75), 4)
  .rename('susceptibility_class')
  .updateMask(analysisMask);


// --- Step 7: Accuracy Assessment with ROC + AUC (ON TESTING SET ONLY) ---
// Attach probabilities back to validation points
var validation = testSet.classify(classifier, 'probability');

// Split by class
var pos = validation.filter(ee.Filter.eq('class', 1)); // actual change
var neg = validation.filter(ee.Filter.eq('class', 0)); // stable

// Sweep thresholds 0–1
var thresholds = ee.List.sequence(0, 1, 0.02);
var rocPts = ee.FeatureCollection(thresholds.map(function(t) {
  t = ee.Number(t);
  var tp = pos.filter(ee.Filter.gte('probability', t)).size();
  var fn = pos.filter(ee.Filter.lt('probability',  t)).size();
  var fp = neg.filter(ee.Filter.gte('probability', t)).size();
  var tn = neg.filter(ee.Filter.lt('probability',  t)).size();
  
  // Protect against division by zero
  var tpr = ee.Algorithms.If(
    ee.Number(tp).add(fn).gt(0),
    ee.Number(tp).divide(ee.Number(tp).add(fn)),
    0
  );
  var fpr = ee.Algorithms.If(
    ee.Number(fp).add(tn).gt(0),
    ee.Number(fp).divide(ee.Number(fp).add(tn)),
    0
  );
  
  return ee.Feature(null, {threshold: t, TPR: tpr, FPR: fpr});
}));

// Sort and integrate for AUC (trapezoidal rule)
var rocSorted = rocPts.sort('FPR');
var fpr = ee.Array(rocSorted.aggregate_array('FPR'));
var tpr = ee.Array(rocSorted.aggregate_array('TPR'));
var dx  = fpr.slice(0, 1).subtract(fpr.slice(0, 0, -1));
var sy  = tpr.slice(0, 1).add(tpr.slice(0, 0, -1));
var auc = dx.multiply(sy).multiply(0.5).reduce('sum', [0]).abs().toList().get(0);

print('Validation AUC:', auc);
print('ROC points (sample):', rocSorted.limit(10));

// --- Step 8: Visualization ---
Map.centerObject(studyArea, 10);

// Probability map
var palette = ['#ffffcc', '#a1dab4', '#41b6c4', '#2c7fb8', '#253494'];
Map.addLayer(prediction.updateMask(analysisMask),
             {min: 0, max: 1, palette: palette}, 'Prediction Probability (2000-2012)');

// Actual change (reference)
Map.addLayer(actualChange.updateMask(actualChange),
             {palette: ['red']}, 'Actual Change (2000-2012)');
             
Map.addLayer(
  susceptibility4,
  {min: 1, max: 4, palette: ['#a1d99b', '#fee08b', '#fdae61', '#d73027']},
  'Prediction Susceptibility (4 Classes)'
);


// --- Step 9: Add Legend ---
var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
  }
});

var legendTitle = ui.Label({
  value: 'Legend',
  style: {fontWeight: 'bold', fontSize: '14px', margin: '0 0 4px 0'}
});
legend.add(legendTitle);

var makeRow = function(color, name) {
  var colorBox = ui.Label({
    style: {
      backgroundColor: color,
      padding: '8px',
      margin: '0 8px 0 0'
    }
  });
  var description = ui.Label({
    value: name,
    style: {margin: '0'}
  });
  return ui.Panel({
    widgets: [colorBox, description],
    layout: ui.Panel.Layout.Flow('horizontal')
  });
};

legend.add(ui.Label('Prediction Susceptibility'));
legend.add(makeRow('#a1d99b', 'Low (0.00–0.25)'));
legend.add(makeRow('#fee08b', 'Moderate (0.25–0.50)'));
legend.add(makeRow('#fdae61', 'High (0.50–0.75)'));
legend.add(makeRow('#d73027', 'Very High (0.75–1.00)'));
legend.add(ui.Label(' ')); // spacer
legend.add(ui.Label('Actual Change'));
legend.add(makeRow('red', 'Green → Urban (2000–2012)'));
Map.add(legend);



// =================================================================================
// Phase 4: Export
// =================================================================================

Export.image.toDrive({
  image: prediction.clip(studyArea),
  description: 'prediction_probability_2000_2012',
  folder: 'GEE_Exports',
  scale: 30,
  region: studyArea.geometry(),
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image: actualChange.clip(studyArea),
  description: 'actual_change_2000_2012',
  folder: 'GEE_Exports',
  scale: 30,
  region: studyArea.geometry(),
  maxPixels: 1e13
});

Export.table.toDrive({
  collection: rocSorted,
  description: 'ROC_points_2000_2012',
  fileFormat: 'CSV'
});
