// Original GEE script: https://code.earthengine.google.com/a2b742eea2a3c97ec34fa5e10deb0609
// Note: link requires Google authentication; full code provided below for archival purposes

// =========================================================
// Phase 3B: Primary Predictive Model Training + Validation (2012-2024)
// Trains the primary susceptibility model and validates it with ROC + AUC
// Part of: BMA LULC change and green space conversion susceptibility study
// =========================================================

// =================================================================================
// Phase 3B: Primary Predictive Model Training Workflow with Validation (AUC + Export)
// =================================================================================

// --- Step 1: Define Inputs and Parameters ---
var lulc_2012_asset = 'projects/ee-agiegizawi/assets/LULC_2012';
var lulc_2024_asset = 'projects/ee-agiegizawi/assets/LULC_2024';

// Points to the separated driver asset (2012)
var drivers_2012_asset = 'projects/agiegizawi/assets/drivers_2012_separated_v1';
var studyArea = ee.FeatureCollection('projects/ee-agiegizawi/assets/ADM_BMA_48S');

var lulc_2012 = ee.Image(lulc_2012_asset);
var lulc_2024 = ee.Image(lulc_2024_asset);
var drivers_2012 = ee.Image(drivers_2012_asset);

// Verify bands in console (Crucial check)
print('Drivers 2012 Bands:', drivers_2012.bandNames());

// Class values
var URBAN_CLASS = 0;
var AGRICULTURE_CLASS = 3;
var FOREST_CLASS = 4;

// --- Step 2: Create Binary Ground Truth Map ---
// We want to learn: What changed from Green (2012) to Urban (2024)?
var greenSpace2012 = lulc_2012.eq(AGRICULTURE_CLASS).or(lulc_2012.eq(FOREST_CLASS));
var urban2024 = lulc_2024.eq(URBAN_CLASS);

var actualChange = greenSpace2012.and(urban2024);
// Stable Green: Was green in 2012 AND did NOT become urban in 2024
var stableGreenSpace = greenSpace2012.and(greenSpace2012.and(urban2024).not());

// 1 = Change, 0 = Stable
var groundTruth = actualChange.unmask(0).where(stableGreenSpace, 0);
var analysisMask = greenSpace2012; // We only care about areas that started as Green
groundTruth = groundTruth.updateMask(analysisMask);

// =================================================================================
// Step 3: Define Band Names and Prepare Training Data
// =================================================================================

// Select the specific bands from the separated driver asset
var PREDICTOR_BANDS = [
  'elevation', 
  'slope', 
  'distance_to_roads', 
  'builtup_density',
  'builtup_proxy',
  'population'
];

// Select the predictors
var drivers_selected = drivers_2012.select(PREDICTOR_BANDS);
var trainingImage = drivers_selected.addBands(groundTruth.rename('class'));

// 1. Stratified sampling (Define 'allSamples' here)
var allSamples = trainingImage.stratifiedSample({
  numPoints: 2000,                // Default number (ignored due to classPoints)
  classBand: 'class',
  region: studyArea.geometry(),
  scale: 30,
  geometries: true,
  classValues: [0, 1],            // 0 = Stable, 1 = Change
  classPoints: [1000, 1000]       // Equal sampling
}).randomColumn('random');        // Add random column for splitting

// 2. Split into Training (70%) and Testing (30%)
var split = 0.7;
var trainSet = allSamples.filter(ee.Filter.lt('random', split));
var testSet  = allSamples.filter(ee.Filter.gte('random', split));

print('Total Samples:', allSamples.size());
print('Training Samples (70%):', trainSet.size());
print('Testing Samples (30%):', testSet.size());

// =================================================================================
// Step 4: Train Random Forest Classifier
// =================================================================================

var primaryClassifier = ee.Classifier.smileRandomForest({numberOfTrees: 100}).train({
  features: trainSet,        
  classProperty: 'class',
  inputProperties: PREDICTOR_BANDS 
});

// =================================================================================
// Step 5: Export Trained Classifier 
// =================================================================================

var assetId = 'projects/ee-agiegizawi/assets/primary_prediction_model_v3_separated'; 

Export.classifier.toAsset({
  classifier: primaryClassifier,
  description: 'Export_Primary_RF_Model_v3_Separated',
  assetId: assetId
});
print('Export task created for classifier:', assetId);

// =================================================================================
// Step 6: AUC Analysis (Validation of the Model)
// =================================================================================

// Train a separate Probability Mode classifier for validation purposes
// (kept independent from the primary classifier above so validation uses
// its own probability-mode training run, consistent with the Phase 3A method)
var probClassifier = ee.Classifier.smileRandomForest({numberOfTrees: 100, seed:42})
  .setOutputMode('PROBABILITY')
  .train({
    features: trainSet,
    classProperty: 'class',
    inputProperties: PREDICTOR_BANDS
  });

// Classify training data to get probabilities
var probTest = testSet.classify(probClassifier, 'probability'); 
// Note: In PROBABILITY mode, the output band is usually named 'classification' or 'probability' depending on GEE version.
// Using default behavior here.

// Split by class 
var pos = probTest.filter(ee.Filter.eq('class', 1)); // actual change
var neg = probTest.filter(ee.Filter.eq('class', 0)); // stable

// Sweep thresholds 0–1, step 0.02 (matches 3A)
var thresholds = ee.List.sequence(0, 1, 0.02);
var rocPts = ee.FeatureCollection(thresholds.map(function(t) {
  t = ee.Number(t);
  var tp = pos.filter(ee.Filter.gte('probability', t)).size();
  var fn = pos.filter(ee.Filter.lt('probability',  t)).size();
  var fp = neg.filter(ee.Filter.gte('probability', t)).size();
  var tn = neg.filter(ee.Filter.lt('probability',  t)).size();

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

// Sort and integrate for AUC (trapezoidal rule) — matches 3A method
var rocSorted = rocPts.sort('FPR');
var fpr = ee.Array(rocSorted.aggregate_array('FPR'));
var tpr = ee.Array(rocSorted.aggregate_array('TPR'));
var dx  = fpr.slice(0, 1).subtract(fpr.slice(0, 0, -1));
var sy  = tpr.slice(0, 1).add(tpr.slice(0, 0, -1));
var auc = dx.multiply(sy).multiply(0.5).reduce('sum', [0]).abs().toList().get(0);

print('Primary Model AUC:', auc);
print('ROC points (sample):', rocSorted.limit(10));

// =================================================================================
// --- Step 7: Generate and Visualize Prediction Probability Map ---
// =================================================================================

var probabilityMap = drivers_selected
  .classify(probClassifier)
  .updateMask(analysisMask);

var probVis = {min: 0, max: 1, palette: ['#006400', '#FFFF00', '#FF0000']}; // Green(Safe) -> Red(High Risk)
var groundTruthVis = {min: 0, max: 1, palette: ['#228B22', '#FF00FF']}; // Green(Stable) -> Magenta(Change)

Map.centerObject(studyArea, 10);
Map.addLayer(lulc_2024.clip(studyArea), {min:0, max:4, palette:['#FF0000', '#C2B280', '#FFFF00', '#006400', '#A0522D']}, 'LULC 2024');
Map.addLayer(groundTruth, groundTruthVis, 'Ground Truth (2012-2024)');
Map.addLayer(probabilityMap, probVis, 'Urban Pressure Probability 2024');

// =================================================================================
// --- Step 8: Legends ---
// =================================================================================

function createLegend(title, palette, names) {
  var legend = ui.Panel({style: {position: 'bottom-left', padding: '8px 15px', border: '1px solid #DDDDDD'}});
  legend.add(ui.Label({value: title, style: {fontWeight: 'bold', fontSize: '18px', margin: '0 0 4px 0'}}));
  
  for (var i = 0; i < names.length; i++) {
    var colorBox = ui.Label({style: {backgroundColor: palette[i], padding: '8px', margin: '0 0 4px 0', border: '1px solid black'}});
    var description = ui.Label({value: names[i], style: {margin: '0 0 4px 6px'}});
    legend.add(ui.Panel({widgets: [colorBox, description], layout: ui.Panel.Layout.Flow('horizontal')}));
  }
  return legend;
}

Map.add(createLegend('Change Classes', ['#228B22', '#FF00FF'], ['Stable Green', 'Converted to Urban']));
