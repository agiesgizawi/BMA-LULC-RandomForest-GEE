// Original GEE script: https://code.earthengine.google.com/c339e392c0012ff4bcd2b2b51c24492a
// Note: link requires Google authentication; full code provided below for archival purposes

// =========================================================
// Susceptibility Projection to 2100 - SSP3 vs SSP5
// Trains RF probability model on 2012-2024 observed change,
// then projects green space conversion susceptibility under
// SSP3 and SSP5 population scenarios for 2100
// Part of: BMA LULC change and green space conversion susceptibility study
// =========================================================

// =================================================================================
// MASTER SCRIPT: Train Probability Model & Predict Future 2100 (SSP3 vs SSP5)
// =================================================================================

// --- Step 1: Define ALL Assets ---

// Study Area & LULC
var studyArea = ee.FeatureCollection('projects/ee-agiegizawi/assets/ADM_BMA_48S');
var lulc_2012 = ee.Image('projects/ee-agiegizawi/assets/LULC_2012');
var lulc_2024 = ee.Image('projects/ee-agiegizawi/assets/LULC_2024');

// Drivers (Separated Versions)
var drivers_2012 = ee.Image('projects/agiegizawi/assets/drivers_2012_separated_v1');
var drivers_2024 = ee.Image('projects/agiegizawi/assets/drivers_2024_separated_v1');

// Future Population Assets for 2100
var pop_assets = {
  '2100_ssp3': 'projects/ee-agiegizawi/assets/FuturePop_SSP3_2100_1km_v0_2',
  '2100_ssp5': 'projects/ee-agiegizawi/assets/FuturePop_SSP5_2100_1km_v0_2'
};

// --- Step 2: Prepare Training Data (The "Teacher") ---

// Define bands to use for prediction
var PREDICTOR_BANDS = [
  'elevation', 'slope', 'distance_to_roads', 
  'builtup_density', 'builtup_proxy', 'population'
];

// 1. Create Ground Truth (What happened 2012->2024?)
var greenSpace2012 = lulc_2012.eq(3).or(lulc_2012.eq(4)); // Ag(3) + Forest(4)
var urban2024 = lulc_2024.eq(0); // Urban(0)
var actualChange = greenSpace2012.and(urban2024); // 1 = Changed
var stableGreen = greenSpace2012.and(urban2024.not()); // 0 = Stable

var groundTruth = actualChange.unmask(0)
    .where(stableGreen, 0)
    .updateMask(greenSpace2012)
    .rename('class');

// 2. Sample Points
var trainingInput = drivers_2012.select(PREDICTOR_BANDS).addBands(groundTruth);
var trainingData = trainingInput.stratifiedSample({
  numPoints: 2000, 
  classBand: 'class', 
  region: studyArea.geometry(), 
  scale: 30, 
  geometries: true,
  classValues: [0, 1], 
  classPoints: [1000, 1000] // Balanced training
});

// --- Step 3: Train Random Forest (PROBABILITY MODE) ---

var classifier = ee.Classifier.smileRandomForest({numberOfTrees: 100})
    .setOutputMode('PROBABILITY') 
    .train({
      features: trainingData,
      classProperty: 'class',
      inputProperties: PREDICTOR_BANDS
    });

print('Model trained successfully in PROBABILITY mode.');

// --- Step 4: Analyze Driver Importance (Optional but recommended) ---

var explanation = classifier.explain();
var importanceDict = ee.Dictionary(explanation.get('importance'));

var importanceFC = ee.FeatureCollection(
  importanceDict.keys().map(function(key){
    return ee.Feature(null, {
      'driver': key,
      'importance': importanceDict.get(key)
    });
  })
);

var importanceChart = ui.Chart.feature.byFeature(
    importanceFC, 'driver', 'importance'
  )
  .setChartType('ColumnChart')
  .setOptions({
    title: 'Driver Importance',
    hAxis: {title: 'Drivers'},
    vAxis: {title: 'Importance Score'},
    legend: {position: 'none'},
    colors: ['#1b9e77']
  });

print(importanceChart);

// --- Step 5: Prepare Future Logic ---

var greenSpaceMask = lulc_2024.eq(3).or(lulc_2024.eq(4));
var projection = lulc_2024.projection();

Map.centerObject(studyArea, 10);
Map.addLayer(lulc_2024.clip(studyArea), {min:0, max:4, palette:['#ff0000','#2b83ba','#fdae61','#ffff00','#008000']}, 'LULC 2024 Base', false);

var studyAreaVis = {color: "ff0000", fillColor: "00000000", width: 2};
Map.addLayer(studyArea.style(studyAreaVis), {}, 'Study Area Boundary');

// --- Step 6: Loop Prediction for Each 2100 Scenario ---

var targetScenarios = ['2100_ssp3', '2100_ssp5'];

targetScenarios.forEach(function(scenario) {
  
  // A. Get Population
  var popAsset = pop_assets[scenario];
  var populationFuture = ee.Image(popAsset)
    .resample('bilinear')
    .reproject({crs: projection, scale: 1000})
    .rename('population')
    .clip(studyArea);

  // B. Build Stack
  var staticBands = [
    'elevation',
    'slope',
    'distance_to_roads', 
    'builtup_density',
    'builtup_proxy'
  ];

  var drivers_future = ee.Image.cat([
    drivers_2024.select(staticBands),
    populationFuture
  ]);

  // C. Predict Probability (0.0 to 1.0)
  var susceptibility = drivers_future
    .classify(classifier)
    .updateMask(greenSpaceMask)
    .rename('susceptibility_' + scenario);

  // D. Visualization (4-color ramp)
  Map.addLayer(susceptibility,
    {min:0, max:1, palette:['green','yellow','orange','red']},
    'Susceptibility ' + scenario
  );

  // E. Export Map Image (GeoTIFF)
  Export.image.toDrive({
    image: susceptibility,
    description: 'Prediction_Susceptibility_' + scenario, 
    folder: 'GEE_Exports',
    scale: 30,
    region: studyArea.geometry(),
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  });
  
  // ============================================================
  // G. Calculate 4-Class Statistics & Export CSV
  // ============================================================
  
  // 1. Reclassify Probability into 4 Integer Classes
  // Class 1 (Low):       0.00 - 0.25
  // Class 2 (Moderate):  0.25 - 0.50
  // Class 3 (High):      0.50 - 0.75
  // Class 4 (Very High): 0.75 - 1.00
  
  var susceptibilityClass = ee.Image(0).byte()
    .where(susceptibility.lte(0.25), 1)
    .where(susceptibility.gt(0.25).and(susceptibility.lte(0.50)), 2)
    .where(susceptibility.gt(0.50).and(susceptibility.lte(0.75)), 3)
    .where(susceptibility.gt(0.75), 4)
    .rename('class_label')
    .updateMask(susceptibility.mask()); 

  // 2. Calculate Area for each class
  var areaImage = ee.Image.pixelArea().addBands(susceptibilityClass);
  
  var classStats = areaImage.reduceRegion({
    reducer: ee.Reducer.sum().group({
      groupField: 1,
      groupName: 'class_label',
    }),
    geometry: studyArea.geometry(),
    scale: 30,
    maxPixels: 1e13
  });
  
  // 3. Convert List to FeatureCollection for Export
  var classStatsList = ee.List(classStats.get('groups'));
  
  var exportFeatures = ee.FeatureCollection(classStatsList.map(function(item) {
    var dict = ee.Dictionary(item);
    var classId = dict.get('class_label');
    var areaSqM = ee.Number(dict.get('sum'));
    var areaHa = areaSqM.divide(10000); 
    
    var label = ee.String('').cat(
       ee.Algorithms.If(ee.Number(classId).eq(1), 'Low',
       ee.Algorithms.If(ee.Number(classId).eq(2), 'Moderate',
       ee.Algorithms.If(ee.Number(classId).eq(3), 'High',
       'Very High')))
    );

    return ee.Feature(null, {
      'Scenario': scenario,
      'Class_ID': classId,
      'Class_Name': label,
      'Area_Hectares': areaHa,
      'Area_SqMeters': areaSqM
    });
  }));

  // 4. Export Table to Drive (CSV)
  Export.table.toDrive({
    collection: exportFeatures,
    description: 'Table_AreaStats_4Classes_' + scenario,
    folder: 'GEE_Exports',
    fileFormat: 'CSV',
    selectors: ['Scenario', 'Class_ID', 'Class_Name', 'Area_Hectares', 'Area_SqMeters']
  });
  
  print('Table export task created for: ' + scenario);

  // ============================================================
  // H. Cross-tabulate Land Cover by Susceptibility Class
  // ============================================================
  
  // 1. Create combined ID band
  var combinedImage = susceptibilityClass.multiply(10).add(lulc_2024).rename('combined');

  // 2. Add pixel area band
  var statsImage = ee.Image.pixelArea().addBands(combinedImage);

  // 3. Reduce region with grouped reducer
  var combinedStats = statsImage.reduceRegion({
    reducer: ee.Reducer.sum().group({
      groupField: 1,   // band index: 0=pixelArea, 1=combined
      groupName: 'group'
    }),
    geometry: studyArea,
    scale: 30,
    maxPixels: 1e13
  });

  // 4. Extract groups list
  var combinedStatsList = ee.List(combinedStats.get('groups'));

  // 5. Decode combined ID back into Susceptibility and LULC
  var detailedExportFeatures = ee.FeatureCollection(combinedStatsList.map(function(item) {
    item = ee.Dictionary(item);

    var combinedId = ee.Number(item.get('group'));
    var areaSqM = ee.Number(item.get('sum'));

    var susId = combinedId.divide(10).floor();
    var lcId = combinedId.mod(10);

    return ee.Feature(null, {
      'Scenario': scenario,
      'Sus_Level': susId,
      'LULC_ID': lcId,
      'Area_Ha': areaSqM.divide(10000)
    });
  }));

  // 6. Export CSV
  Export.table.toDrive({
    collection: detailedExportFeatures,
    description: 'LULC_Cross_Analysis_Final_' + scenario,
    folder: 'GEE_Exports',
    fileFormat: 'CSV'
  });
  
});

// --- Step 7: Add Legend (4 Classes) ---
function createLegend(title, palette, names) {
  var legend = ui.Panel({ style: { position: 'bottom-left', padding: '8px 15px' }});
  legend.add(ui.Label({ value: title, style: { fontWeight: 'bold', fontSize: '14px', margin: '0 0 4px 0' } }));
  for (var i = 0; i < names.length; i++) {
    var colorBox = ui.Label({ style: { backgroundColor: palette[i], padding: '8px', margin: '0 0 4px 0', border: '1px solid black' }});
    var description = ui.Label({ value: names[i], style: { margin: '0 0 4px 6px' }});
    legend.add(ui.Panel({ widgets: [colorBox, description], layout: ui.Panel.Layout.Flow('horizontal') }));
  }
  return legend;
}

// 4 Colors matching Low -> Very High
var susPalette = ['green','yellow','orange','red'];
var susNames = ['Low Risk', 'Moderate Risk', 'High Risk', 'Very High Risk'];
Map.add(createLegend('Susceptibility (2100)', susPalette, susNames));

// --- Extra Analysis: Driver Distributions and Correlation Matrix ---

var changed = trainingData.filter(ee.Filter.eq('class', 1));
var stable = trainingData.filter(ee.Filter.eq('class', 0));

print(ui.Chart.feature.histogram(changed, 'distance_to_roads').setOptions({title: 'Road Dist: Changed Points'}));
print(ui.Chart.feature.histogram(stable, 'distance_to_roads').setOptions({title: 'Road Dist: Stable Points'}));

var bandsOfInterest = [
  'elevation', 'slope', 'distance_to_roads', 
  'builtup_density', 'builtup_proxy', 'population'
];

print('--- Correlation Matrix (r values) ---');
bandsOfInterest.forEach(function(band1) {
  var row = band1 + ': ';
  bandsOfInterest.forEach(function(band2) {
    if (band1 === band2) {
      row += ' 1.00 |'; 
    } else {
      var correlation = trainingData.reduceColumns({
        reducer: ee.Reducer.pearsonsCorrelation(),
        selectors: [band1, band2]
      });
      var r = ee.Number(correlation.get('correlation'));
      var rVal = r.format('%.2f').getInfo();
      row += ' ' + rVal + ' |';
    }
  });
  print(row);
});
print('-----------------------------------------');
