// Original GEE script: https://code.earthengine.google.com/abcc44dd4eaf9805209ddacb8a090b26
// Note: link requires Google authentication; full code provided below for archival purposes


// =========================================================
// Driver Stacking - Bandung Metropolitan Area
// Prepares susceptibility model driver variables (elevation, slope,
// distance to roads, built-up density/proxy, population) for a given year
// Part of: BMA LULC change and green space conversion susceptibility study
// =========================================================

// ============================================================
// 0. Setup & Define Study Area
// ============================================================

var studyArea = ee.FeatureCollection('projects/ee-agiegizawi/assets/ADM_BMA_48S');

// Define study area style
var study_area_vis = {color: "ff0000", fillColor: "00000000", width: 2};
Map.addLayer(studyArea.style(study_area_vis), {}, 'Study Area Boundary');
Map.centerObject(studyArea, 10);

// ============================================================
// 1. SELECT TARGET YEAR
// ============================================================
// Change this to '2000', '2012', or '2024' to switch all inputs automatically
var targetYear = '2012';

// ============================================================
// 2. CONFIGURE INPUTS BASED ON YEAR
// ============================================================
// This maps your Imported Variables to the logic
var inputConfig = {
  '2000': {
    lulc: lulc00,
    density: density00,
    proxy: proxy00,
    pop: pop00,
    roads: roads00 // Uses roads00
  },
  '2012': {
    lulc: lulc12,
    density: density12,
    proxy: proxy12,
    pop: pop12,
    roads: roads00 // Fallback: Uses roads00 for 2012 
  },
  '2024': {
    lulc: lulc24,
    density: density24,
    proxy: proxy24,
    pop: pop24,
    roads: roads24 // Uses roads24 (OSM)
  }
};

// Get the specific data for the chosen year
var data = inputConfig[targetYear];

// Set Projection Reference (from the chosen LULC year)
var projection = data.lulc.projection();
var scale = 30;

// ============================================================
// 3. PREPARE DRIVERS
// ============================================================

// --- A. Built-up Density (Inside City) ---
var densityLayer = data.density
    .rename('builtup_density')
    .clip(studyArea);

// --- B. Built-up Proxy/Distance (Outside City) ---
var proxyLayer = data.proxy
    .rename('builtup_proxy')
    .clip(studyArea);

// --- C. Population (Resample to 1km match) ---
var populationStd = data.pop
    .resample('bilinear')
    .reproject({
        crs: projection,
        scale: 1000
    })
    .rename('population')
    .clip(studyArea);

// --- D. Distance to Roads (Prepared in QGIS) ---
var distanceToRoads = ee.Image(
  'projects/ee-agiegizawi/assets/distance_to_road_2000_30m'
)
  .rename('distance_to_roads')
  .reproject({
    crs: projection,
    scale: scale
  })
  .clip(studyArea);

// --- E. Terrain (Static) ---
var dem = ee.Image("USGS/SRTMGL1_003");
var elevation = dem.select('elevation')
    .reproject(projection)
    .clip(studyArea)
    .rename('elevation');
var slope = ee.Terrain.slope(dem)
    .reproject(projection)
    .clip(studyArea)
    .rename('slope');

// ============================================================
// 4. STACK LAYERS
// ============================================================

var drivers_stacked = ee.Image.cat([
  elevation,
  slope,
  distanceToRoads,
  densityLayer,  // Driver A
  proxyLayer,    // Driver B
  populationStd
]).clip(studyArea);

// Verify in Console
print("Stacked Driver Image (" + targetYear + "):", drivers_stacked);
print("Bands:", drivers_stacked.bandNames());

// ============================================================
// 5. VISUALIZATION
// ============================================================

Map.addLayer(densityLayer, {min: 0, max: 1, palette: ['black', 'blue', 'yellow', 'red']}, 'Built-up Density ' + targetYear);
Map.addLayer(proxyLayer, {min: 0, max: 5000, palette: ['white', 'black']}, 'Built-up Proxy ' + targetYear);

// ============================================================
// 6. EXPORT
// ============================================================

Export.image.toAsset({
  image: drivers_stacked,
  description: 'drivers_' + targetYear + '_separated_v1',
  assetId: 'projects/ee-agiegizawi/assets/drivers_' + targetYear + '_separated_v1',
  scale: scale,
  region: studyArea.geometry(),
  maxPixels: 1e13
});
