# BMA-LULC-RandomForest-GEE

Google Earth Engine scripts for land use/land cover (LULC) classification and green space conversion susceptibility modeling in the Bandung Metropolitan Area (BMA), West Java, Indonesia.

These scripts support the manuscript submitted to *Land Use Policy* (Elsevier): Predicting green space vulnerability to urban development in Indonesia using historical satellite imagery and population projections.

## Study overview

- **Study area:** Bandung Metropolitan Area, West Java, Indonesia
- **Classification method:** Object-Based Image Analysis (OBIA, SNIC segmentation) with Random Forest (100 trees), applied to Landsat 7 (2000, 2012) and Landsat 8 (2024) surface reflectance imagery
- **Land cover classes:** Built-up, Water body, Bare land, Agriculture and other vegetaion, Forest
- **Susceptibility modeling:** Random Forest trained in probability mode on observed 2012→2024 green space-to-urban conversion, validated against 2000→2012 backcasting (Phase 3A) and 2012→2024 (Phase 3B), then projected to 2100 under SSP3 and SSP5 population scenarios
- **Reported accuracy:** Overall accuracy 92–94%, Kappa 0.90–0.92; AUC 0.86 (backcasting validation) and 0.92 (primary model)

## Scripts

Run in the order below — each stage depends on outputs (GEE assets) from the previous one.

| # | Script | Purpose |
|---|--------|---------|
| 1 | `01_LULC_classification_2000.js` | OBIA + Random Forest LULC classification, year 2000 (Landsat 7) |
| 2 | `02_LULC_classification_2012.js` | OBIA + Random Forest LULC classification, year 2012 (Landsat 7) |
| 3 | `03_LULC_classification_2024.js` | OBIA + Random Forest LULC classification, year 2024 (Landsat 8) |
| 4 | `04_driver_stacking.js` | Prepares and stacks susceptibility model driver variables (elevation, slope, distance to roads, built-up density/proxy, population) for a selected year |
| 5 | `05_AUC_validation_phase3A.js` | Backcasting validation: trains a probability-mode RF model on 2000 drivers, predicts 2000→2012 change, reports ROC/AUC (0.86) |
| 6 | `06_AUC_validation_phase3B.js` | Trains the primary susceptibility model on 2012 drivers, predicts 2012→2024 change, reports ROC/AUC (0.92), exports the trained classifier |
| 7 | `07_susceptibility_projection_SSP3_SSP5.js` | Projects green space conversion susceptibility to 2100 under SSP3 and SSP5 population scenarios; exports susceptibility maps and area statistics by class |

## Data sources

- **Imagery:** Landsat 7 Collection 2 Level-2 (`LANDSAT/LE07/C02/T1_L2`), Landsat 8 Collection 2 Level-2 (`LANDSAT/LC08/C02/T1_L2`)
- **Terrain:** SRTM 30m Global (`USGS/SRTMGL1_003`)
- **Population (historical):** [ WorldPop / gridded population dataset used]
- **Population (future, SSP3/SSP5):** Shared Socioeconomic Pathway population projections, 1km resolution, 2100
- **Roads:** [OpenStreetMap / prepared in QGIS]
- **Study area boundary:** Administrative boundary of BMA (`ADM_BMA_48S` asset)

## How to run

1. Open [code.earthengine.google.com](https://code.earthengine.google.com)
2. Copy a script's contents into a new file in the Code Editor
3. Update asset paths (`projects/...`) to point to your own GEE assets, or request access to the original assets
4. Run script in order (01 → 07); each stage's exports feed into the next stage's inputs

Note: scripts reference private GEE assets under the author's account. To reproduce the full pipeline, imagery and intermediate outputs (classified maps, driver stacks) will need to be regenerated or requested from the author.

## Citation

If you use this code, please cite:

> Gizawi, A.S., (2026). Predicting green space vulnerability to urban development in Indonesia using historical satellite imagery and population projections. *Land Use Policy*. [DOI once available]

Code archive: [Zenodo DOI to be added]

## Contact

[Agie Syirban Gizawi/ agie.gizawi@gmail.com/ Forest Engineering Lab, Mie University]

## License

[MIT]
