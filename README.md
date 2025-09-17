# MDS Sample/Cluster Extraction (Shiny)

Interactive Shiny app to manually select sample groups on MDS plots computed from pairwise genetic distance matrices.  
Draw polygons on the plot; points inside (or on the boundary) are selected. Selected Taxon IDs are printed and saved as `groupN.txt` files.

## Features
- MDS plot from TASSEL IBS distance matrices  
- Color samples by pre-defined group (from key file)  
- Freehand polygon selection with live highlight  
- Exports selected sample IDs to `output/groupN.txt`  
- Small test dataset  

## Repository structure
```
├── app.R                      # your Shiny app (this script)
├── example_data_region_bundle.zip
│   ├── data/
│   │   ├── example_IBS.dist.txt
│   │   └── example_key.tsv
│   └── output/               # created/populated by the app
└── README.md
```

## Requirements
- R (tested with version 4.0.2 (2020-06-22))  
- R packages:  
  - data.table  
  - shiny  
  - leaflet  
  - leaflet.extras  
  - sf  
  - R.utils  
  - RColorBrewer  
  - Hmisc (optional)  

## Quick start (with the bundled example)
1. Download `SampleExtraction_MDSplot.R` (app) and `example_data_region_bundle.zip` (example dataset).  
2. Unzip `example_data_region_bundle.zip`. After extraction you should have `./data` and `./output` folders.  
3. Open `SampleExtraction_MDSplot.R` and set the **USER CONFIG** block:  
   - directory of `example_data_region_bundle.zip`  
   - quality filter: `propmis` (proportion missing from taxa summary)  
   - `color_group`: group by which MDS plot should be coloured  
4. Run the code. After `shinyApp(ui, server)` the Shiny app will open.  

## Application
- Pick polygon tool on the left side  
- Draw polygon around samples of interest until closed to select  
- Selected Taxon IDs are highlighted  

## Output
- Each new polygon creates a new file `groupN.txt` in `output/` (`group1.txt`, `group2.txt`, …).  

## Troubleshooting
- **No points appear**  
  - `color_group` doesn’t match any column in the key file (should be e.g. `group` & `group_color`).  
  - `Taxon` values don’t overlap between matrix and key (check spelling/case).  

- **“Cannot correctly read … Expect TASSEL IBS.”**  
  - Ensure the distance file is tab-delimited TASSEL IBS with 5 header lines; no extra delimiters.  

- **sf / GDAL errors on install**  
  - Check installation guide for mac/Linux online.  

- **Nothing exported**  
  - Confirm the polygon encloses points.  
  - Check console messages and output/write permissions.  

## Acknowledgements
- MDS functions adapted from K. Swarts (2019), used with permission
