# Kerachut shoreline: 38 years of beach change at a green turtle nesting beach from free imagery

Code and data for reconstructing shoreline and vegetation line change at **Pantai Kerachut**, Penang National Park, Malaysia (1988 to 2026), using only free satellite imagery. Pantai Kerachut is the main green turtle (*Chelonia mydas*) nesting beach on Penang Island.

**Main finding:** the beach has not eroded as a whole. It has rotated. Since about 2017 the southern waterline has retreated (-1.6 ± 0.3 m/yr, 2018 to 2026) and the northern waterline has advanced (+3.1 ± 0.3 m/yr), while the vegetation line behind the beach has stayed fixed. Dry beach width in the southern third fell from about 54 m (1988 to 2016 mean) to about 30 m (2024 to 2026).

Manuscript in preparation for *Chelonian Conservation and Biology*.

## What is here

```
notebooks/
  01_satellite_check_and_export.ipynb   look at the imagery, pick the clearest scene per year, export to Drive
  02_shoreline_change_best_scenes.ipynb shoreline extraction, co-registration, change rates (37 scenes)
  03_tide_corrected_all_scenes.ipynb    all clear scenes (291), EOT20 tides, beach slope, corrected trends
src/
  kerachut_shoreline.py   shoreline extraction, transects, Landsat/Sentinel-2 co-registration, rates, figures 1-3
  kerachut_tides.py       fixed-transect extraction for many scenes, tide correction, slope estimation, figures 6-9
  kerachut_tides_run.py   Colab helpers: bulk scene download, EOT20 setup, tide-gauge data
scripts/google_earth/
  ge_vegetation_line.py   georeference Google Earth snapshots and trace the vegetation line
  run_vegetation_line.py  example run
  README.md               capture protocol for Google Earth historical imagery
data/
  transects_utm47n.geojson          30 shore-normal transects, 20 m apart, UTM zone 47N (sections use T2 to T26)
  vegetation_line_by_transect.csv   vegetation line distance from baseline, 19 Google Earth dates, 2002 to 2026
  positions_all_scenes_tide.csv     waterline distance per transect for 291 scenes, with modelled tide (m)
figures/                            figures from the analysis (Google Earth backdrops excluded)
docs/
  GEE_tutorial.md                   step-by-step guide to Earth Engine, Colab and this workflow
```

## Method in brief

1. **Imagery:** Landsat 5/7/8/9 (30 m) and Sentinel-2 (10 m) surface reflectance in Google Earth Engine. Scenes are kept if at least 90% cloud-free over a 1.8 km box.
2. **Shoreline:** sub-pixel NDWI contour, with an Otsu threshold between water and sand pixels near the beach (after Vos et al. 2019, CoastSat).
3. **Co-registration:** each scene is fitted to a Sentinel-2 reference along about 66 control transects on the surrounding rocky coast. The fit gives an x/y shift and a threshold bias. This removes a roughly 20 m Landsat to Sentinel-2 offset.
4. **Tide:** EOT20 global tide model (via eo-tides / pyTMD) at each acquisition time. Beach slope is estimated from the imagery by regressing position on tide within years. Results are also checked with tide-filtered scenes and a slope sensitivity test.
5. **Vegetation line:** 19 Google Earth historical images, georeferenced to Sentinel-2 by water-mask matching.
6. **Statistics:** linear regression rates per transect and per section (south T2 to T10, central T11 to T19, north T20 to T26), plus a rotation index (north minus south).

## Open the notebooks in Colab

- [01 Imagery check and export](https://colab.research.google.com/github/nroselnik/kerachut-shoreline/blob/main/notebooks/01_satellite_check_and_export.ipynb)
- [02 Shoreline change, best scene per year](https://colab.research.google.com/github/nroselnik/kerachut-shoreline/blob/main/notebooks/02_shoreline_change_best_scenes.ipynb)
- [03 Tide-corrected analysis, all clear scenes](https://colab.research.google.com/github/nroselnik/kerachut-shoreline/blob/main/notebooks/03_tide_corrected_all_scenes.ipynb)

## Getting started

New to Earth Engine? Start with [docs/GEE_tutorial.md](docs/GEE_tutorial.md).

In every notebook, replace

```python
PROJECT = 'YOUR-GEE-PROJECT-ID'
```

with your own Google Cloud project ID that is registered for Earth Engine. Never commit your own ID to a public repository.

Python packages (Colab installs these in the first cell): see `requirements.txt`.

## Data sources and credits

- Landsat: U.S. Geological Survey. Sentinel-2: Copernicus programme, European Union.
- EOT20: Hart-Davis, M.G. et al. 2021. *Earth System Science Data* 13:3869–3884; data doi:10.17882/79489.
- Penang tide gauge, hourly 1992: University of Hawaii Sea Level Center, via PANGAEA doi:10.1594/PANGAEA.439979.
- Google Earth imagery: © Google and its data providers. It is not redistributed here.

## How to cite

Please cite the paper once published. Until then, cite this repository (see `CITATION.cff`).

## Licence

Code: MIT (see `LICENSE`). Data files in `data/` and figures in `figures/`: CC BY 4.0.

## Contact

Nik Fadzly N. Rosely, School of Biological Sciences, Universiti Sains Malaysia.
