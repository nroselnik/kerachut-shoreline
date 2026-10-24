# Google Earth Engine for beginners: from sign-up to shoreline trends

This guide walks a first-time user through the whole workflow in this repository. No coding experience is assumed beyond running cells in a notebook. It was written for ecologists who want a shoreline history for one beach, but the steps apply to any small coastal site.

Last checked: October 2026. Earth Engine's access rules change from time to time, so check the official pages if something looks different.

---

## 1. What Earth Engine is

Google Earth Engine (GEE) is a free cloud service that holds the full Landsat archive (from the 1980s), Sentinel-2 (from 2015) and hundreds of other datasets. Your code runs on Google's servers, so you never download thousands of images. You only download the small pieces you need.

GEE is a different product from Google Earth (the 3D globe). Google Earth is used in this project only for its historical high-resolution photos; see `scripts/google_earth/README.md`.

## 2. Is it free?

Yes, for noncommercial work such as research, teaching and conservation. Since April 2026, every noncommercial project gets a monthly compute allowance, set by tier:

- **Community tier:** about 150 compute hours a month. No billing account needed. This is plenty for one beach.
- **Contributor tier:** about 1,000 compute hours a month. Needs a billing account attached, but noncommercial Earth Engine use is not charged.
- **Partner tier:** large allowances by application, for climate and biodiversity work.

If you run out, Earth Engine keeps working, but more slowly, until the allowance resets on the first of the next month. Noncommercial status must be re-confirmed once a year.

Official page: https://developers.google.com/earth-engine/guides/noncommercial_tiers

## 3. Sign up and find your project ID

1. Go to https://earthengine.google.com and click **Get started**. Sign in with a Google account. A university account is fine.
2. Create a **Google Cloud project** when asked, and register it for **noncommercial / academic** use.
3. Find your **project ID**. This is the value you paste into the notebooks. It is the *ID*, not the display name. It is all lowercase with hyphens and numbers, for example `ee-yourname` or `my-project-123456`. You can find it in either of two places:
   - **Google Cloud Console** (https://console.cloud.google.com): click the project picker at the top left. The pop-up lists **Name** and **ID** columns. Copy the ID.
   - **Earth Engine Code Editor** (https://code.earthengine.google.com): your profile icon at the top right shows the current project.
4. Check registration: in the Cloud Console go to **Earth Engine > Configuration**. The project should show as registered for noncommercial use, with a tier.

> **Note for future researchers:** every notebook here contains the line `PROJECT = 'YOUR-GEE-PROJECT-ID'`. Replace it with your own project ID. Never share or commit your own project ID in public code.

## 4. Code Editor or Python?

- The **Code Editor** (in the browser) uses JavaScript only.
- **Python** runs in **Google Colab** (https://colab.research.google.com). Colab needs no installation and has the Earth Engine library built in. All notebooks in this repository are written for Colab.

## 5. Running the notebooks

Open a notebook from `notebooks/` in Colab with **File > Open notebook > GitHub**, then paste the repository URL. Run the cells from top to bottom with Shift+Enter.

### Notebook 01: look at the imagery and export the best scene per year

- The first run asks you to sign in to Google (`ee.Authenticate()`). Allow access.
- The **study box** cell shows a map. Check that the red box covers the beach, some sea and the land behind it. If it doesn't, draw your own rectangle with the map's rectangle tool and run the next cell.
- The **inventory** cell lists every Landsat and Sentinel-2 scene, with the percentage of cloud-free pixels over your box.
- The **grid** cell shows the clearest scene for each year. Check it by eye and note any scene with cloud the mask missed. In our case 19 Aug 2013 was excluded, which is the `EXCLUDE` list in `src/kerachut_shoreline.py`.
- The **export** cell sends GeoTIFFs to your Google Drive. Set `EXPORT = True` to run it. Exports run on Google's side, so you can close Colab once they start. Watch progress in the Cloud Console under **Earth Engine > Tasks**.

### Notebook 02: shoreline positions and trends

This notebook reads the exported GeoTIFFs from Drive and writes the results (tables, GeoJSON shorelines and transects, figures) to `results/` in the same Drive folder.

### Notebook 03: tide correction with every clear scene

- Downloads several hundred scenes directly with `computePixels`, skipping the Drive export queue. Run time is about 20 to 30 minutes. If Colab disconnects, run the cell again and it continues where it stopped.
- Downloads the EOT20 tide model (about 2 GB) to Colab's temporary disk, not your Drive.
- Writes everything to `results_tide/`.

## 6. Problems we hit, and fixes

- **`ImageCollection asset 'GOOGLE_CLOUD_SCORE_PLUS/V1/S2_HARMONIZED' not found`**: some projects can't reach the Cloud Score+ collection. The code here uses the Scene Classification band (SCL) that comes with every Sentinel-2 image instead.
- **Several Drive folders with the same name**: Earth Engine creates them when many exports start at once. Move all the GeoTIFFs into one folder before running notebook 02.
- **"Export is off" message**: the export cell is off by default. Set `EXPORT = True`.
- **Files missing in Drive right after "done" in Colab**: Colab only confirms that the exports *started*. Wait for the Tasks page to show **Completed**.
- **Landsat and Sentinel-2 disagree by about 20 m**: this is real misregistration between the two satellite families, not erosion. The co-registration step in `kerachut_shoreline.py` corrects it using the surrounding rocky coast.

## 7. Adapting this to another beach

1. In notebook 01, change `CENTER_LAT`, `CENTER_LON` and `HALF_SIZE`, or draw a box on the map.
2. Choose a clear Sentinel-2 scene as the reference and set `REF` in `src/kerachut_shoreline.py`.
3. In `beach_transects()`, the seaward direction is assumed to be west (`if n[0] > 0: n = -n`). Change this for beaches facing another way.
4. The co-registration needs stable, non-sandy shoreline (rock, seawall) inside the box. Check that `control_transects()` finds at least 40.
5. In `src/kerachut_tides_run.py`, change the tide point (`TIDE_LON`, `TIDE_LAT`) to a spot just offshore of your beach, and change the tide gauge used for validation.

## 8. Data credits

- Landsat: courtesy of the U.S. Geological Survey.
- Sentinel-2: Copernicus programme, European Union.
- EOT20 tide model: Hart-Davis et al. (2021), SEANOE doi:10.17882/79489.
- Penang tide gauge (1992 hourly): University of Hawaii Sea Level Center, via PANGAEA doi:10.1594/PANGAEA.439979.
- Google Earth imagery: © Google and its data providers. It is not included in this repository.
