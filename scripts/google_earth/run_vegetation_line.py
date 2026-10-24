"""
Example: reproduce data/vegetation_line_by_transect.csv from your own Google Earth
snapshots of Pantai Kerachut. Edit the paths at the top.

    python scripts/google_earth/run_vegetation_line.py

Requirements: the snapshots (GE_YYYY-MM-DD.jpg, 800 x 600 px, captured as in
README.md), one Sentinel-2 GeoTIFF exported by notebook 01 close in date to the
georeferencing snapshot, and results from notebook 02 (positions per transect).
"""
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'src')); sys.path.insert(0, HERE)
import kerachut_shoreline as k
import ge_vegetation_line as g

# ---------------------------------------------------------------- edit these
GE_FOLDER = 'ge_snapshots'                                   # GE_YYYY-MM-DD.jpg files
GEOREF_SNAPSHOT = os.path.join(GE_FOLDER, 'GE_2026-01-13.jpg')
GEOREF_S2 = 'Kerachut_GEE/Kerachut_20260105_S2.tif'          # Sentinel-2 scene within a few weeks
S2_SHIFT = (-3.98, 0.34)            # that scene's shift_x_m, shift_y_m from results/coregistration.csv
POSITIONS = 'Kerachut_GEE/results/shoreline_position_by_transect.csv'
CAMERA_LON, CAMERA_LAT = 100.18285, 5.45118973               # camera position used for the snapshots
TRANSECTS = os.path.join(REPO, 'data', 'transects_utm47n.geojson')
OUT = 'vegetation_line_by_transect.csv'

# map labels that could not be hidden (pixel boxes x0, x1, y0, y1); measured on our snapshots
LABELS = [(272, 358, 368, 392), (360, 388, 366, 394), (236, 364, 405, 442), (366, 392, 409, 437)]

# area of the screenshot that shows the map (excludes toolbars and buttons)
valid = np.zeros((600, 800), bool)
valid[130:570, :] = True; valid[505:575, 0:75] = False; valid[515:575, 560:800] = False
# ----------------------------------------------------------------------------

s2 = k.shoreline(GEOREF_S2)
geo = g.georeference(GEOREF_SNAPSHOT, s2, CAMERA_LON, CAMERA_LAT, valid=valid)
print('georeference:', {a: round(float(b), 3) for a, b in geo.items()})
frame = g.Frame(geo, S2_SHIFT)

P = pd.read_csv(POSITIONS, index_col=0)
wm = P[P.sensor == 'S2'].drop(columns='sensor').median()
wm = {int(c): v for c, v in wm.items()}

V = g.extract(GE_FOLDER, TRANSECTS, frame, wm, labels=LABELS)
V.round(1).to_csv(OUT)
print(V.round(0).to_string())
