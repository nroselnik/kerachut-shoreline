"""
Vegetation line from Google Earth historical imagery snapshots.

Workflow (see README.md in this folder for the capture protocol):
1. Capture screenshots of Google Earth web (historical imagery mode) at ONE fixed
   camera position, heading 0, tilt 0, for every date. Name them GE_YYYY-MM-DD.jpg.
2. georeference(): fit pixel scale and offset of the fixed view by matching the
   water mask of one snapshot to a Sentinel-2 scene of a similar date.
3. extract(): along each transect, find the landward edge of bright sand
   (the vegetation line) in every snapshot.

Positions are returned in the same reference frame and baseline distance as
src/kerachut_shoreline.py, so they can be compared directly with waterlines.

Google Earth imagery is (c) Google and its providers and is not redistributed
in this repository; capture your own snapshots following the README.
"""
import glob, json, os
import numpy as np, pandas as pd
from PIL import Image
from scipy import ndimage, optimize
from pyproj import Transformer


# ------------------------------------------------------------------ water mask
def water_mask(img, ui_mask=None):
    """Water = bluish, mid-brightness, smooth texture. img: HxWx3 array."""
    R, G, B = [img[..., i].astype(float) for i in range(3)]
    br = (R + G + B) / 3
    sd = np.sqrt(np.maximum(ndimage.uniform_filter(br ** 2, 9) - ndimage.uniform_filter(br, 9) ** 2, 0))
    w = (B - R > 28) & (br > 75) & (br < 150) & (sd < 12)
    w = ndimage.binary_opening(w, iterations=2)
    lab, n = ndimage.label(w)
    if n == 0: return w
    w = lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)
    w = ndimage.binary_fill_holes(w)
    if ui_mask is not None: w &= ~ui_mask
    return w


# ------------------------------------------------------------------ georeference
def georeference(ge_jpg, s2_shoreline, cam_lon, cam_lat, valid=None,
                 centre_px=(400, 300), s_range=(1.10, 1.26), epsg=32647):
    """Fit metres-per-pixel `s` and offset (dx, dy) of a fixed top-down GE view.
    s2_shoreline: output of kerachut_shoreline.shoreline() for a Sentinel-2 scene
    close in date to the GE snapshot. Returns s, dx, dy, X0, Y0 where
    utm_x = X0 + dx + (col - cx) * s and utm_y = Y0 + dy - (row - cy) * s."""
    img = np.array(Image.open(ge_jpg))
    ge = water_mask(img)
    if valid is None:
        valid = np.ones_like(ge, bool)
    X0, Y0 = Transformer.from_crs(4326, epsg, always_xy=True).transform(cam_lon, cam_lat)
    nd, T, thr = s2_shoreline['ndwi'], s2_shoreline['T'], s2_shoreline['thr']
    yy, xx = np.mgrid[0:ge.shape[0], 0:ge.shape[1]]
    cx, cy = centre_px

    def s2water(p):
        s, dx, dy = p
        ux = X0 + dx + (xx - cx) * s; uy = Y0 + dy - (yy - cy) * s
        c = (ux - T.c) / T.a - 0.5; r = (uy - T.f) / T.e - 0.5
        return ndimage.map_coordinates(nd, [r, c], order=1, mode='nearest') > thr

    def score(p):
        w = s2water(p)
        return -(w & ge & valid).sum() / ((w | ge) & valid).sum()

    best = None
    for s in np.arange(s_range[0], s_range[1], 0.02):
        for dx in range(-40, 41, 8):
            for dy in range(-40, 41, 8):
                v = score((s, dx, dy))
                if best is None or v < best[0]: best = (v, (s, dx, dy))
    res = optimize.minimize(score, best[1], method='Nelder-Mead', options=dict(xatol=0.2, fatol=1e-5))
    s, dx, dy = res.x
    return dict(s=s, dx=dx, dy=dy, X0=X0, Y0=Y0, cx=cx, cy=cy, iou=-res.fun)


# ------------------------------------------------------------------ transects
def load_transects(path):
    fc = json.load(open(path)); out = []
    for f in fc['features']:
        c = np.array(f['geometry']['coordinates'])
        out.append(dict(id=f['properties']['transect'], beach=bool(f['properties']['beach']),
                        o=c[0], n=(c[1] - c[0]) / np.linalg.norm(c[1] - c[0])))
    return out


class Frame:
    """Maps reference-frame UTM coordinates to GE pixel coordinates.
    frame_shift: the co-registration shift (x, y) of the Sentinel-2 scene used for
    georeferencing, relative to the reference scene (from coregistration.csv)."""
    def __init__(self, geo, frame_shift=(0.0, 0.0)):
        self.g = geo; self.sh = np.asarray(frame_shift, float)

    def ref2px(self, xy):
        g = self.g; u = np.atleast_2d(xy) + self.sh
        return np.c_[g['cx'] + (u[:, 0] - g['X0'] - g['dx']) / g['s'],
                     g['cy'] - (u[:, 1] - g['Y0'] - g['dy']) / g['s']]


# ------------------------------------------------------------------ vegetation line
def veg_line(img, t, frame, d_start, step=0.5, search=95, labels=()):
    """Distance from transect origin to the vegetation line.
    Walk landward from d_start (a point on or seaward of the beach); after reaching
    bright sand, the vegetation line is where brightness falls below the midpoint
    between sand and vegetation and stays below for 5 m.
    labels: list of (x0, x1, y0, y1) pixel boxes covered by map labels (masked)."""
    br = img.astype(float).mean(2)
    ds = np.arange(d_start, d_start - search, -step)
    px = frame.ref2px(t['o'] + ds[:, None] * t['n'])
    prof = ndimage.uniform_filter1d(ndimage.map_coordinates(br, [px[:, 1], px[:, 0]], order=1), 7)
    sand = np.percentile(prof[: int(40 / step)], 90)
    veg = np.median(prof[-int(30 / step):])
    if sand - veg < 25: return np.nan, 'low contrast'
    thr = (sand + veg) / 2
    k0 = int(np.argmax(prof > thr)); run = int(5 / step)
    inside = lambda p: any(x0 <= p[0] <= x1 and y0 <= p[1] <= y1 for x0, x1, y0, y1 in labels)
    for k in range(k0, len(prof) - run):
        if np.all(prof[k:k + run] < thr):
            a, b = prof[k - 1], prof[k]
            d = ds[k - 1] - ((a - thr) / (a - b) if a != b else 0) * step
            for dd in (d, d - 5):
                if inside(frame.ref2px(t['o'] + dd * t['n'])[0]): return np.nan, 'label'
            return d, 'ok'
    return np.nan, 'no edge'


def extract(ge_folder, transect_file, frame, waterline_median, labels=(), start_offset=5):
    """Vegetation line for every GE_YYYY-MM-DD.jpg in ge_folder.
    waterline_median: dict transect id -> typical waterline distance (search start)."""
    trs = [t for t in load_transects(transect_file) if t['beach']]
    rows = {}
    for f in sorted(glob.glob(os.path.join(ge_folder, 'GE_*.jpg'))):
        img = np.array(Image.open(f))
        date = os.path.basename(f)[3:13]
        rows[date] = {t['id']: veg_line(img, t, frame, waterline_median[t['id']] + start_offset, labels=labels)[0]
                      for t in trs}
    V = pd.DataFrame(rows).T
    V.index.name = 'date'
    return V
