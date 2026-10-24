"""
Kerachut shoreline change from exported Landsat / Sentinel-2 GeoTIFFs.

Method
1. Shoreline = sub-pixel NDWI contour. Threshold = Otsu between water pixels and
   sand pixels next to the beach (CoastSat-style), per scene.
2. Co-registration: every scene is compared with the reference Sentinel-2 scene
   along ~60 control transects on the rocky coast (assumed stable). A per-scene
   model  offset = shift . normal + delta  gives the geolocation shift (x, y) and
   the threshold bias (delta). Beach positions are corrected with that model.
3. Beach transects every 20 m, normal to the beach, from a landward baseline.
4. Per transect: NSM, EPR, LRR (+95% CI) for the full record and per sensor era.

No tide correction is applied (see notes in the outputs).
"""
import glob, os, re, json
import numpy as np, pandas as pd, rasterio
from scipy import ndimage, stats
from skimage import measure
from skimage.filters import threshold_otsu
from shapely.geometry import LineString, Point, mapping

BANDS = ['blue', 'green', 'red', 'nir', 'swir1']
EXCLUDE = ['Kerachut_20130819_L8']      # cloud the mask missed
REF = 'Kerachut_20180117_S2'            # reference scene (first clear S2)
SPACING, BACK, LENGTH = 20, 150, 350    # transect spacing, baseline setback, length (m)


# ---------------------------------------------------------------- basics
def load(f):
    with rasterio.open(f) as s:
        a = s.read().astype(float); T = s.transform; crs = s.crs; bd = s.bounds
    return dict(zip(BANDS, a)), T, crs, bd

def indices(b):
    e = 1e-9
    ndwi = (b['green'] - b['nir']) / (b['green'] + b['nir'] + e)
    ndvi = (b['nir'] - b['red']) / (b['nir'] + b['red'] + e)
    bright = (b['red'] + b['green'] + b['blue']) / 3
    return ndwi, ndvi, bright

def px2xy(T, rc):
    r, c = rc[:, 0], rc[:, 1]
    return np.c_[T.c + (c + 0.5) * T.a, T.f + (r + 0.5) * T.e]

def biggest(mask):
    lab, n = ndimage.label(mask)
    if n == 0: return mask
    return lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)

def shoreline(f):
    b, T, crs, bd = load(f)
    ndwi, ndvi, br = indices(b)
    sand = (br > 0.10) & (ndvi < 0.25) & (ndwi < 0)
    water = ndwi > threshold_otsu(ndwi)
    near = ndimage.binary_dilation(sand, iterations=3)
    wv, sv = ndwi[near & water], ndwi[sand]
    thr = threshold_otsu(np.r_[wv, sv]) if (len(wv) >= 3 and len(sv) >= 3) else threshold_otsu(ndwi)
    c = max(measure.find_contours(ndwi, thr), key=len)
    return dict(line=LineString(px2xy(T, c)), thr=thr, sand=sand, T=T, crs=crs, bounds=bd,
                ndwi=ndwi, rgb=np.dstack([b['red'], b['green'], b['blue']]))


# ---------------------------------------------------------------- transects
def beach_transects(sand_xy):
    pts = sand_xy - sand_xy.mean(0)
    d = np.linalg.svd(pts, full_matrices=False)[2][0]
    n = np.array([-d[1], d[0]])
    if n[0] > 0: n = -n                              # seaward = west for Kerachut
    if d[1] < 0: d = -d                              # ids run south to north
    proj = pts @ d
    o = sand_xy.mean(0); out = []
    for i, a in enumerate(np.arange(proj.min(), proj.max() + 1e-6, SPACING)):
        p0 = o + a * d - BACK * n
        out.append(dict(id=i + 1, origin=p0, n=n, line=LineString([p0, p0 + LENGTH * n])))
    return out

def dist_on(sl, t, pick='max'):
    inter = sl.intersection(t['line'])
    if inter.is_empty: return np.nan
    g = getattr(inter, 'geoms', [inter])
    ds = [np.dot(np.array(x.coords[0]) - t['origin'], t['n']) for x in g if x.geom_type == 'Point']
    if not ds: return np.nan
    return max(ds) if pick == 'max' else ds[np.argmin(np.abs(ds))]

def control_transects(ref, beach_xy, half=120, edge=60):
    sl, L = ref['line'], ref['line'].length
    c = beach_xy.mean(0); r = np.ptp(beach_xy, 0).max() / 2 + 180
    x0, y0, x1, y1 = ref['bounds']; inv = ~ref['T']; out = []
    for s in np.arange(10, L - 10, 20):
        p = np.array(sl.interpolate(s).coords[0])
        if np.linalg.norm(p - c) < r: continue
        if not (x0 + edge < p[0] < x1 - edge and y0 + edge < p[1] < y1 - edge): continue
        q1 = np.array(sl.interpolate(max(s - 15, 0)).coords[0])
        q2 = np.array(sl.interpolate(min(s + 15, L)).coords[0])
        tv = (q2 - q1) / np.linalg.norm(q2 - q1); n = np.array([tv[1], -tv[0]])
        q = p + 40 * n; col, row = inv * (q[0], q[1])
        if ref['ndwi'][int(row), int(col)] < ref['thr']: n = -n     # point normal to sea
        out.append(dict(origin=p, n=n, line=LineString([p - half * n, p + half * n])))
    return out

def fit_bias(off, ctrl):
    N = np.array([t['n'] for t in ctrl]); ok = np.isfinite(off)
    A = np.c_[N[ok], np.ones(ok.sum())]; y = off[ok]; w = np.ones(len(y))
    for _ in range(10):
        p = np.linalg.lstsq(A * w[:, None], y * w, rcond=None)[0]
        res = y - A @ p; mad = 1.4826 * np.median(np.abs(res)) + 1e-6
        w = 1 / np.maximum(1, np.abs(res) / (2 * mad))
    return p, float(np.sqrt(np.mean(res ** 2))), int(ok.sum())


# ---------------------------------------------------------------- rates
def rates(t_yr, x):
    ok = np.isfinite(x)
    t, y = t_yr[ok], x[ok]
    if len(y) < 3: return dict(n=len(y))
    lr = stats.linregress(t, y)
    ci = stats.t.ppf(0.975, len(y) - 2) * lr.stderr
    i0, i1 = np.argmin(t), np.argmax(t)
    return dict(n=len(y), first=float(t[i0]), last=float(t[i1]),
                NSM_m=float(y[i1] - y[i0]), EPR_m_yr=float((y[i1] - y[i0]) / (t[i1] - t[i0])),
                LRR_m_yr=float(lr.slope), LRR_CI95=float(ci), LRR_r2=float(lr.rvalue ** 2),
                LRR_p=float(lr.pvalue))


# ---------------------------------------------------------------- main
def run(folder, outdir):
    os.makedirs(outdir, exist_ok=True)
    files = sorted(f for f in glob.glob(os.path.join(folder, 'Kerachut_*.tif'))
                   if os.path.basename(f)[:-4] not in EXCLUDE)
    meta = []
    for f in files:
        m = re.search(r'Kerachut_(\d{8})_(\w+)\.tif', os.path.basename(f))
        dt = pd.to_datetime(m.group(1)); meta.append(dict(file=f, date=dt, sensor=m.group(2),
                                     year_dec=dt.year + (dt.dayofyear - 0.5) / 365.25))
    meta = pd.DataFrame(meta)
    S = {r.file: shoreline(r.file) for r in meta.itertuples()}

    ref = S[[f for f in S if REF in f][0]]
    s2 = [f for f in S if f.endswith('_S2.tif')]
    sand_xy = np.vstack([px2xy(S[f]['T'], np.argwhere(biggest(S[f]['sand']))) for f in (s2[0], s2[-1])])
    trs = beach_transects(sand_xy)
    ctrl = control_transects(ref, px2xy(ref['T'], np.argwhere(biggest(ref['sand']))))

    # beach flag: transect crosses the shoreline at >= 60 deg (median over S2 scenes)
    # and has sand within 20 m landward in at least half of the S2 scenes
    def crossing_angle(sl, t):
        dd = dist_on(sl, t)
        if np.isnan(dd): return np.nan
        p = t['origin'] + dd * t['n']; s = sl.project(Point(p))
        q1 = np.array(sl.interpolate(max(s - 15, 0)).coords[0]); q2 = np.array(sl.interpolate(s + 15).coords[0])
        tv = (q2 - q1) / (np.linalg.norm(q2 - q1) + 1e-9)
        return np.degrees(np.arccos(min(1, abs(np.dot(tv, t['n'])))))   # 90 = perpendicular
    for t in trs:
        ang = [crossing_angle(S[f]['line'], t) for f in s2]
        t['angle'] = float(np.nanmedian(ang)) if np.isfinite(ang).any() else np.nan
        hits = 0
        for f in s2:
            dd = dist_on(S[f]['line'], t)
            if np.isnan(dd): continue
            p = t['origin'] + dd * t['n']; inv = ~S[f]['T']; sd = S[f]['sand']
            for k in np.arange(0, 21, 5):                   # sand within 20 m landward
                q = p - k * t['n']; col, row = inv * (q[0], q[1])
                if 0 <= int(row) < sd.shape[0] and 0 <= int(col) < sd.shape[1] and sd[int(row), int(col)]:
                    hits += 1; break
        t['sand_frac'] = hits / len(s2)
        t['beach'] = bool(t['angle'] >= 60 and t['sand_frac'] >= 0.5)

    rows, pos = [], {}
    for r in meta.itertuples():
        sl = S[r.file]['line']
        # offset from the reference shoreline point along each control normal
        off = np.array([dist_on(sl, c, 'nearest') for c in ctrl])
        p, rms, n_ok = fit_bias(off, ctrl)
        rows.append(dict(date=r.date.date(), sensor=r.sensor, threshold=S[r.file]['thr'],
                         shift_x_m=p[0], shift_y_m=p[1], bias_m=p[2], ctrl_rms_m=rms, ctrl_n=n_ok))
        corr = {}
        for t in trs:
            raw = dist_on(sl, t)
            corr[t['id']] = raw - (p[0] * t['n'][0] + p[1] * t['n'][1] + p[2])
        pos[r.file] = corr
    co = pd.DataFrame(rows)

    P = pd.DataFrame(pos).T                               # rows=scenes, cols=transect id
    P.index = meta.set_index('file').loc[P.index, 'date'].dt.date
    P.index.name = 'date'
    P.insert(0, 'sensor', meta['sensor'].values)

    tyr = meta['year_dec'].values
    out = []
    for t in trs:
        x = P[t['id']].values.astype(float)
        base = dict(transect=t['id'], beach=t['beach'], cross_angle_deg=t['angle'], sand_frac_S2=t['sand_frac'],
                    along_m=(t['id'] - 1) * SPACING,
                    lon_lat=None)
        full = rates(tyr, x)
        ls = rates(tyr[meta.sensor != 'S2'], x[(meta.sensor != 'S2').values])
        se = rates(tyr[meta.sensor == 'S2'], x[(meta.sensor == 'S2').values])
        row = {**base, **{f'all_{k}': v for k, v in full.items()},
               **{f'landsat_{k}': v for k, v in ls.items()},
               **{f's2_{k}': v for k, v in se.items()}}
        out.append(row)
    R = pd.DataFrame(out).drop(columns='lon_lat')

    co.to_csv(os.path.join(outdir, 'coregistration.csv'), index=False)
    P.to_csv(os.path.join(outdir, 'shoreline_position_by_transect.csv'))
    R.to_csv(os.path.join(outdir, 'change_rates_by_transect.csv'), index=False)

    crs = ref['crs'].to_string()
    def fc(features): return {'type': 'FeatureCollection',
                              'crs': {'type': 'name', 'properties': {'name': crs}}, 'features': features}
    shift = co.set_index(pd.to_datetime(co['date']))
    feats = []
    for r in meta.itertuples():
        c = shift.loc[r.date]
        xy = np.array(S[r.file]['line'].coords) - np.array([c.shift_x_m, c.shift_y_m])
        feats.append({'type': 'Feature', 'geometry': mapping(LineString(xy)),
                      'properties': {'date': str(r.date.date()), 'sensor': r.sensor,
                                     'note': 'shift-corrected, bias not applied'}})
    json.dump(fc(feats), open(os.path.join(outdir, 'shorelines_utm47n.geojson'), 'w'))
    json.dump(fc([{'type': 'Feature', 'geometry': mapping(t['line']),
                   'properties': {'transect': t['id'], 'beach': bool(t['beach'])}} for t in trs]),
              open(os.path.join(outdir, 'transects_utm47n.geojson'), 'w'))
    return dict(meta=meta, S=S, trs=trs, ctrl=ctrl, co=co, P=P, R=R, ref=ref)


# ---------------------------------------------------------------- figures
def make_figures(o, outdir, backdrop=None):
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    S, meta, trs, P, R, co = o['S'], o['meta'], o['trs'], o['P'], o['R'], o['co']
    bf = backdrop or meta[meta.sensor == 'S2'].file.iloc[-1]
    B = S[bf]; T = B['T']; rgb = np.clip(B['rgb'] / 0.25, 0, 1) ** (1 / 1.3)
    ext = [T.c, T.c + T.a * rgb.shape[1], T.f + T.e * rgb.shape[0], T.f]
    yrs = meta.year_dec.values; norm = colors.Normalize(yrs.min(), yrs.max()); cmap = cm.plasma
    shift = co.set_index(pd.to_datetime(co['date']))
    beach = [t for t in trs if t['beach']]
    allx = np.array([np.array(t['line'].coords) for t in beach]).reshape(-1, 2)
    pad = 60
    xl = (allx[:, 0].min() - pad, allx[:, 0].max() + pad); yl = (allx[:, 1].min() - pad, allx[:, 1].max() + pad)

    # 1 map
    fig, ax = plt.subplots(figsize=(8, 9))
    ax.imshow(rgb, extent=ext, interpolation='nearest')
    for r in meta.sort_values('date').itertuples():
        c = shift.loc[r.date]; xy = np.array(S[r.file]['line'].coords) - [c.shift_x_m, c.shift_y_m]
        ax.plot(xy[:, 0], xy[:, 1], color=cmap(norm(r.year_dec)), lw=1.0 if r.sensor == 'S2' else 0.7,
                alpha=0.9)
    for t in trs:
        x, y = t['line'].xy
        ax.plot(x, y, color='w' if t['beach'] else '0.5', lw=0.5, ls='-' if t['beach'] else ':')
        if t['beach'] and t['id'] % 2 == 1:
            ax.text(x[-1], y[-1], str(t['id']), color='w', fontsize=7, ha='right')
    ax.set_xlim(xl); ax.set_ylim(yl); ax.set_aspect('equal')
    sm = cm.ScalarMappable(norm=norm, cmap=cmap); fig.colorbar(sm, ax=ax, shrink=0.6, label='Year')
    ax.set_title(f'Pantai Kerachut shorelines {int(yrs.min())}-{int(yrs.max())}\n'
                 f'(backdrop {os.path.basename(bf)[9:17]}, UTM 47N)', fontsize=10)
    ax.set_xlabel('Easting (m)'); ax.set_ylabel('Northing (m)')
    plt.tight_layout(); fig.savefig(os.path.join(outdir, 'fig1_shoreline_map.png'), dpi=200); plt.close(fig)

    # 2 rates along the beach
    Rb = R[R.beach]
    fig, ax = plt.subplots(figsize=(10, 4))
    x = Rb.transect.values; w = 0.38
    ax.bar(x - w / 2, Rb.all_LRR_m_yr, w, yerr=Rb.all_LRR_CI95, capsize=2, color='0.35',
           label=f'All years ({int(yrs.min())}-{int(yrs.max())})')
    if 's2_LRR_m_yr' in Rb:
        ax.bar(x + w / 2, Rb.s2_LRR_m_yr, w, yerr=Rb.s2_LRR_CI95, capsize=2, color='tab:orange',
               label='Sentinel-2 era (2018-)')
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xlabel('Transect (south to north, 20 m apart)'); ax.set_ylabel('Shoreline change rate (m/yr)\n+ seaward, - landward')
    ax.legend(fontsize=8); ax.set_xticks(x)
    plt.tight_layout(); fig.savefig(os.path.join(outdir, 'fig2_change_rate_by_transect.png'), dpi=200); plt.close(fig)

    # 3 time series by zone
    ids = Rb.transect.values; k = len(ids) // 3
    zones = {'South': ids[:k], 'Central': ids[k:2 * k], 'North (lake outlet side)': ids[2 * k:]}
    fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    for ax, (zn, zi) in zip(axs, zones.items()):
        y = P[list(zi)].astype(float).mean(1).values
        y = y - np.nanmean(y)
        for sen, mk in [('L5', 'o'), ('L7', 's'), ('L8', '^'), ('L9', 'v'), ('S2', 'D')]:
            m = (meta.sensor == sen).values
            if m.any(): ax.plot(yrs[m], y[m], mk, label=sen, ms=5)
        ok = np.isfinite(y); lr = stats.linregress(yrs[ok], y[ok])
        ax.plot(yrs, lr.intercept + lr.slope * yrs, 'k--', lw=1)
        ax.set_title(f'{zn}: transects {zi[0]}-{zi[-1]}, trend {lr.slope:+.2f} m/yr (r2 {lr.rvalue**2:.2f})', fontsize=9)
        ax.set_ylabel('Position vs mean (m)'); ax.axhline(0, color='0.7', lw=0.5)
    axs[0].legend(ncol=5, fontsize=7); axs[-1].set_xlabel('Year')
    plt.tight_layout(); fig.savefig(os.path.join(outdir, 'fig3_position_timeseries_by_zone.png'), dpi=200); plt.close(fig)
    return zones
