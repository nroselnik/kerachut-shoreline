"""
Kerachut: tide-corrected shoreline change.

Uses the shoreline extraction and co-registration from kerachut_shoreline.py on
every reasonably clear scene, on the SAME transects as the first analysis.

Steps
1. Extract shoreline positions for all scenes on fixed transects; co-register to
   the reference Sentinel-2 scene with rocky-coast control transects; QC.
2. Tide height at each acquisition time (EOT20 via eo-tides / pyTMD).
3. Beach slope per zone from a regression of shoreline position on tide height,
   with one intercept per year to absorb real change (x = a_year + c*tide;
   tan(beta) = -1/c). Bootstrap 95% CI.
4. Three treatments: raw, tide-corrected to MSL (x + tide/tan(beta)),
   tide-filtered (|tide| <= 0.25 m, no slope needed).
5. Rates (LRR) per transect and zone, rotation index (north minus south),
   sensitivity of zone rates to an assumed slope range.
"""
import os, re, json, glob
import numpy as np, pandas as pd
from scipy import stats
from shapely.geometry import LineString
import kerachut_shoreline as k

ZONES = {'South': list(range(2, 11)), 'Central': list(range(11, 20)), 'North': list(range(20, 27))}
FILTER_M = 0.25
SLOPES = [0.04, 0.06, 0.08, 0.10, 0.12, 0.15, 0.20]


# ------------------------------------------------------------------ transects
def load_transects(path, length=k.LENGTH):
    fc = json.load(open(path)); out = []
    for f in fc['features']:
        c = np.array(f['geometry']['coordinates']); o = c[0]
        n = (c[1] - c[0]) / np.linalg.norm(c[1] - c[0])
        out.append(dict(id=f['properties']['transect'], beach=bool(f['properties']['beach']),
                        origin=o, n=n, line=LineString([o, o + length * n])))
    return out


# ------------------------------------------------------------------ extraction
def extract(folder, transect_file, times_csv=None, ref=k.REF, exclude=k.EXCLUDE,
            max_rms=8.0, min_ctrl=40, max_shift=40.0, verbose=True):
    trs = load_transects(transect_file)
    files = sorted(f for f in glob.glob(os.path.join(folder, 'Kerachut_*.tif'))
                   if os.path.basename(f)[:-4] not in exclude)
    rf = [f for f in files if ref in f]
    if not rf: raise FileNotFoundError(f'reference scene {ref} not in {folder}')
    R = k.shoreline(rf[0])
    beach_xy = k.px2xy(R['T'], np.argwhere(k.biggest(R['sand'])))
    ctrl = k.control_transects(R, beach_xy)
    times = None
    if times_csv and os.path.exists(times_csv):
        times = pd.read_csv(times_csv).set_index('file')['datetime_utc']
    rows, pos = [], {}
    for i, f in enumerate(files):
        name = os.path.basename(f)
        m = re.search(r'Kerachut_(\d{8})_(\w+?)\.tif', name)
        if not m: continue
        try:
            S = k.shoreline(f)
        except Exception as e:
            rows.append(dict(file=name, status=f'fail: {e}')); continue
        off = np.array([k.dist_on(S['line'], c, 'nearest') for c in ctrl])
        if np.isfinite(off).sum() < min_ctrl:
            rows.append(dict(file=name, status='few controls', ctrl_n=int(np.isfinite(off).sum()))); continue
        p, rms, n_ok = k.fit_bias(off, ctrl)
        status = 'ok'
        if rms > max_rms: status = 'ctrl rms'
        if np.hypot(p[0], p[1]) > max_shift: status = 'shift'
        dt = pd.to_datetime(times[name]) if (times is not None and name in times) else pd.to_datetime(m.group(1))
        rows.append(dict(file=name, date=dt, sensor=m.group(2), threshold=S['thr'], shift_x_m=p[0],
                         shift_y_m=p[1], bias_m=p[2], ctrl_rms_m=rms, ctrl_n=n_ok, status=status))
        if status == 'ok':
            pos[name] = {t['id']: k.dist_on(S['line'], t) - (p[0] * t['n'][0] + p[1] * t['n'][1] + p[2])
                         for t in trs}
        if verbose and i % 50 == 0: print(f'  {i}/{len(files)} scenes')
    Q = pd.DataFrame(rows)
    P = pd.DataFrame(pos).T
    P.index.name = 'file'
    meta = Q[Q.status == 'ok'].set_index('file')[['date', 'sensor']]
    P = meta.join(P, how='inner').sort_values('date')
    P = robust_clean(P, trs)
    return P, Q, trs


def robust_clean(P, trs, z=4.0):
    """Drop single positions far from the scene's own along-beach pattern and from
    the transect's local temporal neighbours (cloud edges, boats, foam)."""
    ids = [t['id'] for t in trs]
    X = P[ids].astype(float)
    med = X.rolling(9, center=True, min_periods=3).median()
    res = X - med
    mad = 1.4826 * res.abs().median()
    bad = res.abs() > z * mad.clip(lower=3)
    P = P.copy(); P[ids] = X.mask(bad)
    return P


# ------------------------------------------------------------------ tides
def tides_eot20(times_utc, lon, lat, directory):
    """Tide height (m, relative to model MSL) at given UTC times."""
    times = pd.to_datetime(pd.Series(times_utc)).dt.tz_localize(None)
    try:
        from eo_tides.model import model_tides
        df = model_tides(x=lon, y=lat, time=times.values, model='EOT20', directory=directory,
                         parallel=False).reset_index()
        s = df.groupby('time')['tide_height'].first()
        return np.asarray(s.reindex(pd.DatetimeIndex(times.values)).values, float)
    except Exception as e1:
        try:
            import pyTMD
            delta = (times - pd.Timestamp('1992-01-01')).dt.total_seconds().values / 86400.0
            h = pyTMD.compute.tide_elevations(np.full(len(times), lon), np.full(len(times), lat), delta,
                                              DIRECTORY=directory, MODEL='EOT20', EPOCH=(1992, 1, 1, 0, 0, 0),
                                              TYPE='drift', TIME='UTC', EPSG=4326, EXTRAPOLATE=True)
            return np.asarray(h, float)
        except Exception as e2:
            raise RuntimeError(f'eo-tides failed: {e1}\npyTMD failed: {e2}')


# ------------------------------------------------------------------ slope
def zone_mean(P, ids):
    cols = [i for i in ids if i in P.columns]
    X = P[cols].astype(float)
    # require at least half the transects present
    ok = X.notna().sum(1) >= max(2, len(cols) // 2)
    return X.mean(1).where(ok)

def fit_slope(x, tide, year, n_boot=2000, seed=1):
    """x = a_year + c*tide. Returns tan(beta) = -1/c with bootstrap 95% CI."""
    d = pd.DataFrame(dict(x=x, tide=tide, year=year)).dropna()
    d = d[d.groupby('year').x.transform('size') >= 2]
    def c_of(dd):
        xd = dd.x - dd.groupby('year').x.transform('mean')
        td = dd.tide - dd.groupby('year').tide.transform('mean')
        return float((xd * td).sum() / (td * td).sum()), xd, td
    c, xd, td = c_of(d)
    rng = np.random.default_rng(seed); cs = []
    idx = np.arange(len(d))
    for _ in range(n_boot):
        cs.append(c_of(d.iloc[rng.choice(idx, len(idx))])[0])
    cs = np.array(cs)
    r = stats.pearsonr(td, xd)
    out = dict(n=len(d), n_years=d.year.nunique(), c_m_per_m=c, c_lo=np.percentile(cs, 2.5),
               c_hi=np.percentile(cs, 97.5), r=r[0], p=r[1])
    out['tanb'] = -1 / c if c < 0 else np.nan
    out['tanb_lo'] = -1 / out['c_lo'] if out['c_lo'] < 0 else np.nan     # steeper bound
    out['tanb_hi'] = -1 / out['c_hi'] if out['c_hi'] < 0 else np.nan     # gentler bound
    return out, xd, td


# ------------------------------------------------------------------ rates
def lrr(t, y):
    ok = np.isfinite(y)
    if ok.sum() < 4: return dict(n=int(ok.sum()))
    lr = stats.linregress(t[ok], y[ok])
    return dict(n=int(ok.sum()), LRR=lr.slope, CI95=stats.t.ppf(.975, ok.sum() - 2) * lr.stderr,
                r2=lr.rvalue ** 2, p=lr.pvalue)

def decyear(d):
    d = pd.to_datetime(d)
    return (d.dt.year + (d.dt.dayofyear - 0.5) / 365.25).values


def analyse(P, tide, trs, outdir, slope_from='S2era', veg_csv=None):
    os.makedirs(outdir, exist_ok=True)
    P = P.copy(); P['tide_m'] = tide
    P['t'] = decyear(P['date']); P['year'] = pd.to_datetime(P['date']).dt.year
    ids = [t['id'] for t in trs if t['beach']]
    P.to_csv(os.path.join(outdir, 'positions_all_scenes.csv'))

    # --- slope per zone (2013 onward: dense Landsat 8/9 + Sentinel-2) and pooled
    dense = P[P['t'] >= 2013]
    S, fitdata = [], {}
    for z, zi in ZONES.items():
        res, xd, td = fit_slope(zone_mean(dense, zi), dense['tide_m'], dense['year'])
        res['zone'] = z; S.append(res); fitdata[z] = (xd, td, res)
    pooled = pd.concat([zone_mean(dense, zi) - zone_mean(dense, zi).mean() for zi in ZONES.values()])
    res, _, _ = fit_slope(pooled.values, np.tile(dense['tide_m'].values, 3), np.tile(dense['year'].values, 3))
    res['zone'] = 'All (pooled)'; S.append(res)
    S = pd.DataFrame(S).set_index('zone'); S.to_csv(os.path.join(outdir, 'beach_slope_estimates.csv'))
    tanb = {z: (S.loc[z, 'tanb'] if np.isfinite(S.loc[z, 'tanb']) else S.loc['All (pooled)', 'tanb']) for z in ZONES}
    zone_of = {i: z for z, zi in ZONES.items() for i in zi}

    # --- treatments
    C = P.copy()
    for i in ids:
        tb = tanb.get(zone_of.get(i), S.loc['All (pooled)', 'tanb'])
        C[i] = P[i] + P['tide_m'] / tb
    F = P[P['tide_m'].abs() <= FILTER_M]
    treatments = {'raw': P, 'tide_corrected': C, 'tide_filtered': F}

    # --- rates per transect
    rows = []
    for name, D in treatments.items():
        for per, m in [('1988-2026', D['t'] > 0), ('1988-2017', D['t'] < 2018), ('2018-2026', D['t'] >= 2018)]:
            for i in ids:
                r = lrr(D['t'][m].values, D[i][m].values.astype(float))
                rows.append(dict(treatment=name, period=per, transect=i, zone=zone_of.get(i), **r))
    RT = pd.DataFrame(rows); RT.to_csv(os.path.join(outdir, 'rates_by_transect_all_treatments.csv'), index=False)

    # --- zone rates and rotation index
    rows, series = [], {}
    for name, D in treatments.items():
        zm = pd.DataFrame({z: zone_mean(D, zi) for z, zi in ZONES.items()})
        zm['rotation_N_minus_S'] = zm['North'] - zm['South']
        zm['t'] = D['t']; series[name] = zm
        for per, m in [('1988-2026', zm['t'] > 0), ('1988-2017', zm['t'] < 2018), ('2018-2026', zm['t'] >= 2018)]:
            for col in list(ZONES) + ['rotation_N_minus_S']:
                r = lrr(zm['t'][m].values, zm[col][m].values)
                rows.append(dict(treatment=name, period=per, series=col, **r))
    RZ = pd.DataFrame(rows); RZ.to_csv(os.path.join(outdir, 'rates_by_zone_all_treatments.csv'), index=False)

    # --- sensitivity to assumed slope
    rows = []
    for tb in SLOPES:
        D = P.copy()
        for i in ids: D[i] = P[i] + P['tide_m'] / tb
        for per, m in [('1988-2026', D['t'] > 0), ('2018-2026', D['t'] >= 2018)]:
            for z, zi in ZONES.items():
                r = lrr(D['t'][m].values, zone_mean(D, zi)[m].values)
                rows.append(dict(tanb=tb, period=per, zone=z, **r))
            rot = zone_mean(D, ZONES['North']) - zone_mean(D, ZONES['South'])
            r = lrr(D['t'][m].values, rot[m].values)
            rows.append(dict(tanb=tb, period=per, zone='rotation_N_minus_S', **r))
    SE = pd.DataFrame(rows); SE.to_csv(os.path.join(outdir, 'sensitivity_to_slope.csv'), index=False)

    figures(P, C, F, series, fitdata, S, SE, outdir, veg_csv)
    return dict(P=P, C=C, F=F, S=S, RT=RT, RZ=RZ, SE=SE, series=series)


# ------------------------------------------------------------------ figures
def figures(P, C, F, series, fitdata, S, SE, outdir, veg_csv=None):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    col = {'South': '#c0392b', 'Central': '#7f8c8d', 'North': '#2471a3'}

    # fig6: tide vs detrended position
    fig, axs = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    for ax, (z, (xd, td, r)) in zip(axs, fitdata.items()):
        ax.plot(td, xd, '.', color=col[z], alpha=.5, ms=4)
        xx = np.linspace(td.min(), td.max(), 10); ax.plot(xx, r['c_m_per_m'] * xx, 'k-', lw=1)
        tb = r['tanb']
        ax.set_title(f"{z}: tanβ = {tb:.3f} (1:{1/tb:.0f})\nr = {r['r']:.2f}, n = {r['n']}" if np.isfinite(tb)
                     else f"{z}: no tide signal (r = {r['r']:.2f})", fontsize=9)
        ax.set_xlabel('Tide height minus annual mean (m)'); ax.axhline(0, color='.8', lw=.5)
    axs[0].set_ylabel('Waterline position minus annual mean (m)\n+ seaward')
    plt.tight_layout(); fig.savefig(os.path.join(outdir, 'fig6_tide_vs_position.png'), dpi=200); plt.close(fig)

    # fig7: corrected zone time series with vegetation line
    V = None
    if veg_csv and os.path.exists(veg_csv):
        V = pd.read_csv(veg_csv, index_col=0); V.index = pd.to_datetime(V.index)
        V.columns = [int(c) for c in V.columns]
    zc = series['tide_corrected']
    fig, axs = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    for ax, z in zip(axs, ZONES):
        s2 = (C['sensor'] == 'S2').values
        ax.plot(zc['t'][~s2], zc[z][~s2], 'o', ms=3, color=col[z], alpha=.45, label='Waterline, Landsat (tide-corrected)')
        ax.plot(zc['t'][s2], zc[z][s2], 'D', ms=3, color=col[z], alpha=.8, label='Waterline, Sentinel-2 (tide-corrected)')
        y = zc[z].values; ok = np.isfinite(y)
        yr = pd.Series(y, index=np.floor(zc['t'])).groupby(level=0).mean()
        ax.plot(yr.index + 0.5, yr.values, '-', color=col[z], lw=1.5, label='Annual mean')
        if V is not None:
            vz = V[[i for i in ZONES[z] if i in V.columns]].mean(1)
            vt = V.index.year + (V.index.dayofyear - .5) / 365.25
            ax.plot(vt, vz, 's', color='#1e8449', ms=5, label='Vegetation line (Google Earth)')
        ax.set_ylabel('Distance from baseline (m)'); ax.set_title(f'{z} section', fontsize=9)
    axs[0].legend(fontsize=7, ncol=2); axs[-1].set_xlabel('Year')
    plt.tight_layout(); fig.savefig(os.path.join(outdir, 'fig7_corrected_timeseries_with_vegetation.png'), dpi=200); plt.close(fig)

    # fig8: rotation index
    fig, ax = plt.subplots(figsize=(10, 3.8))
    for name, mk, c in [('raw', '.', '0.6'), ('tide_corrected', 'o', 'k')]:
        zs = series[name]
        ax.plot(zs['t'], zs['rotation_N_minus_S'], mk, ms=3, color=c, alpha=.6, label=name.replace('_', ' '))
    zs = series['tide_corrected']
    yr = zs.groupby(np.floor(zs['t']))['rotation_N_minus_S'].mean()
    ax.plot(yr.index + .5, yr.values, '-', color='#8e44ad', lw=2, label='Annual mean (corrected)')
    ax.axhline(0, color='.7', lw=.6)
    ax.set_ylabel('North minus South (m)\n+ north end further seaward'); ax.set_xlabel('Year')
    ax.set_title('Beach rotation index', fontsize=9); ax.legend(fontsize=8)
    plt.tight_layout(); fig.savefig(os.path.join(outdir, 'fig8_rotation_index.png'), dpi=200); plt.close(fig)

    # fig9: sensitivity
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
    for ax, per in zip(axs, ['1988-2026', '2018-2026']):
        for z in list(ZONES) + ['rotation_N_minus_S']:
            d = SE[(SE.period == per) & (SE.zone == z)]
            ax.errorbar(d.tanb, d.LRR, yerr=d.CI95, marker='o', ms=4, capsize=2,
                        color=col.get(z, '#8e44ad'), label=z.replace('_', ' '))
        for z in ZONES:
            if np.isfinite(S.loc[z, 'tanb']): ax.axvline(S.loc[z, 'tanb'], color=col[z], ls=':', lw=1)
        ax.axhline(0, color='.6', lw=.6); ax.set_xlabel('Assumed beach slope tanβ'); ax.set_title(per, fontsize=9)
    axs[0].set_ylabel('Trend (m/yr)'); axs[1].legend(fontsize=7)
    plt.tight_layout(); fig.savefig(os.path.join(outdir, 'fig9_slope_sensitivity.png'), dpi=200); plt.close(fig)


# ------------------------------------------------------------------ gauge check
def validate_with_gauge(gauge_df, lon, lat, directory):
    """gauge_df: columns datetime_utc, level_m. Returns RMSE, r, amplitude ratio."""
    g = gauge_df.dropna().copy()
    g['anom'] = g['level_m'] - g['level_m'].mean()
    m = tides_eot20(g['datetime_utc'], lon, lat, directory)
    m = m - np.nanmean(m)
    ok = np.isfinite(m)
    rmse = float(np.sqrt(np.mean((g['anom'].values[ok] - m[ok]) ** 2)))
    r = float(np.corrcoef(g['anom'].values[ok], m[ok])[0, 1])
    amp = float(np.std(m[ok]) / np.std(g['anom'].values[ok]))
    return dict(n=int(ok.sum()), rmse_m=rmse, r=r, model_to_gauge_std_ratio=amp), m
