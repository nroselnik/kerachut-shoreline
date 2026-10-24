"""
Colab driver for the Kerachut tide-corrected analysis.
Functions: fetch_all_scenes (Earth Engine), setup_eot20 (tide model), gauge_1992 (validation data).
"""
import os, io, glob, zipfile, subprocess
import numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor

X0, Y0 = 630270, 603540          # top-left of the study box, UTM 47N (same grid as the first exports)
TIDE_LON, TIDE_LAT = 100.1790, 5.4505   # point just offshore of the beach


def fetch_all_scenes(ee, outdir, min_clear=0.90, workers=8, start='1987-01-01', end='2026-12-31'):
    """Download every scene with >= min_clear cloud-free pixels over the study box as a
    5-band GeoTIFF (blue, green, red, nir, swir1 reflectance) on a fixed UTM grid.
    Skips files already present, so it can be re-run after an interruption."""
    os.makedirs(outdir, exist_ok=True)
    roi = ee.Geometry.Rectangle([630270, 601740, 632070, 603540], 'EPSG:32647', False)

    def s2_prep(img):
        clear = img.select('SCL').remap([0, 1, 3, 8, 9, 10], [0] * 6, 1).rename('c')
        return img.set('clear_frac', clear.reduceRegion(ee.Reducer.mean(), roi, 10, maxPixels=1e8).get('c'))

    def ls_prep(img):
        clear = img.select('QA_PIXEL').bitwiseAnd(0b11111).eq(0).rename('c')
        return img.set('clear_frac', clear.reduceRegion(ee.Reducer.mean(), roi, 30, maxPixels=1e8).get('c'))

    cols = [('COPERNICUS/S2_SR_HARMONIZED', 'S2', s2_prep),
            ('LANDSAT/LT05/C02/T1_L2', 'L5', ls_prep), ('LANDSAT/LE07/C02/T1_L2', 'L7', ls_prep),
            ('LANDSAT/LC08/C02/T1_L2', 'L8', ls_prep), ('LANDSAT/LC09/C02/T1_L2', 'L9', ls_prep)]
    rows = []
    for cid, sensor, prep in cols:
        c = ee.ImageCollection(cid).filterBounds(roi).filterDate(start, end).map(prep) \
              .filter(ee.Filter.gte('clear_frac', min_clear))
        idx = c.aggregate_array('system:index').getInfo()
        tms = c.aggregate_array('system:time_start').getInfo()
        cf = c.aggregate_array('clear_frac').getInfo()
        for i, t, f in zip(idx, tms, cf):
            rows.append(dict(cid=cid, img_index=i, sensor=sensor, time=pd.to_datetime(t, unit='ms'), clear=f))
        print(f'{sensor}: {len(idx)} scenes >= {min_clear:.0%} clear')
    L = pd.DataFrame(rows)
    L = L[~((L.sensor == 'L7') & (L.time > '2003-05-31'))]          # Landsat 7 SLC-off stripes
    L['file'] = [f"Kerachut_{t:%Y%m%d}_{s}.tif" for t, s in zip(L.time, L.sensor)]
    L = L.sort_values('clear', ascending=False).drop_duplicates('file').sort_values('time')

    def grid(scale, w, h):
        return {'dimensions': {'width': w, 'height': h},
                'affineTransform': {'scaleX': scale, 'shearX': 0, 'translateX': X0,
                                    'shearY': 0, 'scaleY': -scale, 'translateY': Y0},
                'crsCode': 'EPSG:32647'}

    def fetch(r):
        fn = os.path.join(outdir, r.file)
        if os.path.exists(fn) and os.path.getsize(fn) > 1000: return 'skip'
        img = ee.Image(f'{r.cid}/{r.img_index}')
        if r.sensor == 'S2':
            im = img.select(['B2', 'B3', 'B4', 'B8', 'B11']).divide(10000).toFloat(); g = grid(10, 179, 178)
        else:
            b = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5'] if r.sensor in ('L5', 'L7') else \
                ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6']
            im = img.select(b).multiply(0.0000275).add(-0.2).toFloat(); g = grid(30, 60, 60)
        for attempt in range(3):
            try:
                data = ee.data.computePixels({'expression': im, 'fileFormat': 'GEO_TIFF', 'grid': g})
                open(fn, 'wb').write(data); return 'ok'
            except Exception as e:
                err = str(e)
        return 'fail: ' + err[:120]

    with ThreadPoolExecutor(workers) as ex:
        res = list(ex.map(fetch, L.itertuples()))
    L['download'] = res
    L['datetime_utc'] = L.time.dt.strftime('%Y-%m-%d %H:%M:%S')
    L[['file', 'sensor', 'datetime_utc', 'clear', 'download']].to_csv(os.path.join(outdir, '_scene_times.csv'), index=False)
    print(pd.Series([x.split(':')[0] for x in res]).value_counts().to_string())
    return L


def setup_eot20(model_dir='/content/tide_models'):
    """Download EOT20 (about 2 GB, open access, SEANOE doi:10.17882/79489) and arrange it as
    <model_dir>/EOT20/ocean_tides/*.nc for eo-tides / pyTMD."""
    os.makedirs(model_dir, exist_ok=True)
    if glob.glob(os.path.join(model_dir, 'EOT20', 'ocean_tides', '*.nc')):
        return model_dir
    zp = os.path.join(model_dir, 'eot20.zip')
    if not os.path.exists(zp):
        subprocess.run(['wget', '-q', '-O', zp, 'https://www.seanoe.org/data/00683/79489/data/85762.zip'], check=True)
    with zipfile.ZipFile(zp) as z:
        z.extractall(model_dir)
    # nested zips are possible: extract any inner zip that holds ocean tides
    for inner in glob.glob(os.path.join(model_dir, '**', '*.zip'), recursive=True):
        if inner != zp and 'ocean' in os.path.basename(inner).lower():
            with zipfile.ZipFile(inner) as z: z.extractall(os.path.dirname(inner))
    oc = [d for d in glob.glob(os.path.join(model_dir, '**', 'ocean_tides'), recursive=True) if os.path.isdir(d)]
    if not oc: raise RuntimeError('ocean_tides folder not found after unzip; list ' + model_dir)
    target = os.path.join(model_dir, 'EOT20')
    parent = os.path.dirname(oc[0])
    if os.path.abspath(parent) != os.path.abspath(target):
        os.makedirs(target, exist_ok=True)
        link = os.path.join(target, 'ocean_tides')
        if not os.path.exists(link): os.symlink(oc[0], link)
    print('EOT20 ready:', len(glob.glob(os.path.join(target, 'ocean_tides', '*.nc'))), 'constituent files')
    return model_dir


def gauge_1992():
    """Hourly sea level, Penang tide gauge, 1992 (UHSLC research quality data via PANGAEA,
    doi:10.1594/PANGAEA.439979). Returns DataFrame datetime_utc, level_m, or None."""
    import requests
    try:
        txt = requests.get('https://doi.pangaea.de/10.1594/PANGAEA.439979?format=textfile', timeout=60).text
        body = txt.split('*/', 1)[1].strip()
        df = pd.read_csv(io.StringIO(body), sep='\t')
        tcol = [c for c in df.columns if 'date' in c.lower()][0]
        vcol = [c for c in df.columns if c != tcol and pd.api.types.is_numeric_dtype(df[c])][0]
        v = df[vcol].astype(float)
        if '[mm]' in vcol or v.abs().median() > 100: v = v / 1000.0
        out = pd.DataFrame({'datetime_utc': pd.to_datetime(df[tcol]), 'level_m': v})
        print(f'Gauge 1992: {len(out)} hourly values, column "{vcol}"')
        return out
    except Exception as e:
        print('Gauge data not available:', e); return None
