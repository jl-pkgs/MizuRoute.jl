#!/usr/bin/env python3
import argparse, json, math, tarfile
from pathlib import Path
import numpy as np
from netCDF4 import Dataset

# Numerical-Recipes incomplete gamma used by mizuRoute/gamma_func.f90.
def gammln(xx):
    coef = (76.18009172947146, -86.50532032941677, 24.01409824083091,
            -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5)
    x = float(xx)
    tmp = x + 5.5
    tmp = (x + 0.5) * math.log(tmp) - tmp
    s = 1.000000000190015
    for j, c in enumerate(coef, start=1):
        s += c / (x + j)
    return tmp + math.log(2.5066282746310005 * s / x)

def gser(a, x):
    if x == 0.0:
        return 0.0
    ap = a
    summ = 1.0 / a
    delta = summ
    eps = np.finfo(np.float64).eps
    for _ in range(1, 101):
        ap += 1.0
        delta *= x / ap
        summ += delta
        if abs(delta) < abs(summ) * eps:
            break
    return summ * math.exp(-x + a * math.log(x) - gammln(a))

def gcf(a, x):
    if x == 0.0:
        return 1.0
    eps = np.finfo(np.float64).eps
    fpmin = np.finfo(np.float64).tiny / eps
    b = x + 1.0 - a
    c = 1.0 / fpmin
    d = 1.0 / b
    h = d
    for i in range(1, 101):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < fpmin: d = fpmin
        c = b + an / c
        if abs(c) < fpmin: c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) <= eps:
            break
    return math.exp(-x + a * math.log(x) - gammln(a)) * h

def gammp(a, x):
    return gser(a, x) if x < a + 1.0 else 1.0 - gcf(a, x)

def basin_uh(dt, fshape=2.5, tscale=86400.0):
    x = dt / tscale
    cum = gammp(fshape, x)
    if cum > 0.999:
        trial = 1.999
    else:
        lo, hi = 1.0, 1000.0
        trial = 0.5 * (lo + hi)
        for _ in range(100):
            cum = gammp(fshape, dt * trial / tscale)
            if cum < 0.99: lo = trial
            if cum > 0.999: hi = trial
            if 0.99 < cum < 0.999: break
            trial = 0.5 * (lo + hi)
    ntdh = int(math.ceil(trial))
    frac = np.empty(ntdh, dtype=np.float64)
    prev = 0.0
    for j in range(ntdh):
        cum = gammp(fshape, (j + 1) * dt / tscale)
        frac[j] = max(0.0, cum - prev)
        prev = cum
    frac /= frac.sum()
    return frac

p = argparse.ArgumentParser()
p.add_argument('root', help='extracted testCase_cameo_v3.1')
p.add_argument('--out', default='validation/cameo/work')
a = p.parse_args(); root = Path(a.root); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
refdir = out/'reference'; refdir.mkdir(exist_ok=True)
refarc = root/'output'/'ForComparison.tar.gz'
with tarfile.open(refarc, 'r:gz') as tf:
    tf.extractall(refdir)
ref = refdir/'ForComparison'/'case1_1950-1-1.nc'
ntopo = root/'ancillary_data'/'ntopo_nhdplus_cameo_pfaf.nc'
runoff_file = root/'input'/'RUNOFF_case1.nc'

with Dataset(ntopo) as ds:
    area = np.asarray(ds['area'][:], dtype=np.float64)
    hru_id = np.asarray(ds['HRUid'][:], dtype=np.int64)
    hru_seg = np.asarray(ds['hruSegId'][:], dtype=np.int64)
    seg = np.asarray(ds['segId'][:], dtype=np.int64)
    down = np.asarray(ds['downSegId'][:], dtype=np.int64)
    length = np.asarray(ds['length'][:], dtype=np.float64)
    # mizuRoute RPARAM_in uses max(raw slope, min_slope=1e-6).
    slope = np.maximum(np.asarray(ds['slope'][:], dtype=np.float64), 1.0e-6)

idx = {int(x): i for i, x in enumerate(seg)}
hru_to_seg = np.array([idx.get(int(s), -1) for s in hru_seg], dtype=np.int64)
bas = np.zeros(len(seg), dtype=np.float64)
valid_hru = (hru_to_seg >= 0) & (area > 0)
np.add.at(bas, hru_to_seg[valid_hru], area[valid_hru])

indeg = np.zeros(len(seg), dtype=int)
for d in down:
    j = idx.get(int(d))
    if j is not None: indeg[j] += 1
queue = [i for i, x in enumerate(indeg) if x == 0]; order = []; k = 0
while k < len(queue):
    i = queue[k]; k += 1; order.append(i)
    j = idx.get(int(down[i]))
    if j is not None:
        indeg[j] -= 1
        if indeg[j] == 0: queue.append(j)
assert len(order) == len(seg)
total = bas.copy()
for i in order:
    j = idx.get(int(down[i]))
    if j is not None: total[j] += total[i]
active = total > 0.0
print('zero-total-area reaches:', int((~active).sum()))
width = 0.001 * np.sqrt(np.maximum(total, 0.0))
# Geometry is never used to route an inactive reach into an active downstream
# via goodBas, but Julia parameter validation requires a positive width.
width[~active] = 1.0e-6

with Dataset(ref) as ds:
    rid = np.asarray(ds['reachID'][:], dtype=np.int64)
    assert np.array_equal(rid, seg)
    reference_dlay = np.asarray(ds['dlayRunoff'][:], dtype=np.float32)
    irf = np.asarray(ds['IRFroutedRunoff'][:], dtype=np.float32)
    kwt = np.asarray(ds['KWTroutedRunoff'][:], dtype=np.float32)
ntime, nreach = reference_dlay.shape

# Reconstruct the internal double-precision BASIN_QR instead of feeding the
# Float32 history variable dlayRunoff back into Julia. This follows
# get_hru_runoff -> basin2reach -> basinUH/irf_conv for Cameo case1.
with Dataset(runoff_file) as ds:
    input_hru = np.asarray(ds['hru'][:], dtype=np.int64)
    raw = np.asarray(ds['runoff'][:], dtype=np.float64)  # source is Float32, units mm/s
assert raw.shape[0] == ntime
input_ix = {int(x): i for i, x in enumerate(input_hru)}
src = np.array([input_ix.get(int(x), -1) for x in hru_id], dtype=np.int64)
found = src >= 0
if not np.all(found):
    print('network HRUs missing from runoff input:', int((~found).sum()))

inst = np.zeros((ntime, nreach), dtype=np.float64)
for t in range(ntime):
    vals = np.zeros(len(hru_id), dtype=np.float64)
    vals[found] = np.maximum(raw[t, src[found]], 0.0) * 1.0e-3  # mm/s -> m/s
    np.add.at(inst[t], hru_to_seg[valid_hru], vals[valid_hru] * area[valid_hru])

# basin2reach enforces runoffMin=1e-15 m/s for contributing basins and uses
# runoffMin directly for the no-HRU special case.
runoff_min = 1.0e-15
has_basin = bas > 0.0
inst[:, has_basin] = np.maximum(inst[:, has_basin], runoff_min * bas[has_basin])
inst[:, ~has_basin] = runoff_min

frac = basin_uh(86400.0)
print('basin UH:', frac.tolist())
qfuture = np.zeros((nreach, len(frac)), dtype=np.float64)
qlat = np.zeros((ntime, nreach), dtype=np.float64)
for t in range(ntime):
    qfuture += inst[t, :, None] * frac[None, :]
    qlat[t] = qfuture[:, 0]
    qfuture[:, :-1] = qfuture[:, 1:]
    qfuture[:, -1] = 0.0

# History output is Float32, so this comparison diagnoses preparation fidelity
# without contaminating the Julia input with Float32 rounding.
fdiff = qlat - reference_dlay.astype(np.float64)
print('reconstructed dlayRunoff: MAE=', float(np.mean(np.abs(fdiff))),
      'max=', float(np.max(np.abs(fdiff))))

with open(out/'network.csv', 'w') as f:
    f.write('reach_id,downstream_id,length,slope,width,active\n')
    for vals in zip(seg, down, length, slope, width, active.astype(np.int8)):
        f.write(','.join(map(str, vals)) + '\n')
np.ascontiguousarray(qlat, dtype='<f8').tofile(out/'qlat.f64')
np.ascontiguousarray(reference_dlay, dtype='<f4').tofile(out/'reference_dlay.f32')
for name, arr in [('reference_irf', irf), ('reference_kwt', kwt)]:
    np.ascontiguousarray(arr, dtype='<f4').tofile(out/(name+'.f32'))
meta = {'nreach': int(nreach), 'ntime': int(ntime), 'dt': 86400.0,
        'mann_n': 0.01, 'irf_velocity': 1.5, 'irf_diffusivity': 800.0,
        'bankfull_depth': 100000.0, 'floodplain_slope': 1000.0,
        'basin_uh_bins': int(len(frac))}
(out/'meta.json').write_text(json.dumps(meta, indent=2))
print(json.dumps(meta)); print('prepared', out)
