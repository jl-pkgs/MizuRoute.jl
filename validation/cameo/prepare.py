#!/usr/bin/env python3
import argparse, json, tarfile
from pathlib import Path
import numpy as np
from netCDF4 import Dataset

p = argparse.ArgumentParser()
p.add_argument('root', help='extracted testCase_cameo_v3.1')
p.add_argument('--out', default='validation/cameo/work')
a = p.parse_args()
root = Path(a.root)
out = Path(a.out)
out.mkdir(parents=True, exist_ok=True)

# Official NCAR/CESM comparison output is the Fortran system-level oracle.
refdir = out / 'reference'
refdir.mkdir(exist_ok=True)
with tarfile.open(root / 'output' / 'ForComparison.tar.gz', 'r:gz') as tf:
    tf.extractall(refdir)
ref = refdir / 'ForComparison' / 'case1_1950-1-1.nc'
ntopo = root / 'ancillary_data' / 'ntopo_nhdplus_cameo_pfaf.nc'

with Dataset(ntopo) as ds:
    area = np.asarray(ds['area'][:], dtype=np.float64)
    hru_seg = np.asarray(ds['hruSegId'][:], dtype=np.int64)
    seg = np.asarray(ds['segId'][:], dtype=np.int64)
    down = np.asarray(ds['downSegId'][:], dtype=np.int64)
    length = np.asarray(ds['length'][:], dtype=np.float64)
    # RPARAM_in uses max(raw slope, min_slope=1e-6).
    slope = np.maximum(np.asarray(ds['slope'][:], dtype=np.float64), 1.0e-6)

# Reproduce mizuRoute goodBas (>0 total contributing area) only for network
# connectivity semantics. The actual lateral forcing is the official dlayRunoff
# history variable, so basin UH/remapping are deliberately outside this test.
idx = {int(x): i for i, x in enumerate(seg)}
hru_to_seg = np.array([idx.get(int(s), -1) for s in hru_seg], dtype=np.int64)
bas = np.zeros(len(seg), dtype=np.float64)
valid_hru = (hru_to_seg >= 0) & (area > 0)
np.add.at(bas, hru_to_seg[valid_hru], area[valid_hru])

indeg = np.zeros(len(seg), dtype=np.int64)
for d in down:
    j = idx.get(int(d))
    if j is not None:
        indeg[j] += 1
queue = [i for i, x in enumerate(indeg) if x == 0]
order = []
k = 0
while k < len(queue):
    i = queue[k]
    k += 1
    order.append(i)
    j = idx.get(int(down[i]))
    if j is not None:
        indeg[j] -= 1
        if indeg[j] == 0:
            queue.append(j)
assert len(order) == len(seg)

total = bas.copy()
for i in order:
    j = idx.get(int(down[i]))
    if j is not None:
        total[j] += total[i]
active = total > 0.0
print('zero-total-area reaches:', int((~active).sum()))

width = 0.001 * np.sqrt(np.maximum(total, 0.0))
# Inactive reaches are excluded by goodBas. A tiny placeholder width is only
# needed to satisfy ReachParameters validation; it cannot affect active routing.
width[~active] = 1.0e-6

with Dataset(ref) as ds:
    rid = np.asarray(ds['reachID'][:], dtype=np.int64)
    assert np.array_equal(rid, seg)
    dlay = np.asarray(ds['dlayRunoff'][:], dtype=np.float32)
    irf = np.asarray(ds['IRFroutedRunoff'][:], dtype=np.float32)
    kwt = np.asarray(ds['KWTroutedRunoff'][:], dtype=np.float32)
ntime, nreach = dlay.shape

with open(out / 'network.csv', 'w') as f:
    f.write('reach_id,downstream_id,length,slope,width,active\n')
    for vals in zip(seg, down, length, slope, width, active.astype(np.int8)):
        f.write(','.join(map(str, vals)) + '\n')

# Deliberately retain the official Float32 history values. Feeding exactly the
# same reach-local forcing to Julia isolates river-routing differences from
# basin-UH/remapping/version differences upstream of the river network.
np.ascontiguousarray(dlay, dtype='<f4').tofile(out / 'qlat.f32')
for name, arr in [('reference_irf', irf), ('reference_kwt', kwt)]:
    np.ascontiguousarray(arr, dtype='<f4').tofile(out / (name + '.f32'))

meta = {
    'nreach': int(nreach), 'ntime': int(ntime), 'dt': 86400.0,
    'mann_n': 0.01, 'irf_velocity': 1.5, 'irf_diffusivity': 800.0,
    'bankfull_depth': 100000.0, 'floodplain_slope': 1000.0,
    'forcing': 'official ForComparison dlayRunoff (Float32 history output)'
}
(out / 'meta.json').write_text(json.dumps(meta, indent=2))
print(json.dumps(meta))
print('prepared', out)
