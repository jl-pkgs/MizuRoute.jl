#!/usr/bin/env python3
import argparse,json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser(); p.add_argument('--work',default='validation/cameo/work'); a=p.parse_args(); w=Path(a.work)
meta=json.loads((w/'meta.json').read_text()); nr=meta['nreach']; nt=meta['ntime']
rows=[]
for name in ('irf','kwt'):
    ref=np.fromfile(w/f'reference_{name}.f32',dtype='<f4').astype(np.float64).reshape(nt,nr)
    sim=np.fromfile(w/f'julia_{name}.f64',dtype='<f8').reshape((nr,nt),order='F').T
    e=sim-ref; mask=np.isfinite(ref)&np.isfinite(sim); ee=e[mask]; rr=ref[mask]
    rel=np.abs(ee)/np.maximum(np.abs(rr),1e-6)
    flat=np.argmax(np.abs(e)); t,r=np.unravel_index(flat,e.shape)
    row={'method':name,'n':int(mask.sum()),'mae':float(np.mean(np.abs(ee))),'rmse':float(np.sqrt(np.mean(ee**2))),
         'bias':float(np.mean(ee)),'max_abs':float(np.max(np.abs(ee))),'p99_abs':float(np.quantile(np.abs(ee),.99)),
         'median_rel':float(np.median(rel)),'p99_rel':float(np.quantile(rel,.99)),
         'worst_timestep':int(t+1),'worst_reach_index':int(r+1),'reference_at_worst':float(ref[t,r]),'julia_at_worst':float(sim[t,r])}
    rows.append(row)
print(json.dumps(rows,indent=2))
(w/'summary.json').write_text(json.dumps(rows,indent=2))
with open(w/'report.md','w') as f:
    f.write('# Cameo v3.1 reach × timestep regression\n\n')
    f.write('| method | MAE | RMSE | max abs | p99 abs | p99 rel | worst reach index | worst timestep |\n|---|---:|---:|---:|---:|---:|---:|---:|\n')
    for x in rows: f.write(f"| {x['method']} | {x['mae']:.6g} | {x['rmse']:.6g} | {x['max_abs']:.6g} | {x['p99_abs']:.6g} | {x['p99_rel']:.6g} | {x['worst_reach_index']} | {x['worst_timestep']} |\n")