#!/usr/bin/env python3
import argparse, json, tarfile
from pathlib import Path
import numpy as np
from netCDF4 import Dataset

p=argparse.ArgumentParser()
p.add_argument('root', help='extracted testCase_cameo_v3.1')
p.add_argument('--out', default='validation/cameo/work')
a=p.parse_args(); root=Path(a.root); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
refdir=out/'reference'; refdir.mkdir(exist_ok=True)
refarc=root/'output'/'ForComparison.tar.gz'
with tarfile.open(refarc,'r:gz') as tf:
    tf.extractall(refdir)
ref=refdir/'ForComparison'/'case1_1950-1-1.nc'
ntopo=root/'ancillary_data'/'ntopo_nhdplus_cameo_pfaf.nc'
with Dataset(ntopo) as ds:
    area=np.asarray(ds['area'][:],dtype=np.float64)
    hru_seg=np.asarray(ds['hruSegId'][:],dtype=np.int64)
    seg=np.asarray(ds['segId'][:],dtype=np.int64)
    down=np.asarray(ds['downSegId'][:],dtype=np.int64)
    length=np.asarray(ds['length'][:],dtype=np.float64)
    slope=np.asarray(ds['slope'][:],dtype=np.float64)
idx={int(x):i for i,x in enumerate(seg)}
bas=np.zeros(len(seg),dtype=np.float64)
for aa,s in zip(area,hru_seg):
    j=idx.get(int(s))
    if j is not None and aa>0: bas[j]+=aa
indeg=np.zeros(len(seg),dtype=int)
for i,d in enumerate(down):
    j=idx.get(int(d))
    if j is not None: indeg[j]+=1
queue=[i for i,x in enumerate(indeg) if x==0]; order=[]; k=0
while k<len(queue):
    i=queue[k]; k+=1; order.append(i)
    j=idx.get(int(down[i]))
    if j is not None:
        indeg[j]-=1
        if indeg[j]==0: queue.append(j)
assert len(order)==len(seg)
total=bas.copy()
for i in order:
    j=idx.get(int(down[i]))
    if j is not None: total[j]+=total[i]
width=0.001*np.sqrt(np.maximum(total,0.0))
with Dataset(ref) as ds:
    rid=np.asarray(ds['reachID'][:],dtype=np.int64)
    assert np.array_equal(rid,seg)
    qlat=np.asarray(ds['dlayRunoff'][:],dtype=np.float32)
    irf=np.asarray(ds['IRFroutedRunoff'][:],dtype=np.float32)
    kwt=np.asarray(ds['KWTroutedRunoff'][:],dtype=np.float32)
ntime,nreach=qlat.shape
with open(out/'network.csv','w') as f:
    f.write('reach_id,downstream_id,length,slope,width\n')
    for vals in zip(seg,down,length,slope,width): f.write(','.join(map(str,vals))+'\n')
for name,arr in [('qlat',qlat),('reference_irf',irf),('reference_kwt',kwt)]:
    np.ascontiguousarray(arr,dtype='<f4').tofile(out/(name+'.f32'))
meta={'nreach':int(nreach),'ntime':int(ntime),'dt':86400.0,'mann_n':0.01,'irf_velocity':1.5,'irf_diffusivity':800.0,'bankfull_depth':100000.0,'floodplain_slope':1000.0}
(out/'meta.json').write_text(json.dumps(meta,indent=2))
print(json.dumps(meta)); print('prepared',out)