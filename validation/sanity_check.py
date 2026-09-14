from __future__ import annotations
from pathlib import Path
import math
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

def tdma(lower, diag, upper, rhs):
    lower=np.array(lower,float); d=np.array(diag,float); upper=np.array(upper,float); r=np.array(rhs,float)
    n=len(d)
    for i in range(1,n):
        if abs(d[i-1]) <= 1e-14: raise AssertionError('singular TDMA')
        m=lower[i]/d[i-1]
        d[i]-=m*upper[i-1]
        r[i]-=m*r[i-1]
    x=np.empty(n,float); x[-1]=r[-1]/d[-1]
    for i in range(n-2,-1,-1): x[i]=(r[i]-upper[i]*x[i+1])/d[i]
    return x

def solve_ade(qprev,L,dt,qup,ck,dk,qlat=0,wc=1,wd=1):
    qprev=np.asarray(qprev,float); n=len(qprev); assert n>=4
    nx=n-1; dx=L/(nx-1); cd=max(dk,0)*dt/dx**2; ca=max(ck,0)*dt/dx
    lo=np.zeros(n); di=np.zeros(n); up=np.zeros(n); rhs=np.zeros(n)
    di[0]=1; rhs[0]=max(qup,0)
    src=2*dt*max(ck,0)*qlat/L
    for j in range(1,n-1):
        lo[j]=-wc*ca-2*wd*cd; di[j]=2+4*wd*cd; up[j]=wc*ca-2*wd*cd
        rhs[j]=((1-wc)*ca+2*(1-wd)*cd)*qprev[j-1]+(2-4*(1-wd)*cd)*qprev[j]-((1-wc)*ca-2*(1-wd)*cd)*qprev[j+1]+src
    lo[-1]=-1; di[-1]=1; rhs[-1]=qprev[-1]-qprev[-2]
    return np.maximum(tdma(lo,di,up,rhs),0)

def flow_area(y,b,z): return y*(b+z*y)
def water_height(a,b,z): return a/b if z==0 else (-b+math.sqrt(b*b+4*z*a))/(2*z)
def wetted(y,b,z): return b+2*y*math.hypot(1,z)
def manning(y,b,z,S,n):
    A=flow_area(y,b,z); R=A/wetted(y,b,z); return A*R**(2/3)*math.sqrt(S)/n

def flow_depth(q,b,z,S,n):
    lo,hi=0.,1.
    while manning(hi,b,z,S,n)<q: hi*=2
    for _ in range(100):
        mid=(lo+hi)/2
        if manning(mid,b,z,S,n)<q: lo=mid
        else: hi=mid
    return (lo+hi)/2

def bracket_static_check(path: Path):
    text=path.read_text(encoding='utf-8')
    cleaned=[]; in_string=False; triple=False; i=0
    while i < len(text):
        if not in_string and text.startswith('"""',i): in_string=True; triple=True; i+=3; continue
        if in_string and triple and text.startswith('"""',i): in_string=False; triple=False; i+=3; continue
        ch=text[i]
        if not in_string and ch=='"': in_string=True; i+=1; continue
        if in_string and not triple and ch=='"' and (i==0 or text[i-1] != '\\'): in_string=False; i+=1; continue
        if in_string: i+=1; continue
        if ch=='#':
            j=text.find('\n',i); i=len(text) if j<0 else j; continue
        cleaned.append(ch); i+=1
    pairs={')':'(',']':'[','}':'{'}; stack=[]
    for ch in cleaned:
        if ch in '([{': stack.append(ch)
        elif ch in pairs:
            assert stack and stack[-1]==pairs[ch], f'bracket mismatch {path}: {ch}'
            stack.pop()
    assert not stack, f'unclosed bracket {path}: {stack[-5:]}'

def main():
    lines=[]
    module=(ROOT/'src/MizuRoute.jl').read_text()
    for f in sorted((ROOT/'src').rglob('*.jl')): bracket_static_check(f)
    for f in sorted((ROOT/'test').glob('*.jl')): bracket_static_check(f)
    lines.append('PASS: Julia source/test bracket static checks')
    for rel in ['routing/advection_diffusion.jl','routing/irf.jl','routing/kwt.jl','routing/euler_kw.jl','routing/muskingum_cunge.jl','routing/diffusive_wave.jl']:
        assert f'include("{rel}")' in module
    lines.append('PASS: all principal routing kernels included by module')
    for y in [0.05,0.5,2.0,5.0]:
        a=flow_area(y,8.,1.2); y2=water_height(a,8.,1.2); assert abs(y-y2)<1e-12
    q=35.; y=flow_depth(q,12.,1.,0.0015,0.035); assert abs(manning(y,12.,1.,0.0015,0.035)-q)<1e-8
    lines.append('PASS: hydraulic area/depth and Manning depth/discharge round trips')
    qc=np.full(8,3.5); out=solve_ade(qc,5000,300,3.5,1.2,0); assert np.max(np.abs(out-qc))<1e-12
    q0=np.zeros(8); out=solve_ade(q0,5000,300,0,1.2,200); assert np.max(abs(out))<1e-12
    out=solve_ade(q0,5000,300,2.,1.2,200); assert np.all(np.isfinite(out)) and np.all(out>=0) and abs(out[0]-2)<1e-12
    out=solve_ade(q0,5000,300,0.,1.2,200,1.0); assert out.max()>0
    lines.append('PASS: shared ADE constant/zero/boundary/lateral-source invariants')
    X=0.2; Cn=0.55; den=1-X+0.5*Cn
    c0=(-X+0.5*Cn)/den; c1=(X+0.5*Cn)/den; c2=(1-X-0.5*Cn)/den
    assert abs(c0+c1+c2-1)<1e-15
    lines.append('PASS: Muskingum-Cunge coefficient steady-flow conservation')
    down=[2,2,-1]; indeg=[0,0,2]; order=[]; queue=[0,1]
    while queue:
        i=queue.pop(0); order.append(i); d=down[i]
        if d>=0:
            indeg[d]-=1
            if indeg[d]==0: queue.append(d)
    assert order[-1]==2 and len(order)==3
    lines.append('PASS: river-network upstream-to-downstream topological traversal')
    main_typ=(ROOT/'docs/main.typ').read_text()
    for name in ['overview.typ','algorithms.typ','api.typ','validation.typ']:
        assert name in main_typ and (ROOT/'docs'/name).exists()
    lines.append('PASS: Typst document include graph / file presence')
    report='\n'.join(lines)+'\n\nNOTE: These are static and Python numerical mirror checks, not Julia Pkg.test or Typst compilation.\n'
    (ROOT/'validation/python_sanity_report.txt').write_text(report)
    print(report)

if __name__=='__main__': main()
