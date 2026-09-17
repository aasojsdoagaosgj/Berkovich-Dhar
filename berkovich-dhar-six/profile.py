"""Profiles and finite checks for the six remaining cases (diagnostic only)."""
from pathlib import Path
import sys, json, hashlib
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'.local-deps'))
import mpmath as mp
mp.mp.dps=60
BASE=Path(__file__).resolve().parent

def R(p,z):
    if not z: return mp.log(p)
    x=mp.exp(-z)
    return (mp.pi**2*(p-1)/(6*p)+mp.polylog(2,x**p)/p-mp.polylog(2,x))/z

def angles(z):
    x=mp.exp(-z)
    a=mp.atan(mp.sin(2*mp.pi/5)*x/(1-mp.cos(2*mp.pi/5)*x))
    b=mp.atan(mp.sin(4*mp.pi/5)*x/(1-mp.cos(4*mp.pi/5)*x))
    return -mp.pi/10+(a+2*b)/5,-mp.pi/10+(2*a-b)/5

def psi(p,s,a,z):
    if p==3:
        phi=-s*mp.pi/18+s*mp.atan(mp.sqrt(3)*mp.exp(-z)/(2+mp.exp(-z)))/3
        return 2*mp.cos(phi-2*mp.pi*a/3)
    u,v=angles(z)
    return 4*mp.cos(s*u-3*mp.pi*a/5)*mp.cos(s*v+mp.pi*a/5)

def infinite(p,s,M):
    a=[1]+[0]*M
    for j in range(1,M+1):
        if j%p:
            for _ in range(s):
                for k in range(M,j-1,-1):a[k]-=a[k-j]
    return a

def run():
    records=[]
    for p,s in [(3,s) for s in range(5,9)]+[(5,2),(5,3)]:
        row={'p':p,'power':s,'roots':{}}
        for a in ([2] if p==3 else [3,4]):
            if p==3:
                t=mp.tan(mp.pi/6-mp.pi/(2*s))
                z=-mp.log(2*t/(mp.sqrt(3)-t))
            else:
                z=mp.findroot(lambda z:angles(z)[1 if a==3 else 0]+mp.pi/(10*s),(.1,1))
            row['roots'][a]={'tau':str(z),'limit':str(-s*mp.diff(lambda z:R(p,z),z)),
                'psi_derivative':str(mp.diff(lambda z:psi(p,s,a,z),z))}
        vals=infinite(p,s,2000)
        signs=([1,-1,1] if p==3 else ([1,-1,-1,1,1] if s==2 else [1,-1,0,1,0]))
        wrong=[j for j,v in enumerate(vals) if (v*signs[j%p]<0 or (signs[j%p]==0 and v))]
        zeros=[j for j,v in enumerate(vals) if v==0 and signs[j%p]]
        row.update(infinite_max_degree=2000,wrong_signs=wrong,exceptional_zeros=zeros,
            coefficient_sha256=hashlib.sha256(str(vals).encode()).hexdigest())
        if (p,s)==(5,3):
            cumul=0; violations=[]; cz=[]
            for j,v in enumerate(vals):
                cumul+=v
                if j%5 in (2,4):
                    if cumul*(1 if j%5==4 else -1)<0: violations.append(j)
                    if not cumul:cz.append(j)
            row['weak_first_tail_partial_sums']={'violations':violations,'zeros':cz}
        records.append(row)
        print(json.dumps(row),flush=True)
    result={'status':'diagnostic_only','cases':records,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (BASE/'results').mkdir(exist_ok=True)
    (BASE/'results/profiles.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':run()
