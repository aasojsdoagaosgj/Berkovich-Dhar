"""Independent product/EM/transition numerical diagnostics, not proof bounds."""
from pathlib import Path
import sys,json,hashlib,time
from math import factorial,prod
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'.local-deps'))
import mpmath as mp
from profile import R,psi,angles
BASE=Path(__file__).resolve().parent
mp.mp.dps=65

def em(p,s,h,z,order):
    root=mp.exp(2j*mp.pi*h/p)
    total=0
    for b in range(1,p):
        t=mp.mpf(b)/p;r=root**b;x=mp.exp(-z)
        if order==0:
            total+=(t-mp.mpf('.5'))*(mp.log(1-r*x)-mp.log(1-r))
        elif order==1:
            total+=(t*t-t+mp.mpf(1)/6)*z*(r*x/(1-r*x)-r/(1-r))/2
        else:
            total+=(t**3-mp.mpf('1.5')*t*t+t/2)*z*z*(-r*x/(1-r*x)**2+r/(1-r)**2)/6
    return s*total

def amplitudes(p,s,a,z,order):
    total=0
    for h in range(1,p):
        value=mp.exp(-2j*mp.pi*a*h/p+em(p,s,h,z,0))
        if order==1:value*=em(p,s,h,z,1)
        if order==2:value*=em(p,s,h,z,2)+em(p,s,h,z,1)**2/2
        total+=value
    return mp.re(total)

def S1(H,f,z):
    g=mp.diff(H,z,2);h3=mp.diff(H,z,3);h4=mp.diff(H,z,4)
    return -mp.diff(f,z,2)/(2*g)+h3*mp.diff(f,z)/(2*g*g)+f(z)*(h4/(8*g*g)-5*h3*h3/(24*g**3))

def S2(H,f,z):
    g=mp.diff(H,z,2);hd={r:mp.diff(H,z,r) for r in range(3,7)}
    fd={r:mp.diff(f,z,r) for r in range(5)};value=0
    for p in range(5):
      for e3 in range(5):
       for e4 in range(3):
        for e5 in range(2):
         for e6 in range(2):
          es={3:e3,4:e4,5:e5,6:e6}
          if p+sum((r-2)*e for r,e in es.items())!=4:continue
          D=p+sum(r*e for r,e in es.items())
          term=(-1)**(D//2)*prod(range(1,D,2))*fd[p]/(factorial(p)*g**(D//2))
          for r,e in es.items():term*=hd[r]**e/(factorial(r)**e*factorial(e))
          value+=term
    return value

def run():
    start=time.perf_counter();rows=[]
    for p,s in [(3,s) for s in range(5,9)]+[(5,2),(5,3)]:
        H=lambda z:s*R(p,z)
        roots={}
        for a in ([2] if p==3 else [3,4]):
            if p==3:
                t=mp.tan(mp.pi/6-mp.pi/(2*s));z=-mp.log(2*t/(mp.sqrt(3)-t))
            else:z=mp.findroot(lambda z:angles(z)[1 if a==3 else 0]+mp.pi/(10*s),(.1,1))
            f=lambda t:psi(p,s,a,t)
            a1=amplitudes(p,s,a,z,1)+S1(H,f,z)
            beta=mp.diff(H,z,2)*a1/mp.diff(f,z)
            roots[a]=(z,-mp.diff(H,z),beta)
        row={'p':p,'s':s,'roots':{a:{'tau':str(v[0]),'alpha':str(v[1]),'beta':str(v[2])} for a,v in roots.items()},'samples':[]}
        for n in (20,40,80):
            targets={a:[int(mp.floor(al*n*n+be*n-mp.mpf(a)/p))+d for d in (-1,0,1,2)] for a,(t,al,be) in roots.items()}
            cap=max(p*max(ks)+a for a,ks in targets.items());coeff=[1]+[0]*cap
            degree=0
            for part in range(1,p*n):
                if part%p:
                    for _ in range(s):
                        degree=min(cap,degree+part)
                        for k in range(degree,part-1,-1):coeff[k]-=coeff[k-part]
            for a,ks in targets.items():
                f=lambda z:psi(p,s,a,z)
                B=lambda z:amplitudes(p,s,a,z,1)
                samples=[]
                for k in ks:
                    j=p*k+a;al=mp.mpf(j)/(p*n*n)
                    t=mp.findroot(lambda z:-mp.diff(H,z)-al,roots[a][0]);g=mp.diff(H,t,2)
                    norm=mp.exp(n*(H(t)+al*t))/(p*n*mp.sqrt(2*mp.pi*n*g))
                    A1=B(t)+S1(H,f,t)
                    A2=amplitudes(p,s,a,t,2)+S1(H,B,t)+S2(H,f,t)
                    residual=mp.mpf(coeff[j])/norm-f(t)-A1/n-A2/(n*n)
                    assert abs(amplitudes(p,s,a,t,0)-f(t))<mp.mpf('1e-55')
                    samples.append({'k':k,'sign':1 if coeff[j]>0 else -1 if coeff[j]<0 else 0,'n3_residual':mp.nstr(residual*n**3,25)})
                signs=[v['sign'] for v in samples]
                assert 1 in signs and -1 in signs and signs==sorted(signs,reverse=True),(p,s,n,a,signs)
                row['samples'].append({'n':n,'residue':a,'coefficients':samples})
        rows.append(row);print(json.dumps(row),flush=True)
    result={'status':'diagnostic_completed','precision':mp.mp.dps,'cases':rows,'seconds':time.perf_counter()-start,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'dependency_sha256':{'profile.py':hashlib.sha256((BASE/'profile.py').read_bytes()).hexdigest()},
        'scope':'Four coefficients near each transition, n=20,40,80; these are not all-coefficient or interval tests.'}
    (BASE/'results/saddle_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':run()
