"""Exact scalar margins and numerical checks of the modulus-5 endpoint formula."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'.local-deps'))
import mpmath as mp
from profile import R,psi
from saddle_diagnostic import em
from certify_infinite import exp_lower
BASE=Path(__file__).resolve().parent

def run():
    checks={
        'exp5':exp_lower(Q(5),20)>140,
        'C_lower':2*Q(157,50)**2/15>Q(13,10),
        'C_upper':2*Q(22,7)**2/15<Q(4,3),
        'K_second':Q(4*37,139)<Q(11,10),
        'g_lower':2*Q(13,10)-Q(11,10)>=Q(3,2),
        'g_upper':2*Q(4,3)+Q(11,10)<4,
        'arc_decay':Q(16,25)*Q(13,10)-Q(11,20)>Q(1,4),
        'outside_gap':Q(1,25)-Q(6,5*139)>Q(1,40),
        'dual_gap':Q(18,5*139)<Q(1,20),
        'dual_exponent':4*Q(157,50)**2*16/625>1,
    }
    assert all(checks.values()),checks
    mp.mp.dps=80;rows=[]
    for n,tau in [(100,5),(200,5),(100,12),(200,12),(200,30)]:
        z=mp.mpf(tau);w=z/(5*n);x=mp.exp(-z)
        for s in (2,3):
            values={};errors=[]
            for h in (1,2):
                root=mp.exp(2j*mp.pi*h/5)
                logprod=s*sum(mp.log(1-root**(l%5)*mp.exp(-l*w)) for l in range(1,5*n) if l%5)
                values[h]=mp.exp(logprod-n*s*R(5,z))
                values[5-h]=mp.conj(values[h])
                ratio=mp.exp(logprod-n*s*R(5,z)-em(5,s,h,z,0)+s*w/6)
                errors.append(mp.nstr(abs(ratio-1)/(w*x),25))
            amplitudes=[]
            for a in range(5):
                val=mp.re(sum(mp.exp(-2j*mp.pi*a*h/5)*v for h,v in values.items()))
                target=psi(5,s,a,z);weak=s==3 and a in (2,4)
                if weak:val/=x;target/=x
                assert val*target>0,(n,tau,s,a)
                amplitudes.append({'residue':a,'weak_scaled':weak,'value':mp.nstr(val,25),'limit_profile':mp.nstr(target,25),'difference_over_w':mp.nstr((val-target)/w,25)})
            rows.append({'n':n,'tau':tau,'s':s,'root_error_over_wx':errors,'amplitudes':amplitudes})
    deps=['profile.py','saddle_diagnostic.py','certify_infinite.py']
    result={'status':'diagnostic_completed','rational_checks':checks,'samples':rows,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'dependency_sha256':{p:hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in deps},
        'scope':'Exact scalar comparisons and high precision samples; uniform analytical bounds are proved in MOD5_ENDPOINT.md.'}
    (BASE/'results/endpoint_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'rational_checks':checks,'max_root_ratio':max(float(v) for row in rows for v in row['root_error_over_wx']),'samples':len(rows)},indent=2))

if __name__=='__main__':run()
