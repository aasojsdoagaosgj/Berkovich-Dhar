"""Exact finite coefficients and rational margins for the stated Wang bound."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json
BASE=Path(__file__).resolve().parent

def exp_lower(x,n):
    t=s=Q(1)
    for k in range(1,n+1):
        t*=x/k;s+=t
    return s

def coefficients(p,s,M):
    a=[1]+[0]*M
    for j in range(1,M+1):
        if j%p:
            for _ in range(s):
                for k in range(M,j-1,-1):a[k]-=a[k-j]
    return a

def run():
    rows=[]
    for p,s in [(3,s) for s in range(5,9)]+[(5,2),(5,3)]:
        a=coefficients(p,s,2000)
        signs=[1,-1,1] if p==3 else ([1,-1,-1,1,1] if s==2 else [1,-1,0,1,0])
        assert all(v*signs[k%p]>=0 and (signs[k%p]!=0 or v==0) for k,v in enumerate(a))
        zeros=[k for k,v in enumerate(a) if v==0 and signs[k%p]!=0]
        assert zeros==([9] if (p,s)==(5,2) else [])
        rows.append({'p':p,'power':s,'expected_signs':signs,'exceptional_zeros':zeros,
            'max_degree':2000,'coefficient_sha256':hashlib.sha256(str(a).encode()).hexdigest()})
    e31=exp_lower(Q(31),100)
    checks={
        'e_pi_over_five':exp_lower(Q(31,50),8)>Q(50,27),
        'e_pi':exp_lower(Q(31,10),12)>20,
        'first_error_constant':Q(22,7)**7*Q(2,3)<Q(5)**4*8,
        'exponential_sum':Q(20,361)+Q(1350,529)<Q(8,3),
        'second_error_exponent':2+8*Q(22,7)*Q(2,3)<19,
        'third_error_exponent':Q(22,7)*Q(2,3)+8*Q(8,3)<25,
        'constant_error':10*3**19+2*9*5**5*3**25<10**18,
        'bessel_main_prefactor': 4**4*Q(3)**2*Q(1,3)>5**2,
        'phase_margin_mod3':2*(Q(31,180)-Q(31,180)**3/6)>Q(1,3),
        'phase_margin_mod5':Q(7,3)**2>5,
        'X_lower':16*Q(31,10)**2>Q(7)**2*3,
        'X_at_2000':Q(7,5)**2*2000>62**2,
        'power_at_2000':2000**3<300**4,
        'exp31':e31>2*10**13,
        'error_ratio':Q(450000,e31)+Q(9*10**22,e31*e31)<Q(1,100),
    }
    assert all(checks.values()),checks
    result={'status':'passed','cases':rows,'rational_checks':checks,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'Finite signs and elementary bounds; analytic dependencies are in INFINITE_SIGNS.md.'}
    (BASE/'results').mkdir(exist_ok=True)
    (BASE/'results/infinite_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':run()
