"""Finite direct checks of periodic residue bounds; no sampling-based certification."""
import argparse
from math import gcd
from common import metadata,write_result
from mpmath import mp
from parameters import (DIRICHLET_Q,FOURIER_K,ETA_NUMERATOR,ETA_DENOMINATOR,
 NON5_FOURIER_K,NON5_ETA_NUMERATOR,NON5_ETA_DENOMINATOR,THRESHOLD)
from certify_resonant_gap import verify_weight_polynomials


def diagnose(output_dir=None):
    mp.dps=80;Q=mp.mpf(DIRICHLET_Q);rows=[];geometric=[]
    sums={}
    for five in (False,True):
        K=FOURIER_K if five else NON5_FOURIER_K
        eta=mp.mpf(ETA_NUMERATOR)/ETA_DENOMINATOR if five else mp.mpf(NON5_ETA_NUMERATOR)/NON5_ETA_DENOMINATOR
        sums[five]=(K,eta,[mp.exp(-eta*l) for l in range(1,K+1)])
    for b in (1,5,10,149,151,499,505,899,901,1799,1801,5399,5400):
        K,eta,weights=sums[b%5==0]
        S0=sum(weights[l-1]/(l*((l+1)//2)) for l in range(1,K+1))
        S1=sum(weights[l-1]/l**2 for l in range(1,K+1))
        for a in sorted({1,b-1,(b+1)//5,(b-1)//5}):
            if not 0<a<b or gcd(a,b)!=1:continue
            for d in (1,5):
                g=gcd(b,d);B=b//g
                direct=mp.mpf(0)
                for l in range(1,K+1):
                    r=(d*a*l)%b;r=min(r,b-r)
                    if r==0:continue
                    direct+=weights[l-1]/(l*(r-mp.mpf(d)*l/Q))
                q=mp.exp(-eta*B)
                T=S0+2*(1+mp.log(K))*q/(B*(1-q))
                U=mp.pi**2/(3*(1-q))
                correction=mp.mpf(d)/(g*Q*(1-mp.mpf(d)*K/(g*Q)))
                bound=(T+correction*U)/g
                assert direct<bound
                if b%5 and b>2*K:
                    bound=(S1+correction*mp.pi**2/6)/g
                    assert direct<bound
                rows.append({'b':b,'a':a,'d':d,'direct_sum':mp.nstr(direct,45),
                             'bound':mp.nstr(bound,45),'ratio':mp.nstr(direct/bound,45)})
        for mult in (1,5):
            count=mult*THRESHOLD
            z=mp.mpc(mp.mpf('5.5')*K/count,10*mp.pi*K/(b*DIRICHLET_Q))
            if abs(z)>=2*mp.pi:continue
            kernel=1/mp.expm1(z)-1/z
            kappa=mp.mpf('.5')+abs(z)/(12*(1-abs(z)**2/(4*mp.pi**2)))
            assert abs(kernel)<kappa
            geometric.append({'b':b,'M':count,'kernel_modulus':mp.nstr(abs(kernel),45),'bound':mp.nstr(kappa,45)})
    result={**metadata('ordinary multiprecision v0.6 localization diagnostics'),'status':'passed',
            'decimal_precision':80,'residue_samples':rows,'geometric_samples':geometric,
            'exact_weight_polynomials':verify_weight_polynomials(),
            'scope':'Finite direct sums and exact polynomial identities; analytic uniformity is in manuscript.'}
    write_result('localization_diagnostics.json',result,output_dir)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-dir');diagnose(p.parse_args().output_dir)
