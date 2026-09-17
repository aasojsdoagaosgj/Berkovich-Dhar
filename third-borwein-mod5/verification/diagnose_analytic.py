"""Ordinary multiprecision diagnostics. These are NOT interval certificates."""
import argparse
from math import isqrt
from common import write_result, metadata
from mpmath import mp
from verify_exact import rr, convolution


def R(z):
    if z==0:return mp.log(5)
    return (2*mp.pi**2/15+mp.polylog(2,mp.exp(-5*z))/5-mp.polylog(2,mp.exp(-z)))/z


def Rp(t):
    if t==0:return -mp.mpf(1)
    lam=mp.polylog(2,mp.exp(-5*t))/5-mp.polylog(2,mp.exp(-t))
    lamp=mp.log(1-mp.exp(-5*t))-mp.log(1-mp.exp(-t))
    return lamp/t-(2*mp.pi**2/15+lam)/t**2


def Rpp(t):
    return mp.diff(Rp,t)


def saddle(n,k):
    target=mp.mpf(k)/(5*n*n)
    if not 0<target<1:raise ValueError('Saddle requires 0<k<5n^2')
    lo=mp.mpf(0); hi=mp.sqrt((2*mp.pi**2/15)/target)*2
    for _ in range(300):
        mid=(lo+hi)/2
        if -Rp(mid)>target:lo=mid
        else:hi=mid
    return (lo+hi)/2


def lam_root(ell,z):
    zeta=mp.exp(2j*mp.pi/5)
    return sum((mp.mpf(j)/5-mp.mpf('0.5'))*(mp.log(1-zeta**(ell*j)*mp.exp(-z))-mp.log(1-zeta**(ell*j))) for j in range(1,5))


def partitions(limit):
    terms=[]
    for j in range(1,limit+1):
        a=j*(3*j-1)//2;b=j*(3*j+1)//2
        if a>limit:break
        sign=1 if j%2 else -1
        terms.append((a,sign))
        if b<=limit:terms.append((b,sign))
    p=[1]+[0]*limit
    for n in range(1,limit+1):
        value=0
        for g,sign in terms:
            if g>n:break
            value+=sign*p[n-g]
        p[n]=value
        if n%25000==0:print(f'partition recurrence: {n}/{limit}',flush=True)
    return p


def d_from_partitions(k,cumulative):
    value=0
    bound=isqrt(10*k)+3
    for b in range(-bound,bound+1):
        index=k-(b*(3*b-1)//10)
        if index>=0:value-=(-1 if b%2 else 1)*cumulative[index]
    return value


def diagnose(large=False,output_dir=None):
    mp.dps=80
    text=lambda x:mp.nstr(x,45)
    zeta=mp.exp(2j*mp.pi/5)
    rows=[]
    for n in (8,32):
        for tau in ('0','0.5','2','8'):
            for t in ('0','0.1','1.2'):
                z=mp.mpf(tau)-1j*mp.mpf(t)
                for ell in (1,2):
                    exact=sum(mp.log(1-zeta**(ell*k)*mp.exp(-z*k/(5*n))) for k in range(1,5*n+1) if k%5)
                    error=exact-n*R(z)-lam_root(ell,z)
                    bound=mp.mpf(500 if mp.mpf(t)<=mp.mpf('0.1') else 150000)/n
                    assert abs(error)<bound
                    row={'n':n,'tau':tau,'t':t,'ell':ell,'absolute_log_error':text(abs(error)),'bound':text(bound)}
                    if mp.mpf(t)<=mp.mpf('0.1'):
                        real_tau=mp.mpf(tau)
                        integral=mp.mpf(1) if real_tau==0 else -mp.expm1(-real_tau)/real_tau
                        improved=(real_tau**2+mp.mpf('0.01'))*integral/(2*mp.mpf('0.9')**2*n)
                        assert abs(error)<improved
                        row['v04_integral_EM_bound']=text(improved)
                    rows.append(row)
    # Check eta multipliers and transformed products away from prohibitively tiny remainders.
    eta_rows=[]
    A=2*mp.pi**2/75
    for w0 in ('0.03','0.05'):
        for ratio in ('0','0.75'):
            w=mp.mpf(w0)*(1+1j*mp.mpf(ratio))
            for ell in range(1,5):
                q=zeta**ell*mp.exp(-w)
                exact=mp.qp(q)/mp.qp(q**5)
                mult=mp.exp(-1j*mp.pi/5) if ell==1 else mp.exp(1j*mp.pi/5) if ell==4 else mp.mpf(1)
                transformed=zeta**(-pow(ell,-1,5))*mp.exp(-4*mp.pi**2/(25*w))
                denominator=mp.exp(-4*mp.pi**2/(5*w))
                predicted=mult*mp.exp(A/w-w/6)*mp.qp(transformed)/mp.qp(denominator)
                err=abs(exact/predicted-1)
                assert err<mp.mpf('1e-65')
                eta_rows.append({'w0':w0,'imaginary_ratio':ratio,'ell':ell,'relative_error':text(err)})
    max_k=100000 if large else 150
    p=partitions(max_k)
    cumulative=[]; acc=0
    for value in p:acc+=value;cumulative.append(acc)
    g,h=rr(1,150),rr(2,150)
    gh,hh,gg=convolution(g,h,150),convolution(h,h,150),convolution(g,g,150)
    acc=0
    for k in range(151):
        acc+=gh[k]+hh[k]-gg[k]
        assert acc==d_from_partitions(k,cumulative)
    endpoint=[]
    if large:
        for n in (31148,40000,120000,160000,200000):
            K=n-1 if n==31148 else 39200 if n==40000 else n//2
            d=d_from_partitions(K,cumulative)
            for a in (3,4):
                m=5*n+5*K+a;k=m-5*n
                tau=saddle(n,k);w=tau/(5*n)
                logP=n*R(tau)+k*w-w/6-mp.log(5)-mp.log(2*mp.pi*Rpp(tau))/2-mp.mpf('1.5')*mp.log(n)
                ratio=mp.exp(mp.log(d)-logP)
                if n in (31148,40000):
                    budget=15*mp.sqrt(w)+15*mp.exp(-tau)+60*w+50*w**(-3)*mp.exp(-1/(30*w))
                    assert tau>=mp.mpf('5.5') and w<=mp.mpf('.0013') and abs(ratio-1)<budget
                else:
                    budget=5*mp.sqrt(w)+12*mp.exp(-tau)+45*w+50*w**(-3)*mp.exp(-1/(30*w))
                    assert tau>=8 and w<=mp.mpf('0.001') and abs(ratio-1)<budget
                endpoint.append({'n':n,'m':m,'a':a,'K':K,'coefficient':str(-d),'tau':text(tau),'w0':text(w),
                                 'coefficient_over_negative_main':text(ratio),'error_budget':text(budget)})
    rho=-5*Rp(mp.mpf(8))
    data={**metadata('ordinary multiprecision diagnostics'),'decimal_precision':80,'status':'passed',
          'scope':'Finite samples; not uniform bounds and not directed interval arithmetic.',
          'rho':text(rho),'bulk_samples':rows,'eta_samples':eta_rows,
          'partition_formula_checked_through_K':150,'large_endpoint_samples':endpoint,
          'historical_v02_gaussian_budget_threshold':text((1800*mp.exp(8))**2)}
    write_result('analytic_diagnostics.json',data,output_dir)
    return data


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--large',action='store_true')
    parser.add_argument('--output-dir')
    args=parser.parse_args()
    diagnose(args.large,args.output_dir)
