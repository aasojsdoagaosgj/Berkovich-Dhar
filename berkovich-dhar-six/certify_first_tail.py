"""Exact polynomial identities and truncated checks of the positive kernels."""
from pathlib import Path
import hashlib,json
BASE=Path(__file__).resolve().parent

def add(*ps):
    out={}
    for p in ps:
        for k,v in p.items():out[k]=out.get(k,0)+v
    return {k:v for k,v in out.items() if v}
def scale(p,c):return {k:c*v for k,v in p.items() if c*v}
def mul(p,q):
    out={}
    for a,x in p.items():
        for b,y in q.items():
            k=tuple(u+v for u,v in zip(a,b));out[k]=out.get(k,0)+x*y
    return {k:v for k,v in out.items() if v}
def power(p,n):
    out={(0,0,0):1}
    for _ in range(n):out=mul(out,p)
    return out

def run():
    g={(1,0,0):1};h={(0,1,0):1};q={(0,0,1):1};one={(0,0,0):1}
    gh5=mul(g,power(h,5));h6=power(h,6);g3h3=mul(power(g,3),power(h,3))
    A=add(power(g,6),scale(mul(q,gh5),-3))
    B=add(scale(mul(power(g,5),h),-3),scale(mul(q,h6),-1));C=scale(g3h3,5)
    delta=add(power(g,2),scale(mul(g,h),-1),scale(power(h,2),-1))
    P=add(scale(g3h3,5),scale(gh5,-3),scale(h6,-1));u=add(g,scale(h,-1))
    assert add(A,B,C)==add(power(delta,3),mul(add(one,scale(q,-1)),add(scale(gh5,3),h6)))
    assert add(A,B,mul(q,C))==add(power(delta,3),scale(mul(add(one,scale(q,-1)),P),-1))
    assert P==add(h6,scale(mul(u,power(h,5)),12),scale(mul(power(u,2),power(h,4)),15),scale(mul(power(u,3),power(h,3)),5))
    M=80
    def conv(a,b):return [sum(a[j]*b[k-j] for j in range(k+1)) for k in range(M+1)]
    def RR(t):
        out=[0]*(M+1);den=[1]+[0]*M
        for k in range(0,10):
            if k:
                for j in range(k,M+1):den[j]+=den[j-k]
            shift=k*k+(t-1)*k
            if shift>M:break
            for j in range(M+1-shift):out[j+shift]+=den[j]
        return out
    gv,hv=RR(1),RR(2);hh=conv(hv,hv);gh=conv(gv,hv);gg=conv(gv,gv)
    D=[];v=0
    for k in range(M+1):v+=gh[k]+hh[k]-gg[k];D.append(v)
    Dsum=[0]*(M+1)
    for j in range(7):
        shift=2*j*(j+1)
        if shift>M:break
        val=conv(RR(2*j+2),RR(2*j+3))
        for k in range(M+1-shift):Dsum[k+shift]+=val[k]
    assert D==Dsum
    assert all(0<=v<=hh[k] for k,v in enumerate(D))
    result={'status':'passed','polynomial_identities':['F4 numerator','F5 numerator','P in u=g-h'],
        'series_checks_through':M,'positive_kernel_identity':True,'kernel_bounds':True,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'Polynomial identities are exact; infinite coefficientwise inequalities follow from the written positive expansions.'}
    (BASE/'results/first_tail_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':run()
