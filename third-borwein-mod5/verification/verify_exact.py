"""Exact integer series and finite products. No floating-point sign decisions."""
import argparse
import hashlib
import numpy as np
from common import write_result, metadata


def convolution(a,b,limit):
    out=[0]*(limit+1)
    for i,x in enumerate(a[:limit+1]):
        if x:
            for j,y in enumerate(b[:limit-i+1]):
                out[i+j]+=x*y
    return out


def rr(t,limit):
    denominator=[1]+[0]*limit
    out=[0]*(limit+1)
    for j in range(limit+1):
        if j:
            for k in range(j,limit+1): denominator[k]+=denominator[k-j]
        shift=j*j+(t-1)*j
        if shift>limit: break
        for k in range(limit-shift+1): out[k+shift]+=denominator[k]
    return out


def reciprocal_product(residues,limit):
    a=[1]+[0]*limit
    for k in range(1,limit+1):
        if k%5 in residues:
            for m in range(k,limit+1): a[m]+=a[m-k]
    return a


def infinite_coefficients(limit):
    # Euler numerator and partition denominator, independent of RR sums/products.
    parts=[1]+[0]*(limit//5)
    for k in range(1,len(parts)):
        for m in range(k,len(parts)): parts[m]+=parts[m-k]
    out=[0]*(limit+1)
    for b in range(-limit,limit+1):
        g=b*(3*b-1)//2
        if g<=limit:
            sign=-1 if b%2 else 1
            for k in range((limit-g)//5+1): out[g+5*k]+=sign*parts[k]
    return out


def naive_product(n):
    a=[1]
    for k in range(1,5*n+1):
        if k%5:
            old=a
            a=old+[0]*k
            for j,c in enumerate(old): a[j+k]-=c
    return a


def verify(max_n=256,degree=350,output_dir=None):
    if not (1<=max_n<=256): raise ValueError('Resource limit: 1 <= max_n <= 256')
    if not (max_n<=degree<=1000): raise ValueError('Require max_n <= degree <= 1000')
    g,h=rr(1,degree),rr(2,degree)
    assert g==reciprocal_product({1,4},degree)
    assert h==reciprocal_product({2,3},degree)
    checked_t=[1,2,3,4,9]
    for t in checked_t:
        r0,r1,r2=rr(t,degree),rr(t+1,degree),rr(t+2,degree)
        assert all(r0[k]==r1[k]+(r2[k-t] if k>=t else 0) for k in range(degree+1))
        aa,bb=convolution(r1,r1,degree),convolution(r0,r2,degree)
        lhs=[]; running=0
        for k in range(degree+1):
            running+=aa[k]-(bb[k-1] if k else 0)
            lhs.append(running)
        rhs=[0]*(degree+1)
        for j in range(degree+1):
            shift=2*j*(t+j)
            if shift>degree: break
            prod=convolution(rr(t+2*j+1,degree-shift),rr(t+2*j+2,degree-shift),degree-shift)
            for k,v in enumerate(prod):rhs[k+shift]+=v
        assert lhs==rhs
        assert all(v>=0 for v in rhs)
    gh,hh,gg=convolution(g,h,degree),convolution(h,h,degree),convolution(g,g,degree)
    d=[]; running=0
    for k in range(degree+1):
        running+=gh[k]+hh[k]-gg[k]; d.append(running)
    assert d[0]==1 and d[1]==0 and all(v>0 for v in d[2:])
    inf=infinite_coefficients(5*degree+4)
    for k in range(degree+1):
        assert inf[5*k:5*k+5]==[gg[k],-gh[k],-hh[k],0,0]
    rows=[]; count=0
    coefficients=np.array([1],dtype=object)
    for n in range(1,max_n+1):
        for k in range(5*n-4,5*n):
            old=coefficients
            coefficients=np.zeros(len(old)+k,dtype=object)
            coefficients[:len(old)]=old
            coefficients[k:]-=old
        assert len(coefficients)==10*n*n+1
        assert np.array_equal(coefficients,coefficients[::-1])
        for residue in range(5):
            block=coefficients[residue::5]
            assert bool(np.all(block>=0) if residue==0 else np.all(block<=0))
            assert sum(block)==(4 if residue==0 else -1)*5**(n-1)
        for residue in (3,4):
            for k in range(n):assert coefficients[5*n+5*k+residue]==-d[k]
        zeros=set(map(int,np.flatnonzero(coefficients==0)))
        if n>=2:
            left={5*j+a for j in range(n) for a in (3,4)}|{7,5*n+8,5*n+9}
            assert zeros==left|{10*n*n-m for m in left}
            assert len(zeros)==4*n+6
        if n<=10:assert coefficients.tolist()==naive_product(n)
        modular=[]
        if n in {2,8,32,64,128,256}:
            for prime in (1000000007,1000000009):
                value=0
                for c in coefficients[::-1]:value=(2*value+int(c))%prime
                direct=1
                for k in range(1,5*n+1):
                    if k%5:direct=direct*(1-pow(2,k,prime))%prime
                assert value==direct
                modular.append({'prime':prime,'evaluation_at_2':value})
        digest=hashlib.sha256()
        for c in coefficients:digest.update((str(int(c))+'\n').encode('ascii'))
        rows.append({'n':n,'degree':len(coefficients)-1,'zeros':len(zeros),
                     'coefficients_sha256_decimal_lines':digest.hexdigest(),'modular_cross_checks':modular})
        count+=len(coefficients)
        if n%32==0:print(f'exact finite products: n={n}/{max_n}',flush=True)
    data={**metadata('exact integer verification'),'numpy':np.__version__,'status':'passed',
          'series_degree':degree,'RR_identity_t_values':checked_t,'five_dissection_q_degree':5*degree+4,
          'max_n':max_n,'coefficients_checked':count,'d_first':d[:18],
          'scope':'All coefficients for every 1<=n<=max_n; no extrapolation.',
          'finite_rows':rows}
    write_result('exact_certificate.json',data,output_dir)
    return data


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--max-n',type=int,default=256)
    parser.add_argument('--degree',type=int,default=350)
    parser.add_argument('--output-dir')
    args=parser.parse_args()
    verify(args.max_n,args.degree,args.output_dir)
