"""Exact finite-polynomial checks; these do not prove an eventual theorem."""
from pathlib import Path
import json,hashlib,time
BASE=Path(__file__).resolve().parent

def run(maxn=40):
    start=time.perf_counter();rows=[]
    for p,s in [(3,s) for s in range(5,9)]+[(5,2),(5,3)]:
        cap=p*(p-1)*s*maxn*maxn//4
        a=[1]+[0]*cap;degree=0;records=[]
        for n in range(1,maxn+1):
            for part in range(p*(n-1)+1,p*n):
                for _ in range(s):
                    degree=min(cap,degree+part)
                    for k in range(degree,part-1,-1):a[k]-=a[k-part]
            half=p*(p-1)*s*n*n//4
            assert all(a[j]>=0 for j in range(0,half+1,p)),(p,s,n,'residue0')
            for res in ([1] if p==3 else [1,2]):
                assert all(a[j]<=0 for j in range(res,half+1,p)),(p,s,n,res)
            crossings={}
            for res in ([2] if p==3 else [3,4]):
                negative=False;lastpositive=None;firstnegative=None
                for j in range(res,half+1,p):
                    if a[j]>0:
                        assert not negative,(p,s,n,res,j,'reversal')
                        lastpositive=(j-res)//p
                    elif a[j]<0:
                        if not negative:firstnegative=(j-res)//p
                        negative=True
                crossings[res]={'last_positive':lastpositive,'first_negative':firstnegative}
            records.append({'n':n,'lower_half_max_degree':half,'transitions':crossings})
        row={'p':p,'power':s,'all_n_from':1,'all_n_through':maxn,'records':records}
        rows.append(row)
        print(json.dumps({'p':p,'s':s,'n':maxn,'last':records[-1]}),flush=True)
    result={'status':'passed','scope':'Exact integer diagnostics for every lower-half coefficient; reciprocity gives upper half.',
        'cases':rows,'seconds':time.perf_counter()-start,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (BASE/'results/finite_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':run()
