"""Check saved computation provenance and freeze the written dependencies."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json,subprocess,sys,os
BASE=Path(__file__).resolve().parent
ROOT=BASE.parent
JOBS=[('profile.py','profiles.json'),('certify_infinite.py','infinite_certificate.json'),
      ('certify_first_tail.py','first_tail_certificate.json'),
      ('finite_diagnostic.py','finite_diagnostic.json'),
      ('saddle_diagnostic.py','saddle_diagnostic.json'),
      ('endpoint_diagnostic.py','endpoint_diagnostic.json')]
DEPS=['quartic/COMPACT_PROOF.md','quartic/ENDPOINT_PROOF.md','quartic/INFINITE_SIGNS.md',
      'quartic/certify_infinite.py','quartic/results/infinite_certificate.json',
      'quartic/certify_endpoint_constants.py','quartic/results/endpoint_certificate.json',
      'third-borwein-mod5/manuscript.md']
DEPS+=['third-borwein-mod5/verification/'+p for p in ['common.py','parameters.py','certify_localization.py',
      'certify_resonant_gap.py','certify_scalars.py','certify_profiles.py']]

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def execute(path,args=()):
    env=os.environ.copy();env['PYTHONOPTIMIZE']='0'
    r=subprocess.run([sys.executable,'-B',str(path),*args],cwd=ROOT,capture_output=True,text=True,encoding='utf-8',env=env)
    if r.returncode:raise RuntimeError(r.stdout+r.stderr)

def run(replay=False):
    if not __debug__:raise RuntimeError('Do not use -O.')
    if replay:
        execute(ROOT/'quartic/certify_infinite.py')
        execute(ROOT/'quartic/certify_endpoint_constants.py')
    records=[]
    for script,report in JOBS:
        if replay:execute(BASE/script)
        path=BASE/'results'/report;data=json.loads(path.read_text(encoding='utf-8'))
        assert data['source_sha256']==sha(BASE/script),(script,'stale')
        for dep,digest in data.get('dependency_sha256',{}).items():assert sha(BASE/dep)==digest,(script,dep)
        assert data['status'] in ('passed','diagnostic_completed','diagnostic_only')
        records.append({'script':script,'report':report,'status':data['status'],'report_sha256':sha(path)})
        print(script+': '+data['status'],flush=True)
    expected={(3,5),(3,6),(3,7),(3,8),(5,2),(5,3)}
    finite=json.loads((BASE/'results/finite_diagnostic.json').read_text())
    assert {(r['p'],r['power']) for r in finite['cases']}==expected
    assert all([v['n'] for v in r['records']]==list(range(1,41)) for r in finite['cases'])
    saddles=json.loads((BASE/'results/saddle_diagnostic.json').read_text())
    assert {(r['p'],r['s']) for r in saddles['cases']}==expected
    values=[float(v['n3_residual']) for r in saddles['cases'] for sample in r['samples'] for v in sample['coefficients']]
    assert len(values)==96 and max(map(abs,values))<.30
    if replay:execute(ROOT/'third-borwein-mod5/verification/certify_scalars.py',['--output-dir',str(BASE/'results/dependencies')])
    deps_reports={}
    for name in ['localization_certificate.json','resonant_gap_certificate.json','profile_certificate.json','scalar_certificate.json']:
        path=BASE/'results/dependencies'/name
        assert json.loads(path.read_text(encoding='utf-8'))['status']=='passed'
        deps_reports['berkovich-dhar-six/results/dependencies/'+name]=sha(path)
    own=[p for p in BASE.iterdir() if p.suffix in ('.md','.py') and p.name!='LARGE_N_REVIEW_TEXT.md']
    inputs={p.relative_to(ROOT).as_posix().replace('\\','/'):sha(p) for p in own}
    inputs.update({p:sha(ROOT/p) for p in DEPS})
    result={'status':'saved_results_and_provenance_checked','replayed':replay,'checked_at_utc':datetime.now(timezone.utc).isoformat(),
        'jobs':records,'dependency_reports':deps_reports,'input_sha256':dict(sorted(inputs.items())),
        'saddle_sample_count':96,'largest_sample_n3_residual':max(map(abs,values)),
        'limitation':'Hashes and computations do not mechanically verify the written analytic arguments; external report source hashes are frozen here for reproduction.'}
    (BASE/'results/verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print('All six cases and 96 transition samples accounted for.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--replay',action='store_true');run(p.parse_args().replay)
