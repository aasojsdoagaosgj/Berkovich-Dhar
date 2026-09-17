"""Rational checks of the scalar margins in ENDPOINT_PROOF.md."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
from certify_infinite import exp_lower


def run():
    pi2_low,pi2_high=Q(39,4),Q(10)
    a_low,a_high=pi2_low/3,pi2_high/3
    k2=Q(12)*(36+12+2)/399
    k3=Q(12)*(216+108+36+6)*400/(399**2)
    v_upper=Q(2)*Q(7,4)/3/400
    sine_lower=Q(1,6)-Q(1,100)
    checks={
        'exp6_exceeds_400':exp_lower(6,24)>400,
        'k2_bound':k2==Q(600,399),
        'g_lower':8*a_low/3-k2>Q(20,3),
        'g_upper':8*a_high/3<Q(28,3),
        'third_derivative':8*a_high+k3<39,
        'whole_root_arc_decay':4*a_low/(3*(1+Q(9,4)))-k2/2>Q(8,15),
        'real_phase':v_upper<Q(1,100),
        'residue2_amplitude':2*(sine_lower-sine_lower**3/6)>Q(1,4),
        'saddle_tau_square':(4*a_high/3)/(1-Q(1,9))<=5,
    }
    assert all(checks.values())
    result={'status':'passed','checks':checks,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'dependency_sha256':{'certify_infinite.py':hashlib.sha256(Path(__file__).with_name('certify_infinite.py').read_bytes()).hexdigest()},
        'scope':'Scalar margins only; the uniform asymptotic proof is written in ENDPOINT_PROOF.md.'}
    Path(__file__).with_name('results').joinpath('endpoint_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':run()
