"""Independent finite 80-digit checks of grouped roots and B2 EM; not certificates."""
import argparse
from common import metadata, write_result
from mpmath import mp
from diagnose_analytic import R, Rp, lam_root
from parameters import THRESHOLD


def diagnose(output_dir=None):
    mp.dps=80
    fmt=lambda x:mp.nstr(x,45)
    zeta=mp.exp(2j*mp.pi/5)
    beta=mp.mpf('0.39'); d=mp.mpf('0.75'); M=mp.mpf('1.54')
    def bern(j):
        alpha=mp.mpf(j)/5
        return alpha**2-alpha+mp.mpf(1)/6
    def D(ell,z):
        x=mp.exp(-z)
        return z/2*sum(bern(j)*(zeta**(ell*j)*x/(1-zeta**(ell*j)*x)
                       -zeta**(ell*j)/(1-zeta**(ell*j))) for j in range(1,5))
    def C(ell,z):
        x=mp.exp(-z)
        return sum(bern(j)*zeta**(ell*j)/(1-zeta**(ell*j)*x) for j in range(1,5))/2
    def psi(a,z):
        return sum(zeta**(-a*ell)*mp.exp(lam_root(ell,z)) for ell in range(1,5))
    def theta(a,z):
        return sum(zeta**(-a*ell)*mp.exp(lam_root(ell,z))*D(ell,z) for ell in range(1,5))
    em=[]; local=[]; quadratures=[]
    for ts in ('0','0.5','5.5'):
        tau=mp.mpf(ts); x=mp.exp(-tau); Z=mp.sqrt(tau*tau+mp.mpf('0.16'))
        A1=M*mp.mpf('0.8')*x/d
        A2=M/2*(mp.mpf('0.8')*x/d**2+(mp.mpf('0.8')*x/d)**2)
        for n in (8,32,128):
            for tt in ('0','0.4'):
                z=tau-1j*mp.mpf(tt)
                for ell in (1,2):
                    finite=sum(mp.log(1-zeta**(ell*k)*mp.exp(-z*k/(5*n)))
                               for k in range(1,5*n+1) if k%5)
                    error=abs(finite-n*R(z)-lam_root(ell,z)-D(ell,z)/n)
                    discrepancy=abs(D(ell,z)+z/30-z*mp.exp(-z)*C(ell,z))
                    assert error<mp.mpf(3400)/n**2 and discrepancy<mp.mpf('1e-70')
                    em.append({'tau':ts,'t':tt,'n':n,'ell':ell,
                               'B2_log_remainder':fmt(error),'bound':fmt(mp.mpf(3400)/n**2),
                               'correction_identity_discrepancy':fmt(discrepancy)})
        def moment(y,order):
            w=[mp.exp(-tau*y*j) for j in range(5)]; total=sum(w)
            mean=sum(j*w[j] for j in range(5))/total
            return sum((j-mean)**order*w[j] for j in range(5))/total
        V=mp.quad(lambda y:y*y*moment(y,2),[0,1])
        W3=mp.quad(lambda y:y**3*moment(y,2),[0,1])
        W4=mp.quad(lambda y:y**4*moment(y,2),[0,1])
        R3=-mp.quad(lambda y:y**3*moment(y,3),[0,1])
        R0=R(tau); derivative=Rp(tau)
        for a in range(5):
            P=psi(a,tau); p=abs(P)
            linear=mp.diff(lambda t:psi(a,tau-1j*t),0)
            Ta=Z/30*(p+mp.mpf('1.6')*A1)+8*M*Z*x/(25*d)
            Ha=(p*(240*W4/(32*beta**2*V**2)+5*23**2*W3**2/(192*beta**3*V**3))
                +4*(23*A1*W3/(8*beta**2*V**2)+A2/(2*beta*V)))/mp.sqrt(2*beta)
            for n in (1000,THRESHOLD):
                for tt in ('-0.4','-0.2','0.2','0.4'):
                    t=mp.mpf(tt); z=tau-1j*t
                    F=n*(R(z)-R0+1j*derivative*t); G=-n*V*t*t/2
                    C3=1j*n*R3*t**3/6
                    error=abs(psi(a,z)*mp.exp(F)-(P+P*C3+linear*t)*mp.exp(G))
                    budget=(p*(10*n*W4*t**4+mp.mpf(23)**2/72*n*n*W3**2*t**6)
                            +mp.mpf(46)/3*A1*n*W3*t**4+4*A2*t*t)*mp.exp(-beta*n*V*t*t)
                    th=abs(theta(a,z)); assert th<Ta and error<budget
                    local.append({'tau':ts,'residue':a,'n':n,'t':tt,
                                  'grouped_error_over_bound':fmt(error/budget),
                                  'EM_first_over_bound':fmt(th/Ta)})
            if ts in ('0','5.5') and a in (3,4):
                n=THRESHOLD; scale=mp.sqrt(n*V)
                assert 8/scale<mp.mpf('0.4')
                def integrand(u):
                    t=u/scale; z=tau-1j*t
                    F=n*(R(z)-R0+1j*derivative*t)
                    return (psi(a,z)+theta(a,z)/n)*mp.exp(F)-P*mp.exp(-u*u/2)
                error=abs(mp.quad(integrand,[-8,-4,0,4,8])/mp.sqrt(2*mp.pi))
                bound=(Ha+Ta/mp.sqrt(2*beta))/n
                assert error<bound
                quadratures.append({'tau':ts,'residue':a,'n':n,'scaled_interval':'[-8,8]',
                                    'normalized_error':fmt(error),'bound':fmt(bound),
                                    'error_over_phase':fmt(error/p)})
        print(f'grouped diagnostics tau={ts} passed',flush=True)
    result={**metadata('ordinary multiprecision diagnostics of v0.6 grouped lemmas'),
            'status':'passed','decimal_precision':80,
            'scope':'Finite samples and quadrature; not uniform or interval certificates.',
            'EM_samples':em,'grouped_samples':local,'scaled_integral_samples':quadratures}
    write_result('grouped_diagnostics.json',result,output_dir)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output-dir')
    diagnose(parser.parse_args().output_dir)
