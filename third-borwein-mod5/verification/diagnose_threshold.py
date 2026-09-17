"""Finite 80-digit diagnostics of v0.3-v0.4 lemmas; NOT uniform certificates."""
import argparse
from common import metadata, write_result
from mpmath import mp
from diagnose_analytic import R, lam_root
from parameters import THRESHOLD


def diagnose(output_dir=None):
    mp.dps=80
    fmt=lambda value:mp.nstr(value,45)
    geometric=[]
    for M in (1,2,17,257):
        for real in ('0','0.001','0.3'):
            for imag in ('-3.1','-0.1','0','0.1','3.1'):
                d=mp.mpc(real,imag)
                a=M*d
                direct=sum(mp.exp(-d*j) for j in range(1,M+1))
                integral=M if a==0 else -M*mp.expm1(-a)/a
                error=direct-integral
                predicted=mp.mpf(0) if a==0 else -mp.expm1(-a)*(1/mp.expm1(d)-1/d)
                assert abs(d)<mp.mpf('3.2') and abs(error)<2
                assert abs(error-predicted)<mp.mpf('1e-70')
                geometric.append({'M':M,'d_real':real,'d_imag':imag,
                                  'absolute_sum_error':fmt(abs(error)),
                                  'identity_discrepancy':fmt(abs(error-predicted))})

    local=[]
    quadratures=[]
    for tau_string in ('0','0.5','2','8'):
        tau=mp.mpf(tau_string)
        def moment(x,order):
            weights=[mp.exp(-tau*x*j) for j in range(5)]
            denom=sum(weights)
            mean=sum(j*weights[j] for j in range(5))/denom
            if order==1:return mean
            return sum((j-mean)**order*weights[j] for j in range(5))/denom
        V=mp.quad(lambda x:x*x*moment(x,2),[0,1])
        W3=mp.quad(lambda x:x**3*moment(x,2),[0,1])
        W4=mp.quad(lambda x:x**4*moment(x,2),[0,1])
        Rp=-mp.quad(lambda x:x*moment(x,1),[0,1])
        R3=-mp.quad(lambda x:x**3*moment(x,3),[0,1])
        R0=R(tau)
        for n in (10**4,10**6,THRESHOLD):
            a=n*V
            for t_string in ('-0.1','-0.01','0.01','0.1'):
                t=mp.mpf(t_string)
                F=n*(R(tau-1j*t)-R0+1j*Rp*t)
                G=-a*t*t/2
                D=F-G
                C3=1j*n*R3*t**3/6
                assert abs(D)<=a*abs(t)**3
                assert abs(D-C3)<=mp.mpf('1.75')*a*t**4
                assert abs(D)<=n*W3*abs(t)**3
                assert abs(D-C3)<=mp.mpf('1.75')*n*W4*t**4
                assert mp.re(F)<=-mp.mpf('0.49')*a*t*t
                for ell in (1,2,3,4):
                    amplitude=lambda s:mp.exp(lam_root(ell,tau-1j*s))
                    A0=amplitude(0)
                    A1=mp.diff(amplitude,0)
                    At=amplitude(t)
                    assert abs(At-A0-A1*t)<=t*t
                    corrected=At*mp.exp(F)-A0*mp.exp(G)-(A0*C3+A1*t)*mp.exp(G)
                    budget=(mp.mpf('2.75')*a*t**4+a*a*t**6/2+t*t)*mp.exp(-mp.mpf('0.49')*a*t*t)
                    ratio=abs(corrected)/budget
                    assert ratio<1
                    phase_x=mp.exp(-tau)
                    amp1=mp.mpf('1.1')*mp.mpf('0.8')*phase_x/mp.mpf('0.9')
                    amp2=mp.mpf('1.1')/2*(mp.mpf('0.8')*phase_x/mp.mpf('0.9')**2
                            +(mp.mpf('0.8')*phase_x/mp.mpf('0.9'))**2)
                    assert abs(At-A0)<=amp1*abs(t)
                    assert abs(At-A0-A1*t)<=amp2*t*t
                    weighted_budget=(mp.mpf('1.75')*n*W4*t**4+n*n*W3**2*t**6/2
                                      +amp1*n*W3*t**4+amp2*t*t)*mp.exp(-mp.mpf('0.49')*a*t*t)
                    weighted_ratio=abs(corrected)/weighted_budget
                    assert weighted_ratio<1
                    local.append({'tau':tau_string,'n':n,'t':t_string,'ell':ell,
                                  'corrected_integrand_over_bound':fmt(ratio),
                                  'v04_weighted_integrand_over_bound':fmt(weighted_ratio),
                                  'fourth_remainder_over_bound':fmt(abs(D-C3)/(mp.mpf('1.75')*a*t**4))})

        # Actual scaled integrals provide an independent numerical check of
        # the assembled approximation. The interval is inside |t|<=0.1;
        # these finite quadratures are not used to prove the uniform estimate.
        if tau_string in ('0','8'):
            for n in (10**7,2*10**7):
                scale=mp.sqrt(n*V)
                assert 8/scale<mp.mpf('0.1')
                for ell in (1,2):
                    A0=mp.exp(lam_root(ell,tau))
                    def integrand(u):
                        t=u/scale
                        F=n*(R(tau-1j*t)-R0+1j*Rp*t)
                        return mp.exp(lam_root(ell,tau-1j*t)+F)-A0*mp.exp(-u*u/2)
                    error=abs(mp.quad(integrand,[-8,-4,0,4,8])/mp.sqrt(2*mp.pi))
                    phase_x=mp.exp(-tau)
                    amp1=mp.mpf('1.1')*mp.mpf('0.8')*phase_x/mp.mpf('0.9')
                    amp2=mp.mpf('1.1')/2*(mp.mpf('0.8')*phase_x/mp.mpf('0.9')**2
                            +(mp.mpf('0.8')*phase_x/mp.mpf('0.9'))**2)
                    beta=mp.mpf('0.49')
                    weighted_H=(21*W4/(16*beta**2*V**2)+15*W3**2/(16*beta**3*V**3)
                                +3*amp1*W3/(4*beta**2*V**2)+amp2/(2*beta*V))/mp.sqrt(2*beta)
                    assert error<weighted_H/n
                    quadratures.append({'tau':tau_string,'n':n,'ell':ell,
                                        'scaled_interval':'[-8,8]',
                                        'normalized_integral_error':fmt(error),
                                        'n_times_error':fmt(n*error),
                                        'one_root_error_bound':fmt(weighted_H/n)})
    result={**metadata('ordinary multiprecision diagnostics of v0.3-v0.4 lemmas'),
            'decimal_precision':80,'status':'passed',
            'scope':'Finite samples and numerical quadrature, not interval or uniform proofs.',
            'geometric_samples':geometric,'local_remainder_samples':local,
            'scaled_integral_samples':quadratures}
    write_result('threshold_diagnostics.json',result,output_dir)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-dir')
    diagnose(parser.parse_args().output_dir)
