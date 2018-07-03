'''
Compute the critical value for Hill stability given the eccentricities of
two planets.
see Gladman 1993
'''
from imports import *

def compute_mus(m1, m2, M3):
    return rvs.Mearth2kg(m1)/rvs.Msun2kg(M3), rvs.Mearth2kg(m2)/rvs.Msun2kg(M3)

def compute_alpha(mu1, mu2):
    return mu1 + mu2

def compute_delta(a1, a2):
    a1, a2 = float(a1)/a1, float(a2)/a1
    return np.sqrt(a1 + (a2-a1))

def compute_gammas(e1, e2):
    return np.sqrt(1.-e1**2), np.sqrt(1.-e2**2)

def compute_pcrit(mu1, mu2, alpha):
    return 1. + 3**(4./3) * mu1*mu2 / alpha**(4./3)

def compute_p(alpha, mu1, mu2, delta, gamma1, gamma2):
    return 1./alpha**3 * (mu1+mu2/delta**2) * (mu1*gamma1 + mu2*gamma2*delta)**2

def compute_e2_at_crit(m1, m2, M3, a1, a2, e1):
    '''Evaluate the value of e2 for which the Hill stability bifurcation 
    occurs given the other keplerian parameters of the system.'''
    # Get constant values
    mu1, mu2 = compute_mus(m1, m2, M3)
    alpha = compute_alpha(mu1, mu2)
    pcrit = compute_pcrit(mu1, mu2, alpha)
    delta = compute_delta(a1, a2)

    e2s = np.linspace(0,1,1000)
    p_pcrit = np.zeros(e2s.size)
    for i in range(e2s.size):
	gamma1, gamma2 = compute_gammas(e1, e2s[i])
     	p_pcrit[i] = compute_p(alpha, mu1, mu2, delta, gamma1, gamma2) / pcrit

    print e1, p_pcrit.min()
    if p_pcrit.min() > 1 or p_pcrit.max() < 1:
	return np.nan
    else:
	fint = interp1d(p_pcrit, e2s)
	return fint(1)

def compute_pcrit_curve(m1, m2, M3, a1, a2):
    '''Compute the isocurve of p=pcrit (p/pcrit==1) as a function of 
    the eccentricities of each planet.'''
    e1s = np.linspace(0, 1, 1000)
    e2s = np.zeros(e1s.size)
    for i in range(e1s.size):
	e2s[i] = compute_e2_at_crit(m1, m2, M3, a1, a2, e1s[i])
    return e1s, e2s

#if __name__ == '__main__':
#    m1, m2, M3 = 7.6, 7.96, .359
#    a1, a2 = .06, .143
#    ecs, ebs = compute_pcrit_curve(m1, m2, M3, a1, a2)
