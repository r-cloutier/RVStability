import numpy as np

def compute_fraction(numerator, denominator):
    '''
    Compute the fraction and plus/minus "1 sigma" uncertainties 
    given the numerator and denominator. This uses the Wilson score 
    interval correction so it is valid even as p -> 0 or 1 or 
    with small number statistics.
    '''
    N, M = float(numerator), float(denominator)
    p = N / M
    zs = np.array([1,-1])
    outvals = np.zeros(zs.size)  # +1sig, -1sig
    for i in range(zs.size):
	z = float(zs[i])
	F1 = 1./(1+z*z/M)
	F2 = p + .5*z*z/M + z*np.sqrt(p*(1-p)/M + .25*z*z/(M*M))
	outvals[i] = F1*F2 - p
    # Return frac, +1sig, -1sig
    return tuple(np.append(p, outvals))
