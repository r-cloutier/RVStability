'''
Functions to simulate the RVs of a planetary system.
'''
from imports import *
import rebound
import rvmodel

def years2days(t):
    return t*365.25
def days2years(t):
    return t/365.25

global T0
T0 = 2450000.

def initialize_system_parameters():
    '''
    Call this function to initialize the dictionary that will hold the 
    required planet parameters derived from the RV analysis. 

    The output format is for 2 planets but can be modified to include 
    additional planets.
    ''' 
    Pdict = {'b': 0., 'c':0.}
    T0dict = {'b': T0, 'c':T0}
    Kdict = {'b':0., 'c':0.}
    hdict = {'b':0., 'c':0.}
    kdict = {'b':0., 'c':0.}
    return Pdict, T0dict, Kdict, hdict, kdict


def setup_sim(Ms, planettheta, bjd0):
    '''
    Setup a simulation with a central star and planets.

    Parameters
    ----------
    `Ms` : float
	The mass of the central star in solar masses.
    `planettheta` : tuple of dict
	Tuple containing 5 or 6 dictionaries of each planet's 
	1) orbital period in days,
	2) time of inferior conjuction in BJD,
	3) RV semi-amplitude in m/s,
	4) h = sqrt(e)*cos(omega),
	5) k = sqrt(e)*sin(omega).
	If there are 6 dictionaries, the sixth and last is 
	6) orbital inclination in degrees.
    `bjd0` : scalar
	Reference epoch to begin the simulation at.

    Returns
    -------
    `sim` : rebound simulation object
	The rebound.Simulation containing the all injected particles.

    Example
    -------
    >>> Ms, planettheta = .1, initialize_system_parameters()
    >>> sim = setup_sim(Ms, planettheta, 2450000)

    '''
    # Initialize simulation
    sim = rebound.Simulation()
    sim.integrator = "whfast"
    sim.units = ('AU','Msun','yr')

    # Set timestep
    if len(planettheta) == 5:
	Ps, T0s, Ks, hs, ks = planettheta
	incs = {}
 	for p in Ps.keys():
	    incs[p] = 90.
    else:
	Ps, T0s, Ks, hs, ks, incs = planettheta
    hs, ks = np.array(hs.values()), np.array(ks.values())
    eccstmp, omegastmp = hs*hs + ks*ks, np.arctan2(ks, hs)
    eccs, omegas, ind = {}, {}, 0
    for p in Ps.keys():
	eccs[p], omegas[p] = eccstmp[ind], omegastmp[ind]
	ind += 1
    #sim.dt = days2years(np.min(Ps.values())) * 1e-1

    # Add star
    sim.add(m=Ms, hash='star')

    # Add planets
    for p in Ps.keys():
	mp = rvs.kg2Msun(rvs.Mearth2kg(rvs.RV_mp(Ps[p], Ms, Ks[p], ecc=eccs[p])))
        Pp = Ps[p]
        ap = rvs.semimajoraxis(Pp, Ms, 0)
        thetap = 2*np.pi * foldAt(bjd0, Pp, T0s[p])
        eccp, omegap, incp = eccs[p], omegas[p], np.deg2rad(incs[p])
        sim.add(m=mp, a=ap, inc=incp, e=eccp, omega=omegap, theta=thetap, hash=p)

    sim.move_to_com()

    return sim


def integrate_sim(bjd, sim):
    '''
    Integrate the system forward in time with outputs computed over 
    the window function.

    Parameters
    ----------
    `bjd` : array (Nobs,)
	Numpy array of the window function to output measurements at in BJD.
	For dynamical modeling of RVs this should be the window function
	of observations. For dynamical stability it should be much 
	larger.
    `sim` : rebound simulation object
	Output from setup_sim(); the simulation to integrate forward 
	in time.

    Returns
    -------
    `bjd` : array (Nobs,)
	The input window function; time arrray.
    `RVs` : array (Nobs,)
	Numpy array of the stellar RVs in m/s.
    `smas` : array (Nobs, Nplanets,)
	Numpy array of each planets' semimajor axis at the 
	times in bjd in AU.
    `eccs` : array (Nobs, Nplanets,)
        Numpy array of each planets' eccentricity at the
        times in bjd.
    `incs` : array (Nobs, Nplanets,)
        Numpy array of each planets' inclination at the
        times in bjd in degrees.

    '''
    # Get window function for outputs
    times = days2years(bjd-bjd.min())  ##np.linspace(sim.t, sim.t+Tfin_yrs, int(Nout))

    # Get timestep
    nparticles, dts = len(sim.particles.keys()), [np.diff(days2years(times)).min()*1e-1]
    for i in range(1,nparticles):
    	dts.append(sim.particles[i].P * 1e-1)
    sim.dt = np.min(dts)

    # Save output
    nparticles = len(sim.particles.keys())
    RVs = np.zeros(times.size)
    smas, eccs, incs = np.zeros((times.size, nparticles-1)), \
                       np.zeros((times.size, nparticles-1)), \
                       np.zeros((times.size, nparticles-1))
    for i in range(times.size):
        sim.integrate(times[i])
        RVs[i] = -sim.particles['star'].vx
        for j in range(nparticles-1):
            smas[i,j] = sim.particles[j+1].a  # AU
            eccs[i,j] = sim.particles[j+1].e
	    incs[i,j] = np.rad2deg(sim.particles[j+1].inc)
    
    # Convert units (AU/yr -> m/s)
    RVs = rvs.AU2m(RVs) / (365.25*24*60*60)
    
    return bjd, RVs, smas, eccs, incs


def get_keplerians(planettheta, Ms, bjd):
    '''Add keplerians together to compare to dynamical integration.'''
    Ps, T0s, mps, incs, eccs = planettheta
    RVs = np.zeros(bjd.size)
    for p in Ps.keys():
        K = rvs.RV_K(Ps[p], Ms, mps[p], ecc=eccs[p], inc=incs[p])
        h, k = np.sqrt(eccs[p]), 0.
        RVs += rvmodel.get_rv2((Ps[p], T0s[p], K, h, k), bjd)
    return bjd, RVs


# trappist-1 test
if __name__ == '__main__':
    Ps, T0s, mps, incs, eccs = initialize_system_parameters()
    Ps = {'b':1.51087, 'c':2.42183}
    T0s = {'b':2457322.51736, 'c':2457282.80728}
    mps['b'], mps['c'] = .85, 1.38
    incs['b'], incs['c'] = 89.65, 89.67
    planettheta = (Ps, T0s, mps, incs, eccs)
    Ms = .0802
    sim = setup_sim(Ms, planettheta)
    bjd1, rv1, a1, e1 = integrate_sim(sim, 1./12)
    bjd2, rv2 = get_keplerians(planettheta, Ms, bjd1)
