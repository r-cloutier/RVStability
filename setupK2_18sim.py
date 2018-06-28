import rebound
import numpy as np
import rvs
from PyAstronomy.pyasl import foldAt

def setupK2_18sim(Ms, bjd0, Ps, T0s, Ks, hs, ks, incs_deg):
    '''Create a simulation of K2-18 with custom input for the 
    planetary parameters in the form of 1d arrays.'''
    # Initialize
    sim = rebound.Simulation()
    sim.integrator = "whfast"
    sim.units = ('AU','Msun','yr')

    # Get keplerian parameters
    eccs, omegas = hs*hs + ks*ks, np.arctan2(ks, hs)
    mps = rvs.kg2Msun(rvs.Mearth2kg(rvs.RV_Mp(Ps, Ms, Ks, ecc=eccs, inc=incs_deg)))
    smas = rvs.m2AU(rvs.semimajoraxis(Ps, Ms, mps))
    nplanets = Ps.size
    Omegas = np.random.uniform(0, 2*np.pi, nplanets)
    thetas = 2*np.pi * foldAt(bjd0, Ps, T0s)

    # Add star
    sim.add(m=Ms, hash='star')

    # Add planets
    for i in range(nplanets):
        sim.add(m=mps[i], a=smas[i], inc=np.deg2rad(incs_deg[i]), e=eccs[i], 
		omega=omegas[i], Omega=Omegas[i], theta=thetas[i])

    sim.move_to_com()

    return sim


def setupK2_18simv2(bjd0, Ms, eccs, incs_deg, outname, interval_yrs, 
		    random=True, DeltaOmegamax=180):
    '''Create a simulation of K2-18 with custom input for the 
    planetary parameters in the form of 1d arrays.'''
    # Initialize
    sim = rebound.Simulation()
    sim.integrator = "whfast"
    sim.units = ('AU','Msun','yr')
    sim.initSimulationArchive("%s.bin"%outname, interval=interval_yrs)

    # Get keplerian parameters
    Ps, Ks = np.zeros(2), np.zeros(2)
    while np.any(Ps<=0):
    	Ps  = np.array([8.964, 32.939621])
	if random:
	    Ps += np.array([np.random.randn()*8e-3, np.random.randn()*1e-4])
    T0s = np.array([7264.49, 7264.39142])
    if random:
	T0s += 2450000 + np.array([np.random.randn()*.48, np.random.randn()*5.9e-4])
    #while np.any(Ks<=0):
    #	Ks  = np.array([4.69, 3.18])
	#if random:
	 #   Ks += np.array([np.random.randn()*.7, np.random.randn()*.7])
    mps = np.array([7.6, 7.96])
    if random:
    	mps += np.array([np.random.randn()*1.3, np.random.randn()*1.3]) #rvs.kg2Msun(rvs.Mearth2kg(rvs.RV_Mp(Ps, Ms, Ks, ecc=eccs, inc=abs(incs_deg-90))))
    smas = rvs.m2AU(rvs.semimajoraxis(Ps, Ms, mps))
    mps = rvs.kg2Msun(rvs.Mearth2kg(mps))
    nplanets = Ps.size
    omegas = np.random.uniform(-np.pi, np.pi, nplanets)
    Omegas = np.random.uniform(-np.pi, np.pi, nplanets)
    #while abs(np.diff(np.rad2deg(Omegas))) > DeltaOmegamax:
	#Omegas = np.random.uniform(-np.pi, 2*np.pi, nplanets)
    thetas = 2*np.pi * foldAt(bjd0, Ps, T0s) - np.pi

    # Add star
    sim.add(m=Ms, hash='star')

    # Add planets
    for i in range(nplanets):
	#print '\nplanet %i'%(i+1)
	#print 'm, a, i_deg, e, omega_rad, Omega_rad, theta_rad'
	#print mps[i], smas[i], incs_deg[i]-90, eccs[i], omegas[i], Omegas[i], thetas[i]
        sim.add(m=mps[i], a=smas[i], inc=np.deg2rad(incs_deg[i]-90), e=eccs[i],
                omega=omegas[i], Omega=Omegas[i], theta=thetas[i])

    sim.move_to_com()

    RHill = (mps.sum()/(3.*Ms))**(1./3) * smas.mean()
    sim.exit_min_distance = float(RHill)

    return sim
