import rebound
import numpy as np
import rvs
from PyAstronomy.pyasl import foldAt


def setup_simv2(bjd0, Ms, Ps, ePs, T0s, eT0s, mps, emps, eccs, incs_deg,
                outname, interval_yrs, random=True, DeltaOmegamax=180):
    '''Create a simulation of a 2-planet system with custom input for the 
    planetary parameters in the form of 1d arrays.
    example T0s = np.array([7264.49, 7264.39142])
    '''
    # Initialize
    sim = rebound.Simulation()
    sim.integrator = "whfast"
    sim.units = ('AU','Msun','yr')
    sim.initSimulationArchive("%s.bin"%outname, interval=interval_yrs)

    # Get keplerian parameters
    Psin, ePs = np.ascontiguousarray(Ps) + 0, np.ascontiguousarray(ePs)
    assert Ps[0] < Ps[1]
    Ps, Ks = np.zeros(2), np.zeros(2)
    while np.any(Ps<=0):
    	Ps = Psin
	if random:
	    Ps += np.array([np.random.randn()*ePs[0], np.random.randn()ePs[1]])
    if random:
	T0s += 2450000 + np.array([np.random.randn()*eT0s[0],
                                   np.random.randn()*eT0s[1]])
    if random:
    	mps += np.array([np.random.randn()*emps[0], np.random.randn()*emps[1]])
    smas = rvs.m2AU(rvs.semimajoraxis(Ps, Ms, mps))
    mps = rvs.kg2Msun(rvs.Mearth2kg(mps))
    nplanets = Ps.size
    omegas = np.random.uniform(-np.pi, np.pi, nplanets)
    Omegas = np.random.uniform(-np.pi, np.pi, nplanets)
    thetas = 2*np.pi * foldAt(bjd0, Ps, T0s) - np.pi

    # Add star
    sim.add(m=Ms, hash='star')

    # Add planets
    for i in range(nplanets):
        sim.add(m=mps[i], a=smas[i], inc=np.deg2rad(incs_deg[i]-90), e=eccs[i],
                omega=omegas[i], Omega=Omegas[i], theta=thetas[i])

    sim.move_to_com()

    RHill = (mps.sum()/(3.*Ms))**(1./3) * smas.mean()
    sim.exit_min_distance = float(RHill)

    return sim
