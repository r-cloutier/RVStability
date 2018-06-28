from DynamicalSim import *
import numpy as np
import pylab as plt
import rvs, sys, glob, os
from setupK2_18sim import setupK2_18simv2
import cPickle as pickle

#global incb, Rs, Ms
#incb, Rs, Ms = 89.5785, .411, .359


#class K2_18integration
class system_integration:

    def __init__(self, folder, index, index2, Ms, eMs, Rs, bjd0, smac, incs, eincs, Nyrs=1e6, Nout=500):
        '''
        folder = top level directory name to save simulation results
        index = job index 1 for output filename
        index2 = job index 2 for output filename
        Ms = host stellar mass in MSun (eMs is its uncertainty)
        Rs = host stellar radius in RSun
        bjd0 = initial epoch to begin simulation (e.g. the first RV observation epoch)
        smac = planet c semi-major axis in AU
        incs = list of known planet inclinations (2 entries if both planets are transiting) (eincs are uncertainties)
        '''
	self.index, self.index2, self.Nyrs, self.Nout = int(index), int(index2), int(Nyrs), int(Nout)
	self.Rs, self.Ms = float(Rs), 0.
	while self.Ms <= 0:
	    self.Ms = Ms + np.random.randn() * eMs
	self.bjd0 = bjd
        self.incs, self.eincs = np.ascontiguousarray(incs), np.ascontiguousarray(eincs)
        assert self.incs.size == self.eincs.size
        self.incb, self.eincb = float(self.incs[0]), float(self.eincs[0])
        if self.incs > 0:
            self.nplanets = 2
            self.incc, self.eincc = float(self.incs[1]), float(self.eincs[1])
        else:
            self.nplanets = 1
        self.DONE = False	

	# Setup simulation
	self.folder = folder
        self.outname = '%s/SimArchive/archived%.4d_%.4d'%(self.folder, self.index, self.index2)
	sim = setupsim(self, self.outname)
	self.mps,self.sma0,self.ecc0,self.inc0,self.omega0,self.Omega0,self.theta0 = get_initial_parameters(sim)
	self.inc0 += 90
	self.Rhill = get_Rhill_init(self.mps, self.Ms, self.sma0)
	self.pickleobject()

	# Integrate simulation
	print 'Integrating system...'
	self.bjds,self.RVs,self.smas,self.eccs,self.incs,self.dist,self.stable = \
		integrate_sim(np.linspace(0, rvs.yrs2days(self.Nyrs), self.Nout) + self.bjd0, sim)
	self.bs = np.zeros(self.incs.shape)
	for i in range(self.bjds.size):
	    self.bs[i] = rvs.impactparam_inc(rvs.AU2m(self.smas[i])/rvs.Rsun2m(Rs), self.incs[i]+90, ecc=self.eccs[i])
	self.DONE = True
	self.pickleobject()


    def pickleobject(self):
	f = open('%s/pickles/Nbodysimulation%.4d_%.4d'%(self.folder, self.index, self.index2), 'wb')
	pickle.dump(self, f)
	f.close()


#def setupsim_corr(self, outname, N=50):
#    # Sample eccentricities logarithmically
#    eccbs = np.logspace(-3, -1, N)
#    ecccs = np.logspace(-3, -1, N)
#    #eccbs = np.linspace(1e-3, .59, N)
#    #ecccs = np.linspace(1e-3, .43, N)
#    eccs = np.array([ecccs[self.index-1], eccbs[self.index-1]])
#    eccs += np.random.randn(eccs.size) * eccs/10.
#    while np.any(eccs<0) or np.any(eccs>=1):
#	eccs = np.array([ecccs[self.index-1], eccbs[self.index-1]])
#    	eccs += np.random.randn(eccs.size) * eccs/10.
#    incs_deg = np.array([sampleK218cincdeg(0.06, eccs[0]), incb+np.random.randn()*.0084])
#    # How often to save snapshots
#    interval_yrs = self.Nyrs / self.Nout
#    return setupK2_18simv2(self.bjd0, self.Ms, eccs, incs_deg, outname, int(interval_yrs))


def setupsim(self, outname):
    # Sample eccentricities logarithmically
    eccupperlim = 1.
    eccs = np.random.uniform(0, eccupperlim, 2)
    incbtmp = self.incb + np.random.randn() * self.eincb
    if self.nplanets == 1:
        incs_deg = np.array([sampleincdeg_nontransiting(incbtmp, self.smac, self.Rs), incbtmp])
    else:
        incs_deg = np.array([self.incc + np.random.randn() * self.eincc])
    # How often to save snapshots
    interval_yrs = self.Nyrs / self.Nout
    return setupK2_18simv2(self.bjd0, self.Ms, eccs, incs_deg, outname, int(interval_yrs))


def loadpickle(fname):
    f = open(fname, 'rb')
    self = pickle.load(f)
    f.close()
    return self


#def sampleK218cincdeg(smac, eccc):
 #   '''Sample the orbital inclination from its eccentrcitity (<i^2>=<e^2>/4)
 #   ensuring that the planet does not transit. 
 #   smac in AU.'''
 #   b, ind = 0, 0
 #   while abs(b) < 1:
#	incc = incb + np.rad2deg(np.sqrt(eccc**2/4.)) * np.random.choice([-1,1])
#	# Add noise to inc if cannot find solution
#	if ind > 100:
#	    incc += np.random.randn() * ind/100.
 #   	b = rvs.impactparam_inc(smac/rvs.m2AU(rvs.Rsun2m(Rs)), incc)
#	ind += 1
 #   return incc

 
def sampleincdeg_nontransiting(incbb, smac, Rs, sig=1.5):
    '''Sample incc from a Gaussian but reject inclinations that result in a transit.'''
    incc = incbb + np.random.randn()*sig
    b = rvs.impactparam_inc(smac/rvs.m2AU(rvs.Rsun2m(Rs)), incc)
    while abs(b) < 1:
	incc = incbb + np.random.randn()*sig
    	b = rvs.impactparam_inc(smac/rvs.m2AU(rvs.Rsun2m(Rs)), incc)
    return incc


def RHill_from_sim(sim):
    '''Compute the Hill radius from a simulation with 1 star plus 2 planets.'''
    ps = sim.particles
    star, p1, p2 = ps['star'], ps[1], ps[2]
    return ((p1.m+p2.m)/(3.*star.m))**(1./3) * (p1.a+p2.a)/2.   # units of a (e.g. AU)


def get_Rhill_init(mps, Ms, smas):
    '''mps in Msun, Ms in Msun, smas in AU.'''
    return (mps.sum()/(3.*Ms))**(1./3) * smas.sum()/2.


def get_initial_parameters(sim):
    '''Return the initial planetary parameters.'''
    ps = sim.particles
    mps, smas, eccs, incs = np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2)
    omegas, Omegas, thetas = np.zeros(2), np.zeros(2), np.zeros(2)
    for i in range(2):
	mps[i] = ps[i+1].m
        smas[i] = ps[i+1].a
        eccs[i] = ps[i+1].e
        incs[i] = np.rad2deg(ps[i+1].inc)
        omegas[i] = np.rad2deg(ps[i+1].omega)
        Omegas[i] = np.rad2deg(ps[i+1].Omega)
        thetas[i] = np.rad2deg(ps[i+1].theta)
    return mps, smas, eccs, incs, omegas, Omegas, thetas


def do_i_run_this_simulation(fname):
    if not os.path.isfile(fname):
	return True
    else:
	d = loadpickle(fname)
	return False if d.DONE else True


if __name__ == '__main__':
    # Neccentricities == 10 in log space
    # so index in [1,1e2] with 10 simulations per ecc value
    Nyrs, Nsim_per_ecc, Nout = 1e6, 50, 1e3
    folder = 'pickles_uncorr_ALLecc_lin'
    index = int(sys.argv[1])
    for i in range(1,Nsim_per_ecc+1):
	if do_i_run_this_simulation('%s/SimArchive/archived%.4d_%.4d'%(folder, index, i)):
    	    self = system_integration(folder, index, i, Ms, eMs, Rs, bjd0, smac, incs, eincs, Nyrs=Nyrs, Nout=Nout)
