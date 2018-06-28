import numpy as np
import pylab as plt
from integrateK2_18 import *
import rebound
import glob
from scipy.ndimage import gaussian_filter

def plot_evolution(fname, incc_b1=88.1768716, incc_bm1=91.8231284,
		   incb_b1=89.233830, incb_bm1=90.766169, label=False, pltt=True):
    t,ac,ab,ec,eb,ic,ib = np.loadtxt(fname).T

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(311)
    ax.plot(t, ab, 'g-', label='K2-18b')
    ax.plot(t, ac, 'b-', label='K2-18c')
    ax.set_xticklabels('')
    ax.set_title('c = %.3f AU, b = %.3f AU'%(ac[0], ab[0]))
    ax.set_ylabel('Semimajor axis (AU)')
    ax.set_xlim((0,1e6))
    ax.legend(loc='center left')    

    ax = fig.add_subplot(312)
    ax.plot(t, ac, 'b-')
    ax.plot(t, ab, 'g-')
    ax.set_xticklabels('')
    ax.set_title('c = %.3f, b = %.3f'%(ec[0], eb[0]))
    ax.set_ylabel('Eccentricity')
    ax.set_yscale('log')
    ax.set_xlim((0,1e6))

    ax = fig.add_subplot(313)
    ax.plot(t, ic, 'b-', lw=.8)
    ax.plot(t, ib, 'g-', lw=.8)
    ax.fill_between(t, incb_b1, incb_bm1, color='g', alpha=.2)
    ax.fill_between(t, incc_b1, incc_bm1, color='b', alpha=.2)
    ax.set_ylabel('Inclination (deg)'), ax.set_xlabel('Time (Myrs)')
    ax.set_title('c = %.3f deg, b = %.3f deg'%(ic[0], ib[0]))
    ax.set_xlim((0,1))

    if label:
	plt.savefig(fname.replace('pickles','plots').replace('dat','png'))
    if pltt:
    	plt.show()
    plt.close('all')


def plot_archival_evolution(index, nplanets=2):
    sa = rebound.SimulationArchive('pickles/archived%.6d.bin'%index)
    times = np.zeros(0)
    eccs = np.zeros((0,nplanets))
    for i, sim in enumerate(sa):
	times = np.append(times, sim.t)
	eccs = np.insert(eccs, i, np.array([sim.particles[1].e,sim.particles[2].e]),axis=0)
    plt.plot(times, eccs, '-')
    plt.show()

def save_evolution():
    fs = np.array(glob.glob('pickles/archived0*'))
    for i in range(fs.size):
	sa = rebound.SimulationArchive(fs[i])
  	times, smas, eccs, incs = np.zeros(0), np.zeros((0,2)), np.zeros((0,2)), np.zeros((0,2))
	ind = 0
	for j, sim in enumerate(sa):
	    times = np.append(times, sim.t)
	    smas = np.insert(smas, ind, np.array([sim.particles[1].a, sim.particles[2].a]), axis=0)
            eccs = np.insert(eccs, ind, np.array([sim.particles[1].e, sim.particles[2].e]), axis=0)
            incs = np.insert(incs, ind, np.rad2deg([sim.particles[1].inc, sim.particles[2].inc]), axis=0)
	    ind += 1
   	# Save
	fname = fs[i].replace('archived','evolution').replace('bin','dat')
	np.savetxt(fname, np.array([times,smas[:,0],smas[:,1],eccs[:,0],eccs[:,1],incs[:,0],incs[:,1]]).T,
		   delimiter='\t', fmt='%.6f')
	#plot_evolution(fname, label=1, pltt=0)


def plot_Omegab_Omegac_stability(self, bins=10, pltt=True, label=False):
    '''Plot the fraction of stability on a 2d map as a function of
    each planet's Omega.'''
    Omega0, Omega1 = np.linspace(-180, 180, bins), np.linspace(-180, 180, bins)
    stability_frac = np.zeros((bins-1, bins-1))
    for i in range(bins-1):
	for j in range(bins-1):
	    thisbin = (self.Omega0[:,0]>=Omega0[i]) & (self.Omega0[:,0]<=Omega0[i+1]) & \
		      (self.Omega0[:,1]>=Omega1[j]) & (self.Omega0[:,1]<=Omega1[j+1])
	    try:
	    	stability_frac[i,j] = float(self.stable[thisbin].sum()) / self.stable[thisbin].size
	    except ZeroDivisionError:
		stability_frac[i,j] = 0.

    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    im = ax.pcolormesh(Omega0, Omega1, stability_frac, cmap=plt.get_cmap('Greys'))
    cbar = fig.colorbar(im)
    cbar.set_label('Fraction of Stable Systems')   
 
    Omegas = np.linspace(-180,180,10)
    ax.plot(Omegas, Omegas, 'g--', lw=2)
    ax.text(0, 0, '$\Omega_c=\Omega_b$', color='g', horizontalalignment='center', 
	    verticalalignment='bottom')
    offset = 180
    ax.plot(Omegas, Omegas+offset, 'g--', lw=2)
    ax.text(-90, -90+offset, '$\Omega_c=\Omega_b$\n+%i'%offset, color='g', horizontalalignment='center', 
            verticalalignment='bottom')
    ax.plot(Omegas, Omegas-offset, 'g--', lw=2)
    ax.text(90, 90-offset, '$\Omega_c=\Omega_b$\n-%i'%offset, color='g', horizontalalignment='center', 
            verticalalignment='bottom')
    offset = 110
    ax.plot(Omegas, Omegas+offset, 'b--', lw=2)
    ax.text(-45, -45+offset, '$\Omega_c=\Omega_b$\n+%i'%offset, color='b', horizontalalignment='center', 
            verticalalignment='bottom')
    ax.plot(Omegas, Omegas-offset, 'b--', lw=2)
    ax.text(45, 45-offset, '$\Omega_c=\Omega_b$\n-%i'%offset, color='b', horizontalalignment='center', 
            verticalalignment='bottom')


    ax.set_xlabel('$\Omega_b$'), ax.set_ylabel('$\Omega_c$')
    ax.set_xlim((-180,180)), ax.set_ylim((-180,180))
    #ax.invert_xaxis()

    if label:
	plt.savefig('plots/Omega_stability.png')
    if pltt:
	plt.show()
    plt.close('all')


def plot_ecc_ecc_stability(self, bins=20, sig=1, pltt=True):
    # Get stability fraction over grid
    eccs = np.linspace(0, .5, bins)
    Decc = np.diff(eccs)[0] / 2.
    stabfrac = np.zeros((bins-1, bins-1))
    for i in range(bins-1):
	for j in range(bins-1):
	    thisbin = (self.ecc0[:,0]>=eccs[i]) & (self.ecc0[:,0]<=eccs[i+1]) & \
		      (self.ecc0[:,1]>=eccs[j]) & (self.ecc0[:,1]<=eccs[j+1])
	    stabfrac[i,j] = float(self.stable[thisbin].sum()) / self.stable[thisbin].size
	    #print '\ne_c in [%.2f, %.2f]'%(eccs[i], eccs[i+1])
            #print 'e_b in [%.2f, %.2f]'%(eccs[j], eccs[j+1])
	    #print 'StabFrac = %.3f'%stabfrac[i,j]

    # Smooth the map
    stabfrac1 = gaussian_filter(stabfrac, 1)
    stabfrac2 = gaussian_filter(stabfrac, sig)    

    fig = plt.figure(figsize=(6.5,5))
    ax = fig.add_subplot(111)
    im = ax.pcolormesh(eccs, eccs, stabfrac2.T, cmap=plt.get_cmap('hot_r'),
		       vmin=0, vmax=1)
    cbar = fig.colorbar(im)
    cbar.set_label('Stability Fraction')

    ax.contour(eccs[1:]-Decc, eccs[1:]-Decc, stabfrac1.T, levels=[.5,1], colors='w')

    ax.set_xlabel('$e_c$')
    ax.set_ylabel('$e_b$')

    if pltt:
	plt.show()
    plt.close('all')


#if __name__ == '__main__':
#    save_evolution()
