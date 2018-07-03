from integrate_system import *
from scipy.ndimage import gaussian_filter
import matplotlib.gridspec as gridspec
from matplotlib.colors import rgb2hex
from Hillstability_eccentricityspace import *
from uncertainties import unumpy

class integrate_results:

    def __init__(self, folder):

	self.folder = folder
	self.get_results()
	self.pickleobject()


    def get_pickles(self):
	'''Get all the pickles in folder.'''
	self.picklestmp = np.array(glob.glob('%s/pickles/N*'%self.folder))
	self.npicklestmp = self.picklestmp.size


    def get_results(self):
	'''Compile parameters and stability results.'''
	self.get_pickles()
	self.pickles = np.zeros(0, dtype=str)
	# Get initial parameters
	self.Mss = np.zeros(0)
	self.mps = np.zeros((0,2))
	self.sma0 = np.zeros((0,2))
        self.ecc0 = np.zeros((0,2))
        self.inc0 = np.zeros((0,2))
        self.omega0 = np.zeros((0,2))
        self.Omega0 = np.zeros((0,2))
	self.theta0 = np.zeros((0,2))
	self.stable = np.zeros(0, dtype=bool)
	ind = 0
	for i in range(self.npicklestmp):
	    if i % 10 == 0:
		print float(i)/self.npicklestmp
	    d = loadpickle(self.picklestmp[i])
	    if d.DONE:
		self.pickles = np.append(self.pickles, self.picklestmp[i])
		self.Mss = np.append(self.Mss, d.Ms)
		self.mps = np.insert(self.mps, ind, rvs.kg2Mearth(rvs.Msun2kg(d.mps)), axis=0)
                self.sma0 = np.insert(self.sma0, ind, d.sma0, axis=0)
                self.ecc0 = np.insert(self.ecc0, ind, d.ecc0, axis=0)
                self.inc0 = np.insert(self.inc0, ind, d.inc0, axis=0)
                self.omega0 = np.insert(self.omega0, ind, d.omega0, axis=0)
                self.Omega0 = np.insert(self.Omega0, ind, d.Omega0, axis=0)
                self.theta0 = np.insert(self.theta0, ind, d.theta0, axis=0)
		self.stable = np.append(self.stable, stability_check(d))
		ind += 1

	self.npickles = self.pickles.size
	self.Delta_mps = np.diff(self.mps, axis=1).reshape(self.npickles)
	self.Delta_Omegas = np.diff(self.Omega0, axis=1).reshape(self.npickles)
        self.Delta_incs = np.diff(self.inc0, axis=1).reshape(self.npickles)
        self.Delta_eccs = np.diff(self.ecc0, axis=1).reshape(self.npickles)


    def pickleobject(self):
	f = open('%s/integrate_results'%self.folder, 'wb')
	pickle.dump(self, f)
	f.close()



def stability_check(simulation):
    '''Ensure that the simulation is actually stable.'''
    if np.any(simulation.smas < 0):
     	return False
    else:
	return simulation.stable


def plot_one_stability_histograms(folder, stable, arr, bins=10, xlog=False, ylog=False, 
				  xlabel='', label=False, pltt=True):
    '''Plot stacked hisogram of the planet parameter showing the fraction of 
    stable systems.
    e.g. plot_one_stability_histograms(self.stable, self.ecc0[:,0])
    '''
    assert len(arr.shape) == 1
    cols = ['#a63603','#F2F22F']

    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(111)
    bins = np.logspace(np.log10(arr.min()), np.log10(arr.max()), bins) if xlog else bins
    y = ax.hist([arr[stable==1], arr[stable==0]], bins=bins, stacked=1, 
		log=ylog, label=['Stable','Unstable'], color=cols[::-1])
    Nstable, Ntotal, xbins = y[0][0], y[0][1], y[1]
    xbins = xbins[:-1] + np.diff(xbins)[0]/2
    frac_unstable = (Ntotal-Nstable) / Ntotal
    efrac_unstable = np.sqrt(Ntotal-Nstable) / Ntotal
    for i in range(xbins.size):
    	ax.text(xbins[i], Ntotal[i], '%.2f\n$\pm$%.2f'%(frac_unstable[i],efrac_unstable[i]), 
		horizontalalignment='center', verticalalignment='bottom')

    if xlog:
    	ax.set_xscale('log')
    if xlabel != '':
	ax.set_xlabel(xlabel)
    ax.legend(loc='upper right')

    if label and xlabel != '':
	plt.savefig('plots/%s_stabilityhist_%s.png'%(folder,xlabel))
    if pltt:
	plt.show()
    plt.close('all')


def plot_ecc_stability(folder, stable, eccs, m1, m2, M, a1, a2,
		       sig=2, bins=21, pltt=True, label=False):

    eccarr = np.linspace(0,1,bins)
    Decc = np.diff(eccarr)[0] / 2.
    print 'bin width = %.4f'%Decc
    #ecccont = np.arange(0,.6,np.diff(eccarr)[0])
    stabfrac = np.zeros((bins-1, bins-1))
    for i in range(bins-1):
        for j in range(bins-1):
            thisbin = (eccs[:,0] >= eccarr[i]) & (eccs[:,0] < eccarr[i+1]) & \
                      (eccs[:,1] >= eccarr[j]) & (eccs[:,1] < eccarr[j+1])
            stabfrac[i,j] = float(stable[thisbin].sum()) / stable[thisbin].size

    # Smooth the map
    #stabfrac1 = gaussian_filter(stabfrac, 1)
    stabfrac2 = stabfrac#gaussian_filter(stabfrac, sig)

    # Plot 2d stability map
    fig = plt.figure(figsize=(7,6.9))
    gs = gridspec.GridSpec(6,5) # rows, cols
    ax1 = plt.subplot(gs[1:-1,:-1])
    cm = plt.get_cmap('hot')
    im = ax1.pcolormesh(eccarr, eccarr, stabfrac2.T, cmap=cm,
                       	vmin=0, vmax=1)
    cax = plt.axes([.1, .09, .71, .04])
    ticks = np.arange(0,1.01,.1)
    cbar = fig.colorbar(im, cax=cax, orientation='horizontal', ticks=ticks)
    cbar.set_label('Stability Fraction')
    ticklabels = np.array(['%.1f'%i for i in ticks])
    ticklabels[np.arange(1,10,2)] = ''
    cbar.ax.set_xticklabels(ticklabels)
    #CS = ax1.contour(eccarr[1:]-Decc, eccarr[1:]-Decc, stabfrac1.T, levels=[.5], colors='k')
    #ax1.text(.372, .41, '0.5', rotation=-30, verticalalignment='center', horizontalalignment='center')
    #plt.clabel(CS, inline=1, fontsize=12, fmt='%.1f')
    ax1.set_xlabel('$e_c$', labelpad=1)
    ax1.set_ylabel('$e_b$')
    ax1.set_xlim((0,1))
    ax1.set_ylim((0,1))

    # Get analytical predicion
    #ecc1, ecc2 = np.loadtxt('pickles_uncorr_ALLecc_lin/Hillstable_forbubbles.dat').T
    ecc1, ecc2 = compute_pcrit_curve(m1, m2, M, a1, a2)
    ax1.plot(ecc1, ecc2, 'k--', lw=3)

    # Plot MAP ecc from Model 1
    eccval = [.139, .135]
    eccupp = [.305, .197]
    ecclow = [.098, .144]
#    ax1.errorbar(np.array([eccval[0]]), np.array([eccval[1]]), 
#		 xerr=[np.array([ecclow[0]]),np.array([eccupp[0]])], 
#		 yerr=[np.array([ecclow[1]]),np.array([eccupp[1]])], fmt='wo', capsize=5)

    # Plot histogram for c
    ax2 = plt.subplot(gs[0,:-1])
    cols = [rgb2hex(cm(int(f*255))[:3]) for f in [.8,.2]]
    p = ax2.hist([eccs[:,0][stable==1], eccs[:,0][stable==0]], bins=eccarr, stacked=True, 
	     	 label=['Stable','Unstable'], color=cols)
    Nstable, Ntotal, xbins = p[0][0], p[0][1], p[1]
    xbins = xbins[:-1] + np.diff(xbins)[0]/2
    frac_stable = (Nstable) / Ntotal
    #efrac_unstable = np.sqrt(Ntotal-Nstable) / Ntotal
    for i in range(xbins.size):
        ax2.text(xbins[i], Ntotal[i], '%i'%(frac_stable[i]*1e2), fontsize=10,
                 horizontalalignment='center', verticalalignment='bottom')
    ylim2 = ax2.get_ylim()
    ax2.set_yticklabels('')
    ax2.set_xticklabels('')
    ax2.set_ylim((0,ylim2[1]+50))
    ax2.legend(bbox_to_anchor=(1.31,.9), fontsize=13)

    # Plot histogram for b
    ax3 = plt.subplot(gs[1:-1,-1])
    p = ax3.hist([eccs[:,1][stable==1], eccs[:,1][stable==0]], bins=eccarr, stacked=True,
                 label=['Stable','Unstable'], color=cols, orientation='horizontal')
    Nstable, Ntotal, xbins = p[0][0], p[0][1], p[1]
    xbins = xbins[:-1] + np.diff(xbins)[0]/2
    frac_stable = (Nstable) / Ntotal
    #efrac_unstable = np.sqrt(Ntotal-Nstable) / Ntotal
    for i in range(xbins.size):
        ax3.text(Ntotal[i]+10, xbins[i], '%i'%(frac_stable[i]*1e2), fontsize=10,
                 horizontalalignment='left', verticalalignment='center')
    xlim3 = ax3.get_xlim()
    ax3.set_yticklabels('')
    ax3.set_xticklabels('')
    ax3.set_xlim((0,xlim3[1]+90))

    plt.subplots_adjust(bottom=.05, left=.11, top=.97, right=.98, hspace=0, wspace=0)
    if label:
	plt.savefig('plots/ecc_stability.png')
    if pltt:
	plt.show()
    plt.close('all')


#if __name__ == '__main__':
    #folder = 'pickles_uncorr_ALLecc_lin'
    #self = integrate_results(folder)
