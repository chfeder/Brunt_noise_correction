
#from starter1 import *
#from dtools import volavg
from startlaser import *
import volavg
reload(volavg)
from startlaser import *
import yt
import fourier_tools_py3.fourier_filter as Filter

class shell_average():
    def __init__(self,power):
        self.power=power
        ff = Filter.FourierFilter(self.power)
        self.power_1d = np.array([self.power[ff.get_shell(bin)].sum() for bin in range(ff.nx)])
        self.Nzones = np.array([ff.get_shell(bin).sum() for bin in range(ff.nx)])
        self.kd=ff.get_shell_k()
class ft_only():
    def __init__(self,array):
        self.array=array
        self.fft = np.fft.fftn( self.array )
        self.power=self.fft*np.conjugate(self.fft)
        self.power/=self.power.size
        ff = Filter.FourierFilter(self.power)
        self.power_1d = np.array([self.power[ff.get_shell(bin)].sum() for bin in range(ff.nx)])
        self.Nzones = np.array([ff.get_shell(bin).sum() for bin in range(ff.nx)])
        self.kd=ff.get_shell_k()


def get_cubes(sim,frame,do_rho_4=False):
    ds = yt.load("/data/cb1/Projects/P49_EE_BB/%s/DD%04d/data%04d"%(sim, frame, frame))
    print('get cg')
    cg = ds.covering_grid(0, [0.0]*3, [512]*3)

    dds = cg.dds
    rho_full = cg["density"].v
    rho = volavg.volavg(rho_full,rank=3,refine_by=2)
    output = rho_full, rho
    if do_rho_4:
        rho_4 = volavg.volavg(rho,rank=3,refine_by=2)
        output = rho_full, rho, rho_4

    return output

# ADDED.
def get_cubes_cf(void_size, box, pull):
    the_ds = yt.load("/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_%s_bubbles/NIF_hdf5_plt_cnt_0375"%(void_size))
 
    delta = 0.01  
    delta_m = 0.02
    edge_mtop = 0.3675 
    edge_mbot = 0.3275

    if box == 'sx_box':
        the_ledge = the_ds.arr([-delta*pull, edge_mbot, -delta_m],'code_length') 
        the_redge = the_ds.arr([delta*pull, edge_mtop, delta_m],'code_length')  
    if box == 'sz_box':
        the_ledge = the_ds.arr([-delta_m, edge_mbot, -delta*pull],'code_length')
        the_redge = the_ds.arr([delta_m, edge_mtop, delta*pull],'code_length')

    the_center =0.5*the_ledge + 0.5*the_redge
    mindx = the_ds.index.get_smallest_dx()
    size = ((the_redge - the_ledge)/mindx).astype('int') 
    the_level = 0

    cg = the_ds.covering_grid(the_level, the_ledge, size)
    print('getting cg ',pull)

    rho_full = cg["density"].v
    vol_full = cg["cell_volume"].v
    rho = volavg.volavg(rho_full,rank=3,refine_by=1)
    output = rho_full, rho, vol_full

    return output


    
def plot_fft(ftool, outname=None,ax=None):
    savefig=False
    if ax is None:
        savefig=True
        fig,ax=plt.subplots(1,1)
    if ftool.done2:
        ax.plot(ftool.k2d, ftool.power_1d2.real,c='r', label='P2d')
    if ftool.done3:
        ax.plot(ftool.k3d, ftool.power_1d3.real,c='g', label='P3d')

    if ftool.done2 and ftool.done3:
        ax.plot(ftool.k2d,ftool.k2d*ftool.power_1d2.real,c='b', label = 'k P2d')
    if savefig:
        ax.legend(loc=0)

        ax.set(yscale='log',xscale='log')
        fig.savefig(outname)

def plot_fft2(ftool, outname=None):
    fig,ax=plt.subplots(1,3,figsize=(8,4))
    ax[0].imshow( ftool.rho2)
    ax[1].imshow(ftool.rho2p)
    ax[2].plot( ftool.k2d, ftool.power_1d2.real, c='k', label='2d')
    ax[2].plot( ftool.k2dp, ftool.power_1d2p.real, c='r',label='2dp')
    ax[2].set(xscale='log',yscale='log')
    ax[2].legend(loc=0)
    fig.savefig(outname)


def plot_fft3(ftool,apod=None, outname=None):
    fig,ax=plt.subplots(2,2,figsize=(8,4))
    ax[0][0].imshow( ftool.rho2)
    ax[1][0].imshow(ftool.rho2p)
    ax[0][1].imshow( apod.rho2)
    ax[1][1].plot( ftool.k2d, ftool.power_1d2.real, c='k', label='2d')
    ax[1][1].plot( ftool.k2dp, ftool.power_1d2p.real, c='r',label='2dp')
    ax[1][1].plot( apod.k2d, apod.power_1d2.real, c='g', label='apod')
    ax[1][1].set(xscale='log',yscale='log')
    ax[1][1].legend(loc=0)
    fig.savefig(outname)

def plot_brunt(ftool, dim, outname=None, fitrange=None, method='full', ax=None):
    if method=='full':
        sigmas_full(ftool)
    elif method=='norm':
        sigmas_norm(ftool)
    elif method=='range':
        sigmas_range(ftool, fitrange)
    savefig=False
    if ax is None:
        savefig=True
        fig,ax=plt.subplots(1,1)
    if fitrange is None:
        mask = slice(None)
    else:
        mask = slice(fitrange[0],fitrange[1])
    ax.plot(ftool.k2d,          ftool.power_1d2.real,c=[0.5]*3, label='P2d')
    ax.plot(ftool.k2d[mask],    ftool.power_1d2.real[mask],c='r', label='P2d')
    ax.plot(ftool.k3d[mask],    ftool.power_1d3.real[mask],c='g', label='P3d')
    ax.plot(ftool.k2d[mask],    ftool.k2d[mask]*ftool.power_1d2.real[mask],c='b', label = 'k P2d')
    ax.set(xscale='log',yscale='log')
    error = 1-ftool.sigma_Brunt/ftool.sigma_x3d
    text_x=0.05
    text_y=0.3
    dy=0.05
    ax.text(text_x, text_y-1*dy,"%0.2e sigma_x3d"%ftool.sigma_x3d, transform=ax.transAxes)
    ax.text(text_x, text_y-2*dy,"%0.2e sigma_x2d"%ftool.sigma_x2d, transform=ax.transAxes)
    ax.text(text_x, text_y-3*dy,"%0.2e R "%(1./ftool.Rinv), transform=ax.transAxes)
    ax.text(text_x, text_y-4*dy,"%0.2e sigma_B "%ftool.sigma_Brunt,  transform=ax.transAxes)
    ax.text(text_x, text_y-5*dy,"%0.2e  error  "%error,  transform=ax.transAxes)
    ax.text(text_x, text_y-6*dy,"%0.2e  rat  "%(ftool.sigma_Brunt/ftool.sigma_x3d),  transform=ax.transAxes)

    if savefig:
        fig.savefig(outname)


# ADDED.
def sigmas_back(ftool, dim, method='full', fitrange=None):
    if method=='full':
        #ftool.sigmas_full()
        sigmas_full(ftool)
    elif method=='norm':
        sigmas_norm(ftool)
    elif method=='range':
        sigmas_range(ftool, fitrange)

    if dim == 3:
        output = ftool.sigma_k3d
    if dim == 2:
        output = ftool.sigma_x3d
    if dim == 1:
        sigmabrunt = ftool.sigma_Brunt_corr
        sigmax2d = ftool.sigma_x2d
        rinv = ftool.Rinv
        output = sigmabrunt, sigmax2d, rinv   
    return sigmabrunt

def sigmas_full(self):
    #this works.  Not normalized, though.
    self.sigma_x3d = np.sqrt((self.rho**2).sum().real)
    self.sigma_k3d = np.sqrt((self.power_1d3).sum().real)

    self.sigma_x2d = np.sqrt(((self.rho2)**2).sum().real)
    self.sigma_k2d = np.sqrt(self.power_1d2[1:].sum().real)

    self.sigma_k2dk= np.sqrt(( self.k2d*self.power_1d2)[1:].sum())

    self.Rinv = self.sigma_k2dk.real/self.sigma_k2d.real
    self.sigma_Brunt = self.sigma_x2d.real*self.Rinv
    self.R1 =self.sigma_Brunt/self.sigma_x3d
    #just to check that everything works right, do it with the actual 3d power spectrum.
    self.Rinv_actual = self.power_1d3.sum()/self.power_1d2.sum()
    self.sigma_Brunt_actual = self.sigma_x2d*self.Rinv_actual
    #R2 should be 1
    self.R2 = self.sigma_Brunt_actual/self.sigma_x3d

def sigmas_norm(self):
    #works, normalized, only using the fit range
    Nz = self.rho.size
    N2d = self.rho2.size
    self.mean_rho=(self.rho).sum()/Nz
    self.mean_column= self.rho2.sum()/N2d
    #self.sigma_x3d  =np.sqrt( (((self.rho-self.mean_rho)**2)*self.vol).sum().real * 1/self.vol.sum())  #yields the same as below!
    self.sigma_x3d  =np.sqrt( ((self.rho-self.mean_rho)**2).sum().real * 1/Nz)
    self.sigma_k3d  =np.sqrt( (self.power_1d3[1:]).sum().real * 1/Nz )  
    #self.sigma_k3d  =np.sqrt( (self.power_1d3[1:] * 1/Nz).sum().real * 1/Nz )  #yields order of 10-5.

    self.sigma_k2dk =np.sqrt(( self.k2d*self.power_1d2)[1:].sum()/Nz)
    self.sigma_k2dk_corr =np.sqrt(( ((self.k2d/self.k2d[1])*(self.power_1d2/self.power_1d2[3])).sum()/Nz))

    self.sigma_x2d  =np.sqrt(((self.rho2-self.mean_column)**2).sum().real/N2d)
    self.sigma_k2d  =np.sqrt(self.power_1d2[1:].sum().real/N2d)
    self.sigma_k2d_corr  =np.sqrt( (self.power_1d2/self.power_1d2[3]).sum().real/N2d)

    #self.Rinv = ()/((self.power_1d2)[1:].sum()/N2d)
    self.Rinv = (self.sigma_k2dk/self.sigma_k2d).real
    self.Rinv_corr = (self.sigma_k2dk_corr/self.sigma_k2d_corr).real
    #this should give us the right answer.
    #Rinv_actual = (self.power_1d3[1:].sum()/Nz)/(self.power_1d2[1:].sum()/N2d)
    self.Rinv_actual = (self.sigma_k3d/self.sigma_k2d).real

    self.sigma_Brunt = (self.sigma_k2d*self.Rinv).real
    self.sigma_Brunt_corr = (self.sigma_k2d_corr*self.Rinv_corr).real

    self.sigma_Brunt_actual = (self.sigma_k2d*self.Rinv_actual).real

    R1 =self.sigma_Brunt/self.sigma_x3d
    R2 = self.sigma_Brunt_actual/self.sigma_x3d

def sigmas_range(self, fitrange):
    mask = slice(fitrange[0],fitrange[1])
    Nz = self.rho.size
    N2d = self.rho2.size
    self.mean_rho=(self.rho).sum()/Nz
    self.mean_column= self.rho2.sum()/N2d
    self.sigma_x3d  =np.sqrt(((self.rho-self.mean_rho)**2).sum().real/Nz)
    self.sigma_k3d  =np.sqrt((self.power_1d3[mask]).sum().real/Nz)
    self.sigma_k2dk =np.sqrt(( self.k2d*self.power_1d2)[mask].sum()/Nz)
    self.sigma_x2d  =np.sqrt(((self.rho2-self.mean_column)**2).sum().real/N2d)
    self.sigma_k2d  =np.sqrt(self.power_1d2[mask].sum().real/N2d)

    #self.Rinv = ()/((self.power_1d2)[1:].sum()/N2d)
    self.Rinv = (self.sigma_k2dk/self.sigma_k2d).real
    #this should give us the right answer.
    #Rinv_actual = (self.power_1d3[1:].sum()/Nz)/(self.power_1d2[1:].sum()/N2d)
    self.Rinv_actual = (self.sigma_k3d/self.sigma_k2d).real
    self.sigma_Brunt = (self.sigma_x2d*self.Rinv).real
    sigma_Brunt_actual = (self.sigma_k2d*self.Rinv_actual).real

    R1 =self.sigma_Brunt/self.sigma_x3d
    R2 = self.sigma_Brunt_actual/self.sigma_x3d


if 0:
    N2d = ftool.rho2.size
    mean_column = ftool.rho2.sum()/N2d
    sigma_col = ((ftool.rho2-mean_column)**2).sum().real/N2d
    sigma_2d_fft = ftool.power_1d2[mask].sum().real/N2d
    #print("col %0.2e fft %0.2e ratio %0.2e"%(sigma_col, sigma_2d_fft, sigma_col/sigma_2d_fft))
    Rinv = (( ftool.k2d*ftool.power_1d2)[mask].sum()/Nz)/((ftool.power_1d2)[mask].sum()/N2d)
    Rinv_actual = (ftool.power_1d3[mask].sum()/Nz)/(ftool.power_1d2[mask].sum()/N2d)
    sigma_Brunt = (sigma_col*Rinv).real
    sigma_Brunt_actual = (sigma_col*Rinv_actual).real
    #print("Brunt", sigma_Brunt/sigma_rho)
    #print("One", sigma_Brunt_actual/sigma_rho)

    R1 =sigma_Brunt/sigma_rho
    R2 = sigma_Brunt_actual/sigma_rho


class fft_tool():
    def __init__(self,rho,vol=None):
        self.rho=rho
        self.vol=vol
        self.rho2=None
        self.rho2p=None
        self.done2=False
        self.done3=False


    def do3(self):
        #print('3d fourier transform')
        self.fft3 = np.fft.fftn( self.rho )
        self.power=self.fft3*np.conjugate(self.fft3)
        self.power/=self.power.size
        ff = Filter.FourierFilter(self.power)
        self.power_1d3 = np.array([self.power[ff.get_shell(bin)].sum() for bin in range(ff.nx)])
        #self.power_1d3 /= self.rho.size
        self.Nzones = np.array([ff.get_shell(bin).sum() for bin in range(ff.nx)])
        self.k3d=ff.get_shell_k()
        self.done3=True

    def do2(self,projax=0):
        print('2d fourier transform')
        if self.rho2 is None:
            print('MAKE NEW PROJECTION')
            self.rho2=self.rho.sum(axis=projax)
            Nz = self.rho.shape[projax]
            #self.rho2/=Nz   
        '''DEBUG
        plt.clf()
        plt.imshow(self.rho2)
        plt.savefig('density_s2')
        '''
        self.fft2 = np.fft.fftn( self.rho2 )
        self.power2=self.fft2*np.conjugate(self.fft2)
        self.power2/=self.power2.size
        '''DEBUG
        plt.clf()
        plt.imshow(self.power2.real)
        plt.savefig('density_fourier_s2')
        '''
        ff2 = Filter.FourierFilter(self.power2)
        self.power_1d2 = np.array([self.power2[ff2.get_shell(bin)].sum() for bin in range(ff2.nx)])
        #self.power_1d2 /= self.rho2.size
        self.Nzones2 = np.array([ff2.get_shell(bin).sum() for bin in range(ff2.nx)])
        self.k2d=ff2.get_shell_k()
        self.done2=True
        
    def do2_periodic(self,projax=0):
        print('2d fourier transform')
        if self.rho2p is None:
            print('MAKE NEW PROJECTION')
            self.rho2=self.rho.sum(axis=projax)
            Nz = self.rho.shape[projax]
            self.rho2/=Nz
            #self.rho2p = np.concatenate([self.rho2, self.rho2[::-1]], axis=0)
            #self.rho2p = np.concatenate([self.rho2p, self.rho2p[:,::-1]],axis=1)
            self.rho2p = np.concatenate([self.rho2, self.rho2[::-1]], axis=0)
            self.rho2p = np.concatenate([self.rho2p, self.rho2p[:,::-1]],axis=1)
            
            fig,ax=plt.subplots(1,1)
            ax.imshow(self.rho2p)
            fig.savefig('/home/dccollins/PigPen/derp')


        self.fft2p = np.fft.fftn( self.rho2p )
        self.power2p=self.fft2p*np.conjugate(self.fft2p)
        #self.power2/=self.power2.size
        ff2p = Filter.FourierFilter(self.power2p)
        self.power_1d2p = np.array([self.power2p[ff2p.get_shell(bin)].sum() for bin in range(ff2p.nx)])
        self.power_1d2p /= self.rho2p.size
        self.Nzones2p = np.array([ff2p.get_shell(bin).sum() for bin in range(ff2p.nx)])
        self.k2dp=ff2p.get_shell_k()
        self.done2=True

    def apodize1(self,projax=0):
        self.rho2=self.rho.sum(axis=projax)
        self.rho2/=self.rho2.shape[projax]
        shape=np.array(self.rho2.shape)
        baseshape=(1.*shape).astype('int')
        base = np.zeros(baseshape)
        start = baseshape//2-shape//2

        base[start[0]:(start[0]+shape[0]), start[1]:(start[1]+shape[1])] = self.rho2

        if 0:
            #things that don't quite work
            self.rho2 = base
            from scipy.signal import general_gaussian
            from scipy.signal import convolve2d
            window = np.outer(general_gaussian(baseshape[0],6,3),general_gaussian(baseshape[0],6,3))
            #np.roll(window,baseshape[0]//2,axis=0)
            #np.roll(window,baseshape[0]//2,axis=1)
        if 0:
            #actually work ok
            window = np.zeros(baseshape)
            window[0:3,0:3]=1
        if 1:
            #works pretty well.
            x = np.arange(baseshape[0])
            sigma_conv=2
            g = np.exp(-x**2/(2*sigma_conv**2))**6
            window = np.outer(g,g)

        window/=window.sum()
        self.window=window
        
        #self.rho2 = scipy.convolve(base, window)
        #self.rho2 = convolve2d(base, window)
        #convolve the window function with the base
        a = np.fft.fftn(base)
        b = np.fft.fftn(window)
        c = a*b
        self.rho2 = np.fft.ifftn(c)
        q=np.abs(self.rho2.imag).sum()/self.rho2.imag.size
        if q>1e-13:
            print("!!!!!!!!!!!!!!!!!Imaginary",q)
        self.rho2=self.rho2.real

