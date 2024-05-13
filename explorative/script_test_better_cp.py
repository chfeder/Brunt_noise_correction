
from startlaser import *
import yt
import fourier_tools_py3.fourier_filter as Filter
import brunt_tools_cp as bt
reload(bt)
import cfpack as cfp
from cfpack import stop,print 
plt.close('all')

void_size = '12p5mu'
box_vol = 'sx_box'
proj_axis = 0
the_y_3d = []
the_y_b = []
emptylist_b = []
d3dk = 3
d3d = 2
d2d = 1

if 1:
    for pull_num in range(9):  #DEBUG!!
        pull_num = pull_num + 1

        #get cubes; rho_full is straight off disk.  rho is downsampled by 2
        rho_full, rho, vol = bt.get_cubes_cf(void_size, box_vol, pull_num)
        ftool = bt.fft_tool(rho,vol)
        ftool.do3()
        ftool.do2(projax=proj_axis)
        
        #sig_3Dk  = bt.sigmas_back(ftool, d3dk, method='norm') #checked.
        #sig_3D = bt.sigmas_back(ftool, d3d, method='norm') #checked.
        sig_B = bt.sigmas_back(ftool, d2d, method='norm')  #checked.
        #print('sigmas!! 3dk, 3d, B',sig_3Dk, sig_3D, sig_B)

        #the_y_3d = np.append(the_y_3d, sig_3D)
        the_y_b = np.append(the_y_b, sig_B)
if 1:
    the_d = np.linspace(0.02, 0.18, pull_num, endpoint=True)  
    the_l = 0.02 * 2
    the_x = the_d/the_l  #this is the aspect ratio 

    fig,ax=plt.subplots(1,1)
    #ax.scatter(the_x, the_y_3d, color='orange')
    ax.scatter(the_x, the_y_b, color='blue')
    fig.savefig('%s_%s'%(box_vol,void_size))



if 0:
    bt.plot_fft(ftool,"test_fft") #%s/ffts_%s_%04d"%(plotdir,sim,frame))

if 0:
    #outname="%s/fft_full_%s_n%04d"%(plotdir,sim,frame)
    outname= 'test_brunta'
    bt.plot_brunt(ftool, outname=outname,method='full')

    #outname="%s/fft_norm_%s_n%04d"%(plotdir,sim,frame)
    outname= 'test_bruntb'
    bt.plot_brunt(ftool, outname=outname,method='norm')

    #outname="%s/fft_range_%s_n%04d"%(plotdir,sim,frame)
    outname= 'test_bruntc'
    bt.plot_brunt(ftool, outname=outname,method='range', fitrange=[4,25])


if 0:
    #try brunt.  Works to a factor of 2, not bad.
    sigma_rho = (ftool.rho**2).sum().real
    sigma_fft = (ftool.power_1d3).sum().real
    #sigma_fft = (ftool.power).sum().real/ftool.power.size**2
    print("rho %0.2e fft %0.2e ratio %0.2e"%(sigma_rho,sigma_fft,(sigma_fft)/sigma_rho))

    sigma_col = ((ftool.rho2)**2).sum().real
    sigma_2d_fft = ftool.power_1d2.sum().real
    print("col %0.2e fft %0.2e ratio %0.2e"%(sigma_col, sigma_2d_fft, sigma_col/sigma_2d_fft))
    Rinv = ( ftool.k2d*ftool.power_1d2)[1:].sum()/(ftool.power_1d2)[1:].sum()
    Rinv_actual = ftool.power_1d3.sum()/ftool.power_1d2.sum()
    sigma_Brunt = sigma_col*Rinv
    sigma_Brunt_actual = sigma_col*Rinv_actual
    print( sigma_Brunt/sigma_rho)
    print( sigma_Brunt_actual/sigma_rho)
    fig,ax=plt.subplots(1,1)
    Q=ftool.k2d*ftool.power_1d2
    P=ftool.power_1d3
    ax.plot( ftool.k2d, Q)
    ax.plot( ftool.k3d, P)
    ax.set(xscale='log',yscale='log')
    fig.savefig("%s/test"%plotdir)
