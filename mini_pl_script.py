

x=np.linspace(0,50,51)
yone = x**1
yhalf = x**1.5
yonenhalf = x**1.5
ytwo = x**2 
ythree = x**3 
cfp.plot(ythree, x, label='three')
cfp.plot(powspec2D_noise, ks_n, xlog=True, ylog=True, ylim=(1e0, 1e2), save='test_4_three.png')