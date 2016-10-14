import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as py
import math
import numpy as np
#import pylab as py 
#l,clyp2h = np.loadtxt('Cl_yphi.dat',unpack=True,usecols=[0,2])
#z,k = np.loadtxt('kernel.dat',unpack=True,usecols=[0,1])
#l,clyp2h,cl1h,cl2h,cl,cllimber = np.loadtxt('Cl_yphi.dat',unpack=True,usecols=[0,2,4,5,6,7])
l,clyp1h,clyp2h,clyp,cl1h,cl2h,clpp,cllimber,clyy1h,clyy2h,clyy = np.loadtxt('Cl_yphi.dat',unpack=True,usecols=[0,1,2,3,4,5,6,7,8,9,10])
h1 = np.loadtxt('1-halo.dat')
h2 = np.loadtxt('2-halo.dat')
total = np.loadtxt('Total.dat')

py.loglog(l,clyp/np.sqrt(clyy*clpp),label='fig 2')

#print l[33],clyp1h[33]*l[33]**2*(l[33]+1.)/2./math.pi
#exit()
py.xscale('log')

py.yscale('linear')

py.legend(loc=0)

py.xlim(1.e1,1.e4)

py.ylim(0.,1.)

py.savefig('figure_2.pdf')


#py.show()
exit()

fig1 = py.figure()

py.loglog(l,l**4*cl1h,label='lensing (1-halo)')

py.loglog(l,l**4*cl2h,label='lensing (2-halo)')

py.loglog(l,l**4*cllimber,label='lensing (limber approx.)')

py.xlabel(r'$\ell$')

py.ylabel(r'$\ell^4C_\ell$')

py.legend(loc=0)

py.savefig('Cl_phiphi.pdf')

py.close(fig1)

fig1 = py.figure()

py.loglog(l,abs(cllimber-cl2h)/cllimber*100,label='% difference Limber and two halo')

py.hlines(5.,1.,1.e4,color='k',linestyles='dotted')

py.vlines(900.,0.1,1.e2,color='k',linestyles='dotted')

py.xlabel(r'$\ell$')

py.ylabel(r'$|C^{Limber}_\ell-Cl^{2h}_\ell|/C^{Limber}_\ell$')

py.legend(loc=0)

py.savefig('difference.pdf')

py.close(fig1)

fig1 = py.figure()

py.plot(l,l**2*(l+1)*clyp1h/2./math.pi*1.43,label='fiducial (1-halo)')

py.plot(l,l**2*(l+1)*clyp2h/2./math.pi*1.43,label='fiducial (2-halo)')

py.plot(l,l**2*(l+1)*clyp/2./math.pi*1.43,label='fiducial (total)')

py.plot(h1[:,0],h1[:,1]*1.e-11,label='H&S 1-halo')

py.plot(h2[:,0],h2[:,1]*1.e-11,label='H&S 2-halo')

py.plot(total[:,0],total[:,1]*1.e-11,label='H&S total')

py.ylabel(r'$\ell^2(\ell+1)C^{y\phi}_\ell/2\pi$')

py.xlabel(r'$\ell$')

py.xscale('log')

py.yscale('linear')

py.xlim(1.e1,1.e4)

py.ylim(0.e-11,4.e-11)

py.legend(loc=0)

py.savefig('Cl_yphi.pdf')

py.close(fig1)

fig1 = py.figure()

py.loglog(l,clpp,label='lensing (1-halo + 2-halo)')

py.loglog(l,clyp,label='y-phi (1-halo + 2-halo)')

py.loglog(l,clyy,label='y-y (1-halo + 2-halo)')

py.xlabel(r'$\ell$')

py.ylabel(r'$C_\ell$')

py.xlim(1.e1,1.e4)

py.ylim(1.e-25,1.e-9)

py.legend(loc=0)

py.savefig('Cl.pdf')

py.close(fig1)

exit()

