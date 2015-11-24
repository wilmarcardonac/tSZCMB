import math
import numpy as np
import pylab as py 

#l,cl1h,cl2h,cl = np.loadtxt('Cl_yphi.dat',unpack=True,usecols=[0,1,2,3])
l,clyp1h,clyp2h,clyp,cl1h,cl2h,clpp,cllimber = np.loadtxt('Cl_yphi.dat',unpack=True,usecols=[0,1,2,3,4,5,6,7])
h1 = np.loadtxt('1-halo.dat')
h2 = np.loadtxt('2-halo.dat')
total = np.loadtxt('Total.dat')

#py.plot(h1[:,0],h1[:,1]*1.e-11,label='H&S 1-halo')
#py.plot(h2[:,0],h2[:,1]*1.e-11,label='H&S 2-halo')
#py.plot(total[:,0],total[:,1]*1.e-11,label='H&S total')

#py.xscale('log')

#py.yscale('linear')

#py.xlim(1.e1,1.e4)

#py.ylim(0.e-11,4.e-11)

#py.show()
#exit()
#fig1 = py.figure()


#py.loglog(l,l**4*cl1h,label='lensing (1-halo)')

#py.loglog(l,l**4*cl2h,label='lensing (2-halo)')

#py.loglog(l,l**4*cllimber,label='lensing (limber approx.)')

#py.legend(loc=0)

#py.show()
#exit()

#py.loglog(l,clpp,label='lensing (1-halo + 2-halo)')

#py.loglog(l,clyp,label='y-phi (1-halo + 2-halo)')

#py.xlabel(r'$\ell$')

#py.xlim(1.e1,1.e4)

#py.ylim(1.e-25,1.e-9)

#py.ylabel(r'$\ell^2(\ell+1)C^{y\phi}_\ell/2\pi$')

#py.ylabel(r'$ C^{\phi\phi}_{\ell},\, C^{\psi}_{\ell}$')

#py.ylim(1.e5,4.e7)

#py.legend()

#py.savefig('Cl_yphi.pdf')

#py.savefig('Cl_phiphi.pdf')

#py.close(fig1)
#py.show()

#py.plot(l,l**2*(l+1)*clyp1h/2./math.pi,label='fiducial (1-halo)')

#py.plot(l,l**2*(l+1)*clyp2h/2./math.pi,label='fiducial (2-halo)')

#py.plot(l,l**2*(l+1)*clyp/2./math.pi,label='fiducial (total)')

#py.xscale('log')

#py.yscale('linear')

#py.xlim(1.e1,1.e4)

#py.ylim(0.e-11,4.e-11)

#py.legend()

#py.show()

py.loglog(l,abs(cllimber-cl2h)/cllimber*100,label='% difference')

py.hlines(5.,1.,1.e4,color='k',linestyles='dotted')

py.show()

exit()
