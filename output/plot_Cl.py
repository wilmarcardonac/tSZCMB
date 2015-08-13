import math
import numpy as np
import pylab as py 

#l,cl1h,cl2h,cl = np.loadtxt('Cl_yphi.dat',unpack=True,usecols=[0,1,2,3])
l,clyp1h,clyp2h,clyp,cl1h,cl2h,cllimber = np.loadtxt('Cl_yphi.dat',unpack=True,usecols=[0,1,2,3,4,5,7])

fig1 = py.figure()


#py.loglog(l,l**4*cl1h,label='lensing (1-halo)')

py.loglog(l,l**4*cl2h,label='lensing (2-halo)')

py.loglog(l,l**4*cllimber,label='lensing (limber approx.)')

#py.xlabel(r'$\ell$')

#py.ylabel(r'$\ell^2(\ell+1)C^{y\phi}_\ell/2\pi$')

#py.ylabel(r'$ C^{\phi\phi}_{\ell},\, C^{\psi}_{\ell}$')

#py.ylim(1.e5,4.e7)

py.legend()

#py.savefig('Cl_yphi.pdf')

#py.savefig('Cl_phiphi.pdf')

#py.close(fig1)
py.show()

py.loglog(l,l**2*(l+1)*clyp1h/2./math.pi,label='fiducial (1-halo)')

py.loglog(l,l**2*(l+1)*clyp2h/2./math.pi,label='fiducial (2-halo)')

py.loglog(l,l**2*(l+1)*clyp/2./math.pi,label='fiducial (total)')

py.legend()

py.show()

py.loglog(l,abs(cllimber-cl2h)/cllimber*100,label='% difference')

py.show()

exit()
