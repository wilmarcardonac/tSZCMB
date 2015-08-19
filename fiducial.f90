Module fiducial

Implicit none

save 

! Same cosmology as in published version of "Detection of thermal SZ-CMB lensing cross-correlation in Planck  
! nominal mission data" by Hill and Spergel

Real*8,parameter :: c = 2.99792458d8    ! speed of light in m s-1
Real*8,parameter :: G = 6.67384d-11    ! Newtonian constant of gravitation in m3 kg-1 s-2
Real*8,parameter :: M_sun = 1.989d30    ! solar mass in kg
Real*8,parameter :: Mpc = 3.08567758d22    ! 1 Mpc  = 3.08567758d22 m 
Real*8,parameter :: m_e = 9.10938291d-31    ! electron mass in Kg
Real*8,parameter :: sigma_e = 6.652458734d-29    ! Thomson scattering cross section in m^2
Real*8,parameter :: h = 69.7d-2    ! Dimensionless 
Real*8,parameter :: H_0 = 69.7d3    ! Hubble parameter at present time in m s-1 Mpc-1
Real*8,parameter :: Omega_m0 = 0.282d0    ! Matter parameter density at present time
Real*8,parameter :: Omega_L0 = 0.7181d0    ! 1 - Omega_m0 : Dark energy parameter density (LCDM model) at present time in a flat universe
Real*8,parameter :: Omega_b_h2 = 0.02240d0    ! Baryon density today 
Real*8,parameter :: ns = 0.9646d0    ! scalar spectrum power-law index 
Real*8,parameter :: z_dec = 1100d0    ! decoupling red-shift
!Real*8,parameter :: A_s = 2.43d-9    ! dimensionless curvature power spectrum at k_0 = 0.05 Mpc-1
Real*8,parameter :: k_0 = 0.05    ! Mpc-1
Real*8,parameter :: kmin = 1.d-4    ! H_0*1.e-2 # minimum wavenumber to compute angular power spectrum of lensing potential
Real*8,parameter :: kmax = 1.d5    ! #*H_0*1.e-2 # maximum wavenumber to compute angular power spectrum of lensing potential 
Real*8,parameter :: zmax = 1.d1    ! upper red-shift fiducial value in Hill and Spergel paper 1312.4525
Real*8,parameter :: zmin = 5.d-3    ! lower red-shift fiducial value in Hill and Spergel paper 1312.4525 
Real*8,parameter :: sigma8 = 0.817d0     !#0.8285 # RMS matter fluctuations today in linear theory 
Real*8,parameter :: Pi = 3.141592653589793d0
Real*8,parameter :: Mmin = 1.d5/h    ! lower mass fiducial value in Hill and Spergel paper 1312.4525
Real*8,parameter :: Mmax = 5.d15/h    ! upper mass fiducial value in Hill and Spergel paper 1312.4525
Integer*4,parameter :: number_of_k = 1d3    ! size of wavevector array
Integer*4,parameter :: number_of_z = 4d2    ! size of red-shift array
Integer*4,parameter :: number_of_M = 1d3    ! size of virial mass array
Integer*4,parameter :: number_of_l = 2d1    ! size of l array
Integer*4,parameter :: lmax = 1d4    ! highest multipole
Integer*4,parameter :: lmin = 1d1    ! lowest multipole 
Real*8,parameter :: tcmb0 = 2.728d0    ! CMB temperature today in Kelvin
Real*8 :: Normalization    ! Normalization constant for matter power spectrum (Equation (A7) in Eisenstein and Hu)
Real*8,parameter :: DeltaSO = 2.d2    ! Parameter to define Spherical Overdensity   

End Module fiducial
