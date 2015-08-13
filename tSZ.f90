Program tSZ
!#############################################
! Here we load all the modules we need
!#############################################

use fiducial
use arrays
use functions 
use omp_lib

!################################
! We declare variables to be used
!################################

Implicit none
Integer*4 :: index1 ! size of arrays
Real*8 :: wtime
Character(len=15) :: halo_definition

!############
! Assignments
!############

halo_definition = 'virial'    ! halo definition that we are using 

!####################################################################
! Allocate memory : red-shift, multipoles, virial mass,one halo term, 
! two halo term, full Cl's, halo mass function, form factor function,
! lensing potential function, comoving volume per steradian function,
! and linear halo bias function.
!####################################################################

allocate (z(1:number_of_z), M(1:number_of_M), ml(1:number_of_l),k(1:number_of_k),&
Cl1h(1:number_of_l),Cl2h(1:number_of_l),Cl(1:number_of_l),M200d(1:number_of_M,1:number_of_z),&
M200c(1:number_of_M,1:number_of_z),r200c(1:number_of_M,1:number_of_z),&
r200d(1:number_of_M,1:number_of_z),Clphiphi1h(1:number_of_l),Clphiphi2h(1:number_of_l),&
Clphiphi(1:number_of_l),alpha_halo_mass_function(1:number_of_z),&
Clpsilimber(1:number_of_l),stat = status1)

allocate (dndM(1:number_of_M,1:number_of_z),ylMz(1:number_of_l,1:number_of_M,1:number_of_z),&
philMz(1:number_of_l,1:number_of_M,1:number_of_z),d2VdzdO(1:number_of_z),&
bMz(1:number_of_M,1:number_of_z),mbz(1:number_of_z),Scrit(1:number_of_z),stat = status2)

!########################################################
! Filling arrays of red-shift, virial mass and multipoles
!########################################################

! Wavevector array. Units : 1/Mpc
Do index1 = 1, number_of_k    
    k(index1) = 10**(log10(kmin) + real(index1-1)*(log10(kmax) - log10(kmin))/real(number_of_k-1))
End Do

! Red-shift array. Dimensionless.
Do index1 = 1, number_of_z      
    z(index1) = 10**(log10(zmin) + real(index1-1)*(log10(zmax) - log10(zmin))/real(number_of_z-1))
End Do

! Mass array.  Units : Solar mass 
! This will be a virial mass array. Hill and Spergel (1312.4525) claim to use Spherical Overdensity (SO) mass M_{200d}
! (mean background density definition) even though their functions are written in terms of virial mass.
Do index1 = 1, number_of_M    
    M(index1) = 10**(log10(Mmin) + real(index1-1)*(log10(Mmax) - log10(Mmin))/real(number_of_M-1))
End Do

! Multipole array. Dimensionless
Do index1 = 1, number_of_l     
    ml(index1) = int(10**(log10(dble(lmin)) + real(index1-1)*(log10(dble(lmax)) - &
    log10(dble(lmin)))/real(number_of_l-1)),4)
End Do

!####################
! Computing functions 
!####################

!call compute_M_delta_c_from_M_and_z(DeltaSO)

call read_M200dc_r200dc()

wtime = omp_get_wtime()    ! setting starting time 

call compute_normalization()    ! Normalising matter power spectrum with fiducial sigma_8 

!print *, r_delta_d(200.d0,z(1),M(100))/concentration_mass_mean(M(100),z(1))
!stop

!print *, lensing_potential(1,1,1,halo_definition)

!call compute_matter_power_spectrum_at_z(1)
!call compute_transfer_function()

!stop

!call compute_bMz()    ! computing linear halo bias array as a function of mass and red-shift
call read_bMz()

call compute_alpha_halo_mass_function()    ! computing alpha constant in halo mass function in order to fulfill condition 
                                           ! that mean bias of all matter at a fixed red-shift is unity
!call compute_dndM()   ! Computing halo mass function array as a function of mass and red-shift 
call read_dndM()    

call compute_mean_bias_matter()    ! Verify that mean bias of all matter is unity for red-shift array 

call compute_d2VdzdO()   ! Computing comoving volume element per steradian array as a function of red-shift  

call compute_critical_surface_density()    ! Computing critical surface density as a function of red-shift

!call compute_lensing_potential(halo_definition)    ! Computing lensing potential as a function of mass and red-shift
call read_philMz()

call compute_Clphiphi1h()    ! Compute lensing potential angular power spectrum (one halo term)

call compute_Clphiphi2h()    ! Compute lensing potential angular power spectrum (two halo term)

call compute_Clpsilimber()    ! Compute angular power spectrum of lensing potential in the Limber approximation

call compute_form_factor()    ! Compute form factor as a function of mass and red-shift
!call read_ylMz()

call compute_Cl1h()

call compute_Cl2h()
 
call compute_Cl()    ! Compute total angular power spectrum 

call write_Cl()    ! Write angular power spectrum into a text file

print *,omp_get_wtime()-wtime

deallocate (z,M,k,ml,Cl1h,Cl2h,Clphiphi1h,Clphiphi2h,Cl,Clphiphi,&
d2VdzdO,dndM,ylMz,philMz,bMz,mbz,M200c,M200d,r200c,r200d,Scrit,alpha_halo_mass_function,&
Clpsilimber)

End Program tSZ




