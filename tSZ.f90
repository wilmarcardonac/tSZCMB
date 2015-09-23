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
Integer*4 :: index1,q ! size of arrays and auxiliar counter
Real*8 :: wtime
Character(len=15) :: halo_definition

!############
! Assignments
!############

halo_definition = 'virial'    ! halo definition that we are using 

open(20,file='./output/execution_information.txt')

!####################################################################
! Allocate memory : red-shift, multipoles, virial mass,one halo term, 
! two halo term, full Cl's, halo mass function, form factor function,
! lensing potential function, comoving volume per steradian function,
! and linear halo bias function.
!####################################################################

allocate (z(1:number_of_z), M(-1:number_of_M+2), ml(1:number_of_l),k(1:number_of_k),&
Cl1h(1:number_of_l),Cl2h(1:number_of_l),Cl(1:number_of_l),M200d(-1:number_of_M+2,1:number_of_z),&
M200c(-1:number_of_M+2,1:number_of_z),r200c(-1:number_of_M+2,1:number_of_z),dM200ddM(1:number_of_M,1:number_of_z),&
r200d(-1:number_of_M+2,1:number_of_z),Clphiphi1h(1:number_of_l),Clphiphi2h(1:number_of_l),&
Clphiphi(1:number_of_l),alpha_halo_mass_function(1:number_of_z),dM200cdM(1:number_of_M,1:number_of_z),&
Clpsilimber(1:number_of_l),stat = status1)

allocate (dndM(1:number_of_M,1:number_of_z),ylMz(1:number_of_l,1:number_of_M,1:number_of_z),&
philMz(1:number_of_l,1:number_of_M,1:number_of_z),d2VdzdO(1:number_of_z),&
bMz(1:number_of_M,1:number_of_z),mbz(1:number_of_z),Scrit(1:number_of_z),&
sigma_square_M200d(1:number_of_M,1:number_of_z),dsigma_square_M200d(1:number_of_M,1:number_of_z),stat = status2)

If ((status1 .eq. 0) .and. (status2 .eq. 0)) then

    write(20,*) 'Memory allocated succesfully '

End If

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

! Auxiliar index for halo mass function (red-shift greater than 3)
Do q=1,number_of_z

    If (z(q) .gt. 3.d0) then

        indexz_halo_mass_function = q - 1

        exit

    End If

End Do

! Mass array.  Units : Solar mass 
! This will be a virial mass array. Hill and Spergel (1312.4525) claim to use Spherical Overdensity (SO) mass M_{200d}
! (mean background density definition) even though their functions are written in terms of virial mass.
Do index1 = -1, number_of_M+2    
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

wtime = omp_get_wtime()    ! setting starting time 

!write(20,*) 'Executing mass conversion '

!call compute_M_delta_c_from_M_and_z(DeltaSO)

!write(20,*) 'Mass conversion ended '

!print *,omp_get_wtime()-wtime

!stop

call read_M200dc_r200dc()
!print *,dM200ddM(1,1),dM200ddM(number_of_M,number_of_z),M200c(1,1),M200d(2,2)
!stop
call compute_normalization()    ! Normalising matter power spectrum with fiducial sigma_8 
!wtime = omp_get_wtime()    ! setting starting time 
!print *,sigma_squared(M(number_of_M),number_of_z)
!print *,omp_get_wtime()-wtime
!stop
!wtime = omp_get_wtime()    ! setting starting time 
call compute_sigma_square_M200d()
!print *,omp_get_wtime()-wtime
!call read_sigma_square_M200d() 
!stop
call compute_bMz()    ! computing linear halo bias array as a function of mass and red-shift
!call read_bMz()

!wtime = omp_get_wtime()    ! setting starting time 
!print *,indexz_halo_mass_function
!print *,nonnormalised_halo_mass_function(1,indexz_halo_mass_function+1)
!print *,omp_get_wtime()-wtime
!stop
call compute_alpha_halo_mass_function()    ! computing alpha constant in halo mass function in order to fulfill condition 
                                           ! that mean bias of all matter at a fixed red-shift is unity
!print *,nonnormalised_halo_mass_function(number_of_M,indexz_halo_mass_function+1)
!print *,nonnormalised_halo_mass_function(number_of_M,indexz_halo_mass_function+10)
!print *,halo_mass_function(number_of_M,indexz_halo_mass_function+1)
!print *,halo_mass_function(number_of_M,indexz_halo_mass_function+10)
!stop
!write(20,*) 'Computing halo mass function '
call compute_dndM()   ! Computing halo mass function array as a function of mass and red-shift 
!call read_dndM()    
!call write_dndM_at_z(number_of_z-20)
!stop
!write(20,*) 'Computing mean bias of all matter '
!call compute_mean_bias_matter()    ! Verify that mean bias of all matter is unity for red-shift array 

call compute_d2VdzdO()   ! Computing comoving volume element per steradian array as a function of red-shift  

call compute_critical_surface_density()    ! Computing critical surface density as a function of red-shift

!wtime = omp_get_wtime()    ! setting starting time 
!print *,lensing_potential(1,1,10,halo_definition)
!print *,omp_get_wtime()-wtime
!stop

write(20,*) 'Computing lensing potential '
call compute_lensing_potential(halo_definition)    ! Computing lensing potential as a function of mass and red-shift
!call read_philMz()
!print *, pre_Clphiphi(number_of_z-30,number_of_l-2)
!print *,C_l_phiphi_one_halo(5)
!stop

call compute_Clphiphi1h()    ! Compute lensing potential angular power spectrum (one halo term)

call compute_Clphiphi2h()    ! Compute lensing potential angular power spectrum (two halo term)

call compute_Clpsilimber()    ! Compute angular power spectrum of lensing potential in the Limber approximation

write(20,*) 'Computing form factor '
!wtime = omp_get_wtime()
!print *,form_factor(1,1,1)
!print *,omp_get_wtime()-wtime
!stop
call compute_form_factor()    ! Compute form factor as a function of mass and red-shift
!call read_ylMz()
!print *,C_l_yphi_two_halo(number_of_l-2)
!stop
call compute_Cl1h()

call compute_Cl2h()
 
call compute_Cl()    ! Compute total angular power spectrum 

call write_Cl()    ! Write angular power spectrum into a text file

write(20,*) 'Angular power spectra have been written out '

!print *,omp_get_wtime()-wtime

deallocate (z,M,k,ml,Cl1h,Cl2h,Clphiphi1h,Clphiphi2h,Cl,Clphiphi,&
d2VdzdO,dndM,ylMz,philMz,bMz,mbz,M200c,M200d,r200c,r200d,Scrit,alpha_halo_mass_function,&
Clpsilimber)

End Program tSZ




