Program computephilMz
!#############################################
!#############################################
! Here we load all the modules we need
!#############################################
!#############################################
use fiducial
use arrays
use functions 
use omp_lib

!################################
!################################
! We declare variables to be used
!################################
!################################

Implicit none
Integer*4 :: index1 ! size of arrays
Real*8 :: wtime

!############
!############
! Assignments
!############
!############


!####################################################################
! Allocate memory : red-shift, multipoles, virial mass,one halo term, 
! two halo term and full Cl's
!####################################################################

allocate (z(1:number_of_z), M(1:number_of_M), ml(1:number_of_l),&
philMz(1:number_of_l,1:number_of_M,1:number_of_z),M200c(1:number_of_M,1:number_of_z),&
M200d(1:number_of_M,1:number_of_z),r200c(1:number_of_M,1:number_of_z),&
r200d(1:number_of_M,1:number_of_z),stat = status1)

!########################################################
! Filling arrays of red-shift, virial mass and multipoles
!########################################################

Do index1 = 1, number_of_z 
    z(index1) = 10**(log10(zmin) + real(index1-1)*(log10(zmax) - log10(zmin))/real(number_of_z-1))
End Do

Do index1 = 1, number_of_M
    M(index1) = 10**(log10(Mmin) + real(index1-1)*(log10(Mmax) - log10(Mmin))/real(number_of_M-1))/h
End Do

Do index1 = 1, number_of_l
    ml(index1) = int(10**(log10(dble(lmin)) + real(index1-1)*(log10(dble(lmax)) - &
    log10(dble(lmin)))/real(number_of_l-1)),4)
End Do

call read_M200dc_r200dc()

!######################
! Computing form factor
!######################

call compute_lensing_potential()

deallocate (z,M,ml,philMz,M200c,M200d,r200c,r200d)

End Program computephilMz




