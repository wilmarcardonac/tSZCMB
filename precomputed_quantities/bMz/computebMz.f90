Program computebMz
!#############################################
!#############################################
! Here we load all the modules we need
!#############################################
!#############################################

use arrays
use fiducial
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
! two halo term, full Cl's, halo mass function, form factor function,
! lensing potential function, comoving volume per steradian function,
! and linear halo bias function.
!####################################################################

allocate (z(1:number_of_z), M(1:number_of_M), ml(1:number_of_l),&
Cl1h(1:number_of_l),Cl2h(1:number_of_l),Cl(1:number_of_l),stat = status1)

allocate (dndM(1:number_of_M,1:number_of_z),ylMz(1:number_of_l,1:number_of_M,1:number_of_z),&
philMz(1:number_of_l,1:number_of_M,1:number_of_z),d2VdzdO(1:number_of_z),&
bMz(1:number_of_M,1:number_of_z),stat = status2)

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

!####################
! Computing functions 
!####################

call compute_bMz()

deallocate (z,M,ml,Cl1h,Cl2h,Cl,d2VdzdO,dndM,ylMz,philMz,bMz)

End Program computebMz




