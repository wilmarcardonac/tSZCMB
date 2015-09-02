module arrays
    Integer :: status1,status2,status3,status4,status5,status6
    Real*8, allocatable, dimension(:) :: z,M,Cl1h,Cl2h,Cl,d2VdzdO,mbz,k,Scrit,Clphiphi1h,Clphiphi2h,Clphiphi
    Real*8, allocatable, dimension(:) :: Clpsilimber
    Real*8, allocatable, dimension(:) :: alpha_halo_mass_function
    Real*8, allocatable, dimension(:,:) :: dndM,bMz,M200d,M200c,r200d,r200c,dM200ddM,dM200cdM
    Real*8, allocatable, dimension(:,:,:) :: ylMz,philMZ 
    Integer*4,allocatable,dimension(:) :: ml
end module arrays
