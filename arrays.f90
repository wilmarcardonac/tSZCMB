Module arrays

    Integer :: status1,status2,status3,status4,status5,status6
    Integer*4,allocatable,dimension(:) :: ml,ml_functions

    Real*8, allocatable, dimension(:) :: z,M,Cl1h,Cl2h,Cl,d2VdzdO,mbz,k,Clphiphi1h,Clphiphi2h,Clphiphi
    Real*8, allocatable, dimension(:) :: Clpsilimber,z_functions,M_functions,Scrit
    Real*8, allocatable, dimension(:) :: alpha_halo_mass_function,comoving_distance_at_z,angular_diameter_distance_at_z
    Real*8, allocatable, dimension(:,:) :: dndM,bMz,M200d,M200c,r200d,r200c,dM200ddM,dM200cdM,dM200ddM_M_z
    Real*8, allocatable, dimension(:,:) :: bMz_interpolation,dndM_interpolation
    Real*8, allocatable, dimension(:,:,:) :: ylMz,philMZ,philMz_interpolation,ylMz_interpolation

End module arrays
