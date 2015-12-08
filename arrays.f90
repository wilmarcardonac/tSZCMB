Module arrays

    Integer :: status1,status2,status3,status4,status5,status6
    Integer*4,allocatable,dimension(:) :: ml,ml_functions

    Real*8, allocatable, dimension(:) :: z,M,Cl1h,Cl2h,Cl,d2VdzdO,mbz,k,Clphiphi1h,Clphiphi2h,Clphiphi
    Real*8, allocatable, dimension(:) :: Clpsilimber,z_functions,M_functions,Scrit,zlimber
    Real*8, allocatable, dimension(:) :: alpha_halo_mass_function,comoving_distance_at_z,angular_diameter_distance_at_z
    Real*8, allocatable, dimension(:,:) :: dndM,bMz,M200d,M200c,r200d,r200c,dM200ddM,dM200cdM,dM200ddM_M_z,integrand_limber
    Real*8, allocatable, dimension(:,:) :: inte_cl_phiphi_1h,inte_cl_phiphi_2h,sigma_square_M200d,dsigma_square_M200d
    Real*8, allocatable, dimension(:,:,:) :: ylMz,philMZ,inte_pre_cl_phiphi,inte_pre_cl_phiphi_2h

End module arrays
