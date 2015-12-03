Module fiducial

    Implicit none

    save 

    ! SAME COSMOLOGY AS IN PUBLISHED VERSION OF "DETECTION OF THERMAL SZ-CMB LENSING CROSS-CORRELATION IN PLANCK  
    ! NOMINAL MISSION DATA" BY HILL & SPERGEL.

    Integer*4,parameter :: number_of_k = 1000          ! SIZE OF WAVEVECTOR ARRAY 
    Integer*4,parameter :: number_of_z = 320           ! SIZE OF RED-SHIFT ARRAY 
    Integer*4,parameter :: number_of_z_functions = 1000  ! SIZE OF RED-SHIFT ARRAY FOR FUNCTIONS
    Integer*4,parameter :: number_of_z_limber = 1000   ! SIZE OF RED-SHIFT ARRAY FOR LIMBER APPROXIMATION 
    Integer*4,parameter :: number_of_M_log = 500        ! SIZE OF VIRIAL MASS ARRAY
    Integer*4,parameter :: number_of_M_linear = 332 
    Integer*4,parameter :: number_of_M = 1000 !number_of_M_log + number_of_M_linear
    Integer*4,parameter :: number_of_M_functions = 1000 ! SIZE OF VIRIAL MASS ARRAY FOR FUNCTIONS
    Integer*4,parameter :: number_of_l = 5  ! SIZE OF MULTIPOLE ARRAY 
    Integer*4,parameter :: number_of_l_functions = 100 ! SIZE OF MULTIPOLE ARRAY FOR INTERPOLATING FUNCTIONS
    Integer*4,parameter :: lmax = 10000        ! HIGHEST MULTIPOLE
    Integer*4,parameter :: lmin = 1        ! LOWEST MULTIPOLE
    Integer*4,parameter :: UNIT_EXE_FILE = 90 ! UNIT FOR EXECUTION INFORMATION FILE

    Real*8,parameter :: Pi = 3.141592653589793d0
    Real*8,parameter :: c = 2.99792458d8          ! SPEED OF LIGHT IN m s-1
    Real*8,parameter :: G = 6.67384d-11           ! NEWTONIAN CONSTANT OF GRAVITATION IN m3 kg-1 s-2
    Real*8,parameter :: M_sun = 1.989d30          ! SOLAR MASS IN kg
    Real*8,parameter :: Mpc = 3.08567758d22       ! 1 Mpc  = 3.08567758d22 m 
    Real*8,parameter :: m_e = 9.10938291d-31      ! ELECTRON MASS IN Kg
    Real*8,parameter :: sigma_e = 6.652458734d-29 ! THOMSON SCATTERING CROSS SECTION IN m^2

    Real*8,parameter :: h = 69.7d-2                ! DIMENSIONLESS HUBBLE PARAMETER AT PRESENT TIME 
    Real*8,parameter :: H_0 = h*1.d5               ! HUBBLE PARAMETER AT PRESENT TIME IN m s-1 Mpc-1
    Real*8,parameter :: Omega_m0 = 0.282d0         ! MATTER PARAMETER DENSITY AT PRESENT TIME
    Real*8,parameter :: Omega_L0 = 1.d0 - Omega_m0 ! DARK ENERGY PARAMETER DENSITY (LCDM MODEL) AT PRESENT TIME IN A FLAT UNIVERSE
    Real*8,parameter :: Omega_b_h2 = 0.02240d0     ! BARYON DENSITY TODAY
    Real*8,parameter :: ns = 0.9646d0              ! SCALAR SPECTRUM POWER-LAW INDEX
    Real*8,parameter :: z_dec = 1100d0             ! RED-SHIFT AT DECOUPLING
    !Real*8,parameter :: A_s = 2.43d-9             ! DIMENSIONLESS CURVATURE POWER SPECTRUM AT k_0 = 0.05 Mpc-1
    Real*8,parameter :: k_0 = 0.05                 ! PIVOT WAVEVECTOR IN Mpc-1
    Real*8,parameter :: kmin = 1.d-4               ! MINIMUM WAVENUMBER FOR COMPUTATION OF ANGULAR POWER SPECTRUM OF LENSING POTENTIAL
    Real*8,parameter :: kmax = 1.d3                ! MAXIMUM WAVENUMBER FOR COMPUTATION OF ANGULAR POWER SPECTRUM OF LENSING POTENTIAL 
    Real*8,parameter :: zmax = 1.d1                ! UPPER RED-SHIFT IN HILL & SPERGEL PAPER [1312.4525]
    Real*8,parameter :: zmin = 5.d-3               ! LOWER RED-SHIFT IN HILL & SPERGEL PAPER [1312.4525] 
    Real*8,parameter :: sigma8 = 0.817d0           ! RMS OF MATTER FLUCTUATIONS TODAY IN LINEAR THEORY
    Real*8,parameter :: Mmin = 1.d5/h              ! LOWER MASS IN HILL & SPERGEL PAPER [1312.4525]
    Real*8,parameter :: Mmax = 5.d15/h             ! UPPER MASS IN HILL & SPERGEL PAPER [1312.4525]
    Real*8,parameter :: tcmb0 = 2.728d0            ! CMB TEMPERATURE TODAY IN KELVIN
    Real*8 :: Normalization                        ! NORMALIZATION CONSTANT FOR MATTER POWER SPECTRUM (EQUATION (A7) IN EISENSTEIN & HU PAPER)
    Real*8 :: com_dist_at_z_dec                    ! COMOVING DISTANCE AT DECOUPLING
    Real*8,parameter :: DeltaSO = 2.d2             ! PARAMETER TO DEFINE SPHERICAL OVERDENSITY
    Real*8 :: alpha_halo_redshift_3                 ! CONSTANT ALPHA IN HALO MASS FUNCTION AT RED-SHIFT 3 (INTERPOLATION)
    Real*8 :: d_alpha_halo_redshift_3               ! DERIVATIVE OF CONSTANT ALPHA IN HALO MASS FUNCTION AT RED-SHIFT 3 (INTERPOLATION)

    Character(len=*),parameter :: path_to_execution_information = './output/execution_information.txt' ! EXECUTION INFORMATION FILE

    Logical,parameter :: do_mass_conversion = .false. !.true.       ! COMPUTE MASSES FOR OTHER HALO DEFINITIONS IF SET IT TRUE
    Logical,parameter :: compute_functions = .true.                 ! COMPUTE FUNCTIONS HAVING 'number_of_z' AND 'number_of_M' SIZE OF ARRAYS, OTHERWISE INTERPOLATE
    Logical,parameter :: compute_sigma_square = .false. !.true.              ! COMPUTE SIGMA SQUARE AND ITS DERIVATIVE IF SET IT TRUE
    Logical,parameter :: compute_linear_halo_bias = .false.!.true. ! COMPUTE LINEAR HALO BIAS IF SET IT TRUE
    Logical,parameter :: compute_alpha_in_halo_mass_function = .true. ! COMPUTE CONSTANT ALPHA IN HALO MASS FUNCTION IF SET IT TRUE
    Logical,parameter :: compute_halo_mass_function = .true. ! COMPUTE HALO MASS FUNCTION IF SET IT TRUE
    Logical,parameter :: compute_the_lensing_potential = .true. ! COMPUTE LENSING POTENTIAL IF SET IT TRUE
    Logical,parameter :: compute_the_form_factor = .true. !.false. !.true. ! COMPUTE FORM FACTOR IF SET IT TRUE

End Module fiducial
