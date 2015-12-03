Module integrator

  use fgsl
  use, intrinsic :: iso_c_binding
  use fiducial
  use functions
  !use arrays

  Implicit none

contains

  function integrand_sigma_square(k,params) bind(c) ! Equation (6) of 1303.4726. Units of Mass: solar mass. Units : dimensionless

    real(c_double), value :: k
    type(c_ptr), value :: params
    real(c_double) :: integrand_sigma_square

    real(c_double), pointer :: p(:)

    call c_f_pointer(params, p, (/2/)) ! p(1) = R, p(2) = redshift

    integrand_sigma_square =  k**2*window_function(k,p(1))**2*matter_power_spectrum(k,p(2))

  End function integrand_sigma_square

  subroutine sigma_square_at_M_and_z(mass,redshift,output)

    Implicit none

    Integer(fgsl_size_t), parameter :: nmax=1000000

    Real(fgsl_double), target :: fpar(2)
    Real(fgsl_double) :: result, error, output, R, mass, redshift
    Real(fgsl_double) :: lower_limit != kmin !1.0E-3_fgsl_double
    Real(fgsl_double) :: upper_limit != 1.d-1*kmax !1.0E-2_fgsl_double
    Real(fgsl_double),parameter :: absolute_error = 0.0_fgsl_double
    Real(fgsl_double),parameter :: relative_error = 1.0E-5_fgsl_double
    Real(fgsl_double),parameter :: upper_scale = 1.0E4_fgsl_double
    Real(fgsl_double),parameter :: lower_scale = 1.0E-3_fgsl_double

    Integer(fgsl_int) :: status

    Type(c_ptr) :: ptr
    Type(fgsl_function) :: f_obj
    Type(fgsl_integration_workspace) :: wk

    R = (3.d0*mass/4.d0/Pi/mean_density(redshift)*(1.d0 + redshift)**3.d0)**(1.d0/3.d0) ! Units : Mpc. (1+z)**3 factor to use comoving coordinates

    upper_limit = 1/R*upper_scale

    lower_limit = 1/R*lower_scale

    fpar = (/R,redshift/)

    ptr = c_loc(fpar)

    f_obj = fgsl_function_init(integrand_sigma_square, ptr)

    wk = fgsl_integration_workspace_alloc(nmax)

    status = fgsl_integration_qags(f_obj, lower_limit, upper_limit, &
         absolute_error, relative_error, nmax, wk, result, error)

    output = 1.d0/2.d0/Pi**2*result/h**3

    call fgsl_function_free(f_obj)

    call fgsl_integration_workspace_free(wk)

  end subroutine sigma_square_at_M_and_z

  function integrand_dsigma_square(k,params) bind(c) ! Equation (6) of 1303.4726. Units of Mass: solar mass. Units : dimensionless

    real(c_double), value :: k
    type(c_ptr), value :: params
    real(c_double) :: integrand_dsigma_square

    real(c_double), pointer :: p(:)

    call c_f_pointer(params, p, (/2/)) ! p(1) = R, p(2) = redshift

    integrand_dsigma_square =  k**3*window_function(k,p(1))*d_window_function(k,p(1))*matter_power_spectrum(k,p(2))

  End function integrand_dsigma_square

  subroutine dsigma_square_at_M_and_z(mass,redshift,output)

    Implicit none

    Integer(fgsl_size_t), parameter :: nmax=1000000

    Real(fgsl_double), target :: fpar(2)
    Real(fgsl_double) :: result, error, output, R, mass, redshift
    Real(fgsl_double) :: lower_limit != kmin !1.0E-3_fgsl_double
    Real(fgsl_double) :: upper_limit != 1.d-1*kmax !1.0E-2_fgsl_double
    Real(fgsl_double),parameter :: absolute_error = 0.0_fgsl_double
    Real(fgsl_double),parameter :: relative_error = 1.0E-5_fgsl_double
    Real(fgsl_double),parameter :: upper_scale = 1.0E4_fgsl_double
    Real(fgsl_double),parameter :: lower_scale = 1.0E-3_fgsl_double

    Integer(fgsl_int) :: status

    Type(c_ptr) :: ptr
    Type(fgsl_function) :: f_obj
    Type(fgsl_integration_workspace) :: wk

    R = (3.d0*mass/4.d0/Pi/mean_density(redshift)*(1.d0 + redshift)**3.d0)**(1.d0/3.d0) ! Units : Mpc. (1+z)**3 factor to use comoving coordinates

    upper_limit = 1/R*upper_scale

    lower_limit = 1/R*lower_scale

    fpar = (/R,redshift/)

    ptr = c_loc(fpar)

    f_obj = fgsl_function_init(integrand_dsigma_square, ptr)

    wk = fgsl_integration_workspace_alloc(nmax)

    status = fgsl_integration_qags(f_obj, lower_limit, upper_limit, &
         absolute_error, relative_error, nmax, wk, result, error)

    output = 1.d0/Pi**2*result/h**3

    call fgsl_function_free(f_obj)

    call fgsl_integration_workspace_free(wk)

  end subroutine dsigma_square_at_M_and_z

  subroutine compute_sigma_square_M200d()    ! It fills sigma square array in and writes it out to a text file

    use arrays
    use fiducial
    use omp_lib
    Implicit none

    Integer*4 :: indexz,indexM

    open(15,file='./precomputed_quantities/sigma_square_M200d.dat')

    write(15,*) '# Sigma square function (as function of M200d mass and red-shift). Number of lines is ',number_of_z*number_of_M

    write(15,*) '# index_of_M    M200d[solar mass]    index_of_z    red-shift    sigma_square   dsigma_square '

    !$omp Parallel Do Shared(sigma_square_M200d,dsigma_square_M200d)

    Do indexM=1,number_of_M

       Do indexz=1,number_of_z

          call sigma_square_at_M_and_z(M200d(indexM,indexz),z(indexz),sigma_square_M200d(indexM,indexz))

          call dsigma_square_at_M_and_z(M200d(indexM,indexz),z(indexz),dsigma_square_M200d(indexM,indexz))

          write(15,'(i10,es18.10,i5,3es18.10)') indexM,M200d(indexM,indexz),indexz,z(indexz),&
               sigma_square_M200d(indexM,indexz),dsigma_square_M200d(indexM,indexz)

       End Do

    End Do

   !$omp End Parallel Do

    close(15)

  end subroutine compute_sigma_square_M200d

  Function integrand_normalization(k, params) bind(c)

    real(c_double), value :: k
    type(c_ptr), value :: params
    real(c_double) :: integrand_normalization

    real(c_double), pointer :: R

    call c_f_pointer(params, R)

    integrand_normalization =  Delta_squared_non_normalised(k)*(3.d0*(sin(k*R)- k*R*cos(k*R))/k**3/R**3)**2/k 

  End function integrand_normalization

  Subroutine compute_normalization(output)

    Implicit none

    Integer(fgsl_size_t), parameter :: nmax=10000

    Real(fgsl_double), target :: R
    Real(fgsl_double) :: result, error, output
    Real(fgsl_double),parameter :: lower_limit = kmin!1.0E-3_fgsl_double
    Real(fgsl_double),parameter :: upper_limit = kmax!1.0E-2_fgsl_double
    Real(fgsl_double),parameter :: absolute_error = 0.0_fgsl_double
    Real(fgsl_double),parameter :: relative_error = 1.0E-12_fgsl_double

    Integer(fgsl_int) :: status

    Type(c_ptr) :: ptr
    Type(fgsl_function) :: f_obj
    Type(fgsl_integration_workspace) :: wk

    R = 8.0d0/h

    ptr = c_loc(R)

    f_obj = fgsl_function_init(integrand_normalization, ptr)

    wk = fgsl_integration_workspace_alloc(nmax)

    status = fgsl_integration_qags(f_obj, lower_limit, upper_limit, &
         absolute_error, relative_error, nmax, wk, result, error)

    output = sigma8**2/result

    call fgsl_function_free(f_obj)

    call fgsl_integration_workspace_free(wk)

  End subroutine compute_normalization

  Function integrand_comoving_distance(z, params) bind(c)

    real(c_double), value :: z
    type(c_ptr), value :: params
    real(c_double) :: integrand_comoving_distance

    real(c_double), pointer :: RR

    call c_f_pointer(params, RR)

    integrand_comoving_distance =  1.d0/Hubble_parameter(z)

  End function integrand_comoving_distance

  Subroutine comoving_distance_at_redshift(z,output)    !    UNITS: Mpc

    Implicit none

    Integer(fgsl_size_t), parameter :: nmax=10000

    Real(fgsl_double), target :: RR
    Real(fgsl_double) :: result, error, output, z
    Real(fgsl_double),parameter :: lower_limit = 0.0_fgsl_double
    Real(fgsl_double) :: upper_limit !1.0E-2_fgsl_double
    Real(fgsl_double),parameter :: absolute_error = 0.0_fgsl_double
    Real(fgsl_double),parameter :: relative_error = 1.0E-12_fgsl_double

    Integer(fgsl_int) :: status

    Type(c_ptr) :: ptr
    Type(fgsl_function) :: f_obj
    Type(fgsl_integration_workspace) :: wk

    upper_limit = z

    RR = 1.0d0

    ptr = c_loc(RR)

    f_obj = fgsl_function_init(integrand_comoving_distance, ptr)

    wk = fgsl_integration_workspace_alloc(nmax)

    status = fgsl_integration_qags(f_obj, lower_limit, upper_limit, &
         absolute_error, relative_error, nmax, wk, result, error)

    output = c*result

    call fgsl_function_free(f_obj)

    call fgsl_integration_workspace_free(wk)

  End subroutine comoving_distance_at_redshift

  subroutine compute_comoving_and_angular_diameter_distance() ! AND CRITICAL SURFACE DENSITY

    use arrays
    use fiducial
    !use omp_lib
    Implicit none

    Integer*4 :: indexz

    If (compute_functions) then

!!$omp Parallel Do Shared(comoving_distance_at_z,z)

       Do indexz=1,number_of_z

          call comoving_distance_at_redshift(z(indexz),comoving_distance_at_z(indexz))

          angular_diameter_distance_at_z(indexz) = comoving_distance_at_z(indexz)/(1.d0 + z(indexz))

          Scrit(indexz) = c**2*comoving_distance_at_z(indexz)*(1.d0 + &
               z(indexz))/4.d0/Pi/G/comoving_distance_at_z(indexz)/(com_dist_at_z_dec - comoving_distance_at_z(indexz))*Mpc/M_sun

       End Do

!!$omp End Parallel Do

    Else

!!$omp Parallel Do Shared(comoving_distance_at_z,z_functions)

       Do indexz=1,number_of_z_functions

          call comoving_distance_at_redshift(z_functions(indexz),comoving_distance_at_z(indexz))

          angular_diameter_distance_at_z(indexz) = comoving_distance_at_z(indexz)/(1.d0 + z_functions(indexz))

          Scrit(indexz) = c**2*comoving_distance_at_z(indexz)*(1.d0 + &
               z_functions(indexz))/4.d0/Pi/G/comoving_distance_at_z(indexz)/(com_dist_at_z_dec - &
               comoving_distance_at_z(indexz))*Mpc/M_sun

       End Do

!!$omp End Parallel Do

    End If

  end subroutine compute_comoving_and_angular_diameter_distance

  function nonnormalised_halo_mass_function(virial_Mass_index,redshift_index)    ! Equation (8) in "The large-scale bias of dark matter halos: 

    use fiducial    ! numerical calibration and model tests" without \alpha constant. This function corresponds to Equation (C2) in
    use arrays    ! "Toward a halo mass function for precision cosmology: the limits of universality" by Jeremy Tinker et al. Parameters taken 
    use mod_roots

    Implicit none    ! from first line of table 4. Units : 1/(solar mass*Mpc**3)

    Real*8 :: nonnormalised_halo_mass_function,R,sigma,nu,beta,phi,eta,gamma,M200d_M_z
    Real*8 :: f_nu,g_sigma,rdeltad,dsigma
    Real*8,parameter :: delta_c = 1.686d0    ! page 880 in "The large-scale ..."
    Real*8,parameter :: beta0 = 0.589d0
    Real*8,parameter :: phi0 = -0.729d0
    Real*8,parameter :: eta0 = -0.243d0
    Real*8,parameter :: gamma0 = 0.864d0
    Integer*4 :: virial_Mass_index,redshift_index

    If (z(redshift_index) .gt. 3.d0) then

       call compute_r_delta_d_from_M_virial_at_z(M(virial_Mass_index),3.d0,DeltaSO,rdeltad)

       M200d_M_z = M_delta_d_from_M_virial(3.d0,rdeltad,DeltaSO)

       R = (3.d0*M200d_M_z/4.d0/Pi/mean_density(3.d0)*( 1.d0 + 3.d0 )**3.d0)**(1.d0/3.d0)    ! Units : Mpc. (1+z)**3 to use comoving coordinates

       call sigma_square_at_M_and_z(M200d_M_z,3.d0,sigma)

       call dsigma_square_at_M_and_z(M200d_M_z,3.d0,dsigma)

       sigma = sqrt(sigma) !sigma_squared(M200d_M_z,3.d0))    ! Units : dimensionless

       nu = delta_c/sigma           ! page 880 in "The large-scale ... tests"

       beta = beta0*( 1.d0 + 3.d0 )**(0.20d0) 

       phi = phi0*( 1.d0 + 3.d0 )**(-0.08d0)

       eta = eta0*( 1.d0 + 3.d0 )**(0.27d0)

       gamma = gamma0*( 1.d0 + 3.d0 )**(-0.01d0)

       f_nu = (1.d0 + (beta*nu)**(-2.d0*phi))*nu**(2.d0*eta)*dexp(-gamma*nu**2/2.d0)

       g_sigma = nu*f_nu            ! Equation (C2) in "Toward a ... universality"

       nonnormalised_halo_mass_function = -mean_density(3.d0)/(1.d0 + 3.d0 )**3.d0/2.d0/M200d_M_z**2&
            *R/3.d0/sigma**2*dsigma*g_sigma ! dsigma_squared_dR(M200d_M_z,3.d0)*g_sigma

    Else

       R = (3.d0*M200d(virial_Mass_index,redshift_index)/4.d0/Pi/mean_density(z(redshift_index))*(1.d0 + &
            z(redshift_index) )**3.d0)**(1.d0/3.d0)    ! Units : Mpc. (1+z)**3 to use comoving coordinates

       sigma = sqrt(sigma_square_M200d(virial_Mass_index,redshift_index)) !sqrt(sigma_squared(M200d_M_z,redshift))    ! Units : dimensionless

       !sigma = sqrt(sigma_squared(M200d(virial_Mass_index,redshift_index),z(redshift_index)))    ! Units : dimensionless

       dsigma = dsigma_square_M200d(virial_Mass_index,redshift_index)

       nu = delta_c/sigma           ! page 880 in "The large-scale ... tests"

       beta = beta0*( 1.d0 + z(redshift_index) )**(0.20d0) 

       phi = phi0*( 1.d0 + z(redshift_index) )**(-0.08d0)

       eta = eta0*( 1.d0 + z(redshift_index) )**(0.27d0)

       gamma = gamma0*( 1.d0 + z(redshift_index) )**(-0.01d0)

       f_nu = (1.d0 + (beta*nu)**(-2.d0*phi))*nu**(2.d0*eta)*exp(-gamma*nu**2/2.d0)

       g_sigma = nu*f_nu            ! Equation (C2) in "Toward a ... universality"

       nonnormalised_halo_mass_function = -mean_density(z(redshift_index))/(1.d0 + z(redshift_index))**3.d0/&
            2.d0/M200d(virial_Mass_index,redshift_index)**2*&
            R/3.d0/sigma**2*dsigma*g_sigma !dsigma_squared_dR(M200d(virial_Mass_index,redshift_index),z(redshift_index))*g_sigma

    End If

  end function nonnormalised_halo_mass_function

  Function integrand_alpha_halo_mass_function(virial_mass, params) bind(c)

    real(c_double), value :: virial_mass
    type(c_ptr), value :: params
    real(c_double) :: integrand_alpha_halo_mass_function
    Integer(c_int), pointer :: indexz

    call c_f_pointer(params, indexz)

    integrand_alpha_halo_mass_function =  integrand_alpha_halo_mass_function_at_z(virial_mass,indexz)

  End function integrand_alpha_halo_mass_function

  function integrand_alpha_halo_mass_function_at_z(virial_mass,indexz)

    use arrays
    use fiducial

    Implicit none

    Real*8 :: virial_mass,integrand_alpha_halo_mass_function_at_z,dalpha
    Real*8,dimension(number_of_M) :: f

    Integer*4 :: indexz,indexM

    Do indexM=1,number_of_M

       f(indexM) = nonnormalised_halo_mass_function(indexM,indexz)*bMz(indexM,indexz)*&
            M(indexM)/mean_density(z(indexz))*(1.d0 + z(indexz))**3.d0*&
            dM200ddM(indexM,indexz) 

    End Do

    call Interpolate_1D(integrand_alpha_halo_mass_function_at_z,dalpha,virial_mass,M,f)

  End function integrand_alpha_halo_mass_function_at_z

  Subroutine compute_alpha_halo_mass_function_at_z(indexz,output)

    Implicit none

    Integer(fgsl_size_t), parameter :: nmax=10000

    Integer(fgsl_int),target :: pp
    Real(fgsl_double) :: result, error, output
    Real(fgsl_double),parameter :: lower_limit = Mmin!1.0E-3_fgsl_double
    Real(fgsl_double),parameter :: upper_limit = Mmax!1.0E-2_fgsl_double
    Real(fgsl_double),parameter :: absolute_error = 0.0_fgsl_double
    Real(fgsl_double),parameter :: relative_error = 1.0E-5_fgsl_double

    Integer(fgsl_int) :: status
    Integer(fgsl_int) :: indexz

    Type(c_ptr) :: ptr
    Type(fgsl_function) :: f_obj
    Type(fgsl_integration_workspace) :: wk

    pp = indexz

    ptr = c_loc(pp)

    f_obj = fgsl_function_init(integrand_alpha_halo_mass_function, ptr)

    wk = fgsl_integration_workspace_alloc(nmax)

    status = fgsl_integration_qags(f_obj, lower_limit, upper_limit, &
         absolute_error, relative_error, nmax, wk, result, error)

    output = 1.d0/result

    call fgsl_function_free(f_obj)

    call fgsl_integration_workspace_free(wk)

  End subroutine compute_alpha_halo_mass_function_at_z

  subroutine compute_alpha_halo_mass_function()    ! It computes \alpha constant in halo mass function (Equation (8)) of "The large-scale 

    use fiducial    ! bias of dark matter halos: numerical calibration and model tests" for the red-shift array
    use omp_lib
    use arrays
    Implicit none

    Integer*4 :: indexz

    open(15,file='./precomputed_quantities/alpha_halo_mass_function.dat')

    write(15,*) '# Alpha halo mass function (as function of red-shift). Number of lines is ',number_of_z

    write(15,*) '# index_of_z    red-shift    dndM '

    If (compute_functions) then

       !$omp Parallel Do Shared(alpha_halo_mass_function)

       Do indexz=1,number_of_z

          call compute_alpha_halo_mass_function_at_z(indexz,alpha_halo_mass_function(indexz))   ! dimensionless

          write(15,'(i5,2es18.10)') indexz, z(indexz), alpha_halo_mass_function(indexz)

       End Do

       !$omp End Parallel Do

    Else

       print *,'COMPUTE ALPHA HALO MASS FUNCTIONS ROUTINE WORKS ONLY FOR M AND Z ARRAYS OF SIZE ', number_of_M, number_of_z

       stop

    End If

    close(15)

  end subroutine compute_alpha_halo_mass_function

  function halo_mass_function(virial_Mass_index,redshift_index) ! Equation (8) in "The large-scale bias of dark matter halos: numerical calibration and model tests".

    use fiducial                 ! Units of M [M_{200d}] : solar mass. Units : 1/(solar mass*Mpc**3). Halo definition : mean density
    use arrays
    use mod_roots

    Implicit none                 

    Real*8 :: halo_mass_function,R,sigma,nu,beta,phi,eta,gamma,M200d_M_z
    Real*8 :: f_nu,g_sigma,rdeltad,dsigma!,alpha_halo_mass_function_z        
    Real*8,parameter :: delta_c = 1.686d0    ! page 880 in "The large-scale ..."
    Real*8,parameter :: beta0 = 0.589d0
    Real*8,parameter :: phi0 = -0.729d0
    Real*8,parameter :: eta0 = -0.243d0
    Real*8,parameter :: gamma0 = 0.864d0
    Integer*4 :: virial_Mass_index,redshift_index

    If (z(redshift_index) .gt. 3.d0) then

       call compute_r_delta_d_from_M_virial_at_z(M(virial_Mass_index),3.d0,DeltaSO,rdeltad)

       M200d_M_z = M_delta_d_from_M_virial(3.d0,rdeltad,DeltaSO)

       R = (3.d0*M200d_M_z/4.d0/Pi/mean_density(3.d0)*(1.d0 + 3.d0 )**3.d0)**(1.d0/3.d0)    ! Units : Mpc. (1+z)**3 to use comoving coordinates

       call sigma_square_at_M_and_z(M200d_M_z,3.d0,sigma)

       call dsigma_square_at_M_and_z(M200d_M_z,3.d0,dsigma)

       sigma = sqrt(sigma)!sigma_squared(M200d_M_z,3.d0))    ! Units : dimensionless

       nu = delta_c/sigma           ! page 880 in "The large-scale ... tests"

       beta = beta0*( 1.d0 + 3.d0 )**(0.20d0) 

       phi = phi0*( 1.d0 + 3.d0 )**(-0.08d0)

       eta = eta0*( 1.d0 + 3.d0 )**(0.27d0)

       gamma = gamma0*( 1.d0 + 3.d0 )**(-0.01d0)

       f_nu = alpha_halo_redshift_3*(1.d0 + (beta*nu)**(-2.d0*phi))*&
            nu**(2.d0*eta)*exp(-gamma*nu**2/2.d0)

       g_sigma = nu*f_nu ! dimensionless           ! Equation (C2) in "Toward a ... universality"

       halo_mass_function = -mean_density(3.d0)/(1.d0 + 3.d0)**3.d0/2.d0/M200d_M_z**2&    ! (1+z)**3 to use comoving coordinates
            *R/3.d0/sigma**2*dsigma*g_sigma  !dsigma_squared_dR(M200d_M_z,3.d0)*g_sigma

    Else

       R = (3.d0*M200d(virial_Mass_index,redshift_index)/4.d0/Pi/mean_density(z(redshift_index))*(1.d0 + &
            z(redshift_index))**3.d0)**(1.d0/3.d0)    ! Units : Mpc. (1+z)**3 to use comoving coordinates

       sigma = sqrt(sigma_square_M200d(virial_Mass_index,redshift_index))    ! Units : dimensionless

       dsigma = dsigma_square_M200d(virial_Mass_index,redshift_index)

       nu = delta_c/sigma           ! page 880 in "The large-scale ... tests"

       beta = beta0*( 1.d0 + z(redshift_index) )**(0.20d0) 

       phi = phi0*( 1.d0 + z(redshift_index) )**(-0.08d0)

       eta = eta0*( 1.d0 + z(redshift_index) )**(0.27d0)

       gamma = gamma0*( 1.d0 + z(redshift_index) )**(-0.01d0)

       f_nu = alpha_halo_mass_function(redshift_index)*(1.d0 + (beta*nu)**(-2.d0*phi))*nu**(2.d0*eta)*&
            exp(-gamma*nu**2/2.d0)

       g_sigma = nu*f_nu ! dimensionless           ! Equation (C2) in "Toward a ... universality"

       halo_mass_function = -mean_density(z(redshift_index))/(1.d0 + z(redshift_index) )**3.d0/2.d0/&
            M200d(virial_Mass_index,redshift_index)**2*&    ! (1+z)**3 to use comoving coordinates
            R/3.d0/sigma**2*dsigma*g_sigma !dsigma_squared_dR(M200d(virial_Mass_index,redshift_index),z(redshift_index))*g_sigma

    End If

  end function halo_mass_function

  subroutine compute_dndM()    ! It fills halo mass function array in and writes it out to a text file

    use arrays
    use fiducial
    Implicit none

    Integer*4 :: indexz,indexM

    open(15,file='./precomputed_quantities/dndM/dndM.dat')

    write(15,*) '# Halo mass function (as function of virial mass and red-shift). Number of lines is ',&
         number_of_z*number_of_M

    write(15,*) '#  index_of_M    virial_mass[solar mass]    index_of_z    red-shift    dndM '

    If (compute_functions) then

       !$omp Parallel Do Shared(dndM,dM200ddM)

       Do indexM=1,number_of_M

          Do indexz=1,number_of_z

             dndM(indexM,indexz) = halo_mass_function(indexM,indexz)*dM200ddM(indexM,indexz)

             write(15,'(i10,es18.10,i5,2es18.10)') indexM, M(indexM), indexz, z(indexz), dndM(indexM,indexz)

          End Do

       End Do

       !$omp End Parallel Do

    Else

       print *,'HALO MASS FUNCTION COMPUTED ONLY FOR ARRAYS M AND Z OF SIZE ', number_of_M, number_of_z

       stop

    End If

    close(15)

  end subroutine compute_dndM


End Module integrator
