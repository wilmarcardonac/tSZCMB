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
    Real(fgsl_double),parameter :: relative_error = 1.0E-10_fgsl_double
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

End Module integrator
