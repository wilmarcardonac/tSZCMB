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

    use fiducial

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

!!$omp Parallel Do Shared(comoving_distance_at_z,z)

    Do indexz=1,number_of_z

       call comoving_distance_at_redshift(z(indexz),comoving_distance_at_z(indexz))

       angular_diameter_distance_at_z(indexz) = comoving_distance_at_z(indexz)/(1.d0 + z(indexz))

       Scrit(indexz) = c**2*comoving_distance_at_z(indexz)*(1.d0 + &
            z(indexz))/4.d0/Pi/G/comoving_distance_at_z(indexz)/(com_dist_at_z_dec - comoving_distance_at_z(indexz))*Mpc/M_sun

       d2VdzdO(indexz) = c*(1.d0 + z(indexz))**2*angular_diameter_distance_at_z(indexz)**2/Hubble_parameter(z(indexz))

    End Do

!!$omp End Parallel Do

  end subroutine compute_comoving_and_angular_diameter_distance

  function nonnormalised_halo_mass_function_f_nu(nu,beta,phi,eta,gamma)

    Implicit none

    Real*8 :: nonnormalised_halo_mass_function_f_nu,nu,beta,phi,eta,gamma

    nonnormalised_halo_mass_function_f_nu = (1.d0 + (beta*nu)**(-2.d0*phi))*nu**(2.d0*eta)*exp(-gamma*nu**2/2.d0)

  end function nonnormalised_halo_mass_function_f_nu

  function nonnormalised_halo_mass_function(virial_Mass_index,redshift_index)    ! Equation (8) in "The large-scale bias of dark matter halos: 

    use fiducial    ! numerical calibration and model tests" without \alpha constant. This function corresponds to Equation (C2) in
    use arrays    ! "Toward a halo mass function for precision cosmology: the limits of universality" by Jeremy Tinker et al. Parameters taken 
    use mod_roots

    Implicit none    ! from first line of table 4. Units : 1/(solar mass*Mpc**3)

    Real*8 :: nonnormalised_halo_mass_function,R,sigma,nu,M200d_M_z
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

       f_nu = halo_mass_function_f_nu(nu,1.d0,3.d0) 

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

       f_nu = halo_mass_function_f_nu(nu,1.d0,z(redshift_index)) 

       g_sigma = nu*f_nu            ! Equation (C2) in "Toward a ... universality"

       nonnormalised_halo_mass_function = -mean_density(z(redshift_index))/(1.d0 + z(redshift_index))**3.d0/&
            2.d0/M200d(virial_Mass_index,redshift_index)**2*&
            R/3.d0/sigma**2*dsigma*g_sigma !dsigma_squared_dR(M200d(virial_Mass_index,redshift_index),z(redshift_index))*g_sigma

    End If

  end function nonnormalised_halo_mass_function

  Function integrand_alpha_halo_mass_function(nu, params) bind(c)

    real(c_double), value :: nu
    type(c_ptr), value :: params
    real(c_double) :: integrand_alpha_halo_mass_function
    Real(c_double), pointer :: pa(:)

    call c_f_pointer(params, pa,(/2/)) ! pa(1) = alpha, pa(2) = red-shift

    integrand_alpha_halo_mass_function =  linear_halo_bias_b_nu(nu)*halo_mass_function_f_nu(nu,pa(1),pa(2))
!integrand_alpha_halo_mass_function_at_z(virial_mass,indexz)

  End function integrand_alpha_halo_mass_function

!  function integrand_alpha_halo_mass_function_at_z(virial_mass,indexz)

  !  use arrays
 !   use fiducial
!    use omp_lib

!    Implicit none

    !Real*8 :: virial_mass,integrand_alpha_halo_mass_function_at_z,dalpha
   ! Real*8,dimension(number_of_M) :: f

  !  Integer*4 :: indexz,indexM

 !   !$omp Parallel Do Shared(f,bMz,M,z,dM200ddM)

    !Do indexM=1,number_of_M

     !  f(indexM) = nonnormalised_halo_mass_function(indexM,indexz)*bMz(indexM,indexz)*&
      !      M(indexM)/mean_density(z(indexz))*(1.d0 + z(indexz))**3.d0*&
       !     dM200ddM(indexM,indexz) 

!    End Do

  !  !$omp End Parallel Do

   ! call Interpolate_1D(integrand_alpha_halo_mass_function_at_z,dalpha,virial_mass,M,f)

 ! End function integrand_alpha_halo_mass_function_at_z

  Subroutine compute_alpha_halo_mass_function_at_z(redshift,output)

    use fiducial
    use arrays

    Implicit none

    Integer(fgsl_size_t), parameter :: nmax=10000

    Real(fgsl_double),target :: pp(2)
    Real(fgsl_double) :: result, error, output, redshift
    Real(fgsl_double) :: lower_limit != 1.686/sqrt(maxval(sigma_square_M200d))! 1.0E-1_fgsl_double
    Real(fgsl_double) :: upper_limit != 1.13E2_fgsl_double
    Real(fgsl_double),parameter :: absolute_error = 0.0_fgsl_double
    Real(fgsl_double),parameter :: relative_error = 1.0E-8_fgsl_double

    Integer(fgsl_int) :: status

    Type(c_ptr) :: ptr
    Type(fgsl_function) :: f_obj
    Type(fgsl_integration_workspace) :: wk

    lower_limit = 1.686d0/sqrt(maxval(sigma_square_M200d))

    upper_limit = 1.686d0/sqrt(minval(sigma_square_M200d))

    If (redshift .gt. 3.d0) then

       pp = (/1.d0,3.d0/)

    Else

       pp = (/1.d0,redshift/)

    End If
    
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

    !$omp Parallel Do Shared(alpha_halo_mass_function)

    Do indexz=1,number_of_z

       call compute_alpha_halo_mass_function_at_z(z(indexz),alpha_halo_mass_function(indexz))   ! dimensionless

    End Do

    !$omp End Parallel Do

    Do indexz=1,number_of_z

       write(15,'(i5,2es18.10)') indexz, z(indexz), alpha_halo_mass_function(indexz)

    End Do

    close(15)

  end subroutine compute_alpha_halo_mass_function

  function halo_mass_function_f_nu(nu,alpha,z)

    Implicit none

    Real*8 :: halo_mass_function_f_nu,nu,alpha,beta,phi,eta,gamma,z
    Real*8,parameter :: beta0 = 0.589d0    ! page 880 in "The large-scale ..."
    Real*8,parameter :: phi0 = -0.729d0
    Real*8,parameter :: eta0 = -0.243d0
    Real*8,parameter :: gamma0 = 0.864d0

    beta = beta0*( 1.d0 + z )**(0.20d0) 

    phi = phi0*( 1.d0 + z )**(-0.08d0)

    eta = eta0*( 1.d0 + z )**(0.27d0)

    gamma = gamma0*( 1.d0 + z )**(-0.01d0)

    halo_mass_function_f_nu = alpha*(1.d0 + (beta*nu)**(-2.d0*phi))*&
         nu**(2.d0*eta)*exp(-gamma*nu**2/2.d0)

  end function halo_mass_function_f_nu

  function halo_mass_function(virial_Mass_index,redshift_index) ! Equation (8) in "The large-scale bias of dark matter halos: numerical calibration and model tests".

    use fiducial                 ! Units of M [M_{200d}] : solar mass. Units : 1/(solar mass*Mpc**3). Halo definition : mean density
    use arrays
    use mod_roots

    Implicit none                 

    Real*8 :: halo_mass_function,R,sigma,nu,M200d_M_z,alpha
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

       alpha = alpha_halo_redshift_3

       f_nu = halo_mass_function_f_nu(nu,alpha,3.d0)

       g_sigma = nu*f_nu ! dimensionless           ! Equation (C2) in "Toward a ... universality"

       halo_mass_function = -mean_density(3.d0)/(1.d0 + 3.d0)**3.d0/2.d0/M200d_M_z**2&    ! (1+z)**3 to use comoving coordinates
            *R/3.d0/sigma**2*dsigma*g_sigma  !dsigma_squared_dR(M200d_M_z,3.d0)*g_sigma

    Else

       R = (3.d0*M200d(virial_Mass_index,redshift_index)/4.d0/Pi/mean_density(z(redshift_index))*(1.d0 + &
            z(redshift_index))**3.d0)**(1.d0/3.d0)    ! Units : Mpc. (1+z)**3 to use comoving coordinates

       sigma = sqrt(sigma_square_M200d(virial_Mass_index,redshift_index))    ! Units : dimensionless

       dsigma = dsigma_square_M200d(virial_Mass_index,redshift_index)

       nu = delta_c/sigma           ! page 880 in "The large-scale ... tests"

       alpha = alpha_halo_mass_function(redshift_index)

       f_nu = halo_mass_function_f_nu(nu,alpha,z(redshift_index))

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

    !$omp Parallel Do Shared(dndM,dM200ddM)

    Do indexM=1,number_of_M

       Do indexz=1,number_of_z

          dndM(indexM,indexz) = halo_mass_function(indexM,indexz)*dM200ddM(indexM,indexz)

          write(15,'(i10,es18.10,i5,2es18.10)') indexM, M(indexM), indexz, z(indexz), dndM(indexM,indexz)

       End Do

    End Do

    !$omp End Parallel Do

    Do indexM=1,number_of_M

       Do indexz=1,number_of_z

          write(15,'(i10,es18.10,i5,2es18.10)') indexM, M(indexM), indexz, z(indexz), dndM(indexM,indexz)

       End Do

    End Do

    close(15)

  end subroutine compute_dndM

  Function integrand_mean_bias_matter(virial_mass, params) bind(c)

    real(c_double), value :: virial_mass
    type(c_ptr), value :: params
    real(c_double) :: integrand_mean_bias_matter
    Integer(c_int), pointer :: indexz

    call c_f_pointer(params, indexz)

    integrand_mean_bias_matter =  integrand_mean_bias_matter_at_z(virial_mass,indexz)

  End function integrand_mean_bias_matter

  subroutine compute_integrand_mean_bias_matter_at_z()

    use fiducial
    use arrays
    use functions 

    Implicit none

    Integer*4 :: indexM,indexz

    Do indexM=1,number_of_M

       Do indexz=1,number_of_z

          inte_mbz(indexM,indexz) = dndM(indexM,indexz)*bMz(indexM,indexz)*M(indexM)/&
            mean_density(z(indexz))*(1.d0 + z(indexz))**3.d0

       End Do

    End Do

  end subroutine compute_integrand_mean_bias_matter_at_z

  function integrand_mean_bias_matter_at_z(virial_mass,indexz)

    use arrays
    use fiducial

    Implicit none

    Real*8 :: virial_mass,integrand_mean_bias_matter_at_z,dalpha
!    Real*8,dimension(number_of_M) :: f

    Integer*4 :: indexz

!    Do indexM=1,number_of_M

!       f(indexM) = dndM(indexM,indexz)*bMz(indexM,indexz)*M(indexM)/&
 !           mean_density(z(indexz))*(1.d0 + z(indexz))**3.d0

  !  End Do

    call Interpolate_1D(integrand_mean_bias_matter_at_z,dalpha,virial_mass,M,inte_mbz(:,indexz))

  End function integrand_mean_bias_matter_at_z

  Subroutine compute_mean_bias_matter_at_z(indexz,output)

    use fiducial 

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

    f_obj = fgsl_function_init(integrand_mean_bias_matter, ptr)

    wk = fgsl_integration_workspace_alloc(nmax)

    status = fgsl_integration_qags(f_obj, lower_limit, upper_limit, &
         absolute_error, relative_error, nmax, wk, result, error)

    output = result

    call fgsl_function_free(f_obj)

    call fgsl_integration_workspace_free(wk)

  End subroutine compute_mean_bias_matter_at_z

  subroutine compute_mean_bias_matter() ! Dimensionless

    use arrays
    use fiducial
    Implicit none

    Integer*4 :: indexz

    open(15,file='./precomputed_quantities/bMz/mean_bias_at_z.dat')

    write(15,*) '# red-shift        mean bias all matter'

    Do indexz=1,number_of_z

       call compute_mean_bias_matter_at_z(indexz,mbz(indexz))

       write(15,'(2es18.10)') z(indexz),mbz(indexz)

    End Do

    close(15)

  end subroutine compute_mean_bias_matter

  Function integrand_pre_cl_phiphi(virial_mass, params) bind(c)

    real(c_double), value :: virial_mass
    type(c_ptr), value :: params
    real(c_double) :: integrand_pre_cl_phiphi
    Integer(c_int), pointer :: pa(:)

    call c_f_pointer(params, pa,(/2/)) ! pa(1) = indexz, pa(2) = indexl

    integrand_pre_cl_phiphi =  integrand_pre_cl_phiphi_at_z_and_l(virial_mass,pa(1),pa(2))

  End function integrand_pre_cl_phiphi

  subroutine compute_integrand_pre_cl_phiphi_at_z_and_l()

    use arrays
    use fiducial

    Implicit none

    Integer*4 :: indexM,indexz,indexl

    Do indexl=1,number_of_l

       Do indexM=1,number_of_M

          Do indexz=1,number_of_z

             inte_pre_cl_phiphi(indexl,indexM,indexz) = dndM(indexM,indexz)*philMz(indexl,indexM,indexz)**2 ! Units : 1/solar mass/Mpc**3

          End Do

       End Do

    End Do
    
  end subroutine compute_integrand_pre_cl_phiphi_at_z_and_l

  function integrand_pre_cl_phiphi_at_z_and_l(virial_mass,indexz,indexl)

    use arrays
    !use fiducial

    Implicit none

    Real*8 :: virial_mass,integrand_pre_cl_phiphi_at_z_and_l,dalpha
    !Real*8,dimension(number_of_M) :: f

    Integer*4 :: indexz,indexl

    !Do indexM=1,number_of_M

    !   f(indexM) = dndM(indexM,indexz)*philMz(indexl,indexM,indexz)**2 ! Units : 1/solar mass/Mpc**3

    !End Do

    call Interpolate_1D(integrand_pre_cl_phiphi_at_z_and_l,dalpha,virial_mass,M,inte_pre_cl_phiphi(indexl,:,indexz))

  End function integrand_pre_cl_phiphi_at_z_and_l

  Subroutine compute_pre_cl_phiphi_at_z_and_l(indexz,indexl,output)

    use  fiducial 

    Implicit none

    Integer(fgsl_size_t), parameter :: nmax=1000000

    Integer(fgsl_int),target :: pp(2)
    Real(fgsl_double) :: result, error, output
    Real(fgsl_double),parameter :: lower_limit = Mmin!1.0E-3_fgsl_double
    Real(fgsl_double),parameter :: upper_limit = Mmax!1.0E-2_fgsl_double
    Real(fgsl_double),parameter :: absolute_error = 0.0_fgsl_double
    Real(fgsl_double),parameter :: relative_error = 1.0E-5_fgsl_double

    Integer(fgsl_int) :: status
    Integer(fgsl_int) :: indexz,indexl
    Integer(fgsl_int),parameter :: key = 4

    Type(c_ptr) :: ptr
    Type(fgsl_function) :: f_obj
    Type(fgsl_integration_workspace) :: wk

    pp = (/indexz,indexl/)

    ptr = c_loc(pp)

    f_obj = fgsl_function_init(integrand_pre_cl_phiphi, ptr)

    wk = fgsl_integration_workspace_alloc(nmax)

!    status = fgsl_integration_qags(f_obj, lower_limit, upper_limit, &
 !        absolute_error, relative_error, nmax, wk, result, error)
    status = fgsl_integration_qag(f_obj, lower_limit, upper_limit, &
         absolute_error, relative_error, nmax, key, wk, result, error)

    output = result

    call fgsl_function_free(f_obj)

    call fgsl_integration_workspace_free(wk)

  End subroutine compute_pre_cl_phiphi_at_z_and_l

  function integrand_cl_phiphi_one_halo_at_z_and_l(redshift,indexl)

    use arrays
    use fiducial

    Implicit none

    Real*8 :: redshift,integrand_cl_phiphi_one_halo_at_z_and_l,dalpha
    !Real*8,dimension(number_of_z) :: f

    Integer*4 :: indexl !indexl

    !Do indexz=1,number_of_z

    !   call compute_pre_cl_phiphi_at_z_and_l(indexz,indexl,f(indexz)) 

    !   f(indexz) = f(indexz)*d2VdzdO(indexz)

    !End Do

    call Interpolate_1D(integrand_cl_phiphi_one_halo_at_z_and_l,dalpha,redshift,z,inte_cl_phiphi_1h(indexl,:))

  End function integrand_cl_phiphi_one_halo_at_z_and_l

  subroutine compute_integrand_cl_phiphi_one_halo_at_z_and_l()

    use arrays
    use fiducial

    Implicit none

    Integer*4 :: indexl,indexz

    Do indexl=1,number_of_l

       Do indexz=1,number_of_z

          call compute_pre_cl_phiphi_at_z_and_l(indexz,indexl,inte_cl_phiphi_1h(indexl,indexz))

          inte_cl_phiphi_1h(indexl,indexz) = inte_cl_phiphi_1h(indexl,indexz)*d2VdzdO(indexz)

       End Do

    End Do

  end subroutine compute_integrand_cl_phiphi_one_halo_at_z_and_l

  Function integrand_cl_phiphi_one_halo(redshift, params) bind(c)

    real(c_double), value :: redshift
    type(c_ptr), value :: params
    real(c_double) :: integrand_cl_phiphi_one_halo
    Integer(c_int), pointer :: pa

    call c_f_pointer(params, pa) ! pa = indexl

    integrand_cl_phiphi_one_halo =  integrand_cl_phiphi_one_halo_at_z_and_l(redshift,pa)

  End function integrand_cl_phiphi_one_halo

  Subroutine compute_cl_phiphi_one_halo_at_l(indexl,output)

    Implicit none

    Integer(fgsl_size_t), parameter :: nmax=10000

    Integer(fgsl_int),target :: pp
    Real(fgsl_double) :: result, error, output
    Real(fgsl_double),parameter :: lower_limit = zmin!1.0E-3_fgsl_double
    Real(fgsl_double),parameter :: upper_limit = zmax!1.0E-2_fgsl_double
    Real(fgsl_double),parameter :: absolute_error = 0.0_fgsl_double
    Real(fgsl_double),parameter :: relative_error = 1.0E-5_fgsl_double

    Integer(fgsl_int) :: status
    Integer(fgsl_int) :: indexl

    Type(c_ptr) :: ptr
    Type(fgsl_function) :: f_obj
    Type(fgsl_integration_workspace) :: wk

    pp = indexl

    ptr = c_loc(pp)

    f_obj = fgsl_function_init(integrand_cl_phiphi_one_halo, ptr)

    wk = fgsl_integration_workspace_alloc(nmax)

    status = fgsl_integration_qags(f_obj, lower_limit, upper_limit, &
         absolute_error, relative_error, nmax, wk, result, error)

    output = result

    call fgsl_function_free(f_obj)

    call fgsl_integration_workspace_free(wk)

  End subroutine compute_cl_phiphi_one_halo_at_l

  subroutine compute_Clphiphi1h()  ! It fills lensing potential angular power spectrum array. Dimensionless

    use arrays
    use fiducial
    Implicit none

    Integer*4 :: indexl

    Do indexl=1,number_of_l

       call compute_cl_phiphi_one_halo_at_l(indexl,Clphiphi1h(indexl))

    End Do

  end subroutine compute_Clphiphi1h


!  function pre_Clphiphi(indexz,indexl) ! Required function to compute one halo term of lensing potential angular power spectrum. Dimensionless

 !   use fiducial
  !  use arrays
   ! Implicit none

    !Real*8 :: pre_Clphiphi,sum
!    Integer*4 :: indexl,i,indexM,indexz
 !   Integer*4,parameter :: intervals = number_of_M - 1 
  !  Real*8,dimension(number_of_M) :: f

   ! Do indexM = 1, number_of_M

    !   f(indexM) = dndM(indexM,indexz)*philMz(indexl,indexM,indexz)**2 ! Units : 1/solar mass/Mpc**3

!    End Do

 !   open(15,file='./output/preclphiphi.dat')

  !  write(15,*) '# Virial_ Mass[solar mass]  f_in_preclphiphi',z(indexz),ml(indexl)

    !        Do indexM=1,number_of_M

    !            write(15,'(2es18.10)') M(indexM), f(indexM)

    !        End Do

    !        close(15)

!    sum = 0.d0

 !   Do i=1,intervals

  !     sum = (M(i+1)- M(i))/2.d0*( f(i) + f(i+1) ) + sum   ! Units : 1/Mpc**3

   ! End Do

    !pre_Clphiphi = sum     ! Dimensionless                     

  !end function pre_Clphiphi

!  function C_l_phiphi_one_halo(indexl) ! It computes lensing potential angular power spectrum for a single multipole. Dimensionless

 !   use fiducial
  !  use arrays
    !use omp_lib
   ! Implicit none

    !Real*8 :: C_l_phiphi_one_halo,sum
!    Integer*4 :: indexl,indexz
 !   Integer*4,parameter :: number_of_redshift = number_of_z_functions 
  !  Integer*4,parameter :: intervals = number_of_redshift - 1
   ! Real*8,dimension(number_of_redshift):: f

!!$omp Parallel Do Shared(f)

    !Do indexz=1,number_of_redshift

!       f(indexz) = pre_Clphiphi(indexz,indexl)*d2VdzdO(indexz)          ! Dimensionless

 !   End Do

!!$omp End Parallel Do

    !    open(16,file='./output/clphiphi.dat')

    !    Do indexz=1,number_of_z

    !        write(16,'(2es18.10)') z(indexz),f(indexz) 

    !    End Do

    !    close(16)

  !  sum = 0.d0

   ! Do indexz=1,intervals

    !   sum = (z_functions(indexz+1) -  z_functions(indexz))/2.d0*( f(indexz) + f(indexz+1) ) + sum

!    End Do

 !   C_l_phiphi_one_halo = sum

  !end function C_l_phiphi_one_halo

  function integrand_limber_approximation_at_z_and_l(redshift,indexl)

    use arrays
    use fiducial
    use functions

    Implicit none

    Real*8 :: integrand_limber_approximation_at_z_and_l,dalpha,redshift

    Integer*4 :: indexl

!    prefactor = 4.d0/dble(ml(indexl))**2/(dble(ml(indexl))+1.d0)**2*9.d0/4.d0/&
 !        c**4*Hubble_parameter(0.d0)**4*Omega_m(0.d0)**2    !    Units : 

    !Do indexz=1,number_of_z_limber

     !  f(indexz) = (1.d0 + zlimber(indexz))**2*( com_dist_at_z_dec - comoving_distance(zlimber(indexz)) )**2/&
      !      com_dist_at_z_dec**2*c/Hubble_parameter(zlimber(indexz))*&
       !     matter_power_spectrum((dble(ml(indexl)) + 1.d0/2.d0)/comoving_distance(zlimber(indexz)),zlimber(indexz))

!    End Do

    call Interpolate_1D(integrand_limber_approximation_at_z_and_l,dalpha,redshift,zlimber,integrand_limber(:,indexl))

!    integrand_limber_approximation_at_z_and_l = integrand_limber_approximation_at_z_and_l!*prefactor/h**3

  End function integrand_limber_approximation_at_z_and_l

  Function integrand_limber_approximation(redshift, params) bind(c)

    real(c_double), value :: redshift
    type(c_ptr), value :: params
    real(c_double) :: integrand_limber_approximation
    Integer(c_int), pointer :: pa

    call c_f_pointer(params, pa) ! pa = indexl

    integrand_limber_approximation =  integrand_limber_approximation_at_z_and_l(redshift,pa)

  End function integrand_limber_approximation

  Subroutine compute_limber_approximation_at_l(indexl,output)

    use fiducial

    Implicit none

    Integer(fgsl_size_t), parameter :: nmax=10000

    Integer(fgsl_int),target :: pp
    Real(fgsl_double) :: result, error, output
    Real(fgsl_double),parameter :: lower_limit = 1.0E-5_fgsl_double
    Real(fgsl_double),parameter :: upper_limit = z_dec!1.0E-2_fgsl_double
    Real(fgsl_double),parameter :: absolute_error = 0.0_fgsl_double
    Real(fgsl_double),parameter :: relative_error = 1.0E-8_fgsl_double

    Integer(fgsl_int) :: status
    Integer(fgsl_int) :: indexl

    Type(c_ptr) :: ptr
    Type(fgsl_function) :: f_obj
    Type(fgsl_integration_workspace) :: wk

    pp = indexl

    ptr = c_loc(pp)

    f_obj = fgsl_function_init(integrand_limber_approximation, ptr)

    wk = fgsl_integration_workspace_alloc(nmax)

    status = fgsl_integration_qags(f_obj, lower_limit, upper_limit, &
         absolute_error, relative_error, nmax, wk, result, error)

    output = result

    call fgsl_function_free(f_obj)

    call fgsl_integration_workspace_free(wk)

  End subroutine compute_limber_approximation_at_l

  subroutine compute_Clpsilimber()   ! It fills array for angular power spectrum of lensing potential in the Limber approximations. 

    use arrays
    use fiducial
    use omp_lib

    Implicit none

    Integer*4 :: indexl

    !$omp Parallel Do Shared(Clpsilimber)

    Do indexl=1,number_of_l

       call compute_limber_approximation_at_l(indexl,Clpsilimber(indexl))

    End Do

    !$omp End Parallel Do

  end subroutine compute_Clpsilimber


!  function C_l_psi_limber(indexl)    ! It computes the angular power spectrum of the lensing potential in the Limber approximation 

 !   use fiducial    ! for a single multipole as in Lewis and Challinor. Units : Dimensionless
  !  use arrays
   ! Implicit none

!    Real*8 :: C_l_psi_limber,sum,prefactor
 !   Integer*4 :: indexl,indexz
  !  Integer*4,parameter :: intervals = number_of_z_limber - 1 
   ! Real*8,dimension(number_of_z_limber):: f,zlimber
    !Real*8,parameter :: zminlimber = 1.d-5

    ! Red-shift array. Dimensionless.

!    Do indexz = 1, number_of_z_limber      

 !      zlimber(indexz) = 10**(log10(zminlimber) + real(indexz-1)*(log10(z_dec) - log10(zminlimber))/real(number_of_z_limber-1))

  !  End Do

   ! prefactor = 4.d0/dble(ml(indexl))**2/(dble(ml(indexl))+1.d0)**2*9.d0/4.d0/&
    !     c**4*Hubble_parameter(0.d0)**4*Omega_m(0.d0)**2    !    Units : 

!    Do indexz=1,number_of_z_limber

 !      f(indexz) = (1.d0 + zlimber(indexz))**2*( com_dist_at_z_dec - comoving_distance(zlimber(indexz)) )**2/&
  !          com_dist_at_z_dec**2*c/Hubble_parameter(zlimber(indexz))*&
   !         matter_power_spectrum((dble(ml(indexl)) + 1.d0/2.d0)/comoving_distance(zlimber(indexz)),zlimber(indexz))

    !End Do

!    sum = 0.d0

 !   Do indexz=1,intervals

  !     sum = (zlimber(indexz+1)-zlimber(indexz))/2.d0*( f(indexz) + f(indexz+1) ) + sum

   ! End Do

    !C_l_psi_limber = prefactor*sum/h**3    

!  end function C_l_psi_limber

!  function pre_Cl_1(indexz,indexl) ! Units 1/Mpc**3. 

 !   use fiducial
    !use omp_lib
  !  use arrays
   ! Implicit none

!    Real*8 :: pre_Cl_1,sum
 !   Integer*4 :: indexl,indexM,indexz
  !  Integer*4,parameter :: intervals = number_of_M - 1
   ! Real*8,dimension(number_of_M) :: f

!    Do indexM = 1, number_of_M    

 !      f(indexM) = dndM(indexM,indexz)*bMz(indexM,indexz)*philMz(indexl,indexM,indexz) 

  !  End Do

    !    open(16,file='./output/clphiphi.dat')

    !    Do indexM=1,number_of_M

    !        write(16,'(2es18.10)') M(indexM),f(indexM) 

    !    End Do

    !    close(16)

   ! sum = 0.d0

    !Do indexM=1,intervals

!       sum = (M(indexM+1)-M(indexM))/2.d0*( f(indexM) + f(indexM+1) ) + sum

 !   End Do

  !  pre_Cl_1 = sum

!  end function pre_Cl_1

  subroutine compute_integrand_pre_cl_phiphi_2h_at_z_and_l()

    use arrays
    use fiducial

    Implicit none

    Integer*4 :: indexl,indexM,indexz

    Do indexl=1,number_of_l

       Do indexM=1,number_of_M

          Do indexz=1,number_of_z

             inte_pre_cl_phiphi_2h(indexl,indexM,indexz) = dndM(indexM,indexz)*bMz(indexM,indexz)*philMz(indexl,indexM,indexz)  

          End Do

       End Do

    End Do

  end subroutine compute_integrand_pre_cl_phiphi_2h_at_z_and_l

  function integrand_pre_cl_phiphi_2h_at_z_and_l(virial_mass,indexz,indexl)

    use arrays
    !use fiducial

    Implicit none

    Real*8 :: virial_mass,integrand_pre_cl_phiphi_2h_at_z_and_l,dalpha
    !Real*8,dimension(number_of_M) :: f

    Integer*4 :: indexz,indexl

    !Do indexM=1,number_of_M

    !   f(indexM) = dndM(indexM,indexz)*bMz(indexM,indexz)*philMz(indexl,indexM,indexz)  

    !End Do

    call Interpolate_1D(integrand_pre_cl_phiphi_2h_at_z_and_l,dalpha,virial_mass,M,inte_pre_cl_phiphi_2h(indexl,:,indexz))

  End function integrand_pre_cl_phiphi_2h_at_z_and_l

  Function integrand_pre_cl_phiphi_2h(virial_mass, params) bind(c)

    real(c_double), value :: virial_mass
    type(c_ptr), value :: params
    real(c_double) :: integrand_pre_cl_phiphi_2h
    Integer(c_int), pointer :: pa(:)

    call c_f_pointer(params, pa,(/2/)) ! pa(1) = indexz, pa(2) = indexl

    integrand_pre_cl_phiphi_2h =  integrand_pre_cl_phiphi_2h_at_z_and_l(virial_mass,pa(1),pa(2))

  End function integrand_pre_cl_phiphi_2h

  Subroutine compute_pre_cl_phiphi_2h_at_z_and_l(indexz,indexl,output)

    use fiducial 

    Implicit none

    Integer(fgsl_size_t), parameter :: nmax=1000000

    Integer(fgsl_int),target :: pp(2)
    Real(fgsl_double) :: result, error, output
    Real(fgsl_double),parameter :: lower_limit = Mmin!1.0E-3_fgsl_double
    Real(fgsl_double),parameter :: upper_limit = Mmax!1.0E-2_fgsl_double
    Real(fgsl_double),parameter :: absolute_error = 0.0_fgsl_double
    Real(fgsl_double),parameter :: relative_error = 1.0E-5_fgsl_double

    Integer(fgsl_int) :: status
    Integer(fgsl_int) :: indexz,indexl
    Integer(fgsl_int),parameter :: key = 4

    Type(c_ptr) :: ptr
    Type(fgsl_function) :: f_obj
    Type(fgsl_integration_workspace) :: wk

    pp = (/indexz,indexl/)

    ptr = c_loc(pp)

    f_obj = fgsl_function_init(integrand_pre_cl_phiphi_2h, ptr)

    wk = fgsl_integration_workspace_alloc(nmax)

!    status = fgsl_integration_qags(f_obj, lower_limit, upper_limit, &
 !        absolute_error, relative_error, nmax, wk, result, error)

    status = fgsl_integration_qag(f_obj, lower_limit, upper_limit, &
         absolute_error, relative_error, nmax, key, wk, result, error)

    output = result

    call fgsl_function_free(f_obj)

    call fgsl_integration_workspace_free(wk)

  End subroutine compute_pre_cl_phiphi_2h_at_z_and_l

  subroutine compute_integrand_cl_phiphi_two_halo_at_z_and_l()

    use arrays 
    use fiducial

    Implicit none 

    Integer*4 :: indexz,indexl

    Do indexl=1,number_of_l

       Do indexz=1,number_of_z

          call compute_pre_cl_phiphi_2h_at_z_and_l(indexz,indexl,inte_cl_phiphi_2h(indexl,indexz))

          inte_cl_phiphi_2h(indexl,indexz) = inte_cl_phiphi_2h(indexl,indexz)**2*d2VdzdO(indexz)*&
            matter_power_spectrum((dble(ml(indexl))+1.d0/2.d0)/comoving_distance_at_z(indexz),z(indexz))    !  Dimensionless
          
       End Do

    End Do

  end subroutine compute_integrand_cl_phiphi_two_halo_at_z_and_l

  function integrand_cl_phiphi_two_halo_at_z_and_l(redshift,indexl)

    use arrays
    use fiducial
    use functions

    Implicit none

    Real*8 :: redshift,integrand_cl_phiphi_two_halo_at_z_and_l,dalpha
    !Real*8,dimension(number_of_z) :: f

    Integer*4 :: indexl!,indexl

    !Do indexz=1,number_of_z

    !   call compute_pre_cl_phiphi_2h_at_z_and_l(indexz,indexl,f(indexz)) 

    !   f(indexz) = f(indexz)**2*d2VdzdO(indexz)*&
     !       matter_power_spectrum((dble(ml(indexl))+1.d0/2.d0)/comoving_distance_at_z(indexz),z(indexz))    !  Dimensionless

    !End Do

    call Interpolate_1D(integrand_cl_phiphi_two_halo_at_z_and_l,dalpha,redshift,z,inte_cl_phiphi_2h(indexl,:))

  End function integrand_cl_phiphi_two_halo_at_z_and_l

  Function integrand_cl_phiphi_two_halo(redshift, params) bind(c)

    real(c_double), value :: redshift
    type(c_ptr), value :: params
    real(c_double) :: integrand_cl_phiphi_two_halo
    Integer(c_int), pointer :: pa

    call c_f_pointer(params, pa) ! pa = indexl

    integrand_cl_phiphi_two_halo =  integrand_cl_phiphi_two_halo_at_z_and_l(redshift,pa)

  End function integrand_cl_phiphi_two_halo

  Subroutine compute_cl_phiphi_two_halo_at_l(indexl,output)

    Implicit none

    Integer(fgsl_size_t), parameter :: nmax=1000000

    Integer(fgsl_int),target :: pp
    Real(fgsl_double) :: result, error, output
    Real(fgsl_double),parameter :: lower_limit = zmin!1.0E-3_fgsl_double
    Real(fgsl_double),parameter :: upper_limit = zmax!1.0E-2_fgsl_double
    Real(fgsl_double),parameter :: absolute_error = 0.0_fgsl_double
    Real(fgsl_double),parameter :: relative_error = 1.0E-5_fgsl_double

    Integer(fgsl_int) :: status
    Integer(fgsl_int) :: indexl
!    Integer(fgsl_int),parameter :: key = 4

    Type(c_ptr) :: ptr
    Type(fgsl_function) :: f_obj
    Type(fgsl_integration_workspace) :: wk

    pp = indexl

    ptr = c_loc(pp)

    f_obj = fgsl_function_init(integrand_cl_phiphi_two_halo, ptr)

    wk = fgsl_integration_workspace_alloc(nmax)

    status = fgsl_integration_qags(f_obj, lower_limit, upper_limit, &
         absolute_error, relative_error, nmax, wk, result, error)
!    status = fgsl_integration_qag(f_obj, lower_limit, upper_limit, &
 !        absolute_error, relative_error, nmax, key, wk, result, error)


    output = result/h**3

!    print *, output, error, indexl

    call fgsl_function_free(f_obj)

    call fgsl_integration_workspace_free(wk)

  End subroutine compute_cl_phiphi_two_halo_at_l

  subroutine compute_Clphiphi2h()   ! It fills two halo term of lensing potential angular power spectrum array. 

    use arrays
    use fiducial
    Implicit none

    Integer*4 :: indexl

    Do indexl=1,number_of_l

       call compute_cl_phiphi_two_halo_at_l(indexl,Clphiphi2h(indexl))

    End Do

  end subroutine compute_Clphiphi2h

!  function C_l_phiphi_two_halo(indexl)    ! It computes two halo term lensing potential angular power spectrum for a single multipole. 

 !   use fiducial    ! Equation (2.11) in 1312.4525. Dimensionless
  !  use arrays
    !        use omp_lib
   ! Implicit none

!    Real*8 :: C_l_phiphi_two_halo,sum
 !   Integer*4 :: indexl,indexz
  !  Integer*4,parameter :: number_of_redshift = number_of_z!_functions  !1d4
   ! Integer*4,parameter :: intervals = number_of_redshift - 1 
    !Real*8,dimension(number_of_redshift):: f!,redshift

    !      Do indexz=1,number_of_redshift

    !            redshift(indexz) = 10**(log10(zmin) + real(indexz-1)*(log10(zmax) - log10(zmin))/real(number_of_redshift-1))

    !        End Do

    !        !$omp Parallel Do Default(Shared)
!    Do indexz=1,number_of_redshift

 !      f(indexz) = pre_Cl_1(indexz,indexl)**2*d2VdzdO(indexz)*&
  !          matter_power_spectrum((dble(ml(indexl))+1.d0/2.d0)/comoving_distance_at_z(indexz),z(indexz))    !  Dimensionless

   ! End Do
    !        !$omp End Parallel Do

    !    open(16,file='./output/clphiphi.dat')

    !    Do indexz=1,number_of_z

    !        write(16,'(2es18.10)') z(indexz),f(indexz) 

    !    End Do

    !    close(16)

!    sum = 0.d0

 !   Do indexz=1,intervals

  !     sum = (z(indexz+1)-z(indexz))/2.d0*( f(indexz) + f(indexz+1) ) + sum

   ! End Do

    !C_l_phiphi_two_halo = sum/h**3

!  end function C_l_phiphi_two_halo

  function integrand_form_factor_at_M_z_and_l(x,indexM,indexz,indexl)

    use arrays
    use fiducial
    use functions

    Implicit none

    Real*8 :: x,integrand_form_factor_at_M_z_and_l,l_s

    Integer*4 :: indexz,indexM,indexl

    l_s = angular_diameter_distance_at_z(indexz)/r200c(indexM,indexz)         ! dimensionless

    integrand_form_factor_at_M_z_and_l = x**2*sin((dble(ml(indexl))+1.d0/2.d0)*x/l_s)&
         *ICM_electron_pressure(x*r200c(indexM,indexz),indexM,indexz)/(dble(ml(indexl))+1.d0/2.d0)/x*l_s

  End function integrand_form_factor_at_M_z_and_l

  Function integrand_form_factor(x, params) bind(c)

    real(c_double), value :: x
    type(c_ptr), value :: params
    real(c_double) :: integrand_form_factor
    Integer(c_int), pointer :: pa(:)

    call c_f_pointer(params, pa,(/3/)) ! pa(1) = indexM, pa(2) = indexz, pa(3) = indexl

    integrand_form_factor =  integrand_form_factor_at_M_z_and_l(x,pa(1),pa(2),pa(3))

  End function integrand_form_factor

  Subroutine compute_form_factor_at_M_z_and_l(indexM,indexz,indexl,output)

    use arrays
    use functions
    use fiducial

    Implicit none

    Integer(fgsl_size_t), parameter :: nmax=10000

    Integer(fgsl_int),target :: pp(3)

    Real(fgsl_double) :: result, error, output, prefactor, l_s
    Real(fgsl_double),parameter :: lower_limit = 0.0_fgsl_double
    Real(fgsl_double) :: upper_limit !! Mmax!1.0E-2_fgsl_double
    Real(fgsl_double),parameter :: absolute_error = 0.0_fgsl_double
    Real(fgsl_double),parameter :: relative_error = 1.0E-5_fgsl_double
    Real*8,parameter :: betafactor = 1.060d0!1.037d0

    Integer(fgsl_int) :: status
    Integer(fgsl_int) :: indexz,indexl,indexM

    Type(c_ptr) :: ptr
    Type(fgsl_function) :: f_obj
    Type(fgsl_integration_workspace) :: wk

    l_s = angular_diameter_distance_at_z(indexz)/r200c(indexM,indexz)         ! dimensionless

    prefactor = sigma_e*4.d0*Pi*r200c(indexM,indexz)*M_sun/m_e/c**2/l_s**2  ! Units : Mpc*s**2/solar mass

    upper_limit = betafactor*virial_radius(z(indexz),M(indexM))/r200c(indexM,indexz) ! dimensionless

    pp = (/indexM,indexz,indexl/)

    ptr = c_loc(pp)

    f_obj = fgsl_function_init(integrand_form_factor, ptr)

    wk = fgsl_integration_workspace_alloc(nmax)

    status = fgsl_integration_qags(f_obj, lower_limit, upper_limit, &
         absolute_error, relative_error, nmax, wk, result, error)

    output = prefactor*result

    call fgsl_function_free(f_obj)

    call fgsl_integration_workspace_free(wk)

  End subroutine compute_form_factor_at_M_z_and_l

  subroutine compute_form_factor()

    use arrays
    use fiducial
    !use omp_lib
    Implicit none

    Integer*4 :: indexl,indexM,indexz

!!$omp Parallel Do Shared(ylMz)

    Do indexl=1,number_of_l

       Do indexM=1,number_of_M

          Do indexz=1,number_of_z

             call compute_form_factor_at_M_z_and_l(indexM,indexz,indexl,ylMz(indexl,indexM,indexz))

          End Do

       End Do

    End Do

!!$omp End Parallel Do

    open(15,file='./precomputed_quantities/form_factor/form_factor.dat')

    write(15,*) '# Form factor file. Number of lines is ',number_of_l*number_of_z*number_of_M

    write(15,*) '# index_of_l    l    index_of_M    virial_mass[solar mass]    index_of_z    red-shift    y'

    Do indexl=1,number_of_l

       Do indexM=1,number_of_M

          Do indexz=1,number_of_z

             write(15,'(3i10,es18.10,i5,2es18.10)') indexl,ml(indexl),indexM,M(indexM),indexz,z(indexz),&
                  ylMz(indexl,indexM,indexz)

          End Do

       End Do

    End Do

    close(15)

  end subroutine compute_form_factor

!  function pre_Cl(indexz,indexl) ! Dimensionless

 !   use fiducial
  !  use arrays
   ! Implicit none

!    Real*8 :: pre_Cl,sum!,dndM_M_z,ylMz_M_z,philMz_M_z
 !   Integer*4 :: indexl,i,indexM,indexz
  !  Integer*4,parameter :: number_of_virial_Mass = number_of_M_functions 
   ! Integer*4,parameter :: intervals = number_of_virial_Mass - 1 
    !Real*8,dimension(number_of_virial_Mass) :: f!,virial_Mass

!    Do indexM = 1, number_of_virial_Mass    

       !            virial_Mass(indexM) = 10**(log10(Mmin) + real(indexM-1)*(log10(Mmax) - log10(Mmin))/real(number_of_virial_Mass-1))

       !            call Interpolate_2D(dndM_M_z,virial_Mass(indexM),redshift,M_functions(1:number_of_M_functions),&
       !                z_functions(1:number_of_z_functions),dndM(1:number_of_M_functions,1:number_of_z_functions))

       !          call Interpolate_2D(philMz_M_z,virial_Mass(indexM),redshift,M_functions(1:number_of_M_functions),&
       !              z_functions(1:number_of_z_functions),philMz(indexl,1:number_of_M_functions,1:number_of_z_functions))

       !        call Interpolate_2D(ylMz_M_z,virial_Mass(indexM),redshift,M_functions(1:number_of_M_functions),&
       !            z_functions(1:number_of_z_functions),ylMz(indexl,1:number_of_M_functions,1:number_of_z_functions))

 !      f(indexM) = dndM_interpolation(indexM,indexz)*ylMz_interpolation(indexl,indexM,indexz)*&
  !          philMz_interpolation(indexl,indexM,indexz) ! Units : 1/solar mass/Mpc**3

   ! End Do

    !    open(16,file='./output/clphiphi.dat')

    !    Do indexM=1,number_of_M

    !        write(16,'(2es18.10)') M(indexM),f(indexM) 

    !    End Do

    !    close(16)

    !sum = 0.d0

!    Do i=1,intervals

 !      sum = (M_functions(i+1)- M_functions(i))/2.d0*( f(i) + f(i+1) ) + sum 

  !  End Do

   ! pre_Cl = sum       ! dimensionless                   

!  end function pre_Cl

  subroutine compute_integrand_pre_cl_yphi_at_z_and_l()

    use arrays
    use fiducial

    Implicit none

    Integer*4 :: indexl,indexz,indexM

    Do indexl=1,number_of_l

       Do indexM=1,number_of_M

          Do indexz=1,number_of_z

             inte_pre_cl_yphi(indexl,indexM,indexz) = dndM(indexM,indexz)*ylMz(indexl,indexM,indexz)*philMz(indexl,indexM,indexz) ! Units : 1/solar mass/Mpc**3 

          End Do

       End Do

    End Do

  end subroutine compute_integrand_pre_cl_yphi_at_z_and_l

  function integrand_pre_cl_yphi_at_z_and_l(virial_mass,indexz,indexl)

    use arrays
    use fiducial

    Implicit none

    Real*8 :: virial_mass,integrand_pre_cl_yphi_at_z_and_l,dalpha
!    Real*8,dimension(number_of_M) :: f

    Integer*4 :: indexz,indexl!,indexl

!    Do indexM=1,number_of_M

 !      f(indexM) = dndM(indexM,indexz)*ylMz(indexl,indexM,indexz)*philMz(indexl,indexM,indexz) ! Units : 1/solar mass/Mpc**3 

  !  End Do

    call Interpolate_1D(integrand_pre_cl_yphi_at_z_and_l,dalpha,virial_mass,M,inte_pre_cl_yphi(indexl,:,indexz))

  End function integrand_pre_cl_yphi_at_z_and_l

  Function integrand_pre_cl_yphi(virial_mass, params) bind(c)

    real(c_double), value :: virial_mass
    type(c_ptr), value :: params
    real(c_double) :: integrand_pre_cl_yphi
    Integer(c_int), pointer :: pa(:)

    call c_f_pointer(params, pa,(/2/)) ! pa(1) = indexz, pa(2) = indexl

    integrand_pre_cl_yphi =  integrand_pre_cl_yphi_at_z_and_l(virial_mass,pa(1),pa(2))

  End function integrand_pre_cl_yphi

  Subroutine compute_pre_cl_yphi_at_z_and_l(indexz,indexl,output)

    use fiducial 

    Implicit none

    Integer(fgsl_size_t), parameter :: nmax=10000

    Integer(fgsl_int),target :: pp(2)
    Real(fgsl_double) :: result, error, output
    Real(fgsl_double),parameter :: lower_limit = Mmin!1.0E-3_fgsl_double
    Real(fgsl_double),parameter :: upper_limit = Mmax!1.0E-2_fgsl_double
    Real(fgsl_double),parameter :: absolute_error = 0.0_fgsl_double
    Real(fgsl_double),parameter :: relative_error = 1.0E-5_fgsl_double

    Integer(fgsl_int) :: status
    Integer(fgsl_int) :: indexz,indexl

    Type(c_ptr) :: ptr
    Type(fgsl_function) :: f_obj
    Type(fgsl_integration_workspace) :: wk

    pp = (/indexz,indexl/)

    ptr = c_loc(pp)

    f_obj = fgsl_function_init(integrand_pre_cl_yphi, ptr)

    wk = fgsl_integration_workspace_alloc(nmax)

    status = fgsl_integration_qags(f_obj, lower_limit, upper_limit, &
         absolute_error, relative_error, nmax, wk, result, error)

    output = result

    call fgsl_function_free(f_obj)

    call fgsl_integration_workspace_free(wk)

  End subroutine compute_pre_cl_yphi_at_z_and_l

  subroutine compute_integrand_cl_yphi_one_halo_at_z_and_l()

    use arrays
    use fiducial

    Implicit none

    Integer*4 :: indexz, indexl

    Do indexl=1,number_of_l

       Do indexz=1,number_of_z

          call compute_pre_cl_yphi_at_z_and_l(indexz,indexl,inte_cl_yphi_1h(indexl,indexz))

          inte_cl_yphi_1h(indexl,indexz)  = inte_cl_yphi_1h(indexl,indexz)*d2VdzdO(indexz)

       End Do

    End Do

  end subroutine compute_integrand_cl_yphi_one_halo_at_z_and_l

  function integrand_cl_yphi_one_halo_at_z_and_l(redshift,indexl)

    use arrays
    use fiducial

    Implicit none

    Real*8 :: redshift,integrand_cl_yphi_one_halo_at_z_and_l,dalpha
!    Real*8,dimension(number_of_z) :: f

    Integer*4 :: indexl!,indexl

!    Do indexz=1,number_of_z

 !      call compute_pre_cl_yphi_at_z_and_l(indexz,indexl,f(indexz)) 

  !     f(indexz) = f(indexz)*d2VdzdO(indexz)

   ! End Do

    call Interpolate_1D(integrand_cl_yphi_one_halo_at_z_and_l,dalpha,redshift,z,inte_cl_yphi_1h(indexl,:))

  End function integrand_cl_yphi_one_halo_at_z_and_l

  Function integrand_cl_yphi_one_halo(redshift, params) bind(c)

    real(c_double), value :: redshift
    type(c_ptr), value :: params
    real(c_double) :: integrand_cl_yphi_one_halo
    Integer(c_int), pointer :: pa

    call c_f_pointer(params, pa) ! pa = indexl

    integrand_cl_yphi_one_halo =  integrand_cl_yphi_one_halo_at_z_and_l(redshift,pa)

  End function integrand_cl_yphi_one_halo

  Subroutine compute_cl_yphi_one_halo_at_l(indexl,output)

    Implicit none

    Integer(fgsl_size_t), parameter :: nmax=10000

    Integer(fgsl_int),target :: pp
    Real(fgsl_double) :: result, error, output
    Real(fgsl_double),parameter :: lower_limit = zmin!1.0E-3_fgsl_double
    Real(fgsl_double),parameter :: upper_limit = zmax!1.0E-2_fgsl_double
    Real(fgsl_double),parameter :: absolute_error = 0.0_fgsl_double
    Real(fgsl_double),parameter :: relative_error = 1.0E-5_fgsl_double

    Integer(fgsl_int) :: status
    Integer(fgsl_int) :: indexl

    Type(c_ptr) :: ptr
    Type(fgsl_function) :: f_obj
    Type(fgsl_integration_workspace) :: wk

    pp = indexl

    ptr = c_loc(pp)

    f_obj = fgsl_function_init(integrand_cl_yphi_one_halo, ptr)

    wk = fgsl_integration_workspace_alloc(nmax)

    status = fgsl_integration_qags(f_obj, lower_limit, upper_limit, &
         absolute_error, relative_error, nmax, wk, result, error)

    output = result

    call fgsl_function_free(f_obj)

    call fgsl_integration_workspace_free(wk)

  End subroutine compute_cl_yphi_one_halo_at_l

  subroutine compute_Clyphi1h()  ! It fills lensing potential angular power spectrum array. Dimensionless

    use arrays
    use fiducial
    Implicit none

    Integer*4 :: indexl

    Do indexl=1,number_of_l

       call compute_cl_yphi_one_halo_at_l(indexl,Cl1h(indexl))

    End Do

  end subroutine compute_Clyphi1h



  !function C_l_yphi_one_halo(indexl) ! Dimensionless

   ! use fiducial
    !use arrays
!    use omp_lib
 !   Implicit none

  !  Real*8 :: C_l_yphi_one_halo,sum
   ! Integer*4 :: indexl,indexz
    !Integer*4,parameter :: number_of_redshift = number_of_z_functions 
!    Integer*4,parameter :: intervals = number_of_redshift - 1
 !   Real*8,dimension(number_of_redshift):: f

  !  !$omp Parallel Do Shared(f)

   ! Do indexz=1,number_of_redshift

    !   f(indexz) = pre_Cl(indexz,indexl)*d2VdzdO(indexz)          ! Dimensionless

!    End Do

 !   !$omp End Parallel Do

    !    open(16,file='./output/clphiphi.dat')

    !    Do indexz=1,number_of_z

    !        write(16,'(2es18.10)') z(indexz),f(indexz) 

    !    End Do

    !    close(16)

  !  sum = 0.d0

   ! Do indexz=1,intervals

    !   sum = (z_functions(indexz+1) -  z_functions(indexz))/2.d0*( f(indexz) + f(indexz+1) ) + sum

!    End Do

 !   C_l_yphi_one_halo = sum

  !end function C_l_yphi_one_halo

!  subroutine compute_Cl1h()  ! Dimensionless

 !   use arrays
  !  use fiducial
   ! Implicit none

!    Integer*4 :: indexl

 !   If (compute_functions) then

  !     print *, 'ONE HALO TERM OF CROSS-CORRELATION YPHI ONLY IS COMPUTED WHEN INTERPOLATING FUNCTIONS USED'

   !    stop

    !Else

     !  Do indexl=1,number_of_l

      !    Cl1h(indexl) = C_l_yphi_one_halo(indexl)

!       End Do

 !   End If

!  end subroutine compute_Cl1h

!  function form_factor(indexM,indexz,indexl)    ! Form factor. Equation (2.9) in 1312.4525. Units of M: solar mass. 

 !   use fiducial                              ! Units : dimensionless
  !  use arrays
   ! Implicit none

!    Real*8 :: l_s,x_y_min,x_y_max,form_factor,y,prefactor,stepsize,x1,x2,f1,f2!,virial_Mass,redshift,M200c_M_z,r200c_M_z
 !   Integer*4,parameter :: number_of_x = 10000
  !  Integer*4,parameter :: intervals = number_of_x - 1 !  
   ! Integer*4 :: indexx,indexl,indexM,indexz
!    Real*8,dimension(number_of_x) :: x,f
 !   Real*8,parameter :: betafactor = 1.060d0!1.037d0
  !  Integer*4,parameter :: max_iterations = 1000000000
   ! logical :: logscale

    !        call Interpolate_2D(M200c_M_z,virial_Mass,redshift,M(1:number_of_M),z(1:number_of_z),M200c(1:number_of_M,1:number_of_z))

    !        call Interpolate_2D(r200c_M_z,virial_Mass,redshift,M(1:number_of_M),z(1:number_of_z),r200c(1:number_of_M,1:number_of_z))

!    l_s = angular_diameter_distance_at_z(indexz)/r200c(indexM,indexz)         ! dimensionless
    !        l_s = angular_diameter_distance(redshift)/r200c_M_z         ! dimensionless

 !   x_y_max = betafactor*virial_radius(z(indexz),M(indexM))/r200c(indexM,indexz) ! dimensionless
    !        x_y_max = betafactor*virial_radius(redshift,virial_Mass)/r200c_M_z ! dimensionless

  !  If ( x_y_max .le. ( l_s/(dble(ml(indexl))+ 1.d0/2.d0) ) ) then 

   !    logscale = .true.

!    Else

 !      logscale = .false.

  !     stepsize = l_s/(dble(ml(indexl))+ 1.d0/2.d0)*1.d-3

   ! End If

!    prefactor = sigma_e*4.d0*Pi*r200c(indexM,indexz)*M_sun/m_e/c**2/l_s**2  ! Units : Mpc*s**2/solar mass
    !        prefactor = sigma_e*4.d0*Pi*r200c_M_z*M_sun/m_e/c**2/l_s**2  ! Units : Mpc*s**2/solar mass

 !   If (logscale) then

  !     x_y_min = 1.d-5*virial_radius(z(indexz),M(indexM))/r200c(indexM,indexz) ! dimensionless
       !            x_y_min = 1.d-5*virial_radius(redshift,virial_Mass)/r200c_M_z ! dimensionless

   !    Do indexx=1,number_of_x

    !      x(indexx) = 10**(log10(x_y_min) + dble((indexx-1))*(log10(x_y_max) &
     !          - log10(x_y_min))/dble(number_of_x-1)   )

      ! End Do

       !Do indexx=1,number_of_x 

        !  f(indexx) = x(indexx)**2*sin((dble(ml(indexl))+1.d0/2.d0)*x(indexx)/l_s)&
         !      *ICM_electron_pressure(x(indexx)*r200c(indexM,indexz),indexM,indexz)/(dble(ml(indexl))+1.d0/2.d0)/x(indexx)*l_s
          !                f(indexx) = x(indexx)**2*sin((dble(ml(indexl))+1.d0/2.d0)*x(indexx)/l_s)&
          !                *ICM_electron_pressure(x(indexx)*r200c_M_z,virial_Mass,redshift)/(dble(ml(indexl))+1.d0/2.d0)/x(indexx)*l_s

!       End Do

 !      y = 0.d0 

  !     Do indexx=1,intervals

   !       y = (x(indexx+1)-x(indexx))/2.d0*( f(indexx) + f(indexx+1))+ y

    !   End Do

!    Else

 !      x_y_min = 0.d0 ! dimensionless

  !     x1 = 0.d0

   !    f1 = 0.d0

    !   x2 = x1 + stepsize

     !  f2 = x2**2*sin((dble(ml(indexl))+1.d0/2.d0)*x2/l_s)&
      !      *ICM_electron_pressure(x2*r200c(indexM,indexz),indexM,indexz)/(dble(ml(indexl))+1.d0/2.d0)/x2*l_s
       !            f2 = x2**2*sin((dble(ml(indexl))+1.d0/2.d0)*x2/l_s)&
       !            *ICM_electron_pressure(x2*r200c_M_z,virial_Mass,redshift)/(dble(ml(indexl))+1.d0/2.d0)/x2*l_s

       !y = 0.d0

       !y = (x2-x1)/2.d0*( f1 + f2 )+ y

!       Do indexx=1,max_iterations

 !         x1 = x2

  !        f1 = f2

   !       If (x2 .gt. x_y_max) then

    !         x2 = x_y_max

     !        f2 = x2**2*sin((dble(ml(indexl))+1.d0/2.d0)*x2/l_s)&
      !            *ICM_electron_pressure(x2*r200c(indexM,indexz),indexM,indexz)/(dble(ml(indexl))+1.d0/2.d0)/x2*l_s
             !                    f2 = x2**2*sin((dble(ml(indexl))+1.d0/2.d0)*x2/l_s)&
             !                    *ICM_electron_pressure(x2*r200c_M_z,virial_Mass,redshift)/(dble(ml(indexl))+1.d0/2.d0)/x2*l_s

       !   Else If (x2 .eq. x_y_max) then

!             exit

 !         Else

  !           x2 = x1 + stepsize

   !          f2 = x2**2*sin((dble(ml(indexl))+1.d0/2.d0)*x2/l_s)&
    !              *ICM_electron_pressure(x2*r200c(indexM,indexz),indexM,indexz)/(dble(ml(indexl))+1.d0/2.d0)/x2*l_s
             !                    f2 = x2**2*sin((dble(ml(indexl))+1.d0/2.d0)*x2/l_s)&
             !                   *ICM_electron_pressure(x2*r200c_M_z,virial_Mass,redshift)/(dble(ml(indexl))+1.d0/2.d0)/x2*l_s

     !     End IF

      !    y = (x2-x1)/2.d0*( f1 + f2 )+ y

       !   If (indexx .eq. max_iterations) then

        !     print *,'Maximum number of iterations in integral of form factor achieved.'

         !    stop

!          End If

 !      End Do

  !  End If

   ! form_factor = prefactor*y

!  end function form_factor

!  function pre_Cl_2(indexz,indexl) ! Units 1/Mpc**3

 !   use fiducial
  !  use arrays
   ! use omp_lib
    !Implicit none

!    Real*8 :: pre_Cl_2,sum!,redshift,dndM_M_z,bMz_M_z,ylMz_M_z
 !   Integer*4 :: indexl,indexM,indexz
  !  Integer*4,parameter :: number_of_virial_Mass = number_of_M_functions !1d5
   ! Integer*4,parameter :: intervals = number_of_virial_Mass - 1
    !Real*8,dimension(number_of_virial_Mass) :: f!,virial_Mass

!    Do indexM = 1, number_of_virial_Mass    

       !            virial_Mass(indexM) = 10**(log10(Mmin) + real(indexM-1)*(log10(Mmax) - log10(Mmin))/real(number_of_virial_Mass-1))

       !           call Interpolate_2D(dndM_M_z,virial_Mass(indexM),redshift,M_functions(1:number_of_M_functions),&
       !               z_functions(1:number_of_z_functions),dndM(1:number_of_M_functions,1:number_of_z_functions))

       !         call Interpolate_2D(bMz_M_z,virial_Mass(indexM),redshift,M_functions(1:number_of_M_functions),&
       !             z_functions(1:number_of_z_functions),bMz(1:number_of_M_functions,1:number_of_z_functions))

       !       call Interpolate_2D(ylMz_M_z,virial_Mass(indexM),redshift,M_functions(1:number_of_M_functions),&
       !           z_functions(1:number_of_z_functions),ylMz(indexl,1:number_of_M_functions,1:number_of_z_functions))

 !      f(indexM) = dndM_interpolation(indexM,indexz)*bMz_interpolation(indexM,indexz)*&
  !          ylMz_interpolation(indexl,indexM,indexz)

   ! End Do


    !    open(16,file='./output/clphiphi.dat')

    !    Do indexM=1,number_of_M

    !        write(16,'(2es18.10)') M(indexM),f(indexM) 

    !    End Do

    !    close(16)

   ! sum = 0.d0

    !Do indexM=1,intervals

     !  sum = (M_functions(indexM+1)-M_functions(indexM))/2.d0*( f(indexM) + f(indexM+1) ) + sum

!    End Do

 !   pre_Cl_2 = sum

  !end function pre_Cl_2

  subroutine compute_integrand_pre_cl_yphi_2h_at_z_and_l()

    use arrays 
    use fiducial

    Implicit none 

    Integer*4 :: indexl,indexM,indexz

    Do indexl=1,number_of_l

       Do indexM=1,number_of_M

          Do indexz=1,1,number_of_z

             inte_pre_cl_yphi_2h(indexl,indexM,indexz) = dndM(indexM,indexz)*bMz(indexM,indexz)*ylMz(indexl,indexM,indexz) 

          End Do

       End Do

    End Do

  end subroutine compute_integrand_pre_cl_yphi_2h_at_z_and_l

  function integrand_pre_cl_yphi_2h_at_z_and_l(virial_mass,indexz,indexl)

    use arrays
    use fiducial

    Implicit none

    Real*8 :: virial_mass,integrand_pre_cl_yphi_2h_at_z_and_l,dalpha
!    Real*8,dimension(number_of_M) :: f

    Integer*4 :: indexz,indexl!,indexl

 !   Do indexM=1,number_of_M

  !     f(indexM) = dndM(indexM,indexz)*bMz(indexM,indexz)*ylMz(indexl,indexM,indexz) 

   ! End Do

    call Interpolate_1D(integrand_pre_cl_yphi_2h_at_z_and_l,dalpha,virial_mass,M,inte_pre_cl_yphi_2h(indexl,:,indexz))

  End function integrand_pre_cl_yphi_2h_at_z_and_l

  Function integrand_pre_cl_yphi_2h(virial_mass, params) bind(c)

    real(c_double), value :: virial_mass
    type(c_ptr), value :: params
    real(c_double) :: integrand_pre_cl_yphi_2h
    Integer(c_int), pointer :: pa(:)

    call c_f_pointer(params, pa,(/2/)) ! pa(1) = indexz, pa(2) = indexl

    integrand_pre_cl_yphi_2h =  integrand_pre_cl_yphi_2h_at_z_and_l(virial_mass,pa(1),pa(2))

  End function integrand_pre_cl_yphi_2h

  Subroutine compute_pre_cl_yphi_2h_at_z_and_l(indexz,indexl,output)

    use fiducial

    Implicit none

    Integer(fgsl_size_t), parameter :: nmax=10000

    Integer(fgsl_int),target :: pp(2)
    Real(fgsl_double) :: result, error, output
    Real(fgsl_double),parameter :: lower_limit = Mmin!1.0E-3_fgsl_double
    Real(fgsl_double),parameter :: upper_limit = Mmax!1.0E-2_fgsl_double
    Real(fgsl_double),parameter :: absolute_error = 0.0_fgsl_double
    Real(fgsl_double),parameter :: relative_error = 1.0E-5_fgsl_double

    Integer(fgsl_int) :: status
    Integer(fgsl_int) :: indexz,indexl

    Type(c_ptr) :: ptr
    Type(fgsl_function) :: f_obj
    Type(fgsl_integration_workspace) :: wk

    pp = (/indexz,indexl/)

    ptr = c_loc(pp)

    f_obj = fgsl_function_init(integrand_pre_cl_yphi_2h, ptr)

    wk = fgsl_integration_workspace_alloc(nmax)

    status = fgsl_integration_qags(f_obj, lower_limit, upper_limit, &
         absolute_error, relative_error, nmax, wk, result, error)

    output = result

    call fgsl_function_free(f_obj)

    call fgsl_integration_workspace_free(wk)

  End subroutine compute_pre_cl_yphi_2h_at_z_and_l

  subroutine compute_integrand_cl_yphi_two_halo_at_z_and_l()

    use fiducial
    use arrays
    use functions

    Implicit none

    Integer*4 :: indexl,indexz
    Real*8,dimension(number_of_z) :: f1,f2

    Do indexl=1,number_of_l

       Do indexz=1,number_of_z

          call compute_pre_cl_yphi_2h_at_z_and_l(indexz,indexl,f1(indexz)) 

          call compute_pre_cl_phiphi_2h_at_z_and_l(indexz,indexl,f2(indexz))

          inte_cl_yphi_2h(indexl,indexz) = f1(indexz)*f2(indexz)*d2VdzdO(indexz)*&
            matter_power_spectrum((dble(ml(indexl))+1.d0/2.d0)/comoving_distance_at_z(indexz),z(indexz))    !  Dimensionless

       End Do

    End Do

  end subroutine compute_integrand_cl_yphi_two_halo_at_z_and_l
    
  function integrand_cl_yphi_two_halo_at_z_and_l(redshift,indexl)

    use arrays
    use fiducial
!    use functions

    Implicit none

    Real*8 :: redshift,integrand_cl_yphi_two_halo_at_z_and_l,dalpha
!    Real*8,dimension(number_of_z) :: f1,f2,f3

    Integer*4 :: indexl!,indexl

 !   Do indexz=1,number_of_z

  !     call compute_pre_cl_yphi_2h_at_z_and_l(indexz,indexl,f1(indexz)) 

   !    call compute_pre_cl_phiphi_2h_at_z_and_l(indexz,indexl,f2(indexz))

    !   f3(indexz) = f1(indexz)*f2(indexz)*d2VdzdO(indexz)*&
     !       matter_power_spectrum((dble(ml(indexl))+1.d0/2.d0)/comoving_distance_at_z(indexz),z(indexz))    !  Dimensionless

!    End Do

    call Interpolate_1D(integrand_cl_yphi_two_halo_at_z_and_l,dalpha,redshift,z,inte_cl_yphi_2h(indexl,:))

  End function integrand_cl_yphi_two_halo_at_z_and_l

  Function integrand_cl_yphi_two_halo(redshift, params) bind(c)

    real(c_double), value :: redshift
    type(c_ptr), value :: params
    real(c_double) :: integrand_cl_yphi_two_halo
    Integer(c_int), pointer :: pa

    call c_f_pointer(params, pa) ! pa = indexl

    integrand_cl_yphi_two_halo =  integrand_cl_yphi_two_halo_at_z_and_l(redshift,pa)

  End function integrand_cl_yphi_two_halo

  Subroutine compute_cl_yphi_two_halo_at_l(indexl,output)

    Implicit none

    Integer(fgsl_size_t), parameter :: nmax=10000

    Integer(fgsl_int),target :: pp
    Real(fgsl_double) :: result, error, output
    Real(fgsl_double),parameter :: lower_limit = zmin!1.0E-3_fgsl_double
    Real(fgsl_double),parameter :: upper_limit = zmax!1.0E-2_fgsl_double
    Real(fgsl_double),parameter :: absolute_error = 0.0_fgsl_double
    Real(fgsl_double),parameter :: relative_error = 1.0E-5_fgsl_double

    Integer(fgsl_int) :: status
    Integer(fgsl_int) :: indexl

    Type(c_ptr) :: ptr
    Type(fgsl_function) :: f_obj
    Type(fgsl_integration_workspace) :: wk

    pp = indexl

    ptr = c_loc(pp)

    f_obj = fgsl_function_init(integrand_cl_yphi_two_halo, ptr)

    wk = fgsl_integration_workspace_alloc(nmax)

    status = fgsl_integration_qags(f_obj, lower_limit, upper_limit, &
         absolute_error, relative_error, nmax, wk, result, error)

    output = result/h**3

    call fgsl_function_free(f_obj)

    call fgsl_integration_workspace_free(wk)

  End subroutine compute_cl_yphi_two_halo_at_l

  subroutine compute_Clyphi2h()   ! It fills two halo term of lensing potential angular power spectrum array. 

    use arrays
    use fiducial
    Implicit none

    Integer*4 :: indexl

    Do indexl=1,number_of_l

       call compute_cl_yphi_two_halo_at_l(indexl,Cl2h(indexl))

    End Do

  end subroutine compute_Clyphi2h



!  function C_l_yphi_two_halo(indexl) ! Equation (2.11) in 1312.4525

 !   use fiducial                   ! Units : dimensionless
  !  use arrays
   ! use omp_lib
    !Implicit none

!    Real*8 :: C_l_yphi_two_halo,sum
 !   Integer*4 :: indexl,indexz
  !  Integer*4,parameter :: number_of_redshift = number_of_z_functions !1d4
   ! Integer*4,parameter :: intervals = number_of_redshift - 1 
    !Real*8,dimension(number_of_redshift):: f!,redshift

    !        Do indexz=1,number_of_redshift

    !           redshift(indexz) = 10**(log10(zmin) + real(indexz-1)*(log10(zmax) - log10(zmin))/real(number_of_redshift-1))

    !      End Do

    !!$omp Parallel Do Default(Shared)

!    Do indexz=1,number_of_redshift

 !      f(indexz) = pre_Cl_1(indexz,indexl)*pre_Cl_2(indexz,indexl)*d2VdzdO(indexz)*&
  !          matter_power_spectrum((dble(ml(indexl))+1.d0/2.d0)/comoving_distance_at_z(indexz),z_functions(indexz)) ! Units : 1/h**3

   ! End Do

    !!$omp End Parallel Do

    !        open(16,file='./output/clphiphi.dat')

    !        Do indexz=1,number_of_z

    !            write(16,'(2es18.10)') z(indexz),f(indexz) 

    !        End Do

    !        close(16)

    !sum = 0.d0

    !Do indexz=1,intervals

     !  sum = (z_functions(indexz+1)-z_functions(indexz))/2.d0*( f(indexz) + f(indexz+1) ) + sum

!    End Do

 !   C_l_yphi_two_halo = sum/h**3

  !end function C_l_yphi_two_halo

!  subroutine compute_Cl2h()

 !   use arrays
  !  use fiducial
   ! Implicit none

    !Integer*4 :: indexl

!    If (compute_functions) then

 !      print *, 'TWO HALO TERM FOR CROSS-CORRELATIONS YPHI IS COMPUTED ONLY WITH INTERPOLATING FUNCTIONS'

  !     stop

   ! Else

    !   Do indexl=1,number_of_l

     !     Cl2h(indexl) = C_l_yphi_two_halo(indexl)

      ! End Do

 !   End If

  !end subroutine compute_Cl2h

  subroutine compute_Cl()

    use arrays
    use fiducial
    Implicit none

    Integer*4 :: indexl

    Do indexl=1,number_of_l

       Cl(indexl) = Cl1h(indexl) + Cl2h(indexl)

       Clphiphi(indexl) = Clphiphi1h(indexl) + Clphiphi2h(indexl)

    End Do

  end subroutine compute_Cl


End Module integrator
