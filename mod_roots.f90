Module mod_roots

  Use fgsl
  Use fiducial
  Use functions
  Use arrays
  Use, intrinsic :: iso_c_binding

  Implicit none

contains

  function r_delta_c_from_M_virial(x, params) bind(c)

    real(c_double), value :: x
    type(c_ptr), value :: params
    real(c_double) :: r_delta_c_from_M_virial
    
    real(fgsl_double), pointer :: p(:)
    
    call c_f_pointer(params, p, (/3/))

    r_delta_c_from_M_virial = integral_r_delta_minus_total_matter(p(1),p(2),x,p(3))
    
  end function r_delta_c_from_M_virial

  subroutine compute_r_delta_c_from_M_virial_at_z(mass,redshift,delta,rdeltac)

    Implicit none

    Real(fgsl_double), parameter :: eps = 1.0E-15_fgsl_double
    Real(fgsl_double),parameter :: lower_limit = 1.0E-5_fgsl_double
    Real(fgsl_double),parameter :: upper_limit = 1.0E2_fgsl_double

    Integer(fgsl_int), parameter :: itmax = 1000

    Real(fgsl_double), target :: fpar(3)
    Real(fgsl_double) :: root, xlo, xhi,mass,redshift,delta,rdeltac

    Character(kind=fgsl_char,len=fgsl_strmax) :: name

    Integer :: i
    Integer(fgsl_int) :: status

    Type(c_ptr) :: ptr
    Type(fgsl_root_fsolver) :: root_fslv
    Type(fgsl_function) :: func

    fpar = (/mass,redshift,delta/)

    ptr = c_loc(fpar)

    root_fslv = fgsl_root_fsolver_alloc(fgsl_root_fsolver_brent)

    func = fgsl_function_init(r_delta_c_from_M_virial, ptr)

    If (fgsl_well_defined(root_fslv)) then

       status = fgsl_root_fsolver_set(root_fslv, func, lower_limit, upper_limit)

       name = fgsl_root_fsolver_name (root_fslv)

       i = 0

       Do

          i = i + 1

          status = fgsl_root_fsolver_iterate(root_fslv)

          If (status /= fgsl_success .or. i > itmax) then

             write(6, *) 'FAILED TO CONVERGE OR ITERATE'

             exit

          End If

          root = fgsl_root_fsolver_root(root_fslv)

          xlo = fgsl_root_fsolver_x_lower(root_fslv)

          xhi = fgsl_root_fsolver_x_upper(root_fslv)

          status = fgsl_root_test_interval (xlo, xhi, 0.0_fgsl_double, eps)

          If (status == fgsl_success) exit

       End do

    End If

    !write(6, '(''Output for root finding algorithm '',A,'':'')') trim(name)

    !write(6, '(''Number of iterations needed: '',i4)') i

    !write(6, '(''Root of '',F5.2,''*x*x + '',F5.2,''*x +'',F5.2,'' is: '',1PE18.10)') &
    !     fpar, root

    rdeltac = root

    !write(6, '(''Function value at root is: '',1PE18.10)') r_delta_c_from_M_virial(root,ptr)

    call fgsl_root_fsolver_free(root_fslv)

    call fgsl_function_free(func)
      
  end subroutine compute_r_delta_c_from_M_virial_at_z

  function r_delta_d_from_M_virial(x, params) bind(c)

    real(c_double), value :: x
    type(c_ptr), value :: params
    real(c_double) :: r_delta_d_from_M_virial
    
    real(fgsl_double), pointer :: p(:)
    
    call c_f_pointer(params, p, (/3/))

    r_delta_d_from_M_virial = integral_r_delta_d_minus_total_matter(p(1),p(2),x,p(3))
    
  end function r_delta_d_from_M_virial

  subroutine compute_r_delta_d_from_M_virial_at_z(mass,redshift,delta,rdeltad)

    Implicit none

    Real(fgsl_double), parameter :: eps = 1.0E-15_fgsl_double
    Real(fgsl_double),parameter :: lower_limit = 1.0E-5_fgsl_double
    Real(fgsl_double),parameter :: upper_limit = 1.0E2_fgsl_double

    Integer(fgsl_int), parameter :: itmax = 1000

    Real(fgsl_double), target :: fpar(3)
    Real(fgsl_double) :: root, xlo, xhi, mass, redshift, delta, rdeltad

    Character(kind=fgsl_char,len=fgsl_strmax) :: name

    Integer :: i
    Integer(fgsl_int) :: status

    Type(c_ptr) :: ptr
    Type(fgsl_root_fsolver) :: root_fslv
    Type(fgsl_function) :: func

    fpar = (/mass,redshift,delta/)

    ptr = c_loc(fpar)

    root_fslv = fgsl_root_fsolver_alloc(fgsl_root_fsolver_brent)

    func = fgsl_function_init(r_delta_d_from_M_virial, ptr)

    If (fgsl_well_defined(root_fslv)) then

       status = fgsl_root_fsolver_set(root_fslv, func, lower_limit, upper_limit)

       name = fgsl_root_fsolver_name (root_fslv)

       i = 0

       Do

          i = i + 1

          status = fgsl_root_fsolver_iterate(root_fslv)

          If (status /= fgsl_success .or. i > itmax) then

             write(6, *) 'FAILED TO CONVERGE OR ITERATE'

             exit

          End If

          root = fgsl_root_fsolver_root(root_fslv)

          xlo = fgsl_root_fsolver_x_lower(root_fslv)

          xhi = fgsl_root_fsolver_x_upper(root_fslv)

          status = fgsl_root_test_interval (xlo, xhi, 0.0_fgsl_double, eps)

          If (status == fgsl_success) exit

       End do

    End If

    !write(6, '(''Output for root finding algorithm '',A,'':'')') trim(name)

    !write(6, '(''Number of iterations needed: '',i4)') i

    !write(6, '(''Root of '',F5.2,''*x*x + '',F5.2,''*x +'',F5.2,'' is: '',1PE18.10)') &
    !     fpar, root

    rdeltad = root

    !write(6, '(''Function value at root is: '',1PE18.10)') r_delta_c_from_M_virial(root,ptr)

    call fgsl_root_fsolver_free(root_fslv)

    call fgsl_function_free(func)
      
  end subroutine compute_r_delta_d_from_M_virial_at_z

  subroutine compute_M_delta_c_from_M_and_z(delta)

    use fiducial
    use arrays
    use omp_lib
    !use mod_roots

    Implicit none

    Integer*4 :: indexz,indexM
    Real*8,dimension(1:number_of_M,number_of_z) :: rdc,rdd
    Real*8 :: delta

    rdc(:,:) = 0.d0

    rdd(:,:) = 0.d0

    open(15,file='./precomputed_quantities/M_delta_c_d.dat')
    
    write(15,*) '# Mass conversion file. The number of masses is ',number_of_z*number_of_M

    write(15,*) '# index_of_red-shift    red-shift    index_of_M    virial_mass_M[solar mass]    r_delta_c[Mpc]'    

    write(15,*) '# M_delta_c[solar mass]    r_delta_d[Mpc]    M_delta_d[solar mass] '

    !$omp Parallel Do Shared(rdc,rdd)

    Do indexz = 1,number_of_z

       Do indexM = 1,number_of_M

          call compute_r_delta_c_from_M_virial_at_z(M(indexM),z(indexz),delta,rdc(indexM,indexz))

          call compute_r_delta_d_from_M_virial_at_z(M(indexM),z(indexz),delta,rdd(indexM,indexz))

       End Do

    End Do

    !$omp End Parallel Do

    Do indexz = 1,number_of_z

       Do indexM = 1,number_of_M

          write(15,'(i5,es18.10,i10,5es18.10)') indexz, z(indexz), indexM, M(indexM), rdc(indexM,indexz), &
               M_delta_c_from_M_virial(z(indexz),rdc(indexM,indexz),delta), &
               rdd(indexM,indexz),M_delta_d_from_M_virial(z(indexz),rdd(indexM,indexz),delta)

       End Do

    End Do

    close(15)

  end subroutine compute_M_delta_c_from_M_and_z

End module mod_roots

