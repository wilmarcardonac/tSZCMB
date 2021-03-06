Module functions

    Implicit none
     
    contains

    Subroutine Interpolate_1D(function_value,derivative_at_x,point_x,array_x,array_y)
      
        use fgsl
        Implicit none

        Integer(fgsl_int) :: status
        Integer(fgsl_size_t) :: t
        real(fgsl_double) :: point_x, function_value, derivative_at_x, array_x(:), array_y(:)
        type(fgsl_interp_accel) :: acc
        type(fgsl_spline) :: spline

        t = size(array_x)

        acc = fgsl_interp_accel_alloc()

        spline = fgsl_spline_alloc(fgsl_interp_cspline,t)

        status = fgsl_spline_init(spline,array_x,array_y,t)

        function_value = fgsl_spline_eval(spline,point_x,acc)

        derivative_at_x = fgsl_spline_eval_deriv(spline,point_x,acc)

        call fgsl_spline_free(spline)

        call fgsl_interp_accel_free(acc)

    End subroutine Interpolate_1D

    Subroutine Interpolate_2D(function_value,point_x,point_y,array_x,array_y,matrix_function_values)
    
        Implicit none
    
        Real*8 :: function_value, point_x, point_y,flowerleft,flowerright,fupperright,fupperleft,t,u
        Real*8,dimension(:) :: array_x, array_y
        Real*8,dimension(:,:) :: matrix_function_values
        Integer*4 :: q,indexleft,indexright,indexup,indexdown

        indexleft = 1 

        indexright = 1 

        indexup = 1

        indexdown = 1

        Do q=1,size(array_x)

            If ( (point_x .lt. array_x(1)) .or. (point_x .gt. array_x(size(array_x))) ) then

                print *, 'COORDINATE X OF POINT ', point_x,' IS OUT OF RANGE OF GIVEN VALUES. EXTRAPOLATION IS NOT IMPLEMENTED'

                stop

            Else If ( ((array_x(q) .ge. point_x) .and. (array_x(q-1) .le. point_x)) .and. (q .gt. 1) ) then

                indexright = q 

                indexleft = q - 1

                exit

            Else if ( (q .eq. 1) .and. ((array_x(q) .le. point_x) .and. (array_x(q+1) .ge. point_x)) ) then

                indexright = q + 1 

                indexleft = q 

                exit

            End If 

        End Do

        Do q=1,size(array_y)

            If ( (point_y .lt. array_y(1)) .or. (point_y .gt. array_y(size(array_y))) ) then

                print *, 'COORDINATE Y OF POINT ', point_y,' IS OUT OF RANGE OF GIVEN VALUES. EXTRAPOLATION IS NOT IMPLEMENTED'

                stop

            Else If ( ((array_y(q) .ge. point_y) .and. (array_y(q-1) .le. point_y)) .and. (q .gt. 1) ) then

                indexup = q 

                indexdown = q - 1

                exit

            Else if ( (q .eq. 1) .and. ((array_y(q) .le. point_y) .and. (array_y(q+1) .ge. point_y)) ) then

                indexup = q + 1 

                indexdown = q 

                exit

            End If

        End Do

        flowerleft = matrix_function_values(indexleft,indexdown)

        flowerright = matrix_function_values(indexright,indexdown) 

        fupperright = matrix_function_values(indexright,indexup)

        fupperleft = matrix_function_values(indexleft,indexup)

        t = (point_x - array_x(indexleft))/(array_x(indexright) - array_x(indexleft) )

        u = (point_y - array_y(indexdown) )/( array_y(indexup) - array_y(indexdown) )

        function_value = ( 1.d0 - t )*( 1.d0 - u )*flowerleft + t*( 1.d0 - u )*flowerright +&
        t*u*fupperright + ( 1.d0 - t )*u*fupperleft

    End subroutine Interpolate_2D

    function Hubble_parameter(z)    !    Units : m/s/Mpc
    
        use fiducial
        Implicit none

        Real*8 :: Hubble_parameter,z

        Hubble_parameter = H_0*sqrt(Omega_m0*(1.d0 + z)**3 + Omega_L0)

    end function Hubble_parameter

    function comoving_distance(z)    !    UNITS: Mpc

        use fiducial
        use omp_lib
        Implicit none

        Integer*4,parameter :: number_of_x = 100000
        Integer*4,parameter :: intervals = number_of_x - 1
        Integer*4 :: m
        Real*8 :: comoving_distance, z
        Real*8,dimension(number_of_x) :: x,f
        Real*8,parameter :: xmin = 0.d0

        comoving_distance = 0.d0

        Do m=1,number_of_x

            x(m) = xmin + dble(m-1)*z/dble(number_of_x-1)

        End Do

        !$omp Parallel Do Shared(f)    

        Do m=1,number_of_x

            f(m) = 1.d0/Hubble_parameter(x(m))

        End Do

        !$omp End Parallel Do

        Do m = 1,intervals

            comoving_distance = (x(m+1)-x(m))/2.d0*(f(m) + f(m+1)) + comoving_distance

        End Do

        comoving_distance = c*comoving_distance

    end function comoving_distance

    function angular_diameter_distance(redshift)    !    Units : Mpc

        use arrays
        Implicit none

        Real*8 :: angular_diameter_distance,redshift

        angular_diameter_distance = comoving_distance(redshift)/(1.d0 + redshift)

    end function angular_diameter_distance

    function comoving_volume_per_steradian(z)    !    Units : Mpc**3 

        use fiducial
        Implicit none

        Real*8 :: comoving_volume_per_steradian,z

        comoving_volume_per_steradian = c*(1.d0 + z)**2*angular_diameter_distance(z)**2/Hubble_parameter(z)

    end function comoving_volume_per_steradian

    function comoving_volume(z)    !    Units : Mpc**3

        use fiducial
        use omp_lib
        Implicit none

        Integer*4,parameter :: number_of_x = 100000
        Integer*4,parameter :: intervals = number_of_x - 1
        Integer*4 :: m
        Real*8,dimension(number_of_x) :: x,f
        Real*8 :: comoving_volume,z
        Real*8,parameter :: xmin = 0.d0

        Do m=1,number_of_x

            x(m) = xmin + dble(m-1)*z/dble(number_of_x-1)

        End Do

        !$omp Parallel Do Shared(f)    

        Do m=1,number_of_x

            f(m) = comoving_volume_per_steradian(x(m))

        End Do

        !$omp End Parallel Do

        comoving_volume = 0.d0

        Do m = 1,intervals

            comoving_volume = (x(m+1)-x(m))/2.d0*(f(m) + f(m+1)) + comoving_volume

        End Do 

        comoving_volume = 4.d0*Pi*comoving_volume

    end function comoving_volume

    subroutine compute_d2VdzdO()    !    It fills in the vector with comoving volume per steradian

        use arrays
        use fiducial
        use omp_lib
        Implicit none

        Integer*4 :: indexz

        !$omp Parallel Do Shared(d2VdzdO)

        Do indexz=1,number_of_z

           d2VdzdO(indexz) = comoving_volume_per_steradian(z(indexz))

        End Do

        !$omp End Parallel Do

    end subroutine compute_d2VdzdO

    function Omega_m(z)    !    Matter density parameter as seen by an observer at red-shift z

        use fiducial
        Implicit none
        Real*8 :: Omega_m,z

        Omega_m = (Omega_m0*(1.d0 + z)**3)/(Omega_m0*(1.d0 + z)**3 + Omega_L0)

    end function Omega_m

    function Omega_L(z)    !    Dark energy density parameter as seen by an observer at red-shift z

        use fiducial
        Implicit none
        Real*8 :: Omega_L,z

        Omega_L = Omega_L0/(Omega_L0 + Omega_m0*(1.d0 + z)**3)

    end function Omega_L

    function critical_density(z)    !    Units : solar mass/Mpc**3.  

        use fiducial
        Implicit none
        Real*8 :: critical_density,z

        critical_density = 3.d0*Hubble_parameter(z)**2/8.d0/Pi/G 

        critical_density = critical_density*Mpc/M_sun !/h**2  

    end function critical_density

    function mean_density(z)    !    Units : solar mass/Mpc**3. Physical coordinates. Must divide by (1+z)**3 to get equation in comoving coordinates

        Implicit none
        Real*8 :: z,mean_density

        mean_density = Omega_m(z)*critical_density(z)

    end function mean_density

    function virial_radius(z,M)    !    M is the virial mass and must be given in solar mass units. Then the virial radius has units : Mpc

        use fiducial 
        Implicit none
        Real*8 :: z,Delta_cr,M,virial_radius

        Delta_cr = 18.d0*Pi**2 + 82.d0*(Omega_m(z) - 1.d0) - 39.d0*(Omega_m(z) - 1.d0)**2

        virial_radius = (3.d0*M/4.d0/Pi/Delta_cr/critical_density(z))**(1.d0/3.d0)

    end function virial_radius

    function r_delta_c(delta,z,Mdeltac)    !    Radius with respect to critical density at red-shift z. Mdeltac given in solar mass gives units : Mpc

        use fiducial
        Implicit none
        Real*8 :: delta,z,Mdeltac,r_delta_c

        r_delta_c = (3.d0*Mdeltac/4.d0/Pi/delta/critical_density(z))**(1.d0/3.d0)

    end function r_delta_c

    function r_delta_d(delta,z,Mdeltad)    !    Radius with respect to mean density at red-shift z. Mdeltad given in solar mass gives units : Mpc. 

        use fiducial    ! mean_density(z) has been corrected to match comoving coordinates
        Implicit none
        Real*8 :: delta,z,Mdeltad,r_delta_d

        r_delta_d = (3.d0*Mdeltad/4.d0/Pi/delta/mean_density(z)*(1.d0 + z )**3.d0)**(1.d0/3.d0)

    end function r_delta_d

    function critical_surface_density(redshift) ! Units : Solar mass/Mpc**2 

        use arrays
        use fiducial
        Implicit none

        Real*8 :: redshift,critical_surface_density

        critical_surface_density = c**2*com_dist_at_z_dec*(1.d0 + &
             redshift)/4.d0/Pi/G/comoving_distance(redshift)/(com_dist_at_z_dec - comoving_distance(redshift))

        critical_surface_density = critical_surface_density*Mpc/M_sun 

    end function critical_surface_density

    subroutine compute_integrand_limber_approximation()

      use arrays
      use fiducial

      Implicit none

      Real*8 :: prefactor
      Integer*4 :: indexz,indexl

      Do indexz = 1, number_of_z_limber      

         zlimber(indexz) = 10**(log10(zminlimber) + real(indexz-1)*(log10(z_dec) - log10(zminlimber))/real(number_of_z_limber-1))

      End Do

      Do indexl=1,number_of_l

         prefactor = 4.d0/dble(ml(indexl))**2/(dble(ml(indexl))+1.d0)**2*9.d0/4.d0/&
              c**4*Hubble_parameter(0.d0)**4*Omega_m(0.d0)**2    !    Units : 

         Do indexz=1,number_of_z_limber

            integrand_limber(indexz,indexl) = (1.d0 + zlimber(indexz))**2*( com_dist_at_z_dec - &
                 comoving_distance(zlimber(indexz)) )**2/com_dist_at_z_dec**2*c/Hubble_parameter(zlimber(indexz))*&
                 matter_power_spectrum((dble(ml(indexl)) + 1.d0/2.d0)/comoving_distance(zlimber(indexz)),zlimber(indexz))

            integrand_limber(indexz,indexl) = integrand_limber(indexz,indexl)*prefactor/h**3

         End Do
         
      End Do

    end subroutine compute_integrand_limber_approximation

!    subroutine compute_critical_surface_density()    !    It fills in the vector with critical surface density

!        use arrays
!        use fiducial
!        use omp_lib
!        Implicit none

!        Integer*4 :: indexz 

!        If (compute_functions) then

!           !$omp Parallel Do Shared(Scrit,z)

!           Do indexz=1,number_of_z

 !             Scrit(indexz) = critical_surface_density(z(indexz))

  !         End Do

   !        !$omp End Parallel Do

    !    Else

     !      !$omp Parallel Do Shared(Scrit,z_functions)

      !     Do indexz=1,number_of_z_functions

       !       Scrit(indexz) = critical_surface_density(z_functions(indexz))

        !   End Do

         !  !$omp End Parallel Do

!        End If

 !   end subroutine

    function concentration_mass_critical(M,z)    ! Concentration-mass relation for critical density halo definition. It computes equation (4) in

        use fiducial    ! published version of "Dark matter halo concentrations in the Wilkinson Microwave Anisotropy Probe year 5 cosmology"
        Implicit none    ! by Alan R. Duffy et al. Parameters A, B, C from second line in table 1. M_pivot from page L66. Mass must be M_{200c}

        Real*8 :: z,M,concentration_mass_critical !   Units of mass M : solar mass. Dimensionless.   
        Real*8,parameter :: A = 5.71d0    
        Real*8,parameter :: M_pivot = 2.d12/h   !    solar mass
        Real*8,parameter :: B = -0.084d0
        Real*8,parameter :: CC = -0.47d0 
    
        concentration_mass_critical = A*(M/M_pivot)**B*(1.d0 + z)**CC

    end function concentration_mass_critical

    function concentration_mass_virial(M,z) ! concentration-mass relation for virial mass halo definition. It computes equation (4) in

        use fiducial    ! published version of "Dark matter halo concentrations in the Wilkinson Microwave Anisotropy Probe year 5 cosmology" by
        Implicit none    ! Alan R. Duffy et al. Parameters A, B, C from sixth line in table 1. M_pivot from page L66. Mass must be virial mass.

        Real*8 :: z,M,concentration_mass_virial ! Units of the viral mass M : solar mass. Dimensionless.
        Real*8,parameter :: A = 7.85d0                                             
        Real*8,parameter :: M_pivot = 2.d12/h    !    solar mass
        Real*8,parameter :: B = -0.081d0
        Real*8,parameter :: CC = -0.71d0

        concentration_mass_virial = A*(M/M_pivot)**B*(1.d0 + z)**CC

    end function concentration_mass_virial

    function concentration_mass_mean(M,z)    ! concentration-mass relation for mean background density halo definition. It computes equation (4) in

        use fiducial    ! published version of "Dark matter halo concentrations in the Wilkinson Microwave Anisotropy Probe year 5 cosmology" by
        Implicit none    ! Alan R. Duffy et al. Parameters A, B, C from tenth line in table 1. M_pivot from page L66. Mass must be M_{200d}.

        Real*8 :: z,M,concentration_mass_mean ! Units of the mass M : solar mass. Dimensionless. 
        Real*8,parameter :: A = 10.14d0                                             
        Real*8,parameter :: M_pivot = 2.d12/h    !    solar mass/h
        Real*8,parameter :: B = -0.081d0
        Real*8,parameter :: CC = -1.01d0

        concentration_mass_mean = A*(M/M_pivot)**B*(1.d0 + z)**CC

    end function concentration_mass_mean

    subroutine write_c200_at_z(indexz)

        use arrays
        use fiducial
        Implicit none

        Integer*4 :: indexM,indexz
        Integer*4,parameter :: indexr = 100

        open(15,file='./output/c200_at_z.dat')

        write(15,*) '# M200c [solar mass/h]       c200 at redshift ',z(indexz)

        Do indexM =1,number_of_M

            write(15,'(2es18.10)') M200c(indexM,indexz), concentration_mass_virial(M200c(indexM,indexz),z(indexz))

        End Do

        close(15)

    end subroutine write_c200_at_z

    function NFW_density_profile(r,M,z,halo_definition) ! NFW density profile (Equation (1)) from "A universal density profile from  

        use fiducial   ! hierarchical clustering". Units of r : Mpc. Units of M : solar mass
        Implicit none    ! Units of density profile : solar mass/Mpc**3 

        Real*8 :: r,M,z,r_s,delta_c,NFW_density_profile,Delta_cr 
        Character(len=*) :: halo_definition 

        If (halo_definition .eq. 'critical_density') then

            delta_c = DeltaSO*concentration_mass_critical(M,z)**3/3.d0/&
            (log(1.d0 + concentration_mass_critical(M,z)) - &
            concentration_mass_critical(M,z)/(1.d0 + concentration_mass_critical(M,z)))

            r_s = r_delta_c(DeltaSO,z,M)/concentration_mass_critical(M,z)

        Else If (halo_definition .eq. 'mean_background') then 

            delta_c = DeltaSO*concentration_mass_mean(M,z)**3*mean_density(z)/(1.d0 + z)**3.d0/& ! mean_density(z) corrected to match comoving coordinates
            3.d0/critical_density(z)/(log(1.d0 + concentration_mass_mean(M,z)) - &  
            concentration_mass_mean(M,z)/(1.d0 + concentration_mass_mean(M,z)))

            r_s = r_delta_d(DeltaSO,z,M)/concentration_mass_mean(M,z)

        Else If (halo_definition .eq. 'virial') then

            Delta_cr = 18.d0*Pi**2 + 82.d0*(Omega_m(z) - 1.d0) - 39.d0*(Omega_m(z) - 1.d0)**2

            delta_c = Delta_cr*concentration_mass_virial(M,z)**3/3.d0/&
            (log(1.d0 + concentration_mass_virial(M,z)) - &
            concentration_mass_virial(M,z)/(1.d0 + concentration_mass_virial(M,z)))

            r_s = virial_radius(z,M)/concentration_mass_virial(M,z)

        Else

            print *,'Not halo definition implemented with name ', halo_definition

            stop

        End If

        NFW_density_profile = critical_density(z)*delta_c/(r/r_s)/(1.d0 + r/r_s)**2 

    end function NFW_density_profile

    function r2NFW_density_profile(r,M,z,halo_definition)    !    r**2 times NFW density profile. Units : solar mass/Mpc

        use fiducial
        Implicit none                          

        Real*8 :: r,M,z,r2NFW_density_profile 
        Character(len=*) :: halo_definition 

        r2NFW_density_profile = r**2*NFW_density_profile(r,M,z,halo_definition)

    end function r2NFW_density_profile

    function integral_r_delta_c(M,z,rdc)    !    Equation (4) in published version of 1303.4726. Units of M : solar mass. 

        use fiducial    !    Units of rdc : Mpc. Units : solar mass
        Implicit none

        Real*8 :: M,z,rdc,rmin,integral_r_delta_c,Delta_cr,delta_c,r_s

        Delta_cr = 18.d0*Pi**2 + 82.d0*(Omega_m(z) - 1.d0) - 39.d0*(Omega_m(z) - 1.d0)**2

        delta_c = Delta_cr*concentration_mass_virial(M,z)**3/3.d0/&
        (log(1.d0 + concentration_mass_virial(M,z)) - &
        concentration_mass_virial(M,z)/(1.d0 + concentration_mass_virial(M,z)))

        r_s = virial_radius(z,M)/concentration_mass_virial(M,z)

        rmin = rdc/r_s
 
        integral_r_delta_c = 4.d0*Pi*r_s**3*critical_density(z)*delta_c*(-1.d0 + 1.d0/(1.d0 + rmin) + log(1.d0 + rmin)  )

    end function integral_r_delta_c

    function integral_r_delta_minus_total_matter(M,z,rdc,delta)    !    Units of M : solar mass. Units of rdc : Mpc. Units : solar mass

        use fiducial
        Implicit none
        Real*8 :: M,z,rdc,delta,integral_r_delta_minus_total_matter

        integral_r_delta_minus_total_matter = integral_r_delta_c(M,z,rdc) - 4.d0*Pi*rdc**3*critical_density(z)*delta/3.d0

    end function integral_r_delta_minus_total_matter

    subroutine compute_dMdc_dM_and_dMdd_dM()

        use fiducial
        use arrays
        use omp_lib
        Implicit none

        Integer*4 :: indexM,indexz

        open(15,file='./precomputed_quantities/dMdd_dM_dMdc_dM.dat')

        write(15,*) '# Mass derivatives file. The number of data is ',number_of_z*number_of_M

        write(15,*) '# index_of_red-shift    red-shift    index_of_M    dM200ddM    dM200cdM'    

        !$omp Parallel Do Shared(dM200ddM,dM200cdM)

        Do indexz=1,number_of_z

           Do indexM=1,number_of_M

                 dM200ddM(indexM,indexz) = dMdd_dM(indexM,indexz)

                 dM200cdM(indexM,indexz) = dMdc_dM(indexM,indexz)
 
           End Do

        End Do

        !$omp End Parallel Do

        Do indexz=1,number_of_z

           Do indexM=1,number_of_M

              write(15,'(i5,es18.10,i10,3es18.10)') indexz, z(indexz), indexM, M(indexM),&
                   dM200ddM(indexM,indexz),dM200cdM(indexM,indexz)
 
           End Do

        End Do

        close(15)

    end subroutine compute_dMdc_dM_and_dMdd_dM

    subroutine read_dMdc_dM_and_dMdd_dM()

        use arrays
        use fiducial
        Implicit none

        Integer*4 :: index,iz,iM
        Real*8 :: MM,zz,tM200c,tr200d

        open(15,file='./precomputed_quantities/dMdd_dM_dMdc_dM.dat')

        read(15,*)

        read(15,*)

        Do index=1,number_of_z*number_of_M

            read(15,'(i5,es18.10,i10,3es18.10)') iz,zz,iM,MM,tr200d,tM200c

            dM200ddM(iM,iz) = tr200d

            dM200cdM(iM,iz) = tM200c

        End Do

        close(15)

    end subroutine read_dMdc_dM_and_dMdd_dM

    function dMdc_dM(indexM,indexz)

        use fiducial
        use arrays
        Implicit none

        Integer*4 :: indexM,indexz
        Real*8 :: dMdc_dM,M200c_at_M_z

!        Real*8,parameter :: dM = 10**((log10(Mmax) - log10(Mmin))/real(number_of_M-1))

!        dMdc_dM  = (-M200c(indexM+2,indexz) + 8.d0*M200c(indexM+1,indexz) - 8.d0*M200c(indexM-1,indexz) +&
!        M200c(indexM-2,indexz))/12.d0/dM

        call Interpolate_1D(M200c_at_M_z,dMdc_dM,M(indexM),M(:),M200c(:,indexz))

    end function dMdc_dM

    function dMdd_dM(indexM,indexz)

        use fiducial
        use arrays
        Implicit none

        Integer*4 :: indexM,indexz
        Real*8 :: dMdd_dM,M200d_at_M_z
!        Real*8,parameter :: dM = 10**((log10(Mmax) - log10(Mmin))/real(number_of_M-1))

!        dMdd_dM = (-M200d(indexM+2,indexz) + 8.d0*M200d(indexM+1,indexz) - 8.d0*M200d(indexM-1,indexz) +&
!             M200d(indexM-2,indexz))/12.d0/dM

        call Interpolate_1D(M200d_at_M_z,dMdd_dM,M(indexM),M(:),M200d(:,indexz))

    end function dMdd_dM

!    function r_delta_c_from_M_virial(M,z,delta)    !    Units of M : solar mass. Units of r_delta_c_from_M_virial : Mpc

!        use fiducial
!        Implicit none

!        Real*8 :: r_delta_c_from_M_virial,M,z,delta,r1,r2,step,rini,f1,f2,rmid,dr,fmid
!        Integer*4,parameter :: iter = 1000000000
!        Real*8,parameter :: racc = 1.d-10
 !       Integer*4 :: p
  !      logical :: root

   !     rini = (3.d0*M/4.d0/Pi/critical_density(z)/delta)**(1.d0/3.d0) ! Units : Mpc

    !    step = 1.6d-2

     !   r1 = 1.d0*rini  ! Units : Mpc
        
      !  r2 = 1.1d0*rini ! Units : Mpc

       ! f1 = integral_r_delta_minus_total_matter(M,z,r1,delta) ! Units : solar mass

!        f2 = integral_r_delta_minus_total_matter(M,z,r2,delta) ! Units : solar mass

 !       root = .false. 

  !      Do p = 1, iter

   !         If (f1*f2 .lt. 0.d0) then

    !            root = .true.

     !           exit

      !      else if (abs(f1) .lt. abs(f2)) then

       !         r1 = r1 + step*(r1-r2)

        !        f1 = integral_r_delta_minus_total_matter(M,z,r1,delta)

         !   else if (abs(f1) .gt. abs(f2)) then

          !      r2 = r2 + step*(r2-r1)

           !     f2 = integral_r_delta_minus_total_matter(M,z,r2,delta)

            !end if

!        End Do

 !       If (root) then

  !          Do p=1,iter

   !             If (f1 .lt. 0.d0) then 

    !                dr = r2 - r1

     !               rmid = r1 + dr/2.d0

      !              fmid = integral_r_delta_minus_total_matter(M,z,rmid,delta)    

       !             If (fmid .lt. 0.d0) then 

        !                f1 = fmid

         !               r1 = rmid

          !          else if ((fmid .eq. 0.d0) .or. (abs(dr) .lt. racc)) then

           !             exit

            !        else 

             !           f2 = fmid

              !          r2 = rmid

               !     end if

                !else if (f2 .lt. 0.d0) then

!                    dr = r1 - r2

 !                   rmid = r2 + dr/2.d0

  !                  fmid = integral_r_delta_minus_total_matter(M,z,rmid,delta)    

   !                 If (fmid .lt. 0.d0) then 

    !                    f2 = fmid

     !                   r2 = rmid

      !              else if ((fmid .eq. 0.d0) .or. (abs(dr) .lt. racc)) then

       !                 exit

        !            else

         !               f1 = fmid

          !              r1 = rmid

           !         end if

            !    End if

!            End Do

 !           r_delta_c_from_M_virial = rmid

  !      else 

   !        print *, 'NO ROOT FOUND IN "r_delta_c_from_M_virial" '

    !       stop !r_delta_c_from_M_virial = -1.d10

     !   End If

!    end function r_delta_c_from_M_virial

    function M_delta_c_from_M_virial(z,rdc,delta)    !    Units of M and M_delta_c_from_M_virial: solar mass

        use fiducial
        Implicit none

        Real*8 :: rdc,z,delta,M_delta_c_from_M_virial

        M_delta_c_from_M_virial = 4.d0*Pi*rdc**3*critical_density(z)*delta/3.d0

    end function M_delta_c_from_M_virial

    function integral_r_delta_d(M,z,rdd)    !    Similar to 'integral_r_delta_c' above. Units of M : solar mass. Units of rdd : Mpc. Units: solar mass

        use fiducial    !    Units of rdc : Mpc. Units : solar mass
        Implicit none

        Real*8 :: M,z,rdd,rmin,integral_r_delta_d,Delta_cr,delta_c,r_s

        Delta_cr = 18.d0*Pi**2 + 82.d0*(Omega_m(z) - 1.d0) - 39.d0*(Omega_m(z) - 1.d0)**2

        delta_c = Delta_cr*concentration_mass_virial(M,z)**3/3.d0/&
        (log(1.d0 + concentration_mass_virial(M,z)) - &
        concentration_mass_virial(M,z)/(1.d0 + concentration_mass_virial(M,z)))

        r_s = virial_radius(z,M)/concentration_mass_virial(M,z)

        rmin = rdd/r_s
 
        integral_r_delta_d = 4.d0*Pi*r_s**3*critical_density(z)*delta_c*(-1.d0 + 1.d0/(1.d0 + rmin) + log(1.d0 + rmin)  )

    end function integral_r_delta_d

    function integral_r_delta_d_minus_total_matter(M,z,rdd,delta)    !    Units of M : solar mass. Units of rdd : Mpc. Units : solar mass

        use fiducial
        Implicit none

        Real*8 :: M,z,rdd,delta,integral_r_delta_d_minus_total_matter

        integral_r_delta_d_minus_total_matter = integral_r_delta_d(M,z,rdd) - &
        4.d0*Pi*rdd**3*mean_density(z)/(1.d0+z)**3.d0*delta/3.d0

    end function integral_r_delta_d_minus_total_matter

!    function r_delta_d_from_M_virial(M,z,delta)    !    Units of M : solar mass. Units of r_delta_d_from_M_virial : Mpc

 !       use fiducial
  !      Implicit none

   !     Real*8 :: r_delta_d_from_M_virial,M,z,delta,r1,r2,step,rini,f1,f2,rmid,dr,fmid
    !    Integer*4,parameter :: iter = 1000000000
     !   Real*8,parameter :: racc = 1.d-10
      !  Integer*4 :: p
       ! logical :: root 

!        rini = (3.d0*M/4.d0/Pi/mean_density(z)*(1.d0+z)**3.d0/delta)**(1.d0/3.d0) ! Units : Mpc/h
        
 !       step = 1.6d-2

  !      r1 = 1.d0*rini ! Units : Mpc/h

   !     r2 = 1.1d0*rini ! Units : Mpc/h

    !    f1 = integral_r_delta_d_minus_total_matter(M,z,r1,delta)

     !   f2 = integral_r_delta_d_minus_total_matter(M,z,r2,delta)

      !  root = .false.

       ! Do p = 1, iter

        !    If (f1*f2 .lt. 0.d0) then

!                root = .true.

 !               exit

  !          Else If (abs(f1) .lt. abs(f2)) then

   !             r1 = r1 + step*(r1-r2)

    !            f1 = integral_r_delta_d_minus_total_matter(M,z,r1,delta)

     !       Else If (abs(f1) .gt. abs(f2)) then

      !          r2 = r2 + step*(r2-r1)

       !         f2 = integral_r_delta_d_minus_total_matter(M,z,r2,delta)

        !    End if

!        End Do

 !       If (root) then

  !          Do p=1,iter

   !             If (f1 .lt. 0.d0) then 

    !                dr = r2 - r1

     !               rmid = r1 + dr/2.d0

      !              fmid = integral_r_delta_d_minus_total_matter(M,z,rmid,delta)    

       !             If (fmid .lt. 0.d0) then 

        !                f1 = fmid

         !               r1 = rmid

          !          else if ((fmid .eq. 0.d0) .or. (abs(dr) .lt. racc)) then

           !             exit

!                    else 

 !                       f2 = fmid

  !                      r2 = rmid

   !                 end if

    !            else if (f2 .lt. 0.d0) then

     !               dr = r1 - r2

      !              rmid = r2 + dr/2.d0

       !             fmid = integral_r_delta_d_minus_total_matter(M,z,rmid,delta)    

        !            If (fmid .lt. 0.d0) then 

         !               f2 = fmid

          !              r2 = rmid

           !         else if ((fmid .eq. 0.d0) .or. (abs(dr) .lt. racc)) then

            !            exit

             !       else

              !          f1 = fmid

               !         r1 = rmid

!                    end if

 !               End if

  !          End Do

   !         r_delta_d_from_M_virial = rmid

    !    else 

     !      print *,'NO ROOT FOUND IN "r_delta_d_from_M_virial" '

      !     stop!r_delta_d_from_M_virial = -1.d10

       ! End If

!    end function r_delta_d_from_M_virial

    function M_delta_d_from_M_virial(z,rdd,delta)    !    Units of M and M_delta_d_from_M_virial: solar mass

        use fiducial
        Implicit none

        Real*8 :: rdd,z,delta,M_delta_d_from_M_virial

        M_delta_d_from_M_virial = 4.d0*Pi*rdd**3*mean_density(z)/(1.d0+z)**3.d0*delta/3.d0
        
    end function M_delta_d_from_M_virial


    subroutine read_M200dc_r200dc()

        use arrays
        use fiducial
        Implicit none

        Integer*4 :: index,iz,iM
        Real*8 :: MM,zz,tM200c,tM200d,tr200c,tr200d

        open(15,file='./precomputed_quantities/M_delta_c_d.dat')

        read(15,*)

        read(15,*)

        read(15,*)

        Do index=1,number_of_z*number_of_M

            read(15,'(i5,es18.10,i10,5es18.10)') iz,zz,iM,MM,tr200c,tM200c,tr200d,tM200d

            r200c(iM,iz) = tr200c

            M200c(iM,iz) = tM200c

            r200d(iM,iz) = tr200d

            M200d(iM,iz) = tM200d

        End Do

        close(15)

    end subroutine read_M200dc_r200dc

    function ICM_electron_pressure(r,indexM,indexz) ! Units of Mv : solar mass. Units of r : Mpc. Intracluster medium electron pressure from "On 

        use fiducial              ! the cluster physics of Sunyaev-Zel'dovich and x-ray surveys. II. Deconstructing the thermal SZ
        use arrays                ! power spectrum" by Battaglia et al. All the masses and distances in this work are given relative to h.
        Implicit none             ! Units of ICM_electron_pressure : solar mass/(Mpc*s**2)        

        ! Parameters from columns "AGN Feedback \Delta = 200" in Table 1 in published version of 1109.3711
        Real*8,parameter :: A_p = 18.1d0                      ! FOR P_0
        Real*8,parameter :: alpha_pm = 0.154d0                ! FOR P_0
        Real*8,parameter :: alpha_pz = -0.758d0               ! FOR P_0

        Real*8,parameter :: A_xc = 0.497d0                    ! FOR x_c
        Real*8,parameter :: alpha_xcm = -0.00865d0            ! FOR x_c
        Real*8,parameter :: alpha_xcz = 0.731d0               ! FOR x_c

        Real*8,parameter :: A_beta = 4.35d0                   ! FOR beta
        Real*8,parameter :: alpha_betam = 0.0393d0            ! FOR beta
        Real*8,parameter :: alpha_betaz = 0.415d0             ! FOR beta

        Real*8,parameter :: gamma = -0.3d0                    ! Parameters \gamma and \alpha in equation (10) of published version of 1109.3711
        Real*8,parameter :: alpha = 1.0d0
        
        Real*8,parameter :: thermal_gas_pressure_electron_pressure_relation = 1.932d0

        Real*8 :: x,P_0,x_c,r,ICM_electron_pressure,beta,P_th,P200,P_fit

        Integer*4 :: indexM,indexz

        P200 = G*M200c(indexM,indexz)*DeltaSO*critical_density(z(indexz))&
        *Omega_b_h2/h**2/Omega_m0/2.d0/r200c(indexM,indexz)*M_sun/(Mpc)**3 ! Units : solar mass/(Mpc*s**2)

        x = r/r200c(indexM,indexz)        ! Dimensionless 

        P_0 = A_p*(M200c(indexM,indexz)*1.d-14)**alpha_pm*( 1.d0 + z(indexz) )**alpha_pz

        x_c = A_xc*(M200c(indexM,indexz)*1.d-14)**alpha_xcm*( 1.d0 + z(indexz) )**alpha_xcz

        beta = A_beta*(M200c(indexM,indexz)*1.d-14)**alpha_betam*( 1.d0  + z(indexz) )**alpha_betaz

        P_fit = P_0*(x/x_c)**gamma*( 1.d0 + (x/x_c)**alpha )**(-beta) ! Equation (10) in published version of 1109.3711

        P_th = P200*P_fit

        ICM_electron_pressure = P_th/thermal_gas_pressure_electron_pressure_relation

    end function ICM_electron_pressure    

    subroutine write_ICM_electron_pressure_at_z(indexM,indexz)

        use arrays
        use fiducial
        Implicit none

        Integer*4 :: indexM,indexz,index
        
        Real*8 :: P200
        Real*8,dimension(number_of_M) :: r
        Real*8 :: rmin 
        Real*8 :: rmax 
        
        rmin = 1.d-2*r200c(indexM,indexz)

        rmax = 3.d0*r200c(indexM,indexz)

        P200 = G*M200c(indexM,indexz)*DeltaSO*critical_density(z(indexz))&
             *Omega_b_h2/h**2/Omega_m0/2.d0/r200c(indexM,indexz)*M_sun/(Mpc)**3 ! Units : solar mass/(Mpc*s**2)
        
        open(15,file='./output/Pth_at_z.dat')

        write(15,*) '# r/R_200       P_th/P_200*(r/R_200)**3 at redshift and mass',z(indexz),M200c(indexM,indexz)

        Do index=1,number_of_M

           r(index) = 10**(log10(rmin) + real(index-1)*(log10(rmax) - log10(rmin))/real(number_of_M-1))

            write(15,'(2es18.10)') r(index)/r200c(indexM,indexz), &
                 ICM_electron_pressure(r(index),indexM,indexz)/P200*&
                 (r(index)/r200c(indexM,indexz))**3

        End Do

        close(15)

    end subroutine write_ICM_electron_pressure_at_z


    subroutine read_ylMz()

        use arrays
        use fiducial
        Implicit none

        Integer*4 :: index,il,iM,iz,ll
        Real*8 :: MM,zz,tylMz

        open(15,file='./precomputed_quantities/form_factor/form_factor.dat')

        read(15,*)

        read(15,*)

        Do index=1,number_of_l*number_of_z*number_of_M

            read(15,'(3i10,es18.10,i5,2es18.10)') il,ll,iM,MM,iz,zz,tylMz

            ylMz(il,iM,iz) = tylMz

        End Do

        close(15)

    end subroutine read_ylMz

!    subroutine Interpolate_ylMz()

 !     use arrays
  !    use fiducial
      !use omp_lib
   !   Implicit none

    !  Integer*4 :: indexM,indexz,indexl

     ! Do indexl=1,number_of_l

         !!!!$omp Parallel Do Shared(ylMz_interpolation,M_functions,z_functions,M,z,ylMz)

      !   Do indexM=1,number_of_M_functions

       !     Do indexz=1,number_of_z_functions

        !       call Interpolate_2D(ylMz_interpolation(indexl,indexM,indexz),M_functions(indexM),z_functions(indexz),&
         !           M(1:number_of_M),z(1:number_of_z),ylMz(indexl,1:number_of_M,1:number_of_z))

          !  End Do

!         End Do

         !!!!$omp End Parallel Do

 !     End Do

  !  end subroutine Interpolate_ylMz

    function lensing_potential(virial_Mass,redshift,indexl,halo_definition)    ! Lensing potential. Equation (2.10) in 1312.4525. Units of M : solar mass

        use fiducial    ! Units : dimensionless
        use arrays
        Implicit none

        Real*8 :: lensing_potential,r_s,l_s,wavevector,prefactor,Delta_cr,delta_c,virial_Mass,redshift
        Integer*4 :: indexl
        Character(len=*) :: halo_definition 

        If (halo_definition .eq. 'virial') then

            r_s = virial_radius(redshift,virial_Mass)/concentration_mass_virial(virial_Mass,redshift)    ! Mpc

            Delta_cr = 18.d0*Pi**2 + 82.d0*(Omega_m(redshift) - 1.d0) - 39.d0*(Omega_m(redshift) - 1.d0)**2

            delta_c = Delta_cr*concentration_mass_virial(virial_Mass,redshift)**3/3.d0/&
            (log(1.d0 + concentration_mass_virial(virial_Mass,redshift)) - &
            concentration_mass_virial(virial_Mass,redshift)/(1.d0 + concentration_mass_virial(virial_Mass,redshift)))

            l_s = angular_diameter_distance(redshift)/r_s                             ! dimensionless

            prefactor = 8.d0*Pi*critical_density(redshift)*delta_c*r_s/dble(ml(indexl))/dble(ml(indexl)+1.d0)/&
            l_s**2

            wavevector = ( dble(ml(indexl)) + 1.d0/2.d0 )/l_s

        Else

            print *,'Not halo definition implemented with name ', halo_definition

            stop

        End If

        lensing_potential = prefactor*FT_NFW_density_profile(wavevector,virial_Mass,redshift) ! Mpc

    end function lensing_potential

    function FT_NFW_density_profile(k,M,z)    ! It computes the normalized Fourier transform of the NFW density profile

        Implicit none    ! Eq. (9) in "Clustering of submillimetre galaxies in a self-regulated baryon collapse model" by 
                         ! Xia et al. M is the virial mass. The Fourier transform is truncated to the virial radius
        Real*8 :: k,M,z,FT_NFW_density_profile,Si1,Si2,Ci1,Ci2
        Real*8 :: alpha     ! It determines upper limit in Eq. (2.10) of 1312.4525
        Real*8,parameter :: alphafactor = 1.037d0!1.045d0!1.060d0!1.037d0

        alpha = alphafactor*concentration_mass_virial(M,z)

        call trigint((1.d0 + alpha)*k,Ci1,Si1)

        call trigint(k,Ci2,Si2)

        FT_NFW_density_profile = -sin(alpha*k)/(1.d0+alpha)/k + sin(k)*( Si1 - Si2 ) &
        + cos(k)*( Ci1 - Ci2 ) 

    end function FT_NFW_density_profile

    subroutine trigint(x,ci,si)

        use fiducial
        Implicit none

        Real*8 :: x,ci,si,f,g1,x2,y

        x2 = x**2
        y  = 1./x2

        ! ... formulas from wikipedia page http://en.wikipedia.org/wiki/Trigonometric_integral ...

        If (x.le.4.) then

            Si = x*(1. + x2*(-4.54393409816329991e-2 + x2*(1.15457225751016682e-3 + &
            x2*(-1.41018536821330254e-5 + x2*(9.43280809438713025e-8 + &
            x2*(-3.53201978997168357e-10 + x2*(7.08240282274875911e-13 + &
            x2*(-6.05338212010422477e-16)))))))) &
            /(1. + x2*(1.01162145739225565e-2 + x2*(4.99175116169755106e-5 + &
            x2*(1.55654986308745614e-7 + x2*(3.28067571055789734e-10 + &
            x2*(4.5049097575386581e-13 + x2*(3.21107051193712168e-16)))))))

            Ci = 0.577215664901532861 + log(x) + x2*(-0.25 + x2*(7.51851524438898291e-3 + &
            x2*(-1.27528342240267686e-4 + x2*(1.05297363846239184e-6 + &
            x2*(-4.68889508144848019e-9 + x2*(1.06480802891189243e-11 + &
            x2*(-9.93728488857585407e-15))))))) &
            /(1. + x2*(1.1592605689110735e-2 + x2*(6.72126800814254432e-5 + &
            x2*(2.55533277086129636e-7 + x2*(6.97071295760958946e-10 + &
            x2*(1.38536352772778619e-12 + x2*(1.89106054713059759e-15 + &
            x2*(1.39759616731376855e-18))))))))

        Else

            f  = (1. + y*(7.44437068161936700618e2 + y*(1.96396372895146869801e5 + &
            y*(2.37750310125431834034e7 + y*(1.43073403821274636888e9 + &
            y*(4.33736238870432522765e10 + y*(6.40533830574022022911e11 + &
            y*(4.20968180571076940208e12 + y*(1.00795182980368574617e13 + &
            y*(4.94816688199951963482e12 + y*(-4.94701168645415959931e11))))))))))) &
            /(x*(1. + y*(7.46437068161927678031e2 + y*(1.97865247031583951450e5 + &
            y*(2.41535670165126845144e7 + y*(1.47478952192985464958e9 + &
            y*(4.58595115847765779830e10 + y*(7.08501308149515401563e11 + &
            y*(5.06084464593475076774e12 + y*(1.43468549171581016479e13 + &
            y*(1.11535493509914254097e13)))))))))))

            g1 = y*(1. + y*(8.1359520115168615e2 + y*(2.35239181626478200e5 + &
            y*(3.12557570795778731e7 + y*(2.06297595146763354e9 + &
            y*(6.83052205423625007e10 + y*(1.09049528450362786e12 + &
            y*(7.57664583257834349e12 + y*(1.81004487464664575e13 + &
            y*(6.43291613143049485e12 + y*(-1.36517137670871689e12))))))))))) &
            /(1. + y*(8.19595201151451564e2 + y*(2.40036752835578777e5 + &
            y*(3.26026661647090822e7 + y*(2.23355543278099360e9 + &
            y*(7.87465017341829930e10 + y*(1.39866710696414565e12 + &
            y*(1.17164723371736605e13 + y*(4.01839087307656620e13 + &
            y*(3.99653257887490811e13))))))))))

            Si = 0.5*Pi - f*cos(x) - g1*sin(x)

            Ci = f*sin(x) - g1*cos(x)

        End If

        return

    end subroutine trigint

    subroutine compute_lensing_potential(halo_definition)

        use arrays
        use fiducial
        use omp_lib
        Implicit none

        Integer*4 :: indexl,indexM,indexz
        Character(len=*) :: halo_definition 

        open(15,file='./precomputed_quantities/lensing_potential/lensing_potential.dat')

        If (halo_definition .eq. 'critical_density') then

           write(15,*) '#  l  virial mass [solar mass]  red-shift  phi'

           !            !$omp Parallel Do Shared(philMz)        

           !            Do indexl=1,number_of_l

           !                Do indexM=1,number_of_M

           !                    Do indexz=1,number_of_z

           !                        philMz(indexl,indexM,indexz) = lensing_potential(indexM,indexz,indexl,'critical_density')

           !                        write(15,'(i5,3es18.10)') ml(indexl),M(indexM),z(indexz),philMz(indexl,indexM,indexz)

           !                    End Do

           !                End Do

           !            End Do

           !            !$omp End Parallel Do

           close(15)

        Else If (halo_definition .eq. 'mean_background') then 

           write(15,*) '#  l  virial mass [solar mass]  red-shift  phi'

           !            !$omp Parallel Do Shared(philMz)        

           !            Do indexl=1,number_of_l

           !                Do indexM=1,number_of_M

           !                    Do indexz=1,number_of_z

           !                        philMz(indexl,indexM,indexz) = lensing_potential(indexM,indexz,indexl,'mean_background')

           !                        write(15,'(i5,3es18.10)') ml(indexl),M(indexM),z(indexz),philMz(indexl,indexM,indexz)

           !                    End Do

           !                End Do

           !            End Do

           !            !$omp End Parallel Do

           close(15)

        Else If (halo_definition .eq. 'virial') then

           write(15,*) '# Lensing potential file. Number of lines is ',number_of_l*number_of_M*number_of_z

           write(15,*) '#  index_of_l    l   index_of_M    virial_mass[solar mass]    index_of_z    red-shift    phi'

           !$omp Parallel Do Shared(philMz,z,Scrit,M)        

           Do indexl=1,number_of_l

              Do indexM=1,number_of_M

                 Do indexz=1,number_of_z

                    philMz(indexl,indexM,indexz) = lensing_potential(M(indexM),z(indexz),indexl,'virial')/&
                         Scrit(indexz)

                 End Do

              End Do

           End Do

           !$omp End Parallel Do

           Do indexl=1,number_of_l

              Do indexM=1,number_of_M

                 Do indexz=1,number_of_z

                    write(15,'(3i10,es18.10,i5,2es18.10)') indexl,ml(indexl),indexM,M(indexM),indexz,&
                         z(indexz),philMz(indexl,indexM,indexz)

                 End Do

              End Do

           End Do

           close(15)

        End If

    end subroutine compute_lensing_potential

    subroutine read_philMz()

        use arrays
        use fiducial
        Implicit none

        Integer*4 :: ll,q,iz,iM,il
        Real*8 :: MM,zz,tphilMz

        open(15,file='./precomputed_quantities/lensing_potential/lensing_potential.dat')

        read(15,*)

        read(15,*)

        Do q=1,number_of_l*number_of_M*number_of_z

            read(15,'(3i10,es18.10,i5,2es18.10)') il,ll,iM,MM,iz,zz,tphilMz
     
            philMz(il,iM,iz) = tphilMz

        End Do

        close(15)

    end subroutine read_philMz

    !subroutine Interpolate_philMz()

    !  use arrays
    !  use fiducial
      !use omp_lib
    !  Implicit none

     ! Integer*4 :: indexM,indexz,indexl

      !Do indexl=1,number_of_l

         !!!!$omp Parallel Do Shared(philMz_interpolation,M_functions,z_functions,M,z,philMz)

       !  Do indexM=1,number_of_M_functions

        !    Do indexz=1,number_of_z_functions

         !      call Interpolate_2D(philMz_interpolation(indexl,indexM,indexz),M_functions(indexM),z_functions(indexz),&
          !          M(1:number_of_M),z(1:number_of_z),philMz(indexl,1:number_of_M,1:number_of_z))

           ! End Do

         !End Do

         !!!!$omp End Parallel Do

      !End Do

    !end subroutine Interpolate_philMz

    function window_function(k,R) ! Equation (7) of "Cosmology from the Thermal Sunyaev-Zel'dovich Power spectrum: Primordial 
                              ! non-Gaussianity and massive neutrinos [1303.4726]. Units of k : 1/Mpc. Units of R : Mpc. Both k and R 
        Implicit none             ! must by given in comoving or physical coordinates. The window function is dimensionless.

        Real*8 :: k,R,window_function,x

        x = k*R    !    Dimensionless 

        window_function = 3.d0*(sin(x)/x -cos(x))/x**2

    end function window_function

    function d_window_function(k,R) ! Derivative of the window function above. Units of k: 1/Mpc. Units of R : Mpc. Dimensionless. Both

        Implicit none    ! k and R must be given in comoving or physical coordinates.

        Real*8 :: k,R,d_window_function,x

        x = k*R

        d_window_function = 3.d0*(3.d0*cos(x)/x + sin(x)*(1.d0 - 3.d0/x**2))/x**2

    end function d_window_function

    subroutine compute_transfer_function()

        use arrays 
        use fiducial
        Implicit none

        Integer*4 :: indexk

        open(15,file='./output/transfer_function.dat')

        write(15,*) '# k        T(k) '

        Do indexk=1,number_of_k

            write(15,'(2es18.10)') k(indexk),Transfer_function(k(indexk))

        End Do

        close(15)

    end subroutine compute_transfer_function

    function Transfer_function(rk) ! Transfer function from "Baryonic features in the matter transfer function" by  D. Eisenstein and W. Hu.

        use fiducial               ! Units of k : 1/Mpc. Transfer function is dimensionless.
        Implicit none

        Real*8 :: Transfer_function,rk,thet,omegab,obm,b1,b2,zd,ze,rd,re,rke,s,rks,q,a1,a2,ac,bc,f
        Real*8 :: c1,c2,tc,y,g1,ab,bn,ss,bb,tb

        ! constants

        omegab = Omega_b_h2/h**2   

        obm = omegab/Omega_m0

        thet = tcmb0/2.7

        ! Equation 4 - redshift of drag epoch

        b1 = 0.313d0*((Omega_m0*h**2)**(-0.419d0))*(1.d0+0.607d0*((Omega_m0*h**2)**(0.674d0)))

        b2 = 0.238d0*((Omega_m0*h**2)**(0.223d0))

        zd = 1291.d0*(1.d0+b1*((Omega_b_h2)**b2))*((Omega_m0*h**2)**0.251d0)/(1.d0+0.659d0*((Omega_m0*h**2)**(0.828d0)))

        ! Equation 2 - redshift of matter-radiation equality

        ze = 2.50d4*Omega_m0*h**2/thet**4

        ! Value of R=(ratio of baryon-photon momentum density) at drag epoch (Eq. 5)

        rd = 31500.d0*Omega_b_h2/thet**4/zd

        ! Value of R=(ratio of baryon-photon momentum density) at epoch of matter-radiation equality

        re = 31500.d0*Omega_b_h2/thet**4/ze

        ! Equation 3 - scale of ptcle horizon at matter-radiation equality

        rke = 7.46d-2*Omega_m0*h**2/thet**2

        ! Equation 6 - sound horizon at drag epoch

        s = (2.d0/3.d0/rke)*sqrt(6.d0/re)*log((sqrt(1.d0+rd)+sqrt(rd+re))/(1.d0+sqrt(re)))

        ! Equation 7 - silk damping scale

        rks = 1.6d0*((Omega_b_h2)**(0.52d0))*((Omega_m0*h**2)**(0.73d0))*(1.d0+((10.4d0*Omega_m0*h**2)**(-0.95d0)))

        ! Equation 10  - define q

        q = rk/13.41d0/rke

        ! Equations 11 - CDM transfer function fits

        a1 = ((46.9d0*Omega_m0*h**2)**(0.670d0))*(1.d0+((32.1d0*Omega_m0*h**2)**(-0.532d0)))

        a2 = ((12.0d0*Omega_m0*h**2)**(0.424d0))*(1.d0+((45.0d0*Omega_m0*h**2)**(-0.582d0)))

        ac = (a1**(-obm))*(a2**(-(obm)**3))

        ! Equations 12 - CDM transfer function fits

        b1 = 0.944d0/(1.d0+((458.d0*Omega_m0*h**2)**(-0.708d0)))

        b2 = ((0.395d0*Omega_m0*h**2)**(-0.0266d0))

        bc = 1.d0/(1.d0+b1*(((1.d0-obm)**b2)-1.d0))

        ! Equation 18

        f = 1.d0/(1.d0+((rk*s/5.4d0)**4))

        ! Equation 20

        c1 = 14.2d0 + 386.d0/(1.d0+69.9d0*(q**(1.08d0)))

        c2 = 14.2d0/ac + 386.d0/(1.d0+69.9d0*(q**(1.08d0)))

        ! Equation 17 - CDM transfer function

        tc = f*log(exp(1.d0)+1.8d0*bc*q)/(log(exp(1.d0)+1.8d0*bc*q)+c1*q**2)+(1.d0-f)*&
        log(exp(1.d0)+1.8d0*bc*q)/(log(exp(1.d0)+1.8d0*bc*q)+c2*q**2)

        ! Equation 15

        y = (1.d0+ze)/(1.d0+zd)

        g1 = y*(-6.d0*sqrt(1.d0+y)+(2.d0+3.d0*y)*log((sqrt(1.d0+y)+1.d0)/(sqrt(1.d0+y)-1.d0)))

        ! Equation 14

        ab = g1*2.07d0*rke*s/((1.d0+rd)**(0.75d0))

        ! Equation 23

        bn = 8.41d0*((Omega_m0*h**2)**(0.435d0))

        ! Equation 22

        ss = s/((1.d0+((bn/rk/s)**3))**(1.d0/3.d0))

        ! Equation 24

        bb = 0.5d0+(obm) + (3.d0-2.d0*obm)*sqrt(((17.2d0*Omega_m0*h**2)**2)+1.d0)

        ! Equations 19 & 21

        tb = log(exp(1.d0)+1.8d0*q)/(log(exp(1.d0)+1.8d0*q)+c1*q*q)/(1.d0+((rk*s/5.2d0)**2))

        tb=(tb+ab*exp(-((rk/rks)**(1.4d0)))/(1.d0+((bb/rk/s)**3)))*sin(rk*ss)/rk/ss

        ! Equation 8

        Transfer_function = (obm)*tb+(1.-obm)*tc

    end function Transfer_function

    function growth_function(z)    ! Dimensionless. Equation (A4) in "Baryonic features in the matter transfer function" by D. J. Eisenstein and 

        Implicit none    ! W. Hu.

        Real*8 :: growth_function,z,d

        d = Omega_m(z)**(4.d0/7.d0) - Omega_L(z) + (1.d0 + Omega_m(z)/2.d0)*(1.d0 + Omega_L(z)/7.d1)

        growth_function = 5.d0*Omega_m(z)/2.d0/(1.d0 + z)/d

    end function growth_function

    function Delta_squared_non_normalised(k) ! Equation (A1) without constant \delta_H in "Baryonic features in the matter transfer function"

        use fiducial                         ! by Eisenstein and Hu. Units of k : 1/Mpc. Delta_squared_non_normalised is dimensionless
        Implicit none                        ! 

        Real*8 :: Delta_squared_non_normalised,k

        Delta_squared_non_normalised = (c*k/H_0)**(3.d0 + ns)*Transfer_function(k)**2

    end function Delta_squared_non_normalised

!    subroutine compute_normalization() ! Compute normalisation constant in the linear matter power spectrum to agree 

 !       use fiducial                   ! with \sigma_8 in fiducial model 
  !      Implicit none

   !     Real*8 :: sum,x1,x2,f1,f2
    !    Integer*4 :: indexk
     !   Integer*4,parameter :: number_of_k_logscale = 1000
      !  Integer*4,parameter :: intervals = number_of_k_logscale - 1 
       ! Real*8,dimension(number_of_k_logscale) :: f,k
!        Real*8,parameter :: R = 8.d0/h ! Units : Mpc
 !       Real*8,parameter :: kmaxlogscale = 1.d-1
  !      Real*8,parameter :: stepsize = 1.d-3
   !     Integer*4,parameter :: max_iterations = 1000000000

        ! Wavevector array. Units : 1/Mpc
    !    Do indexk = 1, number_of_k_logscale  
  
     !       k(indexk) = 10**(log10(kmin) + real(indexk-1)*(log10(kmaxlogscale) - log10(kmin))/real(number_of_k_logscale-1))

      !  End Do

       ! !$omp Parallel Do Shared(f)

        !Do indexk=1,number_of_k_logscale

!            f(indexk) = Delta_squared_non_normalised(k(indexk))*(3.d0*(sin(k(indexk)*R)-&
 !           k(indexk)*R*cos(k(indexk)*R))/k(indexk)**3/R**3)**2/k(indexk)

  !      End Do

   !     !$omp End Parallel Do

    !    sum = 0.d0

     !   Do indexk=1,intervals

      !      sum = (k(indexk+1) - k(indexk))/2.d0*( f(indexk) + f(indexk+1) ) + sum 

       ! End Do

!        x1 = k(number_of_k_logscale)

 !       f1 = f(number_of_k_logscale)

  !      x2 = x1 + stepsize

   !     f2 = Delta_squared_non_normalised(x2)*(3.d0*(sin(x2*R)- x2*R*cos(x2*R))/x2**3/R**3)**2/x2

    !    sum = (x2-x1)/2.d0*( f1 + f2 )+ sum
        
     !   Do indexk=1,max_iterations

      !      x1 = x2

       !     f1 = f2

        !    If (x2 .gt. kmax) then

!                x2 = kmax

 !               f2 = Delta_squared_non_normalised(x2)*(3.d0*(sin(x2*R)- x2*R*cos(x2*R))/x2**3/R**3)**2/x2

  !          Else If (x2 .eq. kmax) then

   !             exit

    !        Else

     !           x2 = x1 + stepsize

      !          f2 = Delta_squared_non_normalised(x2)*(3.d0*(sin(x2*R)- x2*R*cos(x2*R))/x2**3/R**3)**2/x2
            
       !     End IF

        !    sum = (x2-x1)/2.d0*( f1 + f2 )+ sum
            
         !   If (indexk .eq. max_iterations) then

!                print *,'Maximum number of iterations in integral of Normalization achieved.'

 !               stop

  !          End If

   !     End Do

    !    Normalization = sigma8**2/sum    !    

!    end subroutine compute_normalization

    function matter_power_spectrum(k,z)    !    P(k,z) in Equation (3.15) of A. Lewis and A. Challinor, Physics Reports 429 (2006) 1-65

        use fiducial    !    P(k,z) is normalized with fiducial \sigma_8. It corresponds to P(k) in equation (A1) of "Baryonic features in the 
        Implicit none

        Real*8 :: matter_power_spectrum,k,z    ! matter transfer function" with red-shift evolution added. Units of k : 1/Mpc. Units : (Mpc/h)**3. 

        matter_power_spectrum = Normalization*2.d0*Pi**2*Delta_squared_non_normalised(k)/&
        k**3*growth_function(z)**2/growth_function(0.d0)**2*h**3

    end function matter_power_spectrum

    function power_spectrum_weyl_potential(k,z)    !    It computes power spectrum of the Weyl potential, equation (3.15) in published version 

        use fiducial    !    of review by Lewis and Challinor. There is an additional factor to correct for comoving wavevector. Units of k : 1/Mpc. 
        Implicit none    !    Units : Dimensionless

        Real*8 :: power_spectrum_weyl_potential,k,z 
    
        power_spectrum_weyl_potential = 9.d0*Omega_m(z)**2*Hubble_parameter(z)**4*&
        matter_power_spectrum(k,z)/k/8.d0/Pi**2/c**4

        power_spectrum_weyl_potential = power_spectrum_weyl_potential/h**3

    end function power_spectrum_weyl_potential

    subroutine compute_matter_power_spectrum_at_z(indexz)

        use arrays 
        use fiducial
        Implicit none

        Integer*4 :: indexk,indexz

        open(15,file='./output/P_of_k.dat')

        write(15,*) '# k        P(k) at red-shift', z(indexz)

        Do indexk=1,number_of_k

            write(15,'(2es18.10)') k(indexk),matter_power_spectrum(k(indexk),z(indexz))

        End Do

        close(15)

    end subroutine compute_matter_power_spectrum_at_z

    function sigma_squared(Mass,redshift) ! Equation (6) of 1303.4726. Units of Mass: solar mass. Units : dimensionless

        use fiducial
        Implicit none

        Real*8 :: sigma_squared,R,sum,prefactor,Mass,x1,x2,f1,f2,redshift
        Integer*4 :: indexk
        Integer*4,parameter :: number_of_k_logscale = 100
        Integer*4,parameter :: intervals = number_of_k_logscale - 1 
        Real*8,dimension(number_of_k) :: f,kf
        Real*8 :: stepsize 
        Integer*4,parameter :: max_iterations = 1000000000 

        prefactor = 1.d0/2.d0/Pi**2

        R = (3.d0*Mass/4.d0/Pi/mean_density(redshift)*(1.d0 + redshift)**3.d0)**(1.d0/3.d0) ! Units : Mpc. (1+z)**3 factor to use comoving coordinates

        stepsize = (kmax - 1.d0/R)*1.d-7

        ! Wavevector array. Units : 1/Mpc
        Do indexk = 1, number_of_k_logscale  
  
            kf(indexk) = 10**(log10(kmin) + real(indexk-1)*(log10(1.d0/R) - log10(kmin))/real(number_of_k_logscale-1))

        End Do

        Do indexk=1,number_of_k_logscale

            f(indexk) = kf(indexk)**2*window_function(kf(indexk),R)**2*matter_power_spectrum(kf(indexk),redshift)

        End Do

        sum = 0.d0

        Do indexk=1,intervals

            sum = (kf(indexk+1) - kf(indexk))/2.d0*( f(indexk) + f(indexk+1) ) + sum 

        End Do

        x1 = kf(number_of_k_logscale)

        f1 = f(number_of_k_logscale)

        x2 = x1 + stepsize

        f2 = x2**2*window_function(x2,R)**2*matter_power_spectrum(x2,redshift)

        sum = (x2-x1)/2.d0*(f1 + f2) + sum

        Do indexk=1,max_iterations

            x1 = x2

            f1 = f2

            If (x2 .gt. kmax) then

                x2 = kmax

                f2 = x2**2*window_function(x2,R)**2*matter_power_spectrum(x2,redshift)

            Else If (x2 .eq. kmax) then

                exit

            Else

                x2 = x1 + stepsize

                f2 = x2**2*window_function(x2,R)**2*matter_power_spectrum(x2,redshift)
            
            End IF

            sum = (x2-x1)/2.d0*( f1 + f2 ) + sum
            
            If (indexk .eq. max_iterations) then

                print *,'Maximum number of iterations in integral of sigma_squared achieved.'

                stop

            End If

        End Do
        
        sigma_squared = prefactor*sum/h**3

    end function sigma_squared

    function dsigma_squared_dR(Mass,redshift) ! Derivative w.r.t R of Eq. (6) in 1303.4726. Units of Mass: solar mass. Units : 1/Mpc

        use fiducial
        Implicit none

        Real*8 :: dsigma_squared_dR,R,sum,prefactor,Mass,x1,x2,f1,f2,redshift
        Integer*4 :: indexk
        Integer*4,parameter :: number_of_k_logscale = 100
        Integer*4,parameter :: intervals = number_of_k_logscale - 1 
        Real*8,dimension(number_of_k_logscale) :: f,kf
        Real*8 :: stepsize 
        Integer*4,parameter :: max_iterations = 1000000000

        prefactor = 1.d0/Pi**2

        R = (3.d0*Mass/4.d0/Pi/mean_density(redshift)*(1.d0 + redshift)**3.d0)**(1.d0/3.d0) ! Units : Mpc. (1+z)**3 factor to use comoving coordinates

        stepsize = (kmax - 1.d-1/R)*1.d-3

        ! Wavevector array. Units : 1/Mpc
        Do indexk = 1, number_of_k_logscale    

            kf(indexk) = 10**(log10(kmin) + real(indexk-1)*(log10(1.d-1/R) - log10(kmin))/real(number_of_k_logscale-1))

        End Do

        Do indexk=1,number_of_k_logscale

            f(indexk) = kf(indexk)**3*window_function(kf(indexk),R)*d_window_function(kf(indexk),R)&
            *matter_power_spectrum(kf(indexk),redshift)

        End Do

        sum = 0.d0

        Do indexk=1,intervals

            sum = (kf(indexk+1) - kf(indexk))/2.d0*( f(indexk) + f(indexk+1) ) + sum 

        End Do

        x1 = kf(number_of_k_logscale)

        f1 = f(number_of_k_logscale)

        x2 = x1 + stepsize

        f2 = x2**3*window_function(x2,R)*d_window_function(x2,R)*matter_power_spectrum(x2,redshift)

        sum = (x2-x1)/2.d0*(f1 + f2) + sum

        Do indexk=1,max_iterations

            x1 = x2

            f1 = f2

            If (x2 .gt. kmax) then

                x2 = kmax

                f2 = x2**3*window_function(x2,R)*d_window_function(x2,R)*matter_power_spectrum(x2,redshift)

            Else If (x2 .eq. kmax) then

                exit

            Else

                x2 = x1 + stepsize

                f2 = x2**3*window_function(x2,R)*d_window_function(x2,R)*matter_power_spectrum(x2,redshift)
            
            End IF

            sum = (x2-x1)/2.d0*( f1 + f2 ) + sum
            
            If (indexk .eq. max_iterations) then

                print *,'Maximum number of iterations in integral of dsigma_squared achieved.'

                stop

            End If

        End Do

        dsigma_squared_dR = prefactor*sum/h**3

    end function dsigma_squared_dR



!    subroutine compute_dM200ddM_M_z()

 !     use fiducial
  !    use arrays
   !   use omp_lib
    !  Implicit none

     ! Integer*4 :: indexM,indexz

!      If (compute_functions) then

 !        print *, 'NEED TO IMPLEMENT DERIVATIVES ROUTINE FOR MASSES'
         
  !       stop

   !   Else

    !     !$omp Parallel Do Shared(dM200ddM_M_z,M_functions,z_functions,M,z,dM200ddM)

     !    Do indexz=1,number_of_z_functions 

      !      Do indexM=1,number_of_M_functions

       !        call Interpolate_2D(dM200ddM_M_z(indexM,indexz),M_functions(indexM),z_functions(indexz),&
        !            M(1:number_of_M),z(1:number_of_z),dM200ddM(1:number_of_M,1:number_of_z))

         !   End Do

!         End Do

 !        !$omp End Parallel Do

  !    End If

   ! end subroutine compute_dM200ddM_M_z


    subroutine read_alpha_halo_mass_function()

        use arrays
        use fiducial
        Implicit none

        Integer*4 :: index,iz
        Real*8 :: zz,ts2

        open(15,file='./precomputed_quantities/alpha_halo_mass_function.dat')

        read(15,*)

        read(15,*)

        Do index=1,number_of_z

            read(15,'(i5,2es18.10)') iz,zz,ts2

            alpha_halo_mass_function(iz) = ts2

        End Do

        close(15)

    end subroutine read_alpha_halo_mass_function

    subroutine read_sigma_square_M200d()

        use arrays
        use fiducial
        Implicit none

        Integer*4 :: index,iM,iz
        Real*8 :: MM,zz,ts2,tds2

        open(15,file='./precomputed_quantities/sigma_square_M200d.dat')

        read(15,*)

        read(15,*)

        Do index=1,number_of_M*number_of_z

            read(15,'(i10,es18.10,i5,3es18.10)') iM,MM,iz,zz,ts2,tds2

            sigma_square_M200d(iM,iz) = ts2

            dsigma_square_M200d(iM,iz) = tds2 

        End Do

        close(15)

    end subroutine read_sigma_square_M200d

    subroutine read_dndM()

        use arrays
        use fiducial
        Implicit none

        Integer*4 :: index,iM,iz
        Real*8 :: MM,zz,tdndM

        open(15,file='./precomputed_quantities/dndM/dndM.dat')

        read(15,*)

        read(15,*)

        Do index=1,number_of_M*number_of_z

            read(15,'(i10,es18.10,i5,2es18.10)') iM,MM,iz,zz,tdndM

            dndM(iM,iz) = tdndM

        End Do

        close(15)

    end subroutine read_dndM

!    subroutine Interpolate_dndM()

 !     use arrays
  !    use fiducial
      !use omp_lib
   !   Implicit none

    !  Integer*4 :: indexM,indexz

      !!!!$omp Parallel Do Shared(dndM_interpolation,M_functions,z_functions,M,z,dndM)

     ! Do indexM=1,number_of_M_functions

      !   Do indexz=1,number_of_z_functions

       !     call Interpolate_2D(dndM_interpolation(indexM,indexz),M_functions(indexM),z_functions(indexz),&
        !         M(1:number_of_M),z(1:number_of_z),dndM(1:number_of_M,1:number_of_z))

         !End Do

!      End Do

      !!!!$omp End Parallel Do

 !   end subroutine Interpolate_dndM

    subroutine write_dndM_at_z(indexz)

        use arrays
        use fiducial
        Implicit none

        Integer*4 :: indexM,indexz

        open(15,file='./output/dndM_at_z.dat')

        write(15,*) '# Virial_ Mass[solar mass]  dndM_at_redshift',z(indexz)

        Do indexM=1,number_of_M

            write(15,'(2es18.10)') M(indexM), M(indexM)**2*dndM(indexM,indexz)/mean_density(z(indexz))

        End Do

        close(15)

    end subroutine write_dndM_at_z



!    function linear_halo_bias(virial_Mass,redshift)    ! From table 2 and equation (6) of "The large-scale bias of dark matter halos:  

 !       use arrays    ! numerical calibration and model tests". It computes the linear halo bias. (See paper for details about parameters). 
  !      use fiducial  ! Mass given must be M_{200d}. Dimensionless. 
   !     Implicit none    

    !    Real*8 :: linear_halo_bias,sigma,nu,virial_Mass,redshift,M200d_M_z
     !   Real*8,parameter :: y = log10(DeltaSO)
      !  Real*8,parameter :: A = 1.d0 + 2.4d-1*y*dexp(-(4.d0/y)**4)
       ! Real*8,parameter :: aa = 4.4d-1*y - 8.8d-1
       ! Real*8,parameter :: B = 1.83d-1
       ! Real*8,parameter :: bb = 1.5d0 
       ! Real*8,parameter :: CC = 1.9d-2 + 1.07d-1*y + 1.9d-1*dexp(-(4.d0/y)**4)
       ! Real*8,parameter :: ccc = 2.4d0
       ! Real*8,parameter :: delta_c = 1.686d0
        
       ! call Interpolate_2D(M200d_M_z,virial_Mass,redshift,M(1:number_of_M),z(1:number_of_z),M200d(1:number_of_M,1:number_of_z))

       ! sigma = sqrt(sigma_squared(M200d_M_z,redshift))    ! M must be the SO mass M_{200d} 
!        sigma = sqrt(sigma_square_M200d(indexM,indexz))    ! M must be the SO mass M_{200d} 

       ! nu = delta_c/sigma

       ! linear_halo_bias = 1.d0 - A*nu**aa/(nu**aa + delta_c**aa) + B*nu**b + CC*nu**ccc

    !end function linear_halo_bias

    function linear_halo_bias(indexM_virial,index_redshift)    ! From table 2 and equation (6) of "The large-scale bias of dark matter halos:  

        use arrays    ! numerical calibration and model tests". It computes the linear halo bias. (See paper for details about parameters). 
        use fiducial  ! Mass given must be M_{200d}. Dimensionless. 

        Implicit none    

        Real*8 :: linear_halo_bias,nu,sigma
        Real*8,parameter :: delta_c = 1.686d0
        Integer*4 :: indexM_virial,index_redshift
        
!        sigma = sqrt(sigma_squared(M200d(indexM_virial,index_redshift),z(index_redshift)))    ! M must be the SO mass M_{200d} 
        sigma = sqrt(sigma_square_M200d(indexM_virial,index_redshift))    ! M must be the SO mass M_{200d} 

        nu = delta_c/sigma

        linear_halo_bias = linear_halo_bias_b_nu(nu)

    end function linear_halo_bias

    function linear_halo_bias_b_nu(nu)

      use fiducial

      Implicit none

      Real*8 :: linear_halo_bias_b_nu,nu
      Real*8,parameter :: y = log10(DeltaSO)
      Real*8,parameter :: A = 1.d0 + 2.4d-1*y*exp(-(4.d0/y)**4)
      Real*8,parameter :: aa = 4.4d-1*y - 8.8d-1
      Real*8,parameter :: B = 1.83d-1
      Real*8,parameter :: bb = 1.5d0 
      Real*8,parameter :: CC = 1.9d-2 + 1.07d-1*y + 1.9d-1*exp(-(4.d0/y)**4)
      Real*8,parameter :: ccc = 2.4d0
      Real*8,parameter :: delta_c = 1.686d0

      linear_halo_bias_b_nu = 1.d0 - A*nu**aa/(nu**aa + delta_c**aa) + B*nu**b + CC*nu**ccc

    end function linear_halo_bias_b_nu

    subroutine compute_bMz()

      use arrays
      use fiducial
      use omp_lib

      Implicit none

      Integer*4 :: indexz,indexM

      open(15,file='./precomputed_quantities/bMz/bMz.dat')

      write(15,*) '# Linear bias file. Number of lines is ',number_of_M*number_of_z

      write(15,*) '# index_of_M    SO_virial_mass[solar mass]   index_of_z    red-shift  b '

      !$omp Parallel Do Shared(bMz)

      Do indexM=1,number_of_M

         Do indexz=1,number_of_z

            bMz(indexM,indexz) = linear_halo_bias(indexM,indexz)

         End Do

      End Do

      !$omp End Parallel Do

      Do indexM=1,number_of_M

         Do indexz=1,number_of_z

            write(15,'(i10,es18.10,i5,2es18.10)') indexM,M(indexM),indexz,z(indexz),bMz(indexM,indexz)

         End Do

      End Do

      close(15)

    end subroutine compute_bMz


    subroutine read_bMz()

        use arrays
        use fiducial
        Implicit none

        Integer*4 :: iz,iM,q
        Real*8 :: MM,zz,tbMz

        open(15,file='./precomputed_quantities/bMz/bMz.dat')
        
        read(15,*)

        read(15,*)

        Do q=1,number_of_M*number_of_z

            read(15,'(i10,es18.10,i5,2es18.10)') iM,MM,iz,zz,tbMz

            bMz(iM,iz) = tbMz

        End Do

        close(15)

    end subroutine read_bMz

!    subroutine Interpolate_bMz()

 !     use arrays
  !    use fiducial
   !   !use omp_lib
    !  Implicit none

     ! Integer*4 :: indexM,indexz

      !!!!$omp Parallel Do Shared(bMz_interpolation,M_functions,z_functions,M,z,bMz)

      !Do indexM=1,number_of_M_functions

       !  Do indexz=1,number_of_z_functions

        !    call Interpolate_2D(bMz_interpolation(indexM,indexz),M_functions(indexM),z_functions(indexz),&
         !        M(1:number_of_M),z(1:number_of_z),bMz(1:number_of_M,1:number_of_z))

         !End Do

      !End Do

      !!!!!$omp End Parallel Do

    !end subroutine Interpolate_bMz

    subroutine write_bMz_at_z(indexz)

        use arrays
        use fiducial
        Implicit none

        Integer*4 :: indexM,indexz

        open(15,file='./output/bMz_at_z.dat')

        write(15,*) '# Virial Mass [solar mass/h]  bMz at redshift'

        Do indexM=1,number_of_M

            write(15,'(2es18.10)') M(indexM), bMz(indexM,indexz)

        End Do

        close(15)

    end subroutine write_bMz_at_z

!    subroutine write_bnu_at_z(indexz)

!        use arrays
!        use fiducial
!        Implicit none

!        Integer*4 :: indexM,indexz
!        Real*8 :: sigma,delta_c,nu

!        open(15,file='./output/bnu_at_z.dat')
        
!        write(15,*) '# Virial Mass [solar mass/h]  bnu at redshift'

!        Do indexM=1,number_of_M

!            sigma = sqrt(sigma_squared(M(indexM),indexz))

            !    sigma = sqrt(sigma_squared(M200d(indexM,indexz),indexz))

!            delta_c = 1.686d0

!            nu = delta_c/sigma

!            write(15,'(2es18.10)') nu, bMz(indexM,indexz)

!        End Do

!        close(15)

!    end subroutine write_bnu_at_z



    subroutine write_Cl()

        use arrays
        use fiducial
        Implicit none

        Integer*4 :: indexl

        open(15,file='./output/Cl_yphi.dat')

        write(15,*) '#   l   Clyphi1h    Clyphi2h    Clyphi    Clphiphi1h    Clphiphi2h    Clphiphi    Clpsilimber'&
             //trim('Clyy1h    Clyy2h    Clyy')//''

        Do indexl=1,number_of_l

            write(15,'(i5,10es18.10)') ml(indexl), Cl1h(indexl), Cl2h(indexl), Cl(indexl),Clphiphi1h(indexl),&
                 Clphiphi2h(indexl),Clphiphi(indexl),Clpsilimber(indexl),Clyy1h(indexl),Clyy2h(indexl),Clyy(indexl)

        End Do

        close(15)

    end subroutine write_Cl


End module functions

