module functions

    Implicit none
     
contains

function Hubble_parameter(z)    !    Units : m/s/Mpc
    use fiducial
    Implicit none
    Real*8 :: Hubble_parameter,z

    Hubble_parameter = H_0*sqrt(Omega_m0*(1.d0 + z)**3 + Omega_L0)

end function Hubble_parameter

function comoving_distance(z)    !    Units : Mpc
    use fiducial
    Implicit none
    Real*8 :: comoving_distance,xmin,z
    Integer*4,parameter :: number_of_x = 1d5
    Integer*4,parameter :: intervals = number_of_x - 1
    Real*8,dimension(number_of_x) :: x,f
    Integer*4 :: m

    xmin = 1.d-10

    Do m=1,number_of_x

        x(m) = 10**(log10(xmin) + dble((m-1))*(log10(z) - log10(xmin))/dble(number_of_x-1) )

    End Do

    !$omp Parallel Do Shared(f)    

    Do m=1,number_of_x

        f(m) = 1.d0/Hubble_parameter(x(m))

    End Do

    !$omp End Parallel Do

    comoving_distance = 0.d0

    Do m = 1,intervals

        comoving_distance = (x(m+1)-x(m))/2.d0*(f(m) + f(m+1)) + comoving_distance

    End Do 

    comoving_distance = c*comoving_distance

end function comoving_distance

function angular_diameter_distance(z)    !    Units : Mpc
    Implicit none
    Real*8 :: angular_diameter_distance,z

    angular_diameter_distance = comoving_distance(z)/(1.d0 + z)

end function angular_diameter_distance

function comoving_volume_per_steradian(z)    !    Units : Mpc**3 
    use fiducial
    Implicit none
    Real*8 :: comoving_volume_per_steradian,z

    comoving_volume_per_steradian = c*(1.d0 + z)**2*angular_diameter_distance(z)**2/Hubble_parameter(z)

end function comoving_volume_per_steradian

subroutine compute_d2VdzdO()    !    It fills in the vector with comoving volume per steradian
    use arrays
    use fiducial
    Implicit none
    Integer*4 :: indexz

    !$omp Parallel Do Shared(d2VdzdO)

    Do indexz=1,number_of_z

        d2VdzdO(indexz) = comoving_volume_per_steradian(z(indexz))

    End Do

    !$omp End Parallel Do

end subroutine

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

function critical_surface_density(z) ! Units : Solar mass/Mpc**2 
    use fiducial
    Implicit none
    Real*8 :: z,critical_surface_density,czdec,cz

    czdec = comoving_distance(z_dec) ! Units Mpc

    cz = comoving_distance(z)        ! Units Mpc

    critical_surface_density = c**2*czdec*(1.d0 + z)/4.d0/Pi/G/cz/(czdec-cz)

    critical_surface_density = critical_surface_density*Mpc/M_sun 

end function critical_surface_density

subroutine compute_critical_surface_density()    !    It fills in the vector with critical surface density
    use arrays
    use fiducial
    Implicit none
    Integer*4 :: indexz 

    !$omp Parallel Do Shared(Scrit)

    Do indexz=1,number_of_z

        Scrit(indexz) = critical_surface_density(z(indexz))

    End Do

    !$omp End Parallel Do

end subroutine

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
    Integer*4,parameter :: indexr = 1d2
    Real*8 :: rmin,rmax
    Real*8,dimension(indexr) :: r

    print *, M200c(indexM,indexz)

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
    use omp_lib
    Implicit none
    Real*8 :: M,z,rdc,rmin,integral_r_delta_c
    Integer*4,parameter :: number_of_r = 1d5
    Integer*4,parameter :: intervals = number_of_r - 1
    Integer*4 :: n
    Real*8,dimension(number_of_r) :: r,f 

    rmin = 1.d-10   
 
    Do n=1,number_of_r

        r(n) = 10**(log10(rmin) + dble((n-1))*(log10(rdc) - log10(rmin))/dble(number_of_r-1) )

    End Do

    !$omp Parallel Do Shared(f)

    Do n=1,number_of_r

        f(n) = r2NFW_density_profile(r(n),M,z,'virial')

    End Do

    !$omp End Parallel Do

    integral_r_delta_c = 0.d0

    Do n = 1, intervals

        integral_r_delta_c = (r(n+1)-r(n))/2.d0*(f(n)+ f(n+1)) + integral_r_delta_c

    End Do

    integral_r_delta_c = integral_r_delta_c*4.d0*Pi

end function integral_r_delta_c

function integral_r_delta_minus_total_matter(M,z,rdc,delta)    !    Units of M : solar mass. Units of rdc : Mpc. Units : solar mass
    use fiducial
    Implicit none
    Real*8 :: M,z,rdc,delta,integral_r_delta_minus_total_matter

    integral_r_delta_minus_total_matter = integral_r_delta_c(M,z,rdc) - 4.d0*Pi*rdc**3*critical_density(z)*delta/3.d0

end function integral_r_delta_minus_total_matter

function r_delta_c_from_M_virial(M,z,delta)    !    Units of M : solar mass. Units of r_delta_c_from_M_virial : Mpc
    use fiducial
    Implicit none
    Real*8 :: r_delta_c_from_M_virial,M,z,delta,r1,r2,step,rini,f1,f2,rmid,dr,fmid
    Integer*4,parameter :: iter = 1d3
    Real*8,parameter :: racc = 1.d-12
    Integer*4 :: p
    logical :: root

    rini = (3.d0*M/4.d0/Pi/critical_density(z)/delta)**(1.d0/3.d0) ! Units : Mpc

    step = 1.6d0

    r1 = 1.d0*rini  ! Units : Mpc

    r2 = 1.1d0*rini ! Units : Mpc

    f1 = integral_r_delta_minus_total_matter(M,z,r1,delta) ! Units : solar mass

    f2 = integral_r_delta_minus_total_matter(M,z,r2,delta) ! Units : solar mass

    Do p = 1, iter

        If (f1*f2 .lt. 0.d0) then

            root = .true.

            exit

        else if (abs(f1) .lt. abs(f2)) then

            r1 = r1 + step*(r1-r2)

            f1 = integral_r_delta_minus_total_matter(M,z,r1,delta)

        else if (abs(f1) .gt. abs(f2)) then

            r2 = r2 + step*(r2-r1)

            f2 = integral_r_delta_minus_total_matter(M,z,r2,delta)

        end if

    End Do

    If (root) then

        Do p=1,iter

            If (f1 .lt. 0.d0) then 

                dr = r2 - r1

                rmid = r1 + dr/2.d0

                fmid = integral_r_delta_minus_total_matter(M,z,rmid,delta)    

                If (fmid .lt. 0.d0) then 

                    f1 = fmid

                    r1 = rmid

                else if ((fmid .eq. 0.d0) .or. (abs(dr) .lt. racc)) then

                    exit

                else 

                    f2 = fmid

                    r2 = rmid

                end if

            else if (f2 .lt. 0.d0) then

                dr = r1 - r2

                rmid = r2 + dr/2.d0

                fmid = integral_r_delta_minus_total_matter(M,z,rmid,delta)    

                If (fmid .lt. 0.d0) then 

                    f2 = fmid

                    r2 = rmid

                else if ((fmid .eq. 0.d0) .or. (abs(dr) .lt. racc)) then

                    exit

                else

                    f1 = fmid

                    r1 = rmid

                end if

            End if

        End Do

        r_delta_c_from_M_virial = rmid

    else 

        print *,'Maximum number of iterations achieved without finding root for function '

        stop 

    End If

end function r_delta_c_from_M_virial

function M_delta_c_from_M_virial(z,rdc,delta)    !    Units of M and M_delta_c_from_M_virial: solar mass
    use fiducial
    Implicit none
    Real*8 :: rdc,z,delta,M_delta_c_from_M_virial

    M_delta_c_from_M_virial = 4.d0*Pi*rdc**3*critical_density(z)*delta/3.d0

end function M_delta_c_from_M_virial


function integral_r_delta_d(M,z,rdd)    !    Similar to 'integral_r_delta_c' above. Units of M : solar mass. Units of rdd : Mpc. Units: solar mass
    use fiducial
    Implicit none
    Real*8 :: M,z,rdd,rmin,integral_r_delta_d
    Integer*4,parameter :: number_of_r = 1d5
    Integer*4,parameter :: intervals = number_of_r - 1
    Integer*4 :: n
    Real*8,dimension(intervals+1) :: r,f

    rmin = 1.d-10
    
    Do n=1,number_of_r

        r(n) = 10**(log10(rmin) + dble((n-1))*(log10(rdd) - log10(rmin))/dble(number_of_r-1) )

    End Do

    !$omp Parallel Do Shared(f)

    Do n=1,number_of_r

        f(n) = r2NFW_density_profile(r(n),M,z,'virial')

    End Do

    !$omp End Parallel Do

    integral_r_delta_d = 0.d0

    Do n = 1, intervals

        integral_r_delta_d = (r(n+1)-r(n))/2.d0*(f(n) + f(n+1)) + integral_r_delta_d

    End Do

    integral_r_delta_d = integral_r_delta_d*4.d0*Pi

end function integral_r_delta_d

function integral_r_delta_d_minus_total_matter(M,z,rdd,delta)    !    Units of M : solar mass. Units of rdd : Mpc. Units : solar mass
    use fiducial
    Implicit none
    Real*8 :: M,z,rdd,delta,integral_r_delta_d_minus_total_matter

    integral_r_delta_d_minus_total_matter = integral_r_delta_d(M,z,rdd) - 4.d0*Pi*rdd**3*mean_density(z)/(1.d0+z)**3.d0*delta/3.d0

end function integral_r_delta_d_minus_total_matter

function r_delta_d_from_M_virial(M,z,delta)    !    Units of M : solar mass. Units of r_delta_d_from_M_virial : Mpc
    use fiducial
    Implicit none
    Real*8 :: r_delta_d_from_M_virial,M,z,delta,r1,r2,step,rini,f1,f2,rmid,dr,fmid
    Integer*4,parameter :: iter = 1d3
    Real*8,parameter :: racc = 1.d-12
    Integer*4 :: p
    logical :: root 

    rini = (3.d0*M/4.d0/Pi/mean_density(z)*(1.d0+z)**3.d0/delta)**(1.d0/3.d0) ! Units : Mpc/h

    step = 1.6d0

    r1 = 1.d0*rini ! Units : Mpc/h

    r2 = 1.1d0*rini ! Units : Mpc/h

    f1 = integral_r_delta_d_minus_total_matter(M,z,r1,delta)

    f2 = integral_r_delta_d_minus_total_matter(M,z,r2,delta)

    Do p = 1, iter

        If (f1*f2 .lt. 0.d0) then

            root = .true.

            exit

        else if (abs(f1) .lt. abs(f2)) then

            r1 = r1 + step*(r1-r2)

            f1 = integral_r_delta_d_minus_total_matter(M,z,r1,delta)

        else if (abs(f1) .gt. abs(f2)) then

            r2 = r2 + step*(r2-r1)

            f2 = integral_r_delta_d_minus_total_matter(M,z,r2,delta)

        end if

    End Do

    If (root) then

        Do p=1,iter

            If (f1 .lt. 0.d0) then 

                dr = r2 - r1

                rmid = r1 + dr/2.d0

                fmid = integral_r_delta_d_minus_total_matter(M,z,rmid,delta)    

                If (fmid .lt. 0.d0) then 

                    f1 = fmid

                    r1 = rmid

                else if ((fmid .eq. 0.d0) .or. (abs(dr) .lt. racc)) then

                    exit

                else 

                    f2 = fmid

                    r2 = rmid

                end if

            else if (f2 .lt. 0.d0) then

                dr = r1 - r2

                rmid = r2 + dr/2.d0

                fmid = integral_r_delta_d_minus_total_matter(M,z,rmid,delta)    

                If (fmid .lt. 0.d0) then 

                    f2 = fmid

                    r2 = rmid

                else if ((fmid .eq. 0.d0) .or. (abs(dr) .lt. racc)) then

                    exit

                else

                    f1 = fmid

                    r1 = rmid

                end if

            End if

        End Do

        r_delta_d_from_M_virial = rmid

    else 

        print *,'Maximum number of iterations achieved without finding root for function '

        stop 

    End If

end function r_delta_d_from_M_virial

function M_delta_d_from_M_virial(z,rdd,delta)    !    Units of M and M_delta_d_from_M_virial: solar mass
    use fiducial
    Implicit none
    Real*8 :: rdd,z,delta,M_delta_d_from_M_virial

    M_delta_d_from_M_virial = 4.d0*Pi*rdd**3*mean_density(z)/(1.d0+z)**3.d0*delta/3.d0

end function M_delta_d_from_M_virial

subroutine compute_M_delta_c_from_M_and_z(delta)
    use fiducial
    use arrays
    Implicit none
    Integer*4 :: indexz,indexM
    Real*8,dimension(number_of_M,number_of_z) :: rdc,rdd
    Real*8 :: delta

    rdc(:,:) = 0.d0

    rdd(:,:) = 0.d0

    open(15,file='./precomputed_quantities/M_delta_c_d.dat')

    write(15,*) '# red-shift z    virial mass M [solar mass]    r_delta_c [Mpc]    M_delta_c [solar mass]'

    write(15,*) '# r_delta_d [Mpc]    M_delta_d [solar mass]'

    Do indexz = 1,number_of_z

        Do indexM = 1,number_of_M

            rdc(indexM,indexz) = r_delta_c_from_M_virial(M(indexM),z(indexz),delta)

            rdd(indexM,indexz) = r_delta_d_from_M_virial(M(indexM),z(indexz),delta)

        End Do

    End Do

    Do indexz = 1,number_of_z

        Do indexM = 1,number_of_M

            write(15,'(6es18.10)') z(indexz), M(indexM), rdc(indexM,indexz), &
            M_delta_c_from_M_virial(z(indexz),rdc(indexM,indexz),delta), &
            rdd(indexM,indexz),M_delta_d_from_M_virial(z(indexz),rdd(indexM,indexz),delta)

        End Do

    End Do

    close(15)

end subroutine compute_M_delta_c_from_M_and_z

subroutine read_M200dc_r200dc()
    use arrays
    use fiducial
    Implicit none
    Integer*4 :: indexz,indexM
    Real*8 :: MM,zz,tM200c,tM200d,tr200c,tr200d

    open(15,file='./precomputed_quantities/M_delta_c_d.dat')

    read(15,*)

    read(15,*)

    Do indexz=1,number_of_z

        Do indexM=1,number_of_M

            read(15,'(6es18.10)') zz,MM,tr200c,tM200c,tr200d,tM200d

            If ( ( abs( M(indexM) - MM) .lt. 1.d5 ) .and. ( abs( z(indexz) - zz ) .lt. 1.d-8 ) ) then

                r200c(indexM,indexz) = tr200c

                M200c(indexM,indexz) = tM200c

                r200d(indexM,indexz) = tr200d

                M200d(indexM,indexz) = tM200d

            else 

                print *, 'Problem reading the file '

                print *, zz,z(indexz),MM,M(indexM)

                stop

            End If

        End Do

    End Do

    close(15)

end subroutine read_M200dc_r200dc

function ICM_electron_pressure(r,indexM,indexz) ! Units of Mv : solar mass. Units of r : Mpc. Intracluster medium electron pressure from "On 
    use fiducial              ! the cluster physics of Sunyaev-Zel'dovich and x-ray surveys. II. Deconstructing the thermal SZ
    use arrays                ! power spectrum" by Battaglia et al. All the masses and distances in this work are given relative to h.
    Implicit none             ! Units of ICM_electron_pressure : solar mass/(Mpc*s**2)        
    Real*8 :: x,A_p,alpha_pm,alpha_pz,P_0,A_xc,alpha_xcm
    Real*8 :: alpha_xcz,x_c,A_beta,r,ICM_electron_pressure
    Real*8 :: alpha_betam,alpha_betaz,beta,gamma,alpha,P_th,P_e,P200,P_fit
    Integer*4 :: indexM,indexz

    P200 = G*M200c(indexM,indexz)*DeltaSO*critical_density(z(indexz))&
    *Omega_b_h2/h**2/Omega_m0/2.d0/r200c(indexM,indexz)*M_sun/(Mpc)**3 ! Units : solar mass/(Mpc*s**2)

    x = r/r200c(indexM,indexz)        ! Dimensionless 
                                      ! Parameters from columns "AGN Feedback \Delta = 200" in Table 1 in published version of 1109.3711
    A_p = 18.1d0                      ! For P_0

    alpha_pm = 0.154d0

    alpha_pz = -0.758d0

    P_0 = A_p*(M200c(indexM,indexz)*1.d-14)**alpha_pm*(1+z(indexz))**alpha_pz

    A_xc = 0.497d0                    ! For x_c

    alpha_xcm = -0.00865d0

    alpha_xcz = 0.731d0

    x_c = A_xc*(M200c(indexM,indexz)*1.d-14)**alpha_xcm*(1+z(indexz))**alpha_xcz

    A_beta = 4.35d0                   ! For beta

    alpha_betam = 0.0393d0

    alpha_betaz = 0.415d0

    beta = A_beta*(M200c(indexM,indexz)*1.d-14)**alpha_betam*(1+z(indexz))**alpha_betaz

    gamma = -0.3d0                    ! Parameters \gamma and \alpha in equation (10) of published version of 1109.3711

    alpha = 1.0d0

    P_fit = P_0*(x/x_c)**gamma*(1.d0+(x/x_c)**alpha)**(-beta) ! Equation (10) in published version of 1109.3711

    P_th = P200*P_fit

    P_e = P_th/1.932d0

    ICM_electron_pressure = P_e

end function ICM_electron_pressure    

subroutine write_ICM_electron_pressure_at_z(indexM,indexz)
    use arrays
    use fiducial
    Implicit none
    Integer*4 :: indexM,indexz,index
    Integer*4,parameter :: indexr = 1d2
    Real*8 :: rmin,rmax
    Real*8,dimension(indexr) :: r

    rmin = 1.d-2*r200c(indexM,indexz)

    rmax = 3.d0*r200c(indexM,indexz)

    Do index = 1, indexr 

        r(index) = 10**(log10(rmin) + real(index-1)*(log10(rmax) - log10(rmin))/real(indexr-1))

    End Do

    open(15,file='./output/Pth_at_z.dat')

    write(15,*) '# r/R_200       P_th at redshift and mass',z(indexz),M200c(indexM,indexz)

    Do index =1,indexr

        write(15,'(2es18.10)') r(index)/r200c(indexM,indexz), ICM_electron_pressure(r(index),indexM,indexz)

    End Do

    close(15)

end subroutine write_ICM_electron_pressure_at_z


function form_factor(indexM,indexz,indexl)    ! Form factor. Equation (2.9) in 1312.4525. Units of M: solar mass. 
    use fiducial                              ! Units : dimensionless
    use arrays
    Implicit none
    Real*8 :: l_s,x_y_min,x_y_max,form_factor,y,prefactor
    Integer*4,parameter :: number_of_x = 1d5
    Integer*4,parameter :: intervals = number_of_x - 1 ! It gives error less than 1% 
    Integer*4 :: indexx,indexM,indexz,indexl
    Real*8,dimension(number_of_x) :: x,f

    x_y_min = 1.d-10*virial_radius(z(indexz),M(indexM))/r200c(indexM,indexz) ! dimensionless

    l_s = angular_diameter_distance(z(indexz))/r200c(indexM,indexz)         ! dimensionless

    x_y_max = 1.05d0*virial_radius(z(indexz),M(indexM))/r200c(indexM,indexz) ! dimensionless

    prefactor = sigma_e*4.d0*Pi*r200c(indexM,indexz)*M_sun/m_e/c**2/l_s**2  ! Units : Mpc*s**2/solar mass

    Do indexx=1,number_of_x

        x(indexx) = 10**(log10(x_y_min) + dble((indexx-1))*(log10(x_y_max) &
        - log10(x_y_min))/dble(number_of_x-1)   )

    End Do

    Do indexx=1,number_of_x 

        f(indexx) = x(indexx)**2*sin((dble(ml(indexl))+1.d0/2.d0)*x(indexx)/l_s)&
        *ICM_electron_pressure(x(indexx)*r200c(indexM,indexz),indexM,indexz)/(dble(ml(indexl))+1.d0/2.d0)/x(indexx)*l_s

    End Do

    y = 0.d0 

    Do indexx=1,intervals

        y = (x(indexx+1)-x(indexx))/2.d0*( f(indexx) + f(indexx+1))+ y

    End Do

    form_factor = prefactor*y

end function form_factor

subroutine compute_form_factor()
    use arrays
    use fiducial
    use omp_lib
    Implicit none
    Integer*4 :: indexl,indexM,indexz

    !$omp Parallel Do Shared(ylMz)

    Do indexl=1,number_of_l

        Do indexM=1,number_of_M

            Do indexz=1,number_of_z

                ylMz(indexl,indexM,indexz) = form_factor(indexM,indexz,indexl)

            End Do

        End Do

    End Do

    !$omp End Parallel Do

    open(15,file='./precomputed_quantities/form_factor/form_factor.dat')

    write(15,*) '# l  virial mass[solar mass]  red-shift  y'

    Do indexl=1,number_of_l

        Do indexM=1,number_of_M

            Do indexz=1,number_of_z

                write(15,'(i5,3es18.10)') ml(indexl),M(indexM),z(indexz),ylMz(indexl,indexM,indexz)

            End Do

        End Do

    End Do

    close(15)

end subroutine compute_form_factor


subroutine read_ylMz()
    use arrays
    use fiducial
    Implicit none
    Integer*4 :: indexz,indexM,indexl
    Real*8 :: ll,MM,zz,tylMz

    open(15,file='./precomputed_quantities/form_factor/form_factor.dat')

    read(15,*)

    Do indexl=1,number_of_l

        Do indexM=1,number_of_M

            Do indexz=1,number_of_z

                read(15,'(i5,3es18.10)') ll,MM,zz,tylMz

                If ( ( ml(indexl) .eq. ll ) .and. ( ( abs( M(indexM) - MM) .lt. 1.d5 ) &
                .and. ( abs( z(indexz) - zz ) .lt. 1.d-8 ) ) ) then

                    ylMz(indexl,indexM,indexz) = tylMz

                else 

                    print *, 'Problem reading the file '

                    stop

                End If

            End Do

        End Do

    End Do

    close(15)

end subroutine read_ylMz

function lensing_potential(indexM,indexz,indexl,halo_definition)    ! Lensing potential. Equation (2.10) in 1312.4525. Units of M : solar mass
    use fiducial    ! Units : dimensionless
    use arrays
    use omp_lib
    Implicit none
    Real*8 :: lensing_potential,r_s,l_s,phi_l,prefactor,xmin,xmax
    Integer*4 :: n,points,i,indexM,indexz,indexl
    Integer*4,parameter :: number_of_x = 1d5
    Integer*4,parameter :: intervals = number_of_x - 1  ! It gives less than 1% error
    Real*8,dimension(number_of_x) :: x,f    
    Character(len=*) :: halo_definition 

    If (halo_definition .eq. 'critical_density') then

        r_s = r_delta_c(DeltaSO,z(indexz),M200c(indexM,indexz))/concentration_mass_critical(M200c(indexM,indexz),z(indexz))    ! Mpc

        xmin = 1.d-3*r_delta_c(DeltaSO,z(indexz),M200c(indexM,indexz))/r_s    ! dimensionless

        xmax = 1.5d0*r_delta_c(DeltaSO,z(indexz),M200c(indexM,indexz))/r_s    ! dimensionless


    Else If (halo_definition .eq. 'mean_background') then 

        r_s = r_delta_d(DeltaSO,z(indexz),M200d(indexM,indexz))/concentration_mass_mean(M200d(indexM,indexz),z(indexz))    ! Mpc. 
                                                                                                     
        xmin = 1.d-5*r_delta_d(DeltaSO,z(indexz),M200d(indexM,indexz))/r_s    ! dimensionless                  

!        xmax = 1.05150d0*r_delta_d(DeltaSO,z(indexz),M200d(indexM,indexz))/r_s    ! dimensionless

        xmax = 1.65d0*virial_radius(z(indexz),M(indexM))/r_s    ! dimensionless

    Else If (halo_definition .eq. 'virial') then

        r_s = virial_radius(z(indexz),M(indexM))/concentration_mass_virial(M(indexM),z(indexz))    ! Mpc

        xmin = 1.d-10*virial_radius(z(indexz),M(indexM))/r_s    ! dimensionless

        xmax = 1.05d0*virial_radius(z(indexz),M(indexM))/r_s    ! dimensionless

    Else

         print *,'Not halo definition implemented with name ', halo_definition

         stop

    End If

    l_s = angular_diameter_distance(z(indexz))/r_s                             ! dimensionless

    prefactor = 2.d0*4.d0*Pi*r_s/dble(ml(indexl))/(dble(ml(indexl))+1.d0)/l_s**2 ! Mpc

    Do n=1,number_of_x

        x(n) = 10**(log10(xmin) + dble((n-1))*(log10(xmax) - log10(xmin))/dble(number_of_x-1) ) ! dimensionless

    End Do

    !$omp Parallel Do Shared(f)    

    Do n=1,number_of_x

        If (halo_definition .eq. 'critical_density') then

            f(n) = x(n)**2*sin((dble(ml(indexl))+1.d0/2.d0)*x(n)/l_s)*&
            NFW_density_profile(x(n)*r_s,M200c(indexM,indexz),z(indexz),'critical_density')/&
            (dble(ml(indexl))+1.d0/2.d0)/x(n)*l_s/Scrit(indexz)    ! 1/Mpc
        
        Else If (halo_definition .eq. 'mean_background') then 

            f(n) = x(n)**2*sin((dble(ml(indexl))+1.d0/2.d0)*x(n)/l_s)*&
            NFW_density_profile(x(n)*r_s,M200d(indexM,indexz),z(indexz),'mean_background')/&
            (dble(ml(indexl))+1.d0/2.d0)/x(n)*l_s/Scrit(indexz)    ! 1/Mpc
        
        Else If (halo_definition .eq. 'virial') then

            f(n) = x(n)**2*sin((dble(ml(indexl))+1.d0/2.d0)*x(n)/l_s)*&
            NFW_density_profile(x(n)*r_s,M(indexM),z(indexz),'virial')/&
            (dble(ml(indexl))+1.d0/2.d0)/x(n)*l_s/Scrit(indexz)    ! 1/Mpc
        
        End If

    End Do

    !$omp End Parallel Do

!    open(15,file='./output/array_lensing.dat')
!    Do n=1,number_of_x

!        write(15,*) x(n), f(n)

!    End Do

!    close(15)
    phi_l = 0.d0

    Do i=1,intervals

        phi_l = (x(i+1)-x(i))/2.d0*(f(i) + f(i+1)) + phi_l

    End Do

    lensing_potential = prefactor*phi_l

end function lensing_potential

subroutine compute_lensing_potential(halo_definition)
    use arrays
    use fiducial
    Implicit none
    Integer*4 :: indexl,indexM,indexz
    Character(len=*) :: halo_definition 

    open(15,file='./precomputed_quantities/lensing_potential/lensing_potential.dat')

    If (halo_definition .eq. 'critical_density') then

        write(15,*) '#  l  virial mass [solar mass]  red-shift  phi'

        Do indexl=1,number_of_l

            Do indexM=1,number_of_M

                Do indexz=1,number_of_z

                    philMz(indexl,indexM,indexz) = lensing_potential(indexM,indexz,indexl,'critical_density')

                    write(15,'(i5,3es18.10)') ml(indexl),M(indexM),z(indexz),philMz(indexl,indexM,indexz)

                End Do

            End Do

        End Do

        close(15)

        
    Else If (halo_definition .eq. 'mean_background') then 

        write(15,*) '#  l  virial mass [solar mass]  red-shift  phi'

        Do indexl=1,number_of_l

            Do indexM=1,number_of_M

                Do indexz=1,number_of_z

                    philMz(indexl,indexM,indexz) = lensing_potential(indexM,indexz,indexl,'mean_background')

                    write(15,'(i5,3es18.10)') ml(indexl),M(indexM),z(indexz),philMz(indexl,indexM,indexz)

                End Do

            End Do

        End Do

        close(15)

        
    Else If (halo_definition .eq. 'virial') then

        write(15,*) '#  l  virial mass[solar mass]  red-shift  phi'

        Do indexl=1,number_of_l

            Do indexM=1,number_of_M

                Do indexz=1,number_of_z

                    philMz(indexl,indexM,indexz) = lensing_potential(indexM,indexz,indexl,'virial')

                    write(15,'(i5,3es18.10)') ml(indexl),M(indexM),z(indexz),philMz(indexl,indexM,indexz)

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
    Integer*4 :: indexz,indexM,indexl,ll
    Real*8 :: MM,zz,tphilMz

    open(15,file='./precomputed_quantities/lensing_potential/lensing_potential.dat')

    read(15,*)

    Do indexl=1,number_of_l

        Do indexM=1,number_of_M

            Do indexz=1,number_of_z

                read(15,'(i5,3es18.10)') ll,MM,zz,tphilMz

                If ( ( ml(indexl) .eq. ll ) .and. ( ( abs( M(indexM) - MM) .lt. 1.d5 ) &
                .and. ( abs( z(indexz) - zz ) .lt. 1.d-8 ) ) ) then

                    philMz(indexl,indexM,indexz) = tphilMz

                else 

                    print *, 'Problem reading the file '

                    stop

                End If

            End Do

        End Do

    End Do

    close(15)

end subroutine read_philMz

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

subroutine compute_normalization() ! Compute normalisation constant in the linear matter power spectrum to agree 
    use fiducial                   ! with \sigma_8 in fiducial model 
    use arrays
    Implicit none
    Real*8 :: sum
    Integer*4 :: indexk,indexz
    Integer*4,parameter :: intervals = number_of_k - 1 
    Real*8,dimension(number_of_k) :: f 
    Real*8,parameter :: R = 8.d0/h ! Units : Mpc

    !$omp Parallel Do Shared(f)

    Do indexk=1,number_of_k

        f(indexk) = Delta_squared_non_normalised(k(indexk))*(3.d0*(sin(k(indexk)*R)-&
        k(indexk)*R*cos(k(indexk)*R))/k(indexk)**3/R**3)**2/k(indexk)

    End Do

    !$omp End Parallel Do

    sum = 0.d0

    Do indexk=1,intervals

        sum = (k(indexk+1) - k(indexk))/2.d0*( f(indexk) + f(indexk+1) ) + sum 

    End Do

    Normalization = sigma8**2/sum    !    

end subroutine compute_normalization

function matter_power_spectrum(k,z)    !    P(k,z) in Equation (3.15) of A. Lewis and A. Challinor, Physics Reports 429 (2006) 1-65
    use fiducial    !    P(k,z) is normalized with fiducial \sigma_8. It corresponds to P(k) in equation (A1) of "Baryonic features in the 
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


function sigma_squared(Mass,indexz) ! Equation (6) of 1303.4726. Units of Mass: solar mass. Units : dimensionless
    use fiducial
    use arrays
    Implicit none
    Real*8 :: sigma_squared,R,sum,prefactor,Mass
    Integer*4 :: indexk,indexz
    Integer*4,parameter :: intervals = number_of_k - 1 
    Real*8,dimension(number_of_k) :: f 

    prefactor = 1.d0/2.d0/Pi**2

    R = (3.d0*Mass/4.d0/Pi/mean_density(z(indexz))*(1.d0 + z(indexz))**3.d0)**(1.d0/3.d0) ! Units : Mpc. (1+z)**3 factor to use comoving coordinates

    !$omp Parallel Do Shared(f)

    Do indexk=1,number_of_k

        f(indexk) = k(indexk)**2*window_function(k(indexk),R)**2*matter_power_spectrum(k(indexk),z(indexz))

    End Do

    !$omp End Parallel Do

    sum = 0.d0

    Do indexk=1,intervals

        sum = (k(indexk+1) - k(indexk))/2.d0*( f(indexk) + f(indexk+1) ) + sum 

    End Do

    sigma_squared = prefactor*sum/h**3

end function sigma_squared

function dsigma_squared_dR(Mass,indexz) ! Derivative w.r.t R of Eq. (6) in 1303.4726. Units of Mass: solar mass. Units : 1/Mpc
    use fiducial
    use arrays 
    Implicit none
    Real*8 :: dsigma_squared_dR,R,sum,prefactor,Mass
    Integer*4 :: indexk,indexz
    Integer*4,parameter :: intervals = number_of_k - 1 
    Real*8,dimension(number_of_k) :: f 

    prefactor = 1.d0/Pi**2

    R = (3.d0*Mass/4.d0/Pi/mean_density(z(indexz))*(1.d0 + z(indexz))**3.d0)**(1.d0/3.d0) ! Units : Mpc. (1+z)**3 factor to use comoving coordinates

    !$omp Parallel Do Shared(f)

    Do indexk=1,number_of_k

        f(indexk) = k(indexk)**3*window_function(k(indexk),R)*d_window_function(k(indexk),R)&
        *matter_power_spectrum(k(indexk),z(indexz))

    End Do

    !$omp End Parallel Do

    sum = 0.d0

    Do indexk=1,intervals

        sum = (k(indexk+1)-k(indexk))/2.d0*( f(indexk) + f(indexk+1) ) + sum 

    End Do

    dsigma_squared_dR = prefactor*sum/h**3

end function dsigma_squared_dR

function halo_mass_function(indexM,indexz) ! Equation (8) in "The large-scale bias of dark matter halos: numerical calibration and model tests".
    use fiducial                 ! Units of M [M_{200d}] : solar mass. Units : 1/(solar mass*Mpc**3). Halo definition : mean density
    use arrays                    
    Implicit none                 
    Real*8 :: halo_mass_function,R,sigma,nu,beta,phi,eta,gamma
    Real*8 :: f_nu,g_sigma        
    Integer*4 :: indexM,indexz
    Real*8,parameter :: delta_c = 1.686d0    ! page 880 in "The large-scale ..."
    Real*8,parameter :: beta0 = 0.589d0
    Real*8,parameter :: phi0 = -0.729d0
    Real*8,parameter :: eta0 = -0.243d0
    Real*8,parameter :: gamma0 = 0.864d0

    R = (3.d0*M200d(indexM,indexz)/4.d0/Pi/mean_density(z(indexz))*(1.d0 + z(indexz))**3.d0)**(1.d0/3.d0)    ! Units : Mpc. (1+z)**3 to use comoving coordinates

    sigma = sqrt(sigma_squared(M200d(indexM,indexz),indexz))    ! Units : dimensionless

    nu = delta_c/sigma           ! page 880 in "The large-scale ... tests"

    beta = beta0*(1+z(indexz))**(0.20d0) 

    phi = phi0*(1+z(indexz))**(-0.08d0)

    eta = eta0*(1+z(indexz))**(0.27d0)

    gamma = gamma0*(1+z(indexz))**(-0.01d0)

    f_nu = alpha_halo_mass_function(indexz)*(1.d0 + (beta*nu)**(-2.d0*phi))*nu**(2.d0*eta)*dexp(-gamma*nu**2/2.d0)

    g_sigma = nu*f_nu ! dimensionless           ! Equation (C2) in "Toward a ... universality"

    halo_mass_function = -mean_density(z(indexz))/(1.d0 + z(indexz))**3.d0/2.d0/M200d(indexM,indexz)**2&    ! (1+z)**3 to use comoving coordinates
    *R/3.d0/sigma_squared(M200d(indexM,indexz),indexz)*dsigma_squared_dR(M200d(indexM,indexz),indexz)*g_sigma

end function halo_mass_function

function nonnormalised_halo_mass_function(indexM,indexz)    ! Equation (8) in "The large-scale bias of dark matter halos: 
    use fiducial    ! numerical calibration and model tests" without \alpha constant. This function corresponds to Equation (C2) in
    use arrays    ! "Toward a halo mass function for precision cosmology: the limits of universality" by Jeremy Tinker et al. Parameters taken 
    Implicit none    ! from first line of table 4. Units : 1/(solar mass*Mpc**3)
    Real*8 :: nonnormalised_halo_mass_function,R,sigma,nu,beta,phi,eta,gamma
    Real*8 :: f_nu,g_sigma
    Integer*4 :: indexM,indexz
    Real*8,parameter :: delta_c = 1.686d0    ! page 880 in "The large-scale ..."
    Real*8,parameter :: beta0 = 0.589d0
    Real*8,parameter :: phi0 = -0.729d0
    Real*8,parameter :: eta0 = -0.243d0
    Real*8,parameter :: gamma0 = 0.864d0

    R = (3.d0*M200d(indexM,indexz)/4.d0/Pi/mean_density(z(indexz))*(1.d0 + z(indexz))**3.d0)**(1.d0/3.d0)    ! Units : Mpc. (1+z)**3 to use comoving coordinates

    sigma = sqrt(sigma_squared(M200d(indexM,indexz),indexz))    ! Units : dimensionless

    nu = delta_c/sigma           ! page 880 in "The large-scale ... tests"

    beta = beta0*(1+z(indexz))**(0.20d0) 

    phi = phi0*(1+z(indexz))**(-0.08d0)

    eta = eta0*(1+z(indexz))**(0.27d0)

    gamma = gamma0*(1+z(indexz))**(-0.01d0)

    f_nu = (1.d0 + (beta*nu)**(-2.d0*phi))*nu**(2.d0*eta)*dexp(-gamma*nu**2/2.d0)

    g_sigma = nu*f_nu            ! Equation (C2) in "Toward a ... universality"

    nonnormalised_halo_mass_function = -mean_density(z(indexz))/(1.d0 + z(indexz))**3.d0/2.d0/M200d(indexM,indexz)**2&
    *R/3.d0/sigma_squared(M200d(indexM,indexz),indexz)*dsigma_squared_dR(M200d(indexM,indexz),indexz)*g_sigma

end function nonnormalised_halo_mass_function

subroutine compute_alpha_halo_mass_function()    ! It computes \alpha constant in halo mass function (Equation (8)) of "The large-scale 
    use fiducial    ! bias of dark matter halos: numerical calibration and model tests" for the red-shift array
    use omp_lib
    use arrays
    Implicit none
    Real*8 :: sum
    Integer*4 :: indexM,indexz
    Integer*4,parameter :: intervals = number_of_M - 1
    Real*8,dimension(number_of_M) :: f

    Do indexz=1,number_of_z

        !$omp Parallel Do Shared(f)

        Do indexM=1,number_of_M

            f(indexM) = nonnormalised_halo_mass_function(indexM,indexz)*bMz(indexM,indexz)*&
            M(indexM)/mean_density(z(indexz))*(1.d0 + z(indexz))**3.d0

        End Do

        !$omp End Parallel Do

        sum = 0.d0

        Do indexM=1,intervals

            sum = (M(indexM+1)-M(indexM))/2.d0*( f(indexM) + f(indexM+1) ) + sum

        End Do

        alpha_halo_mass_function(indexz) = 1.d0/sum    ! dimensionless

    End Do

end subroutine compute_alpha_halo_mass_function

subroutine compute_dndM()    ! It fills halo mass function array in and writes it out to a text file
    use arrays
    use fiducial
    Implicit none
    Integer*4 :: indexz,indexM

    open(15,file='./precomputed_quantities/dndM/dndM.dat')

    write(15,*) '#  virial mass [solar mass]    red-shift        dndM '

    Do indexM=1,number_of_M

        Do indexz=1,number_of_z

            dndM(indexM,indexz) = halo_mass_function(indexM,indexz)

            write(15,'(3es18.10)') M(indexM), z(indexz), dndM(indexM,indexz)

        End Do

    End Do

    close(15)

end subroutine compute_dndM


subroutine read_dndM()
    use arrays
    use fiducial
    Implicit none
    Integer*4 :: indexz,indexM
    Real*8 :: MM,zz,tdndM

    open(15,file='./precomputed_quantities/dndM/dndM.dat')

    read(15,*)

    Do indexM=1,number_of_M

        Do indexz=1,number_of_z

            read(15,'(3es18.10)') MM,zz,tdndM

            If ( ( abs( M(indexM) - MM) .lt. 1.d5 ) .and. ( abs( z(indexz) - zz ) .lt. 1.d-8 ) ) then

                dndM(indexM,indexz) = tdndM

            else 

                print *, 'Problem reading the file '

                stop

            End If

        End Do

    End Do

    close(15)

end subroutine read_dndM

subroutine write_dndM_at_z(indexz)
    use arrays
    use fiducial
    Implicit none
    Integer*4 :: indexM,indexz

    open(15,file='./output/dndM_at_z.dat')

    write(15,*) '# Virial Mass [solar mass/h]  dndM at redshift'

    Do indexM=1,number_of_M

        write(15,'(2es18.10)') M(indexM), dndM(indexM,indexz)

    End Do

    close(15)

end subroutine write_dndM_at_z

function pre_Clphiphi(indexz,indexl) ! Required function to compute one halo term of lensing potential angular power spectrum. Dimensionless
    use fiducial
    use arrays
    Implicit none
    Real*8 :: pre_Clphiphi,sum
    Integer*4 :: indexl,i,indexM,indexz
    Integer*4,parameter :: intervals = number_of_M - 1 
    Real*8,dimension(number_of_M) :: f

    Do indexM=1,number_of_M

        f(indexM) = dndM(indexM,indexz)*philMz(indexl,indexM,indexz)**2    ! Units : 1/solar mass/Mpc**3

    End Do

    sum = 0.d0

    Do i=1,intervals

        sum = (M(i+1)- M(i))/2.d0*( f(i) + f(i+1) ) + sum   ! Units : 1/Mpc**3

    End Do

    pre_Clphiphi = d2VdzdO(indexz)*sum     ! Dimensionless                     

end function pre_Clphiphi

function C_l_phiphi_one_halo(indexl) ! It computes lensing potential angular power spectrum for a single multipole. Dimensionless
    use fiducial
    use arrays
    Implicit none
    Real*8 :: C_l_phiphi_one_halo,sum
    Integer*4 :: indexl,i,indexz
    Integer*4,parameter :: intervals = number_of_z - 1
    Real*8,dimension(number_of_z):: f

    Do indexz=1,number_of_z

        f(indexz) = pre_Clphiphi(indexz,indexl)          ! Dimensionless

    End Do

    sum = 0.d0

    Do indexz=1,intervals

        sum = (z(indexz+1) -  z(indexz))/2.d0*( f(indexz) + f(indexz+1) ) + sum

    End Do

    C_l_phiphi_one_halo = sum

end function C_l_phiphi_one_halo

subroutine compute_Clphiphi1h()  ! It fills lensing potential angular power spectrum array. Dimensionless
    use arrays
    use fiducial
    Implicit none
    Integer*4 :: indexl

    !$omp Parallel Do Shared(Cl1h)

    Do indexl=1,number_of_l

        Clphiphi1h(indexl) = C_l_phiphi_one_halo(indexl)

    End Do

    !$omp End Parallel Do

end subroutine compute_Clphiphi1h

function pre_Cl(indexz,indexl) ! Dimensionless
    use fiducial
    use arrays
    Implicit none
    Real*8 :: pre_Cl,sum
    Integer*4 :: indexl,i,indexM,indexz
    Integer*4,parameter :: intervals = number_of_M - 1 
    Real*8,dimension(number_of_M) :: f

    Do indexM=1,number_of_M

        f(indexM) = dndM(indexM,indexz)*ylMz(indexl,indexM,indexz)*philMz(indexl,indexM,indexz) ! Units : 1/solar mass/Mpc**3

    End Do

    sum = 0.d0

    Do i=1,intervals

        sum = (M(i+1)- M(i))/2.d0*( f(i) + f(i+1) ) + sum 

    End Do

    pre_Cl = d2VdzdO(indexz)*sum       ! dimensionless                   

end function pre_Cl

function C_l_yphi_one_halo(indexl) ! Dimensionless
    use fiducial
    use arrays
    Implicit none
    Real*8 :: C_l_yphi_one_halo,sum
    Integer*4 :: indexl,i,indexz
    Integer*4,parameter :: intervals = number_of_z - 1
    Real*8,dimension(number_of_z):: f

    Do indexz=1,number_of_z

        f(indexz) = pre_Cl(indexz,indexl)          ! Dimensionless

    End Do

    sum = 0.d0

    Do indexz=1,intervals

        sum = (z(indexz+1) -  z(indexz))/2.d0*( f(indexz) + f(indexz+1) ) + sum

    End Do

    C_l_yphi_one_halo = sum

end function C_l_yphi_one_halo

subroutine compute_Cl1h()  ! Dimensionless
    use arrays
    use fiducial
    Implicit none
    Integer*4 :: indexl

    !$omp Parallel Do Shared(Cl1h)

    Do indexl=1,number_of_l

        Cl1h(indexl) = C_l_yphi_one_halo(indexl)

    End Do

    !$omp End Parallel Do

end subroutine compute_Cl1h

function linear_halo_bias(indexM,indexz)    ! From table 2 and equation (6) of "The large-scale bias of dark matter halos:  
    use arrays    ! numerical calibration and model tests". It computes the linear halo bias. (See paper for details about parameters). 
    use fiducial  ! Mass given must be M_{200d}. Dimensionless. 
    Implicit none    
    Real*8 :: linear_halo_bias,sigma,nu
    Real*8,parameter :: y = log10(DeltaSO)
    Real*8,parameter :: A = 1.d0 + 2.4d-1*y*dexp(-(4.d0/y)**4)
    Real*8,parameter :: aa = 4.4d-1*y - 8.8d-1
    Real*8,parameter :: B = 1.83d-1
    Real*8,parameter :: bb = 1.5d0 
    Real*8,parameter :: CC = 1.9d-2 + 1.07d-1*y + 1.9d-1*dexp(-(4.d0/y)**4)
    Real*8,parameter :: ccc = 2.4d0
    Real*8,parameter :: delta_c = 1.686d0
    Integer*4 :: indexM,indexz
        
!    sigma = sqrt(sigma_squared(M(indexM),indexz))    ! M must be the SO mass M_{200d} 
    sigma = sqrt(sigma_squared(M200d(indexM,indexz),indexz))    ! M must be the SO mass M_{200d} 

    nu = delta_c/sigma

    linear_halo_bias = 1.d0 - A*nu**aa/(nu**aa + delta_c**aa) + B*nu**b + CC*nu**ccc

end function linear_halo_bias

subroutine compute_bMz()
    use arrays
    use fiducial
    Implicit none
    Integer*4 :: indexz,indexM

    open(15,file='./precomputed_quantities/bMz/bMz.dat')

    write(15,*) '#  SO virial mass [solar mass]  red-shift  b '

    Do indexM=1,number_of_M

        Do indexz=1,number_of_z

            bMz(indexM,indexz) = linear_halo_bias(indexM,indexz)

            write(15,'(3es18.10)') M(indexM), z(indexz), bMz(indexM,indexz)

        End Do

    End Do

    close(15)

end subroutine compute_bMz

subroutine compute_mean_bias_matter() ! Dimensionless
    use arrays
    use fiducial
    Implicit none
    Integer*4 :: indexz,indexM

    open(15,file='./precomputed_quantities/bMz/mean_bias_at_z.dat')

    write(15,*) '# red-shift        mean bias all matter'

    Do indexz=1,number_of_z

        mbz(indexz) = pre_mbz(indexz)

        write(15,'(2es18.10)') z(indexz),mbz(indexz)

    End Do

    close(15)

end subroutine compute_mean_bias_matter

function pre_mbz(indexz)    !    It computes mean bias of all matter at a given red-shift. Dimensionless.  
    use fiducial
    use omp_lib
    use arrays
    Implicit none
    Real*8 :: pre_mbz,sum
    Integer*4 :: indexM,indexz
    Integer*4,parameter :: intervals = number_of_M - 1
    Real*8,dimension(number_of_M) :: f

    Do indexM=1,number_of_M

        f(indexM) = dndM(indexM,indexz)*bMz(indexM,indexz)*M(indexM)/mean_density(z(indexz))*(1.d0 + z(indexz))**3.d0

    End Do

    sum = 0.d0

    Do indexM=1,intervals

        sum = (M(indexM+1)-M(indexM))/2.d0*( f(indexM) + f(indexM+1) ) + sum

    End Do

    pre_mbz = sum

end function pre_mbz


subroutine read_bMz()
    use arrays
    use fiducial
    Implicit none
    Integer*4 :: indexz,indexM
    Real*8 :: MM,zz,tbMz

    open(15,file='./precomputed_quantities/bMz/bMz.dat')

    read(15,*)

    Do indexM=1,number_of_M

        Do indexz=1,number_of_z

            read(15,'(3es18.10)') MM,zz,tbMz

            If ( ( abs( M(indexM) - MM) .lt. 1.d5 ) .and. ( abs( z(indexz) - zz ) .lt. 1.d-8 ) ) then

                bMz(indexM,indexz) = tbMz

            else 

                print *, 'Problem reading the file '

                stop

            End If

        End Do

    End Do

    close(15)

end subroutine read_bMz

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

subroutine write_bnu_at_z(indexz)
    use arrays
    use fiducial
    Implicit none
    Integer*4 :: indexM,indexz
    Real*8 :: sigma,delta_c,nu

    open(15,file='./output/bnu_at_z.dat')

    write(15,*) '# Virial Mass [solar mass/h]  bnu at redshift'

    Do indexM=1,number_of_M

        sigma = sqrt(sigma_squared(M(indexM),indexz))

    !    sigma = sqrt(sigma_squared(M200d(indexM,indexz),indexz))

        delta_c = 1.686d0

        nu = delta_c/sigma

        write(15,'(2es18.10)') nu, bMz(indexM,indexz)

    End Do

    close(15)

end subroutine write_bnu_at_z

function C_l_psi_limber(indexl)    ! It computes the angular power spectrum of the lensing potential in the Limber approximation 
    use fiducial    ! for a single multipole as in Lewis and Challinor. Units : Dimensionless
    use arrays
    Implicit none
    Real*8 :: C_l_psi_limber,sum,prefactor
    Integer*4 :: indexl,indexz
    Integer*4,parameter :: intervals = number_of_z - 1 
    Real*8,dimension(number_of_z):: f,zlimber
    Real*8,parameter :: zminlimber = 1.d-5

    ! Red-shift array. Dimensionless.

    Do indexz = 1, number_of_z      

        zlimber(indexz) = 10**(log10(zminlimber) + real(indexz-1)*(log10(z_dec) - log10(zminlimber))/real(number_of_z-1))

    End Do

!    prefactor = 8.d0*Pi**2/dble(ml(indexl))**3

    prefactor = 4.d0/dble(ml(indexl))**2/(dble(ml(indexl))+1.d0)**2*9.d0/4.d0/&
    c**4*Hubble_parameter(0.d0)**4*Omega_m(0.d0)**2    !    Units : 

    Do indexz=1,number_of_z

!        f(indexz) = comoving_distance(zlimber(indexz))*c/Hubble_parameter(zlimber(indexz))*&
!        power_spectrum_weyl_potential(dble(ml(indexl))/comoving_distance(zlimber(indexz)),zlimber(indexz))*&
!        ( comoving_distance(z_dec) - comoving_distance(zlimber(indexz)) )**2/&
!        comoving_distance(z_dec)**2/comoving_distance(zlimber(indexz))**2/(1.d0 + zlimber(indexz))**4

        f(indexz) = (1.d0 + zlimber(indexz))**2*( comoving_distance(z_dec) - comoving_distance(zlimber(indexz)) )**2/&
        comoving_distance(z_dec)**2*c/Hubble_parameter(zlimber(indexz))*&
        matter_power_spectrum((dble(ml(indexl)) + 1.d0/2.d0)/comoving_distance(zlimber(indexz)),zlimber(indexz))

    End Do

    sum = 0.d0

    Do indexz=1,intervals

        sum = (zlimber(indexz+1)-zlimber(indexz))/2.d0*( f(indexz) + f(indexz+1) ) + sum

    End Do

    C_l_psi_limber = prefactor*sum/h**3    

end function C_l_psi_limber
 
subroutine compute_Clpsilimber()   ! It fills array for angular power spectrum of lensing potential in the Limber approximations. 
    use arrays
    use fiducial
    Implicit none
    Integer*4 :: indexl

    !$omp Parallel Do Shared(Clpsilimber)

    Do indexl=1,number_of_l

        Clpsilimber(indexl) = C_l_psi_limber(indexl)

    End Do

    !$omp End Parallel Do

end subroutine compute_Clpsilimber


function C_l_phiphi_two_halo(indexl)    ! It computes two halo term lensing potential angular power spectrum for a single multipole. 
    use fiducial    ! Equation (2.11) in 1312.4525. Dimensionless
    use arrays
    Implicit none
    Real*8 :: C_l_phiphi_two_halo,sum
    Integer*4 :: indexl,indexz
    Integer*4,parameter :: intervals = number_of_z - 1 
    Real*8,dimension(number_of_z):: f

    Do indexz=1,number_of_z

        f(indexz) = pre_Cl_1(indexz,indexl)**2*d2VdzdO(indexz)*&
        matter_power_spectrum((dble(ml(indexl))+1.d0/2.d0)/comoving_distance(z(indexz)),z(indexz))    !  Dimensionless

    End Do

    sum = 0.d0

    Do indexz=1,intervals

        sum = (z(indexz+1)-z(indexz))/2.d0*( f(indexz) + f(indexz+1) ) + sum

    End Do

    C_l_phiphi_two_halo = sum/h**3

end function C_l_phiphi_two_halo

subroutine compute_Clphiphi2h()   ! It fills two halo term of lensing potential angular power spectrum array. 
    use arrays
    use fiducial
    Implicit none
    Integer*4 :: indexl

    !$omp Parallel Do Shared(Cl2h)

    Do indexl=1,number_of_l

        Clphiphi2h(indexl) = C_l_phiphi_two_halo(indexl)

    End Do

    !$omp End Parallel Do

end subroutine compute_Clphiphi2h


function pre_Cl_1(indexz,indexl) ! Units 1/Mpc**3. 
    use fiducial
    use omp_lib
    use arrays
    Implicit none
    Real*8 :: pre_Cl_1,sum
    Integer*4 :: indexl,indexM,indexz
    Integer*4,parameter :: intervals = number_of_M - 1
    Real*8,dimension(number_of_M) :: f

    Do indexM=1,number_of_M

        f(indexM) = dndM(indexM,indexz)*bMz(indexM,indexz)*philMz(indexl,indexM,indexz) 

    End Do

    sum = 0.d0

    Do indexM=1,intervals

        sum = (M(indexM+1)-M(indexM))/2.d0*( f(indexM) + f(indexM+1) ) + sum

    End Do

    pre_Cl_1 = sum

end function pre_Cl_1

function pre_Cl_2(indexz,indexl) ! Units 1/Mpc**3
    use fiducial
    use arrays
    use omp_lib
    Implicit none
    Real*8 :: pre_Cl_2,sum
    Integer*4 :: indexl,indexM,indexz
    Integer*4,parameter :: intervals = number_of_M - 1
    Real*8,dimension(number_of_M) :: f

    Do indexM=1,number_of_M

        f(indexM) = dndM(indexM,indexz)*ylMz(indexl,indexM,indexz)*bMz(indexM,indexz)

    End Do

    sum = 0.d0

    Do indexM=1,intervals

        sum = (M(indexM+1)-M(indexM))/2.d0*( f(indexM) + f(indexM+1) ) + sum

    End Do

    pre_Cl_2 = sum

end function pre_Cl_2

function C_l_yphi_two_halo(indexl) ! Equation (2.11) in 1312.4525
    use fiducial                   ! Units : dimensionless
    use arrays
    Implicit none
    Real*8 :: C_l_yphi_two_halo,sum
    Integer*4 :: indexl,indexz
    Integer*4,parameter :: intervals = number_of_z - 1 
    Real*8,dimension(number_of_z):: f

    Do indexz=1,number_of_z

        f(indexz) = pre_Cl_1(indexz,indexl)*pre_Cl_2(indexz,indexl)*d2VdzdO(indexz)*&
        matter_power_spectrum((dble(ml(indexl))+1.d0/2.d0)/comoving_distance(z(indexz)),z(indexz)) ! Units : 1/h**3

    End Do

    sum = 0.d0

    Do indexz=1,intervals

        sum = (z(indexz+1)-z(indexz))/2.d0*( f(indexz) + f(indexz+1) ) + sum

    End Do

    C_l_yphi_two_halo = sum/h**3

end function C_l_yphi_two_halo

subroutine compute_Cl2h()
    use arrays
    use fiducial
    use omp_lib
    Implicit none
    Integer*4 :: indexl

    !$omp Parallel Do Shared(Cl2h)

    Do indexl=1,number_of_l

        Cl2h(indexl) = C_l_yphi_two_halo(indexl)

    End Do

    !$omp End Parallel Do

end subroutine compute_Cl2h

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

subroutine write_Cl()
    use arrays
    use fiducial
    Implicit none
    Integer*4 :: indexl

    open(15,file='./output/Cl_yphi.dat')

    write(15,*) '#   l   Clyphi1h    Clyphi2h    Clyphi    Clphiphi1h    Clphiphi2h    Clphiphi    Clpsilimber '

    Do indexl=1,number_of_l

        write(15,'(i5,7es18.10)') ml(indexl), Cl1h(indexl), Cl2h(indexl), Cl(indexl),Clphiphi1h(indexl),&
                                  Clphiphi2h(indexl),Clphiphi(indexl),Clpsilimber(indexl)

    End Do

    close(15)

end subroutine write_Cl


end module functions
