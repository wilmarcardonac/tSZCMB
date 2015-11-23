Program tSZ

    ! LOAD NEEDED MODULES
    use fiducial
    use arrays
    use functions 
    use omp_lib
    use fgsl

    ! DECLARATION AND INITIALIZATION OF VARIABLES
    Implicit none
    Integer*4 :: index1 ! COUNTERS
    Real*8 :: wtime       ! MEASURES TIME
    Character(len=15),parameter :: halo_definition = 'virial'

    ! OPEN FILE TO STORE EXECUTION INFORMATION 
    open(20,file=path_to_execution_information)

    ! ALLOCATING MEMORY FOR : RED-SHIFT, VIRIAL MASS, MULTIPOLES, WAVEVECTOR, ONE- AND TWO-HALO Y-tSZ CROSS CORRELATION AND THEIR SUM, 
    ! MEAN DENSITY MASS, CRITICAL DENSITY MASS, CRITICAL DENSITY RADIUS, DERIVATIVE OF MEAN DENSITY MASS w.r.t VIRIAL MASS, 
    ! MEAN DENSITY RADIUS, ONE- AND TWO-HALO LENSING POTENTIAL AUTO-CORRELATION AND THEIR SUM, NORMALIZATION CONSTANT FOR HALO MASS FUNCTION,
    ! DERIVATIVE OF CRITICAL MASS DENSITY w.r.t VIRIAL MASS, ANGULAR POWER SPECTRUM OF LENSING POTENTIAL IN THE LIMBER APPROXIMATION, 
    ! HALO MASS FUNCTION, FORM FACTOR, LENSING POTENTIAL, COMOVING VOLUME, LINEAR BIAS, MEAN BIAS OF ALL MATTER, CRITICAL SURFACE,
    ! SIGMA SQUARE FOR MEAN DENSITY MASS, DERIVATIVE OF SIGMA SQUARE FOR MEAN DENSITY MASS
    allocate (z(1:number_of_z), M(-1:number_of_M+2), ml(1:number_of_l),k(1:number_of_k),&
    Cl1h(1:number_of_l),Cl2h(1:number_of_l),Cl(1:number_of_l),M200d(-1:number_of_M+2,1:number_of_z),&
    M200c(-1:number_of_M+2,1:number_of_z),r200c(-1:number_of_M+2,1:number_of_z),dM200ddM(1:number_of_M,1:number_of_z),&
    r200d(-1:number_of_M+2,1:number_of_z),Clphiphi1h(1:number_of_l),Clphiphi2h(1:number_of_l),&
    Clphiphi(1:number_of_l),alpha_halo_mass_function(1:number_of_z),dM200cdM(1:number_of_M,1:number_of_z),&
    Clpsilimber(1:number_of_l),dndM(1:number_of_M,1:number_of_z),ylMz(1:number_of_l,1:number_of_M,1:number_of_z),&
    philMz(1:number_of_l,1:number_of_M,1:number_of_z),d2VdzdO(1:number_of_z),&
    bMz(1:number_of_M,1:number_of_z),mbz(1:number_of_z),comoving_distance_at_z(1:number_of_z),stat = status1)

    If (status1 .eq. 0) then
       
        write(20,*) 'MEMORY ALLOCATED SUCCESSFULLY'

    Else

        write(20,*) 'PROBLEM ALLOCATING MEMORY'

        stop

    End If

    ! FILLING ARRAYS 

    ! WAVEVECTOR ARRAY. UNITS : 1/Mpc
    Do index1 = 1, number_of_k    

        k(index1) = 10**(log10(kmin) + real(index1-1)*(log10(kmax) - log10(kmin))/real(number_of_k-1))

    End Do

    ! RED-SHIFT ARRAY.
    Do index1 = 1, number_of_z      

        z(index1) = 10**(log10(zmin) + real(index1-1)*(log10(zmax) - log10(zmin))/real(number_of_z-1))

    End Do

    ! VIRIAL MASS ARRAY. UNITS: Solar mass 
    Do index1 = -1, number_of_M+2    

        M(index1) = 10**(log10(Mmin) + real(index1-1)*(log10(Mmax) - log10(Mmin))/real(number_of_M-1))

    End Do

    ! MULTIPOLE ARRAY.
    Do index1 = 1, number_of_l     

        ml(index1) = int(10**(log10(dble(lmin)) + real(index1-1)*(log10(dble(lmax)) - &
        log10(dble(lmin)))/real(number_of_l-1)),4)

    End Do

    ! COMPUTATION STARTS

    ! COMPUTE COMOVING DISTANCE ARRAY
    call compute_comoving_distance()

    ! COMPUTE COMOVING DISTANCE AT DECOUPLING
    com_dist_at_z_dec = comoving_distance(z_dec)

    ! REDO COMPUTATIONS ONLY IF MASS CONVERSION REQUIRED
    If (do_mass_conversion) then

        !wtime = omp_get_wtime() ! SETTING STARTING TIME OF MASS CONVERSION

        !write(20,*) 'EXECUTING MASS CONVERSION'

        !call compute_M_delta_c_from_M_and_z(DeltaSO)

        !write(20,*) 'MASS CONVERSION ENDED'

        !write(20,*) 'MASS CONVERSION FOR RED-SHIFT ARRAY OF SIZE ', size(z),' AND VIRIAL MASS ARRAY OF SIZE ', size(M),&
        !'TOOK ', (omp_get_wtime()-wtime)/3.6d3, 'HOURS'

        ! READING MASS CONVERSION FILE AND COMPUTING DERIVATIVES OF MASS CONVERSION
        call read_M200dc_r200dc()

        write(20,*) 'MASS CONVERSION FILE READ AND DERIVATIVES OF MASS CONVERSION COMPUTED'

        call compute_normalization() ! IN MATTER POWER SPECTRUM TO FIT FIDUCIAL SIGMA_8

        write(20,*) 'NORMALIZATION OF MATTER POWER SPECTRUM TO MATCH SIGMA_8 (',sigma8,') WAS COMPUTED'

        write(20,*) 'COMPUTING LINEAR HALO BIAS'

        call compute_bMz() ! LINEAR HALO BIAS AS A FUNCTION OF MASS AND RED-SHIFT

        write(20,*) 'COMPUTING NORMALIZATION OF HALO MASS FUNCTION'

        call compute_alpha_halo_mass_function() ! NORMALIZE HALO MASS FUNCTION TO FULLFILL CONDITION THAT MEAN BIAS OF ALL MATTER
                                                ! AT A FIXED RED-SHIFT IS UNITY
        
        write(20,*) 'COMPUTING HALO MASS FUNCTION'

        call compute_dndM() ! COMPUTING HALO MASS FUNCTION AS A FUNCTION OF MASS AND RED-SHIFT 

        write(20,*) 'COMPUTING COMOVING VOLUME ELEMENT PER STERADIAN'

        call compute_d2VdzdO()   ! COMPUTES COMOVING VOLUME ELEMENT PER STERADIAN  

        wtime = omp_get_wtime() ! SETTING STARTING TIME OF LENSING POTENTIAL COMPUTATION

        write(20,*) 'COMPUTING LENSING POTENTIAL'

        call compute_lensing_potential(halo_definition)

        write(20,*) 'LENSING POTENTIAL COMPUTATION FOR RED-SHIFT ARRAY OF SIZE ', size(z),', VIRIAL MASS ARRAY OF SIZE ',&
        size(M),', AND MULTIPOLES ARRAY OF SIZE ', size(ml),' TOOK ', (omp_get_wtime()-wtime)/3.6d3, 'HOURS'

        write(20,*) 'COMPUTING ANGULAR POWER SPECTRUM OF LENSING POTENTIAL'

        call compute_Clphiphi1h() ! ONE HALO TERM

        call compute_Clphiphi2h() ! TWO HALO TERM

        call compute_Clpsilimber()! LIMBER APPROXIMATION 

        call compute_Cl()   ! TOTAL 

        call write_Cl()     ! OF ALL COMPUTATIONS 
        stop
        wtime = omp_get_wtime() ! SETTING STARTING TIME OF FORM FACTOR COMPUTATION

        write(20,*) 'COMPUTING FORM FACTOR'

        call compute_form_factor()

        write(20,*) 'FORM FACTOR COMPUTATION FOR RED-SHIFT ARRAY OF SIZE ', size(z),', VIRIAL MASS ARRAY OF SIZE ', size(M),&
        ', AND MULTIPOLES ARRAY OF SIZE ', size(ml),' TOOK ', (omp_get_wtime()-wtime)/3.6d3, 'HOURS'

        write(20,*) 'COMPUTING ANGULAR POWER SPECTRUM OF Y-tSZ CROSS-CORRELATION'

        call compute_Cl1h() ! ONE HALO TERM

        call compute_Cl2h() ! TWO HALO TERM
 
        call compute_Cl()   ! TOTAL 

        call write_Cl()     ! OF ALL COMPUTATIONS 

        write(20,*) 'COMPUTATION ENDED'

    Else

        deallocate(mbz,M200c,M200d,r200c,r200d,dM200ddM,dM200cdM)

        write(20,*) 'USING MASS DATA FILE PREVIOUSLY COMPUTED'

        ! READING MASS CONVERSION FILE AND COMPUTING DERIVATIVES OF MASS CONVERSION
!        call read_M200dc_r200dc()

        write(20,*) 'MASS CONVERSION FILE READ AND DERIVATIVES OF MASS CONVERSION COMPUTED'

        call compute_normalization() ! IN MATTER POWER SPECTRUM TO FIT FIDUCIAL SIGMA_8

        write(20,*) 'NORMALIZATION OF MATTER POWER SPECTRUM TO MATCH SIGMA_8 (',sigma8,') WAS COMPUTED'

!        call read_bMz() ! LINEAR HALO BIAS AS A FUNCTION OF MASS AND RED-SHIFT
        call compute_bMz() ! LINEAR HALO BIAS AS A FUNCTION OF MASS AND RED-SHIFT

!        call read_alpha_halo_mass_function() ! NORMALIZE HALO MASS FUNCTION TO FULLFILL CONDITION THAT MEAN BIAS OF ALL MATTER
                                                ! AT A FIXED RED-SHIFT IS UNITY
        call compute_alpha_halo_mass_function() ! NORMALIZE HALO MASS FUNCTION TO FULLFILL CONDITION THAT MEAN BIAS OF ALL MATTER
                                                ! AT A FIXED RED-SHIFT IS UNITY


!        call read_dndM() ! READING HALO MASS FUNCTION    
        call compute_dndM() ! COMPUTING HALO MASS FUNCTION AS A FUNCTION OF MASS AND RED-SHIFT 

        write(20,*) 'COMPUTING COMOVING VOLUME ELEMENT PER STERADIAN'

        call compute_d2VdzdO()   ! COMPUTES COMOVING VOLUME ELEMENT PER STERADIAN  

        write(20,*) 'COMPUTING CRITICAL SURFACE DENSITY'

        call compute_lensing_potential(halo_definition)
!        call read_philMz() ! READS LENSING POTENTIAL

        write(20,*) 'COMPUTING ANGULAR POWER SPECTRUM OF LENSING POTENTIAL'

        call compute_Clphiphi1h() ! ONE HALO TERM

        call compute_Clphiphi2h() ! TWO HALO TERM

        call compute_Clpsilimber()! LIMBER APPROXIMATION 

        call compute_Cl()   ! TOTAL 

        call write_Cl()     ! OF ALL COMPUTATIONS 

        stop

        call read_ylMz() ! READS FORM FACTOR

        write(20,*) 'COMPUTING ANGULAR POWER SPECTRUM OF Y-tSZ CROSS-CORRELATION'

        call compute_Cl1h() ! ONE HALO TERM

        call compute_Cl2h() ! TWO HALO TERM
 
        call compute_Cl()   ! TOTAL 

        call write_Cl()     ! OF ALL COMPUTATIONS 

        write(20,*) 'COMPUTATION ENDED'

    End If

    ! DEALLOCATING MEMORY
    deallocate (z,M,k,ml,Cl1h,Cl2h,Clphiphi1h,Clphiphi2h,Cl,Clphiphi,&
    d2VdzdO,dndM,ylMz,philMz,bMz,alpha_halo_mass_function,&
    Clpsilimber,comoving_distance_at_z)

    ! CLOSE EXECUTION INFORMATION FILE
    close(20)

End Program tSZ




