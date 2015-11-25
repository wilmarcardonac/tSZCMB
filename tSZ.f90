Program tSZ

    use fiducial   ! CONTAINS PARAMETERS 
    use arrays     ! CONTAINS ARRAYS 
    use functions  ! CONTAINS ALL THE FUNCTIONS
    use omp_lib
    use fgsl

    ! DECLARATION AND INITIALIZATION OF VARIABLES
    Implicit none
    Integer*4 :: index1,index2                                       ! COUNTER
    Real*8 :: wtime,dndM_M_z                            ! STORES TIME OF EXECUTION
    Character(len=15),parameter :: halo_definition = 'virial' ! HALO DEFINITION USED IN THE COMPUTATIONS
    Real*8,dimension(10) :: zdebug,Mdebug
    Real*8,dimension(10,10) :: M200d_M_z,bMz_M_z

    com_dist_at_z_dec = comoving_distance(z_dec) ! COMPUTE COMOVING DISTANCE AT DECOUPLING

    open(UNIT_EXE_FILE,file=path_to_execution_information)     ! OPEN FILE TO STORE EXECUTION INFORMATION 

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
    Clphiphi(1:number_of_l),alpha_halo_mass_function(1:number_of_z_functions),dM200cdM(1:number_of_M,1:number_of_z),&
    Clpsilimber(1:number_of_l),dndM(1:number_of_M_functions,1:number_of_z_functions),&
    ylMz(1:number_of_l,1:number_of_M_functions,1:number_of_z_functions),&
    philMz(1:number_of_l,1:number_of_M_functions,1:number_of_z_functions),d2VdzdO(1:number_of_z_com_dist),&
    bMz(1:number_of_M_functions,1:number_of_z_functions),mbz(1:number_of_z_functions),&
    comoving_distance_at_z(1:number_of_z_com_dist),&
    z_com_dist(1:number_of_z_com_dist),z_functions(1:number_of_z_functions),&
    M_functions(1:number_of_M_functions),stat = status1)

    If (status1 .eq. 0) then
       
        write(UNIT_EXE_FILE,*) 'MEMORY ALLOCATED SUCCESSFULLY'

    Else

        write(UNIT_EXE_FILE,*) 'PROBLEM ALLOCATING MEMORY'

        stop

    End If

    Do index1 = 1, number_of_k ! FILLS WAVEVECTOR ARRAY. UNITS : 1/Mpc

        k(index1) = 10**(log10(kmin) + real(index1-1)*(log10(kmax) - log10(kmin))/real(number_of_k-1))

    End Do


    Do index1 = 1, number_of_z ! FILLS RED-SHIFT ARRAY.     

        z(index1) = 10**(log10(zmin) + real(index1-1)*(log10(zmax) - log10(zmin))/real(number_of_z-1))

    End Do
    
    Do index1 = -1, number_of_M+2 ! FILLS VIRIAL MASS ARRAY. UNITS: Solar mass    

        M(index1) = 10**(log10(Mmin) + real(index1-1)*(log10(Mmax) - log10(Mmin))/real(number_of_M-1))

    End Do

    Do index1 = 1, number_of_z_functions ! FILLS RED-SHIFT ARRAY.     

        z_functions(index1) = 10**(log10(zmin) + real(index1-1)*(log10(zmax) - log10(zmin))/real(number_of_z_functions-1))

    End Do
    
    Do index1 = 1, number_of_M_functions ! FILLS VIRIAL MASS ARRAY. UNITS: Solar mass    

        M_functions(index1) = 10**(log10(Mmin) + real(index1-1)*(log10(Mmax) - log10(Mmin))/real(number_of_M_functions-1))

    End Do

    Do index1 = 1, number_of_l ! FILLS MULTIPOLE ARRAY.    

        ml(index1) = int(10**(log10(dble(lmin)) + real(index1-1)*(log10(dble(lmax)) - &
        log10(dble(lmin)))/real(number_of_l-1)),4)

    End Do

    ! COMPUTATION STARTS

    call compute_comoving_distance() ! COMPUTE COMOVING DISTANCE ARRAY

    If (do_mass_conversion) then ! COMPUTE MASSES

        wtime = omp_get_wtime() ! SETTING STARTING TIME OF MASS CONVERSION

        write(UNIT_EXE_FILE,*) 'EXECUTING MASS CONVERSION'

        call compute_M_delta_c_from_M_and_z(DeltaSO)

        write(UNIT_EXE_FILE,*) 'MASS CONVERSION ENDED'

        write(UNIT_EXE_FILE,*) 'MASS CONVERSION FOR RED-SHIFT ARRAY OF SIZE ', size(z),' AND VIRIAL MASS ARRAY OF SIZE ', size(M),&
        'TOOK ', (omp_get_wtime()-wtime)/3.6d3, 'HOURS'
        
        call read_M200dc_r200dc() ! READING MASS CONVERSION FILE AND COMPUTING DERIVATIVES OF MASS CONVERSION

        write(UNIT_EXE_FILE,*) 'MASS CONVERSION FILE READ AND DERIVATIVES OF MASS CONVERSION COMPUTED'

     Else

        write(UNIT_EXE_FILE,*) 'USING MASS DATA FILE PREVIOUSLY COMPUTED'

        call read_M200dc_r200dc() ! READING MASS CONVERSION FILE AND COMPUTING DERIVATIVES OF MASS CONVERSION

        write(UNIT_EXE_FILE,*) 'MASS CONVERSION FILE READ AND DERIVATIVES OF MASS CONVERSION COMPUTED'

     End If


     Do index1 = 1, 10 ! FILLS RED-SHIFT ARRAY.     

        zdebug(index1) = 10**(log10(zmin) + real(index1-1)*(log10(zmax) - log10(zmin))/real(10-1))

     End Do
    
     Do index1 = 1, 10 ! FILLS VIRIAL MASS ARRAY. UNITS: Solar mass    
        
        Mdebug(index1) = 10**(log10(Mmin) + real(index1-1)*(log10(Mmax) - log10(Mmin))/real(10-1))
        
     End Do

     print *, r_delta_d_from_M_virial(Mdebug(10),zdebug(10),DeltaSO)

     stop

     Do index1=1,10

        Do index2=1,10

           call Interpolate_2D(M200d_M_z(index1,index2),Mdebug(index1),zdebug(index2),M(1:number_of_M),&
                z(1:number_of_z),M200d(1:number_of_M,1:number_of_z))

           print *, Mdebug(index1), M_delta_d_from_M_virial(zdebug(index2),&
                r_delta_d_from_M_virial(Mdebug(index1),zdebug(index2),DeltaSO),DeltaSO),&
                zdebug(index2),(M200d_M_z(index1,index2)-M_delta_d_from_M_virial(zdebug(index2),&
                r_delta_d_from_M_virial(Mdebug(index1),zdebug(index2),DeltaSO),DeltaSO))/M_delta_d_from_M_virial(zdebug(index2),&
                r_delta_d_from_M_virial(Mdebug(index1),zdebug(index2),DeltaSO),DeltaSO)*100.d0,&
                r_delta_d_from_M_virial(Mdebug(index1),zdebug(index2),DeltaSO)

        End Do

     End Do
     
     stop


     call compute_normalization() ! IN MATTER POWER SPECTRUM TO FIT FIDUCIAL SIGMA_8

     write(UNIT_EXE_FILE,*) 'NORMALIZATION OF MATTER POWER SPECTRUM TO MATCH SIGMA_8 (',sigma8,') WAS COMPUTED'

     If (compute_linear_halo_bias) then

        write(UNIT_EXE_FILE,*) 'COMPUTING LINEAR HALO BIAS'

        call compute_bMz() ! LINEAR HALO BIAS AS A FUNCTION OF MASS AND RED-SHIFT

     Else

        write(UNIT_EXE_FILE,*) 'READING LINEAR HALO BIAS'

        call read_bMz() ! LINEAR HALO BIAS AS A FUNCTION OF MASS AND RED-SHIFT

     End If

     Do index1=1,10

        Do index2=1,10

            call Interpolate_2D(bMz_M_z(index1,index2),Mdebug(index1),zdebug(index2),M_functions(1:number_of_M_functions),&
                 z_functions(1:number_of_z_functions),bMz(1:number_of_M_functions,1:number_of_z_functions))

            print *, Mdebug(index1),zdebug(index2),(bMz_M_z(index1,index2)- linear_halo_bias(Mdebug(index1),zdebug(index2)))

        End Do

     End Do
     
     stop


     If (compute_alpha_in_halo_mass_function) then

        write(UNIT_EXE_FILE,*) 'COMPUTING NORMALIZATION OF HALO MASS FUNCTION'

        call compute_alpha_halo_mass_function() ! NORMALIZE HALO MASS FUNCTION TO FULLFILL CONDITION THAT MEAN BIAS OF ALL MATTER
                                                ! AT A FIXED RED-SHIFT IS UNITY

     Else 

        write(UNIT_EXE_FILE,*) 'READING NORMALIZATION OF HALO MASS FUNCTION'

        call read_alpha_halo_mass_function() ! NORMALIZE HALO MASS FUNCTION TO FULLFILL CONDITION THAT MEAN BIAS OF ALL MATTER
                                             ! AT A FIXED RED-SHIFT IS UNITY

     End If
     
     If (compute_halo_mass_function) then

        write(UNIT_EXE_FILE,*) 'COMPUTING HALO MASS FUNCTION'

        call compute_dndM() ! COMPUTING HALO MASS FUNCTION AS A FUNCTION OF MASS AND RED-SHIFT 

     Else

        write(UNIT_EXE_FILE,*) 'READING HALO MASS FUNCTION'

        call read_dndM() ! READING HALO MASS FUNCTION    

     End If

     stop

     call Interpolate_2D(dndM_M_z,(M(10)+M(100))/2.d0,(z(10)+z(100))/2.d0,M_functions(1:number_of_M_functions),&
          z_functions(1:number_of_z_functions),dndM(1:number_of_M_functions,1:number_of_z_functions))

     print *, dndM_M_z, halo_mass_function((M(10)+M(100))/2.d0,(z(10)+z(100))/2.d0)

     stop

     write(UNIT_EXE_FILE,*) 'COMPUTING COMOVING VOLUME ELEMENT PER STERADIAN'

     call compute_d2VdzdO()   ! COMPUTES COMOVING VOLUME ELEMENT PER STERADIAN  

     If (compute_the_lensing_potential) then 

        wtime = omp_get_wtime() ! SETTING STARTING TIME OF LENSING POTENTIAL COMPUTATION

        write(UNIT_EXE_FILE,*) 'COMPUTING LENSING POTENTIAL'

        call compute_lensing_potential(halo_definition)

        write(UNIT_EXE_FILE,*) 'LENSING POTENTIAL COMPUTATION FOR RED-SHIFT ARRAY OF SIZE ', size(z),&
             ', VIRIAL MASS ARRAY OF SIZE ',size(M),', AND MULTIPOLES ARRAY OF SIZE ', size(ml),&
             ' TOOK ', (omp_get_wtime()-wtime)/3.6d3, 'HOURS'

     Else

        write(UNIT_EXE_FILE,*) 'READING LENSING POTENTIAL'

        call read_philMz() ! READS LENSING POTENTIAL

     End If

     write(UNIT_EXE_FILE,*) 'COMPUTING ANGULAR POWER SPECTRUM OF LENSING POTENTIAL'

     !call compute_Clphiphi1h() ! ONE HALO TERM

     call compute_Clphiphi2h() ! TWO HALO TERM

     call compute_Clpsilimber()! LIMBER APPROXIMATION 

     call compute_Cl()   ! TOTAL 

     call write_Cl()     ! OF ALL COMPUTATIONS 

     stop

     If (compute_the_form_factor) then

        wtime = omp_get_wtime() ! SETTING STARTING TIME OF FORM FACTOR COMPUTATION

        write(UNIT_EXE_FILE,*) 'COMPUTING FORM FACTOR'

        call compute_form_factor()

        write(UNIT_EXE_FILE,*) 'FORM FACTOR COMPUTATION FOR RED-SHIFT ARRAY OF SIZE ', size(z),&
             ', VIRIAL MASS ARRAY OF SIZE ', size(M),', AND MULTIPOLES ARRAY OF SIZE ', size(ml),&
             ' TOOK ', (omp_get_wtime()-wtime)/3.6d3, 'HOURS'

     Else

        write(UNIT_EXE_FILE,*) 'READING FORM FACTOR'

        call read_ylMz() ! READS FORM FACTOR

     End If
        
     write(UNIT_EXE_FILE,*) 'COMPUTING ANGULAR POWER SPECTRUM OF Y-tSZ CROSS-CORRELATION'

     call compute_Cl1h() ! ONE HALO TERM

     call compute_Cl2h() ! TWO HALO TERM
 
     call compute_Cl()   ! TOTAL 

     call write_Cl()     ! OF ALL COMPUTATIONS 

     write(UNIT_EXE_FILE,*) 'COMPUTATION ENDED'

    ! DEALLOCATING MEMORY
    deallocate (z,M,k,ml,Cl1h,Cl2h,Clphiphi1h,Clphiphi2h,Cl,Clphiphi,&
    d2VdzdO,dndM,ylMz,philMz,bMz,alpha_halo_mass_function,&
    Clpsilimber,comoving_distance_at_z,mbz,M200c,M200d,r200c,r200d,dM200ddM,dM200cdM,&
    z_com_dist,z_functions,M_functions)

    ! CLOSE EXECUTION INFORMATION FILE
    close(UNIT_EXE_FILE)

End Program tSZ




