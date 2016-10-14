Program tSZ

    use fiducial   ! CONTAINS PARAMETERS 
    use arrays     ! CONTAINS ARRAYS 
    use functions  ! CONTAINS ALL THE FUNCTIONS
    use mod_roots
    use omp_lib
    use fgsl
    use integrator

    ! DECLARATION AND INITIALIZATION OF VARIABLES
    Implicit none
    Integer*4 :: index1                                       ! COUNTER
    Real*8 :: wtime,dalpha                          ! STORES TIME OF EXECUTION
    Character(len=15),parameter :: halo_definition = 'virial' ! HALO DEFINITION USED IN THE COMPUTATIONS

    call comoving_distance_at_redshift(z_dec,com_dist_at_z_dec) ! COMPUTE COMOVING DISTANCE AT DECOUPLING

    open(UNIT_EXE_FILE,file=path_to_execution_information)     ! OPEN FILE TO STORE EXECUTION INFORMATION 

    ! ALLOCATING MEMORY FOR : RED-SHIFT, VIRIAL MASS, MULTIPOLES, WAVEVECTOR, ONE- AND TWO-HALO Y-tSZ CROSS CORRELATION AND THEIR SUM, 
    ! MEAN DENSITY MASS, CRITICAL DENSITY MASS, CRITICAL DENSITY RADIUS, DERIVATIVE OF MEAN DENSITY MASS w.r.t VIRIAL MASS, 
    ! MEAN DENSITY RADIUS, ONE- AND TWO-HALO LENSING POTENTIAL AUTO-CORRELATION AND THEIR SUM, NORMALIZATION CONSTANT FOR HALO MASS FUNCTION,
    ! DERIVATIVE OF CRITICAL MASS DENSITY w.r.t VIRIAL MASS, ANGULAR POWER SPECTRUM OF LENSING POTENTIAL IN THE LIMBER APPROXIMATION, 
    ! HALO MASS FUNCTION, FORM FACTOR, LENSING POTENTIAL, COMOVING VOLUME, LINEAR BIAS, MEAN BIAS OF ALL MATTER, CRITICAL SURFACE,
    ! SIGMA SQUARE FOR MEAN DENSITY MASS, DERIVATIVE OF SIGMA SQUARE FOR MEAN DENSITY MASS
    allocate (z(1:number_of_z), M(1:number_of_M), ml(1:number_of_l),Clyy2h(number_of_l),Clyy(number_of_l),&
    Cl1h(1:number_of_l),Cl2h(1:number_of_l),Cl(1:number_of_l),M200d(1:number_of_M,1:number_of_z),&
    M200c(1:number_of_M,1:number_of_z),r200c(1:number_of_M,1:number_of_z),dM200ddM(1:number_of_M,1:number_of_z),&
    r200d(1:number_of_M,1:number_of_z),Clphiphi1h(1:number_of_l),Clphiphi2h(1:number_of_l),&
    Clphiphi(1:number_of_l),alpha_halo_mass_function(1:number_of_z),dM200cdM(1:number_of_M,1:number_of_z),&
    Clpsilimber(1:number_of_l),dndM(1:number_of_M,1:number_of_z),&
    ylMz(1:number_of_l,1:number_of_M,1:number_of_z),inte_cl_yy_1h(number_of_l,number_of_z),&
    philMz(1:number_of_l,1:number_of_M,1:number_of_z),d2VdzdO(1:number_of_z),&
    bMz(1:number_of_M,1:number_of_z),mbz(1:number_of_z),inte_mbz(number_of_M,number_of_z),&
    comoving_distance_at_z(1:number_of_z),Scrit(1:number_of_z),Clyy1h(number_of_l),&
    angular_diameter_distance_at_z(1:number_of_z),inte_cl_yy_2h(number_of_l,number_of_z),&
    dM200ddM_M_z(1:number_of_M,1:number_of_z),inte_dndM(number_of_M,number_of_z),&
    sigma_square_M200d(1:number_of_M,1:number_of_z),f4(number_of_l,number_of_z),&
    dsigma_square_M200d(1:number_of_M,1:number_of_z),f5(number_of_l,number_of_z),&
    zlimber(number_of_z_limber),integrand_limber(number_of_z_limber,number_of_l),&
    inte_pre_cl_phiphi(number_of_l,number_of_M,number_of_z),inte_cl_phiphi_1h(number_of_l,number_of_z),&
    inte_pre_cl_phiphi_2h(number_of_l,number_of_M,number_of_z),inte_cl_yphi_1h(number_of_l,number_of_z),&
    inte_cl_phiphi_2h(number_of_l,number_of_z),inte_pre_cl_yphi(number_of_l,number_of_M,number_of_z),&
    inte_pre_cl_yphi_2h(number_of_l,number_of_M,number_of_z),inte_cl_yphi_2h(number_of_l,number_of_z),stat = status1)

    If (status1 .eq. 0) then
       
        write(UNIT_EXE_FILE,*) 'MEMORY ALLOCATED SUCCESSFULLY'

    Else

        write(UNIT_EXE_FILE,*) 'PROBLEM ALLOCATING MEMORY'

        stop

    End If

    Do  index1 = 1,number_of_z!+1  ! FILLS RED-SHIFT ARRAY.     

       z(index1) = 10**(log10(zmin) + real(index1-1)*(log10(zmax) - log10(zmin))/real(number_of_z-1))

!       If ( index1 .le. 50 ) then
          
 !         z(index1) = zmin + dble(index1-1)*1.d-4

  !     Else If ( index1 .le. 140 ) then 
          
   !       z(index1) = z(50) + dble(index1-50)*1.d-3

    !   Else If ( index1 .le. 230 ) then

     !     z(index1) = z(140) + dble(index1-140)*1.d-2

      ! Else If ( index1 .le. 320) then

       !   z(index1) = z(230) + dble(index1-230)*1.d-1

!       End If

    End Do

    Do index1 = 1, number_of_M ! FILLS VIRIAL MASS ARRAY. UNITS: Solar mass    
       
       M(index1) = 10**(log10(Mmin) + real(index1-1)*(log10(Mmax) - log10(Mmin))/real(number_of_M-1))

    End Do

    Do index1 = 1, number_of_l ! FILLS MULTIPOLE ARRAY.    

        ml(index1) = int(10**(log10(dble(lmin)) + real(index1-1)*(log10(dble(lmax)) - &
        log10(dble(lmin)))/real(number_of_l-1)),4)

    End Do

    ! COMPUTATION STARTS

    call compute_comoving_and_angular_diameter_distance() ! COMPUTE COMOVING DISTANCE, ANGULAR DIAMETER DISTANCE, CRITICAL SURFACE DENSITY,
                                                          ! AND COMOVING VOLUME PER STERADIAN

    If (do_mass_conversion) then ! COMPUTE MASSES

        wtime = omp_get_wtime() ! SETTING STARTING TIME OF MASS CONVERSION

        write(UNIT_EXE_FILE,*) 'EXECUTING MASS CONVERSION'

        call compute_M_delta_c_from_M_and_z(DeltaSO)

        write(UNIT_EXE_FILE,*) 'MASS CONVERSION ENDED'

        write(UNIT_EXE_FILE,*) 'MASS CONVERSION FOR RED-SHIFT ARRAY OF SIZE ', size(z),' AND VIRIAL MASS ARRAY OF SIZE ', size(M),&
        'TOOK ', (omp_get_wtime()-wtime)/3.6d3, 'HOURS'
        
        call read_M200dc_r200dc() ! READING MASS CONVERSION FILE 

        call compute_dMdc_dM_and_dMdd_dM() ! COMPUTING DERIVATIVES OF MASS CONVERSION

        write(UNIT_EXE_FILE,*) 'MASS CONVERSION FILE READ AND DERIVATIVES OF MASS CONVERSION COMPUTED'

     Else

        write(UNIT_EXE_FILE,*) 'USING MASS DATA FILE PREVIOUSLY COMPUTED'

        call read_M200dc_r200dc() ! READING MASS CONVERSION FILE 

        call read_dMdc_dM_and_dMdd_dM() ! COMPUTING DERIVATIVES OF MASS CONVERSION

        write(UNIT_EXE_FILE,*) 'MASS CONVERSION FILE READ AND DERIVATIVES OF MASS CONVERSION COMPUTED'

     End If

     call compute_normalization(Normalization) ! IN MATTER POWER SPECTRUM TO FIT FIDUCIAL SIGMA_8

     write(UNIT_EXE_FILE,*) 'NORMALIZATION OF MATTER POWER SPECTRUM TO MATCH SIGMA_8 (',sigma8,') WAS COMPUTED'

     If (compute_sigma_square) then

        call compute_sigma_square_M200d()

     Else

        call read_sigma_square_M200d()

     End If

     If (compute_linear_halo_bias) then

        write(UNIT_EXE_FILE,*) 'COMPUTING LINEAR HALO BIAS'

        call compute_bMz() ! LINEAR HALO BIAS AS A FUNCTION OF MASS AND RED-SHIFT

     Else

        write(UNIT_EXE_FILE,*) 'READING LINEAR HALO BIAS'

        call read_bMz() ! LINEAR HALO BIAS AS A FUNCTION OF MASS AND RED-SHIFT

     End If

     If (compute_alpha_in_halo_mass_function) then

        write(UNIT_EXE_FILE,*) 'COMPUTING NORMALIZATION OF HALO MASS FUNCTION'

        call compute_integrand_alpha_halo_mass_function()

        call compute_alpha_halo_mass_function() ! NORMALIZE HALO MASS FUNCTION TO FULLFILL CONDITION THAT MEAN BIAS OF ALL MATTER
                                                ! AT A FIXED RED-SHIFT IS UNITY

     Else 

        write(UNIT_EXE_FILE,*) 'READING NORMALIZATION OF HALO MASS FUNCTION'

        call read_alpha_halo_mass_function() ! NORMALIZE HALO MASS FUNCTION TO FULLFILL CONDITION THAT MEAN BIAS OF ALL MATTER
                                             ! AT A FIXED RED-SHIFT IS UNITY

     End If

     call Interpolate_1D(alpha_halo_redshift_3,dalpha,3.d0,z,alpha_halo_mass_function)

     If (compute_halo_mass_function) then

        write(UNIT_EXE_FILE,*) 'COMPUTING HALO MASS FUNCTION'

        call compute_dndM() ! COMPUTING HALO MASS FUNCTION AS A FUNCTION OF MASS AND RED-SHIFT 

     Else

        write(UNIT_EXE_FILE,*) 'READING HALO MASS FUNCTION'

        call read_dndM() ! READING HALO MASS FUNCTION    

     End If

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

     write(UNIT_EXE_FILE,*) 'COMPUTING ANGULAR POWER SPECTRUM OF LENSING POTENTIAL'

     call compute_integrand_limber_approximation()

     call compute_Clpsilimber()! LIMBER APPROXIMATION 

     call compute_integrand_pre_cl_phiphi_at_z_and_l()

     call compute_integrand_cl_phiphi_one_halo_at_z_and_l()

     call compute_integrand_pre_cl_phiphi_2h_at_z_and_l()

     call compute_integrand_cl_phiphi_two_halo_at_z_and_l()

     call compute_Clphiphi1h() ! ONE HALO TERM

     call compute_Clphiphi2h() ! TWO HALO TERM

     write(UNIT_EXE_FILE,*) 'COMPUTING ANGULAR POWER SPECTRUM OF Y-tSZ CROSS-CORRELATION'

     call compute_integrand_pre_cl_yphi_at_z_and_l()

     call compute_integrand_pre_cl_yphi_2h_at_z_and_l() 

     call compute_integrand_cl_yphi_one_halo_at_z_and_l()

     call compute_integrand_cl_yphi_two_halo_at_z_and_l()

     call compute_Clyphi1h() ! ONE HALO TERM

     call compute_Clyphi2h() ! TWO HALO TERM

     call compute_integrand_cl_yy_one_halo_at_z_and_l()

     call compute_integrand_cl_yy_two_halo_at_z_and_l()

     call compute_Clyy1h()

     call compute_Clyy2h()

     call compute_Cl()   ! TOTAL 

     call write_Cl()     ! OF ALL COMPUTATIONS 

     write(UNIT_EXE_FILE,*) 'COMPUTATION ENDED'

    ! DEALLOCATING MEMORY
    deallocate (z,M,ml,Cl1h,Cl2h,Clphiphi1h,Clphiphi2h,Cl,Clphiphi,&
    d2VdzdO,dndM,ylMz,philMz,bMz,alpha_halo_mass_function,&
    Clpsilimber,comoving_distance_at_z,mbz,M200c,M200d,r200c,r200d,dM200ddM,dM200cdM,&
    Scrit,dM200ddM_M_z)

    ! CLOSE EXECUTION INFORMATION FILE
    close(UNIT_EXE_FILE)

End Program tSZ




