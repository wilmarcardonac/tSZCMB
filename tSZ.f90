Program tSZ

    use fiducial   ! CONTAINS PARAMETERS 
    use arrays     ! CONTAINS ARRAYS 
    use functions  ! CONTAINS ALL THE FUNCTIONS
    use mod_roots
    use omp_lib
    use fgsl

    ! DECLARATION AND INITIALIZATION OF VARIABLES
    Implicit none
    Integer*4 :: index1!,index2,index3                                       ! COUNTER
    Real*8 :: wtime!,hola                          ! STORES TIME OF EXECUTION
    Character(len=15),parameter :: halo_definition = 'virial' ! HALO DEFINITION USED IN THE COMPUTATIONS
    !Integer*4,parameter :: hola1=10
    !Integer*4,parameter :: hola2=100
    !Real*8,dimension(hola1) :: zdebug
    !Real*8,dimension(hola2) :: Mdebug
!    Real*8,dimension(hola2,hola1) :: dndM_M_z

    com_dist_at_z_dec = comoving_distance(z_dec) ! COMPUTE COMOVING DISTANCE AT DECOUPLING

    open(UNIT_EXE_FILE,file=path_to_execution_information)     ! OPEN FILE TO STORE EXECUTION INFORMATION 

    If ((number_of_M .lt. 100) .or. (number_of_M_functions .lt. 100)) then

       write(UNIT_EXE_FILE,*) 'MASS ARRAYS TOO SMALL. THEY MUST BE AT LEAST SIZE 100'

       stop

    End If

    ! ALLOCATING MEMORY FOR : RED-SHIFT, VIRIAL MASS, MULTIPOLES, WAVEVECTOR, ONE- AND TWO-HALO Y-tSZ CROSS CORRELATION AND THEIR SUM, 
    ! MEAN DENSITY MASS, CRITICAL DENSITY MASS, CRITICAL DENSITY RADIUS, DERIVATIVE OF MEAN DENSITY MASS w.r.t VIRIAL MASS, 
    ! MEAN DENSITY RADIUS, ONE- AND TWO-HALO LENSING POTENTIAL AUTO-CORRELATION AND THEIR SUM, NORMALIZATION CONSTANT FOR HALO MASS FUNCTION,
    ! DERIVATIVE OF CRITICAL MASS DENSITY w.r.t VIRIAL MASS, ANGULAR POWER SPECTRUM OF LENSING POTENTIAL IN THE LIMBER APPROXIMATION, 
    ! HALO MASS FUNCTION, FORM FACTOR, LENSING POTENTIAL, COMOVING VOLUME, LINEAR BIAS, MEAN BIAS OF ALL MATTER, CRITICAL SURFACE,
    ! SIGMA SQUARE FOR MEAN DENSITY MASS, DERIVATIVE OF SIGMA SQUARE FOR MEAN DENSITY MASS
    allocate (z(1:number_of_z), M(1:number_of_M), ml(1:number_of_l),k(1:number_of_k),&
    Cl1h(1:number_of_l),Cl2h(1:number_of_l),Cl(1:number_of_l),M200d(1:number_of_M,1:number_of_z),&
    M200c(1:number_of_M,1:number_of_z),r200c(1:number_of_M,1:number_of_z),dM200ddM(1:number_of_M,1:number_of_z),&
    r200d(1:number_of_M,1:number_of_z),Clphiphi1h(1:number_of_l),Clphiphi2h(1:number_of_l),&
    Clphiphi(1:number_of_l),alpha_halo_mass_function(1:number_of_z_functions),dM200cdM(1:number_of_M,1:number_of_z),&
    Clpsilimber(1:number_of_l),dndM(1:number_of_M_functions,1:number_of_z_functions),&
    ylMz(1:number_of_l,1:number_of_M_functions,1:number_of_z_functions),&
    philMz(1:number_of_l,1:number_of_M_functions,1:number_of_z_functions),d2VdzdO(1:number_of_z_functions),&
    bMz(1:number_of_M,1:number_of_z),mbz(1:number_of_z_functions),&
    comoving_distance_at_z(1:number_of_z_functions),z_functions(1:number_of_z_functions),Scrit(1:number_of_z_functions),&
    M_functions(1:number_of_M_functions),angular_diameter_distance_at_z(1:number_of_z_functions),&
    dM200ddM_M_z(1:number_of_M_functions,1:number_of_z_functions),stat = status1)

    If (status1 .eq. 0) then
       
        write(UNIT_EXE_FILE,*) 'MEMORY ALLOCATED SUCCESSFULLY'

    Else

        write(UNIT_EXE_FILE,*) 'PROBLEM ALLOCATING MEMORY'

        stop

    End If

    Do index1 = 1, number_of_k ! FILLS WAVEVECTOR ARRAY. UNITS : 1/Mpc

        k(index1) = 10**(log10(kmin) + real(index1-1)*(log10(kmax) - log10(kmin))/real(number_of_k-1))

    End Do

    Do  index1 = 1,number_of_z+1  ! FILLS RED-SHIFT ARRAY.     

       If ( index1 .le. 50 ) then
          
          z(index1) = zmin + dble(index1-1)*1.d-4

       Else If ( index1 .le. 140 ) then 
          
          z(index1) = z(50) + dble(index1-50)*1.d-3

       Else If ( index1 .le. 230 ) then

          z(index1) = z(140) + dble(index1-140)*1.d-2

       Else If ( index1 .le. 320) then

          z(index1) = z(230) + dble(index1-230)*1.d-1

       End If

        !z(index1) = 10**(log10(zmin) + real(index1-1)*(log10(zmax) - log10(zmin))/real(number_of_z-1))

    End Do

    Do index1 = 1, number_of_M ! FILLS VIRIAL MASS ARRAY. UNITS: Solar mass    
       
!       If (index1 .le. number_of_M_log) then

       M(index1) = 10**(log10(Mmin) + real(index1-1)*(log10(Mmax) - log10(Mmin))/real(number_of_M-1))

!       Else

!       If ( index1 .le. 90 ) then
          
 !         M(index1) = 1.d5 + dble(index1 - 1)*1.d4

  !     Else If ( index1 .le. 180 ) then 
          
   !       M(index1) = M(90) + dble(index1 - 90)*1.d5

    !   Else If ( index1 .le. 270 ) then

     !     M(index1) = M(180) + dble(index1 - 180)*1.d6

      ! Else If ( index1 .le. 360 ) then

       !   M(index1) = M(270) + dble(index1 - 270)*1.d7

       !Else If ( index1 .le. 450 ) then

!          M(index1) = M(360) + dble(index1 - 360)*1.d8

 !      Else If ( index1 .le. 540 ) then

  !        M(index1) = M(450) + dble(index1 - 450)*1.d9

   !    Else If ( index1 .le. 630 ) then

    !      M(index1) = M(540) + dble(index1 - 540)*1.d10

     !  Else If ( index1 .le. 720 ) then

      !    M(index1) = M(630) + dble(index1 - 630)*1.d11

!       Else If ( index1 .le. 810 ) then

 !         M(index1) = M(720) + dble(index1 - 720)*1.d12

  !     Else If ( index1 .le. 900 ) then

   !       M(index1) = M(810) + dble(index1 - 810)*1.d13

    !   Else 

     !     M(index1) = M(900) + dble(index1 - 900)*1.d14

!       Else 

!          M(index1) = M(number_of_M_log+270) + dble(index1 - number_of_M_log - 270)*1.d14

      ! End If

    !End If

    End Do

!    print *, M_delta_d_from_M_virial(3.d0,r_delta_d_from_M_virial(M(10),3.d0,DeltaSO),DeltaSO)

!    stop

!    Do index1 = -1, number_of_M+2 ! FILLS VIRIAL MASS ARRAY. UNITS: Solar mass    

!       M(index1) = M(index1)/h

!       If (M(index1) .lt. M(index1-1)) then

!          print *, 'FOUND', index1

!          stop

!       End If

!    End Do

!    print *, M(-1),M(0),M(1),Mmin,M(number_of_M),Mmax

    !stop

    Do index1 = 1, number_of_l ! FILLS MULTIPOLE ARRAY.    

        ml(index1) = int(10**(log10(dble(lmin)) + real(index1-1)*(log10(dble(lmax)) - &
        log10(dble(lmin)))/real(number_of_l-1)),4)

    End Do

!    If (compute_functions) then 

!       continue

!    Else

!       Do index1 = 1, number_of_z_functions ! FILLS RED-SHIFT ARRAY.     

!          z_functions(index1) = 10**(log10(zmin) + real(index1-1)*(log10(zmax) - log10(zmin))/real(number_of_z_functions-1))

!       End Do
    
!       Do index1 = 1, number_of_M_functions ! FILLS VIRIAL MASS ARRAY. UNITS: Solar mass    

!          M_functions(index1) = 10**(log10(Mmin) + real(index1-1)*(log10(Mmax) - log10(Mmin))/real(number_of_M_functions-1))

!       End Do

!       Do index1 = 1, number_of_l_functions ! FILLS MULTIPOLE ARRAY.    

!          ml_functions(index1) = int(10**(log10(dble(lmin)) + real(index1-1)*(log10(dble(lmax)) - &
!               log10(dble(lmin)))/real(number_of_l_functions-1)),4)

!       End Do

!    End If 

    ! COMPUTATION STARTS

    call compute_comoving_distance() ! COMPUTE COMOVING DISTANCE ARRAY

    call compute_angular_diameter_distance() ! COMPUTE ANGULAR DIAMETER DISTANCE ARRAY

    call compute_critical_surface_density()

    If (do_mass_conversion) then ! COMPUTE MASSES

        wtime = omp_get_wtime() ! SETTING STARTING TIME OF MASS CONVERSION

        write(UNIT_EXE_FILE,*) 'EXECUTING MASS CONVERSION'

        call compute_M_delta_c_from_M_and_z(DeltaSO)

        write(UNIT_EXE_FILE,*) 'MASS CONVERSION ENDED'

        write(UNIT_EXE_FILE,*) 'MASS CONVERSION FOR RED-SHIFT ARRAY OF SIZE ', size(z),' AND VIRIAL MASS ARRAY OF SIZE ', size(M),&
        'TOOK ', (omp_get_wtime()-wtime)/3.6d3, 'HOURS'
        
        !call read_M200dc_r200dc() ! READING MASS CONVERSION FILE AND COMPUTING DERIVATIVES OF MASS CONVERSION

        call compute_dMdc_dM_and_dMdd_dM()

        write(UNIT_EXE_FILE,*) 'MASS CONVERSION FILE READ AND DERIVATIVES OF MASS CONVERSION COMPUTED'

     Else

        write(UNIT_EXE_FILE,*) 'USING MASS DATA FILE PREVIOUSLY COMPUTED'

        call read_M200dc_r200dc() ! READING MASS CONVERSION FILE AND COMPUTING DERIVATIVES OF MASS CONVERSION

        call compute_dMdc_dM_and_dMdd_dM()

        write(UNIT_EXE_FILE,*) 'MASS CONVERSION FILE READ AND DERIVATIVES OF MASS CONVERSION COMPUTED'

     End If

     print *, (3.d0*M200d(400,320)/4.d0/Pi/mean_density(z(320))*(1.d0 + z(320))**3.d0)**(1.d0/3.d0)!sigma_squared(M200d(400,320),z(320))

     stop
     !!!!!!!!!!!!!!!!!!!!!!!!!DEBUGGING 
     !Do index1 = 1, hola1 ! FILLS RED-SHIFT ARRAY.     

     !   zdebug(index1) = 10**(log10(zmin) + real(index1-1)*(log10(zmax) - log10(zmin))/real(hola1-1))

     !End Do
    
     !Do index1 = 1, hola2 ! FILLS VIRIAL MASS ARRAY. UNITS: Solar mass    
        
     !   Mdebug(index1) = 10**(log10(Mmin) + real(index1-1)*(log10(Mmax) - log10(Mmin))/real(hola2-1))
        
     !End Do
     !!END DEBUGGING

     call compute_normalization() ! IN MATTER POWER SPECTRUM TO FIT FIDUCIAL SIGMA_8

     write(UNIT_EXE_FILE,*) 'NORMALIZATION OF MATTER POWER SPECTRUM TO MATCH SIGMA_8 (',sigma8,') WAS COMPUTED'

     If (compute_linear_halo_bias) then

        write(UNIT_EXE_FILE,*) 'COMPUTING LINEAR HALO BIAS'

        call compute_bMz() ! LINEAR HALO BIAS AS A FUNCTION OF MASS AND RED-SHIFT

     Else

        write(UNIT_EXE_FILE,*) 'READING LINEAR HALO BIAS'

        call read_bMz() ! LINEAR HALO BIAS AS A FUNCTION OF MASS AND RED-SHIFT

        If (compute_functions) then

           continue

        Else

           allocate(bMz_interpolation(1:number_of_M_functions,1:number_of_z_functions),&
                stat=status2)

           If (status2 .eq. 0) then

              call Interpolate_bMz()

              deallocate(bMz)

           Else

              write(UNIT_EXE_FILE,*) 'PROBLEM ALLOCATING MEMORY FOR "bMz_interpolation" '

              stop

           End If

        End If

     End If

     !call compute_dM200ddM_M_z()

     If (compute_alpha_in_halo_mass_function) then

        write(UNIT_EXE_FILE,*) 'COMPUTING NORMALIZATION OF HALO MASS FUNCTION'

        call compute_alpha_halo_mass_function() ! NORMALIZE HALO MASS FUNCTION TO FULLFILL CONDITION THAT MEAN BIAS OF ALL MATTER
                                                ! AT A FIXED RED-SHIFT IS UNITY

     Else 

        write(UNIT_EXE_FILE,*) 'READING NORMALIZATION OF HALO MASS FUNCTION'

        call read_alpha_halo_mass_function() ! NORMALIZE HALO MASS FUNCTION TO FULLFILL CONDITION THAT MEAN BIAS OF ALL MATTER
                                             ! AT A FIXED RED-SHIFT IS UNITY

     End If
     
     call Interpolate_1D(alpha_halo_redshift_3,d_alpha_halo_redshift_3,3.d0,z,alpha_halo_mass_function)

     If (compute_halo_mass_function) then

        write(UNIT_EXE_FILE,*) 'COMPUTING HALO MASS FUNCTION'

        call compute_dndM() ! COMPUTING HALO MASS FUNCTION AS A FUNCTION OF MASS AND RED-SHIFT 

     Else

        write(UNIT_EXE_FILE,*) 'READING HALO MASS FUNCTION'

        call read_dndM() ! READING HALO MASS FUNCTION    

        If (compute_functions) then

           continue

        Else

           allocate(dndM_interpolation(1:number_of_M_functions,1:number_of_z_functions),&
                stat=status2)

           If (status2 .eq. 0) then

              call Interpolate_dndM()

              deallocate(dndM)

           Else

              write(UNIT_EXE_FILE,*) 'PROBLEM ALLOCATING MEMORY FOR "dndM_interpolation" '

              stop

           End If

        End If

     End If

     call compute_mean_bias_matter()

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

        If (compute_functions) then

           continue

        Else

           allocate(philMz_interpolation(1:number_of_l,1:number_of_M_functions,1:number_of_z_functions),&
                stat=status2)

           If (status2 .eq. 0) then

              call Interpolate_philMz()

              deallocate(philMz)

           Else

              write(UNIT_EXE_FILE,*) 'PROBLEM ALLOCATING MEMORY FOR "philMz_interpolation" '

              stop

           End If

        End If

     End If

!     Do index3=1,number_of_l

!        Do index1=1,hola2

!           Do index2=1,hola1

!              call Interpolate_2D(dndM_M_z(index1,index2),Mdebug(index1),zdebug(index2),M_functions(1:number_of_M_functions),&
!                   z_functions(1:number_of_z_functions),philMz(index3,1:number_of_M_functions,1:number_of_z_functions))

!              hola = (dndM_M_z(index1,index2)- halo_mass_function(Mdebug(index1),zdebug(index2)))/&
!                   halo_mass_function(Mdebug(index1),zdebug(index2))*100.d0

!              If (hola .gt. 1.d-5) then

!                 print *, Mdebug(index1),zdebug(index2),(dndM_M_z(index1,index2)- lensing_potential(Mdebug(index1),zdebug(index2),&
!                      index3,halo_definition))

!              End If

!           End Do

!        End Do

!     End Do
     
!     stop

!     print *, pre_Clphiphi(320,3)

!     stop

     write(UNIT_EXE_FILE,*) 'COMPUTING ANGULAR POWER SPECTRUM OF LENSING POTENTIAL'

!     call compute_Clphiphi1h() ! ONE HALO TERM

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

        If (compute_functions) then

           continue

        Else

           allocate(ylMz_interpolation(1:number_of_l,1:number_of_M_functions,1:number_of_z_functions),&
                stat=status2)

           If (status2 .eq. 0) then

              call Interpolate_ylMz()

              deallocate(ylMz)

           Else

              write(UNIT_EXE_FILE,*) 'PROBLEM ALLOCATING MEMORY FOR "ylMz_interpolation" '

              stop

           End If

        End If

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
    z_functions,M_functions,Scrit,dM200ddM_M_z)

    ! CLOSE EXECUTION INFORMATION FILE
    close(UNIT_EXE_FILE)

End Program tSZ




