program main
  !! @note  Main function of KORC. @endnote
  !! The main program contains the calls to the main functions and subroutines. 
  !! Also, it contains the variables that control
  !! the behavior of the core of KORC and all other external/optional modules.

  use korc_types
  use korc_units
  use korc_hpc
  use korc_HDF5
  use korc_fields
  use korc_ppusher
  use korc_interp
  use korc_collisions
  use korc_initialize
  use korc_finalize
  use korc_profiles
  use korc_input
#ifdef FIO
  use korc_fio
#endif
  
  implicit none

  TYPE(KORC_PARAMS) :: params
  !! Contains the parameters that control the core of KORC: 
  !! time steping, output list, etc.
  TYPE(SPECIES), DIMENSION(:), ALLOCATABLE :: spp
  !! Contains the initial parameters of each species, which 
  !! can be different electrons with different
  !! distribution functions.
  TYPE(FIELDS) :: F
  !! F: Contains the parameters of the analytical magnetic 
  !! and electric fields, or in the case of using 
  !! external fields it contains the data used in the interpolations. 
  !!See [[korc_fields(module)]] for details.
  TYPE(PROFILES) :: P
  !! P: Contains the parameters of the analytical plasma profiles, 
  !! or in the case of using external 
  !! fields it contains the data used in the interpolations. 
  !! See [[korc_profiles(module)]] for details.
  INTEGER(ip) :: it 
  !! Time iteration
  INTEGER 				:: mpierr
    
  call initialize_communications(params)
  !!<h2>Order of KORC operations</h2>
  !!
  !!<h3>Communication and Timing</h3>
  !! <h4>1\. Parallel Communications</h4>
  !!
  !! Subroutine [[initialize_communications]] in [[korc_hpc]] that 
  !! initializes MPI and OpenMP communications.

  if (params%mpi_params%rank .EQ. 0) then
     flush(output_unit_write)
  end if
  
  call timing_KORC(params)
  !! <h4>2\. Timers</h4>
  !!
  !! Subroutine [[timing_KORC]] in [[korc_hpc]] that times the 
  !! execution of any parallel sections of KORC.
  
  ! * * * INITIALIZATION STAGE * * *!

  if (params%mpi_params%rank .EQ. 0) then
     flush(output_unit_write)
  end if
  
  call initialize_HDF5()
  !!<h3>Initialization</h3>
  !!
  !! <h4>1\. HDF5</h4>
  !!
  !! Subroutine [[initialize_HDF5]] in [[korc_HDF5]] that initializes
  !! HDF5 library. 

  if (params%mpi_params%rank .EQ. 0) then
     flush(output_unit_write)
  end if
  
  call initialize_korc_parameters(params)
  !! <h4>2\. Initialize korc parameters</h4>
  !!
  !! Subroutine [[initialize_korc_parameters]] in [[korc_initialize]] that 
  !! initializes paths and KORC parameters through [[load_korc_params]]
  !! on MPI processes.

  if (params%mpi_params%rank .EQ. 0) then
     flush(output_unit_write)
  end if
  
  call initialize_fields(params,F)
  !! <h4>3\. Initialize fields</h4>
  !!
  !! Subroutine [[initialize_fields]] in [[korc_fields]] that initializes 
  !! parameters of the EM fields, either analytically or from an external HDF5
  !! file. Reads in &amp;analytical_fields_params and 
  !! &amp;externalPlasmaModel namelists from input file.

  if (params%mpi_params%rank .EQ. 0) then
     flush(output_unit_write)
  end if

  call initialize_profiles(params,P,F)
  !! <h4>4\. Initialize Profiles</h4>
  !! 
  !! Subroutine [[initialize_profiles]] in [[korc_profiles]] that initializes 
  !! parameters of the plasma profiles, either analytically or from an
  !! external HDF5
  !! file. Reads in &amp;plasmaProfiles namelist from input file.
  !! Only initialized if collisions (params%collisions==T) are

  if (params%mpi_params%rank .EQ. 0) then
     flush(output_unit_write)
  end if
  
  call initialize_particles(params,F,P,spp) ! Initialize particles
  !! <h4>5\. Initialize Particle Velocity Phase Space</h4>
  !! 
  !! Subroutine [[initialize_particles]] in [[korc_initialize]] that 
  !! initializes particle parameters from &amplasma_species namelist, 
  !! allocates arrays for individual particles, including location, velocity, 
  !! local EM fields and plasma profiles, etc., and 
  !! calls [[initial_energy_pitch_dist]] to assign particles' energy and pitch
  !! angle according to the chosen distribution.

!  write(output_unit_write,'("init eta: ",E17.10)') spp(1)%vars%eta

  
#ifdef FIO
  if (TRIM(params%field_model) .eq. 'M3D_C1') then

     if (params%mpi_params%rank .EQ. 0) then
        write(output_unit_write,*) "* * * * INITIALIZING M3D-C1 INTERFACE * * * *"
     endif
     
     call initialize_m3d_c1(params, F, P, spp,.true.)

     if (params%mpi_params%rank .EQ. 0) then
        write(output_unit_write,*) "* * * * * * * * * * * * * * * * * * * * * * *"
     endif

  elseif (TRIM(params%field_model) .eq. 'NIMROD') THEN

     if (params%mpi_params%rank .EQ. 0) then
        write(output_unit_write,*) "* * * * INITIALIZING NIMROD INTERFACE * * * *"
     endif
     
     call initialize_nimrod(params, F, P, spp,.true.)

     if (params%mpi_params%rank .EQ. 0) then
        write(output_unit_write,*) "* * * * * * * * * * * * * * * * * * * * * * *"
     endif
     
  endif
#endif  

  if (params%mpi_params%rank .EQ. 0) then
     flush(output_unit_write)
  end if
  
  call compute_charcs_plasma_params(params,spp,F)
  !! <h4>9\. Compute Characteristic Plasma Parameters</h4>
  !!
  !! Subroutine [[compute_charcs_plasma_params]] in [[korc_units]] calculates
  !! the characteristic plasma parameters params%cpp that are used for normalizations.
  !! Also finds the maximum non-relativistic and relativistic cyclotron frequencies
  !! to be used for setting the timstep for the time-evolution algorithms.

  call initialize_collision_params(params,spp,P,F,.true.)
  !! <h4>6\. Initialize Collision Parameters</h4>
  !!
  !! Subroutine [[initialize_collision_params]] in [[korc_collisions]] that
  !! initializes collision parameters for the SS (single-species) and MS
  !! (multiple-species) data types, reading in namefiles from the KORC input file.
  !! MS reads in namelist &CollisionParamsMultipleSpecies while SS reads in
  !! namelist &CollisionParamsSingleSpecies. 
  
  call define_time_step(params,F)
  !! <h4>10\. Define Time Step</h4>
  !!
  !! Subroutine [[define_time_step]] in [[korc_initialize]] either loads
  !! time-stepping parameters for a restart, or defines new parameters based
  !! on a maximum timestep
  !! set by the inverse of the relativistic cyclotron frequency.

  
  if (params%mpi_params%rank .EQ. 0) then
     flush(output_unit_write)
  end if
  
  call initialize_particle_pusher(params)
  !! <h4>11\. Initialize Particle Pusher</h4>    
  
  if (params%SC_E) then
     call define_SC_time_step(params,F)
  end if
     
  call normalize_variables(params,spp,F,P)
  !! <h4>12\. Normalize Variables</h4>
  !!
  !! Subroutine [[normalize_variables]] in [[korc_units]] normalizes 
  !! variables consistent with characteristic plasma parameters 
  !! calculated in [[compute_charcs_plasma_params]].


  
  call normalize_collisions_params(params)
  !! <h4>13\. Normalize Collision Parameters </h4>
  !!
  !! Subroutine [[normalize_collisions_params]] in [[korc_collisions]] that
  !! normalizes collision parameters for the SS (single-species) and MS
  !! (multiple-species) data types.

  
  call define_collisions_time_step(params,F,.true.)
  !! <h4>14\. Define Collision Time Step</h4>
  !!
  !! Subroutine [[define_collisions_time_step]] in [[korc_collisions]] that
  !! sets subcycling iteration number for collisions based off of the collision
  !! frequency model used.

  ! *** *** *** *** *** ***   *** *** *** *** *** *** ***
  ! *** BEYOND THIS POINT VARIABLES ARE DIMENSIONLESS ***
  ! *** *** *** *** *** ***   *** *** *** *** *** *** ***

  if (params%mpi_params%rank .EQ. 0) then
     flush(output_unit_write)
  end if

#ifdef PSPLINE
  call initialize_fields_interpolant(params,F)
  !! <h4>15\. Initialize Fields Interpolant</h4>
  !!
  !! Subroutine [[initialize_fields_interpolant]] in [[korc_interp]] calls
  !! EZspline
  !! subroutines EZspline_init for memory allocation and boundary condition
  !! setup
  !! and EZspline_setup to compute the necessary cubic coefficients needed
  !! for subsequent
  !! field interpolations. The magnetic field can be defined in terms of an
  !! axisymmetric
  !! scalar flux function, axisymmetric field, or 3D field, while the
  !! electric field
  !! can be defined as an axisymmetric or 3D field.

  if (params%mpi_params%rank .EQ. 0) then
     flush(output_unit_write)
  end if
  
  call initialize_profiles_interpolant(params,P)
#endif
  !! <h4>16\. Initialize Profiles Interpolant</h4>
  !!
  !! Subroutine [[initialize_profiles_interpolant]] in [[korc_interp]]
  !! calls EZspline
  !! subroutines EZlinear_init for axisymmetric (flux-surface quantities) or
  !! EZspline_init for 3D profiles for memory allocation and boundary
  !! condition setup
  !! and EZspline_setup to compute the necessary cubic coefficients needed
  !! for subsequent
  !! field interpolations. 
  !! Only initialized if collisions (params%collisions==T) are present for
  !! ne, Te, Zeff

  if (params%mpi_params%rank .EQ. 0) then
     flush(output_unit_write)
  end if
  
  if (params%mpi_params%rank .EQ. 0) then
     write(output_unit_write,'("* * * * INITIALIZING INITIAL CONDITIONS * * * *",/)')
     flush(output_unit_write)
  end if
  call set_up_particles_ic(params,F,spp,P)
  
  if (params%mpi_params%rank .EQ. 0) then
     write(output_unit_write,'("* * * * * * * * * * * * * * * * * * * * * * * *",/)')
     flush(output_unit_write)
  end if

  !write(6,*) 'V',spp(1)%vars%V
  
!  write(output_unit_write,'("post ic eta: ",E17.10)') spp(1)%vars%eta
  
  !! <h4>17\. Set Particle Initial Conditions</h4>  
  !!
  !! Subroutine [[set_up_particles_ic]] in [[korc_initialize]] calls
  !! subroutines to prescribe initial conditions or load them 
  !! from file for a restart. Initial spatial values are prescribed with 
  !! [[intitial_spatial_distribution]] in [[korc_spatial_distribution]] and 
  !! initial velocity values are prescribed with [[initial_gyro_distribution]]
  !! in [[korc_velocity_distribution]].

!  if (minval(spp(1)%vars%Y(:,1)).lt.1._rp/params%cpp%length) stop 'error with init'
  
  ! * * * INITIALIZATION STAGE * * *

  
  if (params%mpi_params%rank .EQ. 0) then
     flush(output_unit_write)
  end if

!  write(output_unit_write,'("GC init eta: ",E17.10)') spp(1)%vars%eta

  if (.NOT.(params%restart.OR.params%proceed.or.params%reinit)) then
     if (params%orbit_model(1:2).eq.'FO') then

        call FO_init(params,F,spp,.true.,.false.)

     else if (params%orbit_model(1:2).eq.'GC') then

        call GC_init(params,F,spp)

     end if

     if (params%SC_E) then

        if (params%field_model(1:1).eq.'A') then
           call init_SC_E1D(params,F,spp(1))
        else if (params%field_model(1:1).eq.'E') then
           call init_SC_E1D_FS(params,F,spp(1))
        end if
           
     end if

  else

     call get_fields(params,spp(1)%vars,F)

     if (params%SC_E) then


        if (params%field_model(1:1).eq.'A') then
           call reinit_SC_E1D(params,F)
        else if (params%field_model(1:1).eq.'E') then
           call reinit_SC_E1D_FS(params,F)
        end if

        
     
     end if
     
  end if

  !write(6,*) 'V',spp(1)%vars%V
  !write(6,*) 'eta',spp(1)%vars%eta
  
!  write(6,*) '1Y_R',spp(1)%vars%Y(1:4,1)*params%cpp%length
  
  ! * * * SAVING INITIAL CONDITION AND VARIOUS SIMULATION PARAMETERS * * * !
  
  call save_simulation_parameters(params,spp,F,P)

  call save_collision_params(params)
  !! <h4>18\. Save Simulation and Collision Parameters</h4>  
  !!
  !! Subroutines [[save_simulation_parameters]] in [[korc_HDF5]] and
  !! [[save_collision_params]] in [[korc_collisions]] call
  !! subroutines to save simulation and collision parameters.
  
  
  if (.NOT.(params%restart.OR.params%proceed.or.params%reinit)) then
     
     call save_simulation_outputs(params,spp,F) ! Save initial condition
     call save_restart_variables(params,spp,F)

  end if
  
  
  ! * * * SAVING INITIAL CONDITION AND VARIOUS SIMULATION PARAMETERS * * * !

!  write(output_unit_write,'("pre ppusher loop eta: ",E17.10)') spp(1)%vars%eta

  call timing_KORC(params)
  if (params%mpi_params%rank .EQ. 0) then
     flush(output_unit_write)
  end if

  if (params%orbit_model(1:2).eq.'FO'.and.((params%field_model(1:3).eq.'ANA') &
       .or.(params%field_model(1:3).eq.'UNI'))) then
     call FO_init(params,F,spp,.false.,.true.)
     ! Initial half-time particle push
     
     do it=params%ito,params%t_steps,params%t_skip
        call adv_FOeqn_top(params,F,P,spp)

        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip,rp)*params%dt        
        params%it = it-1_ip+params%t_skip

        call save_simulation_outputs(params,spp,F)
        call save_restart_variables(params,spp,F)
     end do
  end if

#ifdef PSPLINE
  if (params%orbit_model(1:2).eq.'FO'.and.params%field_model(1:3).eq.'EXT' &
       .and..not.((params%field_model(10:13).eq.'MARS').or. &
       (params%field_model(10:14).eq.'AORSA'))) then
     call FO_init(params,F,spp,.false.,.true.)
     ! Initial half-time particle push
     
     do it=params%ito,params%t_steps,params%t_skip
        call adv_FOinterp_top(params,F,P,spp)
        
        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip,rp)*params%dt        
        params%it = it-1_ip+params%t_skip

        call save_simulation_outputs(params,spp,F)
        call save_restart_variables(params,spp,F)
     end do
  end if
#endif

#ifdef FIO
  if (params%orbit_model(1:2).eq.'FO'.and. &
       (TRIM(params%field_model).eq.'M3D_C1'.or. &
       TRIM(params%field_model).eq.'NIMROD').and. &
       .not.F%ReInterp_2x1t) then
     call FO_init(params,F,spp,.false.,.true.)
     ! Initial half-time particle push
     
     do it=params%ito,params%t_steps,params%t_skip
        call adv_FOfio_top(params,F,P,spp)
        
        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip,rp)*params%dt        
        params%it = it-1_ip+params%t_skip

        call save_simulation_outputs(params,spp,F)
        call save_restart_variables(params,spp,F)
     end do
  end if

  if (params%orbit_model(1:2).eq.'FO'.and. &
       (params%field_model.eq.'M3D_C1'.or. &
       TRIM(params%field_model).eq.'NIMROD') &
       .and.F%ReInterp_2x1t) then
     call FO_init(params,F,spp,.false.,.true.)
     ! Initial half-time particle push
     
     do it=F%ind0_2x1t,params%time_slice

        !write(6,*) it,F%ind0_2x1t
        
        if (it.gt.F%ind0_2x1t) then
           if (params%field_model.eq.'M3D_C1') then
              call initialize_m3d_c1(params, F, P, spp,.false.)
           else
              call initialize_nimrod(params, F, P, spp,.false.)
           end if
           if (params%collisions) then
              if (params%field_model.eq.'M3D_C1') then
                 call initialize_m3d_c1_imp(params,F,P, &
                      params%num_impurity_species,.false.)
              end if
           end if
        end if

        if (params%mpi_params%rank .EQ. 0) then
           write(output_unit_write,*) 'tskip',params%t_skip
           flush(output_unit_write)
        end if

        call adv_FOfio_top(params,F,P,spp)
               
        params%it = params%it+params%t_skip
        params%time = params%init_time &
             +REAL(params%it,rp)*params%dt 
        
        call save_simulation_outputs(params,spp,F)
        call save_restart_variables(params,spp,F)

        F%ind_2x1t=F%ind_2x1t+1_ip
        if (params%mpi_params%rank .EQ. 0) then
           write(output_unit_write,*) 'KORC time ',params%time*params%cpp%time
           flush(output_unit_write)
        end if
              
     end do
     
  end if
#endif

#ifdef PSPLINE
  if (params%orbit_model(1:2).eq.'FO'.and. &
       params%field_model(10:13).eq.'MARS') then
     call FO_init(params,F,spp,.false.,.true.)
     ! Initial half-time particle push
     
     do it=params%ito,params%t_steps,params%t_skip
        call adv_FOinterp_mars_top(params,F,P,spp)
        
        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip,rp)*params%dt        
        params%it = it-1_ip+params%t_skip

        call save_simulation_outputs(params,spp,F)
        call save_restart_variables(params,spp,F)
     end do
  end if

  if (params%orbit_model(1:2).eq.'FO'.and. &
       params%field_model(10:14).eq.'AORSA') then
     call FO_init(params,F,spp,.false.,.true.)
     ! Initial half-time particle push

     do it=params%ito,params%t_steps,params%t_skip
        call adv_FOinterp_aorsa_top(params,F,P,spp)

        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip,rp)*params%dt        
        params%it = it-1_ip+params%t_skip

        call save_simulation_outputs(params,spp,F)
        call save_restart_variables(params,spp,F)
     end do
  end if
#endif
  
  if (params%orbit_model(1:2).eq.'GC'.and.params%field_eval.eq.'eqn'.and..not.params%field_model.eq.'M3D_C1') then
     do it=params%ito,params%t_steps,params%t_skip*params%t_it_SC
        call adv_GCeqn_top(params,F,P,spp)
        
        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip*params%t_it_SC,rp)*params%dt        
        params%it = it-1_ip+params%t_skip*params%t_it_SC

        call save_simulation_outputs(params,spp,F)        
        call save_restart_variables(params,spp,F)

        if (params%mpi_params%rank .EQ. 0) then
           flush(output_unit_write)
        end if
        
     end do
  end if

#ifdef PSPLINE
  if (params%orbit_model(1:2).eq.'GC'.and.params%field_eval.eq.'interp'.and. &
       F%axisymmetric_fields.and.params%field_model(10:12).eq.'PSI'.and. &
       params%SC_E.and..not.params%field_model.eq.'M3D_C1') then
     do it=params%ito,params%t_steps,params%t_skip
        call adv_GCinterp_psi_top_FS(params,spp,P,F)
        
        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip,rp)*params%dt        
        params%it = it-1_ip+params%t_skip

        
        call save_simulation_outputs(params,spp,F)
        call save_restart_variables(params,spp,F)
     end do
  end if

  if (params%orbit_model(1:2).eq.'GC'.and.params%field_eval.eq.'interp'.and. &
       F%axisymmetric_fields.and.(params%field_model(10:12).eq.'PSI'.OR. &
       params%field_model(12:14).eq.'PSI').and. &
       (.not.params%SC_E).and.(.not.F%Dim2x1t).and..not.params%field_model.eq.'M3D_C1') then
     
     do it=params%ito,params%t_steps,params%t_skip
        call adv_GCinterp_psi_top(params,spp,P,F)
        
        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip,rp)*params%dt        
        params%it = it-1_ip+params%t_skip
        
        call save_simulation_outputs(params,spp,F)
        call save_restart_variables(params,spp,F)
     end do
  end if

  if (params%orbit_model(1:2).eq.'GC'.and.params%field_eval.eq.'interp'.and. &
       F%axisymmetric_fields.and.(params%field_model(10:12).eq.'PSI'.OR. &
       params%field_model(12:14).eq.'PSI').and. &
       (.not.params%SC_E).and.F%Dim2x1t.and.(.not.F%ReInterp_2x1t).and..not.params%field_model.eq.'M3D_C1') then
     
     do it=params%ito,params%t_steps,params%t_skip
        call adv_GCinterp_psi2x1t_top(params,spp,P,F)

        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip,rp)*params%dt        
        params%it = it-1_ip+params%t_skip

        call save_simulation_outputs(params,spp,F)
        call save_restart_variables(params,spp,F)
        
     end do
  end if
  
  if (params%orbit_model(1:2).eq.'GC'.and. &
       params%field_eval.eq.'interp'.and. &
       F%axisymmetric_fields.and. &
       (params%field_model(10:12).eq.'PSI'.OR. &
       params%field_model(12:14).eq.'PSI').and. &
       (.not.params%SC_E).and. &
       F%Dim2x1t.and.F%ReInterp_2x1t.and. &
       .not.params%field_model.eq.'M3D_C1') then

     if (params%mpi_params%rank .EQ. 0) then
        write(output_unit_write,*) 'initial 2x1t_ind time',F%X%PHI(F%ind_2x1t)*params%cpp%time
        flush(output_unit_write)  
     end if
        
     do it=params%ito,params%t_steps,params%t_skip
        call adv_GCinterp_psiwE_top(params,spp,P,F)

        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip,rp)*params%dt        
        params%it = it-1_ip+params%t_skip

        call save_simulation_outputs(params,spp,F)

        F%ind_2x1t=F%ind_2x1t+1_ip
        if (params%mpi_params%rank .EQ. 0) then
           write(output_unit_write,*) 'KORC time',params%time*params%cpp%time
           write(output_unit_write,*) '2x1t_ind time',F%X%PHI(F%ind_2x1t)*params%cpp%time
        end if
        call initialize_fields_interpolant(params,F)

        if (params%LargeCollisions) then
           call initialize_collision_params(params,spp,P,F,.false.)
           call define_collisions_time_step(params,F,.false.)
        end if

        call save_restart_variables(params,spp,F)
        
        if (params%mpi_params%rank .EQ. 0) then
           flush(output_unit_write)  
        end if
        
     end do
  end if
  

  if (params%orbit_model(1:2).eq.'GC'.and.params%field_eval.eq.'interp'.and. &
       F%axisymmetric_fields.and.F%dBfield.and..not.params%field_model.eq.'M3D_C1') then
     do it=params%ito,params%t_steps,params%t_skip
        call adv_GCinterp_2DBdB_top(params,spp,P,F)
        
        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip,rp)*params%dt        
        params%it = it-1_ip+params%t_skip

        call save_simulation_outputs(params,spp,F)
        call save_restart_variables(params,spp,F)
     end do
  end if
  
  if (params%orbit_model(1:2).eq.'GC'.and.params%field_eval.eq.'interp'.and. &
       F%axisymmetric_fields.and.(params%field_model(10:12).eq.'2DB'.or. &
       params%field_model(12:13).eq.'2D').and..not.(F%dBfield).and..not.params%field_model.eq.'M3D_C1') then
     do it=params%ito,params%t_steps,params%t_skip
        call adv_GCinterp_B2D_top(params,spp,P,F)
        
        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip,rp)*params%dt        
        params%it = it-1_ip+params%t_skip

        call save_simulation_outputs(params,spp,F)
        call save_restart_variables(params,spp,F)
     end do
  end if
  


  if (params%orbit_model(1:2).eq.'GC'.and.params%field_eval.eq.'interp'.and. &
         .not.(F%axisymmetric_fields).and.(F%dBfield).and. &
         (params%field_model(10:14).eq.'3DBdB').and..not.params%field_model.eq.'M3D_C1') then
     do it=params%ito,params%t_steps,params%t_skip
        call adv_GCinterp_3DBdB_top(params,spp,P,F)
        
        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip,rp)*params%dt        
        params%it = it-1_ip+params%t_skip

        call save_simulation_outputs(params,spp,F)
        call save_restart_variables(params,spp,F)
     end do
  end if

  if (params%orbit_model(1:2).eq.'GC'.and.params%field_eval.eq.'interp'.and. &
         .not.(F%axisymmetric_fields).and.(F%dBfield).and. &
         .not.(params%field_model(10:14).eq.'3DBdB').and..not.params%field_model.eq.'M3D_C1') then
     do it=params%ito,params%t_steps,params%t_skip
        call adv_GCinterp_3DBdB1_top(params,spp,P,F)
        
        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip,rp)*params%dt        
        params%it = it-1_ip+params%t_skip

        call save_simulation_outputs(params,spp,F)
        call save_restart_variables(params,spp,F)
     end do
  end if

  if (params%orbit_model(1:2).eq.'GC'.and.params%field_eval.eq.'interp'.and. &
       .not.(F%axisymmetric_fields).and..not.(F%dBfield).and..not.params%field_model.eq.'M3D_C1') then
     do it=params%ito,params%t_steps,params%t_skip
        call adv_GCinterp_B_top(params,spp,P,F)
        
        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip,rp)*params%dt        
        params%it = it-1_ip+params%t_skip

        call save_simulation_outputs(params,spp,F)
        call save_restart_variables(params,spp,F)
     end do
  end if
#endif

#ifdef FIO
  if (params%orbit_model(1:2).eq.'GC'.and.params%field_model.eq.'M3D_C1'.and. &
       .not.F%ReInterp_2x1t) then
     do it=params%ito,params%t_steps,params%t_skip
        call adv_GCinterp_fio_top(params,spp,P,F)
        
        params%time = params%init_time &
             +REAL(it-1_ip+params%t_skip*params%t_it_SC,rp)*params%dt        
        params%it = it-1_ip+params%t_skip*params%t_it_SC

        call save_simulation_outputs(params,spp,F)
        call save_restart_variables(params,spp,F)

        if (params%mpi_params%rank .EQ. 0) then
           flush(output_unit_write)
        end if
        
     end do
     
  end if

  if (params%orbit_model(1:2).eq.'GC'.and.params%field_model.eq.'M3D_C1'.and. &
       F%ReInterp_2x1t) then

     do it=F%ind0_2x1t,params%time_slice

!        write(6,*) it,F%ind0_2x1t
        
        if (it.gt.F%ind0_2x1t) then
           call initialize_m3d_c1(params, F, P, spp,.false.)
           if (params%collisions.or.params%radiation) then
              call initialize_m3d_c1_imp(params,F,P, &
                   params%num_impurity_species,.false.)
           end if
        end if
        
        if (params%mpi_params%rank .EQ. 0) then
           write(output_unit_write,*) 'tskip',params%t_skip
           flush(output_unit_write)
        end if
           
        call adv_GCinterp_fio_top(params,spp,P,F)
               
        params%it = params%it+params%t_skip
        params%time = params%init_time &
             +REAL(params%it,rp)*params%dt 
        
        call save_simulation_outputs(params,spp,F)
        call save_restart_variables(params,spp,F)

        !comment out for debugging only
        F%ind_2x1t=F%ind_2x1t+1_ip


        if (params%mpi_params%rank .EQ. 0) then
           write(output_unit_write,*) 'KORC time ',params%time*params%cpp%time
           flush(output_unit_write)
        end if
              
     end do
  end if
#endif
  
  call timing_KORC(params)

  ! * * * FINALIZING SIMULATION * * * 
  call finalize_HDF5()

#ifdef PSPLINE
  call finalize_interpolants(params)
#endif
  
#ifdef FIO
  if (TRIM(params%field_model) .eq. 'M3D_C1'.or. &
      TRIM(params%field_model) .eq. 'NIMROD') then
     call finalize_FIO(params,F,P)
  end if
#endif
  
  ! DEALLOCATION OF VARIABLES
  call deallocate_variables(params,F,P,spp)

  
  call deallocate_collisions_params(params)

  
  call finalize_communications(params)
  ! * * * FINALIZING SIMULATION * * *

  if (params%mpi_params%rank .EQ. 0) then
     write(output_unit_write,'("KORC ran successfully!")')
     close(output_unit_write)
  end if
  
end program main

