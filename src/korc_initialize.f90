module korc_initialize
  !! @note Module with subroutines to load simulation parameters 
  !! and to define the time step in the simulation.@endnote
  use korc_types
  use korc_constants
  use korc_hpc
  use korc_HDF5
  use korc_fields
  use korc_rnd_numbers
  use korc_spatial_distribution
  use korc_velocity_distribution
  use korc_coords
  use korc_input

  IMPLICIT NONE


  PRIVATE :: set_paths,&
       load_korc_params
  PUBLIC :: initialize_korc_parameters,&
       initialize_particles,&
       define_time_step

CONTAINS

  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !
  ! ** SUBROUTINES FOR INITIALIZING KORC PARAMETERS ** !
  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !


  subroutine load_korc_params(params)
    !! @note Subroutine that loads the simulation parameters from the 
    !! file specified in params\%path_to_inputs @endnote
    TYPE (KORC_PARAMS), INTENT(INOUT) 	:: params
    !! Core KORC simulation parameters.
    !LOGICAL 				:: restart
    !! Flag to indicate if the simulations restarts (restart=T) or not
    !! (restart=F). Restart simulation that exited before simulation_time
    !! reached.
    !LOGICAL 				:: proceed
    !! Flag to indicate if the simulations proceeds (proceed=T) or not
    !! (proceed=F). Append simulation results after previous simulation_time
    !! reached.
    !LOGICAL  :: reinit
    !! Flag to begin a new simulation, reinitializing from restart file state
    !REAL(rp) 				:: simulation_time
    !! Total simulation time in seconds.
    !REAL(rp) 				:: snapshot_frequency
    !! Time between snapshots in time of the simulation.
    !REAL(rp) 				:: restart_overwrite_frequency
    !! Time between overwrites of restart file in time of the simulation.
    !REAL(rp) 				:: dt
    !! Time step in the simulation as a fraction of the relativistic 
    !! electron gyro-period @f$\tau_e = 2\pi\gamma m_e/eB_0@f$
    !REAL(rp) 				:: minimum_particle_energy
    !! Minimum allowed relativistic factor @f$\gamma@f$ of simulated electrons.
    !LOGICAL 				:: radiation
    !! Flag to indicate if synchrotron radiation losses are included
    !! (radiation=T) or not (radiation=F).
    !LOGICAL 				:: collisions
    !! Flag to indicate if collisionsare included (collisions=T) or not
    !! (collisions=F).
    !CHARACTER(MAX_STRING_LENGTH) 		:: GC_rad_model
    !CHARACTER(MAX_STRING_LENGTH) 		:: collisions_model
    !! String with the name of the collisions model to be used in the simulation.
    !CHARACTER(MAX_STRING_LENGTH) 		:: bound_electron_model
    !CHARACTER(MAX_STRING_LENGTH) 		:: profile_model
    !! String with the name of the model for the plasma profiles.
    !CHARACTER(MAX_STRING_LENGTH) 		:: field_model
    !! String with the name of the model for the field profiles.
    !CHARACTER(MAX_STRING_LENGTH) 		:: magnetic_field_filename
    !! String with the name of the model for the fields and plasma profiles.
    !CHARACTER(MAX_STRING_LENGTH) 		:: outputs_list
    !! List of electron variables to include in the outputs.
    !INTEGER 				:: num_species
    !! Number of different populations of simulated relativistic electrons
    !! in KORC.
    INTEGER 				:: imax
    !! Auxiliary variable used to parse the output_list
    INTEGER 				:: imin
    !! Auxiliary variable used to parse the output_list
    INTEGER 				:: ii
    !! Iterator used to parse the output_list
    INTEGER 				:: jj
    !! Iterator used to parse the output_list
    INTEGER 				:: num_outputs
    !! Auxiliary variable used to parse the output_list
    INTEGER, DIMENSION(2) 		:: indices
    !! Auxiliary variable used to parse the output_list
    !LOGICAL 				:: HDF5_error_handling
    !! Flag for HDF5 error handling
    !LOGICAL 		:: FO_GC_compare
    !CHARACTER(MAX_STRING_LENGTH) 		:: orbit_model
    !! String with the name of the orbit model ('FO' or 'GC').
    !CHARACTER(MAX_STRING_LENGTH) :: field_eval
    !! String with the name of the field evaluation method for
    !! analytical fields ('interp' or 'eqn')
    !LOGICAL 				:: FokPlan
    !! Flag to decouple spatial-dependence of evolution
    !LOGICAL :: SameRandSeed
    !LOGICAL :: SC_E
    !LOGICAL :: SC_E_add
    !INTEGER                           :: time_slice
    !REAL(rp)                          :: rmax,rmin,zmax,zmin
    !INTEGER                           :: pchunk

    !NAMELIST /input_parameters/ restart,field_model,magnetic_field_filename, &
    !     simulation_time,snapshot_frequency,dt,num_species,radiation, &
    !     collisions,collisions_model,outputs_list,minimum_particle_energy, &
    !     HDF5_error_handling,orbit_model,field_eval,proceed,profile_model, &
    !     restart_overwrite_frequency,FokPlan,GC_rad_model, &
    !     bound_electron_model,FO_GC_compare,SameRandSeed,SC_E,reinit, &
    !     SC_E_add,time_slice,rmax,rmin,zmax,zmin,pchunk

    !open(unit=default_unit_open,file=TRIM(params%path_to_inputs), &
    !     status='OLD',form='formatted')
    !read(default_unit_open,nml=input_parameters)
    !close(default_unit_open)

    params%restart = restart
    params%proceed = proceed
    params%reinit  = reinit

    params%load_balance = load_balance
    
    params%simulation_time = simulation_time
    params%snapshot_frequency = snapshot_frequency
    params%restart_overwrite_frequency=restart_overwrite_frequency
    params%dt = dt

    params%num_species = num_species
    params%profile_model = TRIM(profile_model)
    params%field_model = TRIM(field_model)
    params%magnetic_field_filename = TRIM(magnetic_field_filename)
    params%time_slice = time_slice
    params%rmax = rmax
    params%rmin = rmin
    params%zmax = zmax
    params%zmin = zmin
    params%minimum_particle_energy = minimum_particle_energy*C_E
    params%minimum_particle_g = 1.0_rp + params%minimum_particle_energy/ &
         (C_ME*C_C**2) ! Minimum value of relativistic gamma factor
    params%radiation = radiation
    params%collisions = collisions
    params%LargeCollisions = LargeCollisions
    params%collisions_model = TRIM(collisions_model)
    params%bound_electron_model = TRIM(bound_electron_model)
    params%GC_rad_model = TRIM(GC_rad_model)

    if (HDF5_error_handling) then
       params%HDF5_error_handling = 1_idef
    else
       params%HDF5_error_handling = 0_idef
    end if

    params%orbit_model = orbit_model
    params%FO_GC_compare = FO_GC_compare
    params%field_eval = field_eval

    params%GC_coords=.FALSE.

    params%FokPlan=FokPlan

    params%SameRandSeed = SameRandSeed

    params%SC_E=SC_E
    params%SC_E_add=SC_E_add

    params%pchunk=pchunk

    ! Loading list of output parameters (parsing)
    imin = SCAN(outputs_list,'{')
    imax = SCAN(outputs_list,'}')

    ii = 1_idef
    jj = 1_idef
    num_outputs = 1_idef
    do while (ii.NE.0)
       ii = SCAN(outputs_list(jj:),",")
       if (ii.NE.0) then
          jj = jj + ii
          num_outputs = num_outputs + 1_idef
       end if
    end do

    ALLOCATE(params%outputs_list(num_outputs))

    if (num_outputs.GT.1_idef) then
       indices = 0_idef
       indices(2) = SCAN(outputs_list,",")
       params%outputs_list(1) = TRIM(outputs_list(imin+1_idef:indices(2)-1_idef))
       indices(1) = indices(1) + indices(2) + 1_idef
       do ii=2_idef,num_outputs
          indices(2) = SCAN(outputs_list(indices(1):),",")
          if (indices(2).EQ.0_idef) then
             params%outputs_list(ii) = TRIM(outputs_list(indices(1):imax-1_idef))
          else
             params%outputs_list(ii) = TRIM(outputs_list(indices(1):indices(1)+indices(2)-2_idef))
             indices(1) = indices(1) + indices(2)
          end if
       end do
    else
       params%outputs_list(1) = TRIM(outputs_list(imin+1_idef:imax-1_idef))
    end if

    if (params%mpi_params%rank .EQ. 0) then
       write(output_unit_write,'(/,"* * * * * SIMULATION PARAMETERS * * * * *")')
       write(output_unit_write,'("Restarting simulation: ",L1)') params%restart
       write(output_unit_write,'("Continuing simulation: ",L1)') params%proceed
       write(output_unit_write,'("Number of electron populations: ",I16)') params%num_species
       write(output_unit_write,*) 'Orbit model: ',TRIM(params%orbit_model)
       write(output_unit_write,*) 'Magnetic field model: ',TRIM(params%field_model)
       write(output_unit_write,*) 'Magnetic field evaluation: ',TRIM(params%field_eval)
       if (TRIM(params%field_model).EQ.'EXTERNAL') then
          write(output_unit_write,*) 'Magnetic field file: ',TRIM(params%magnetic_field_filename)
       end if

       write(output_unit_write,'("Radiation losses included: ",L1)') params%radiation
       if (params%radiation.and.(params%orbit_model(1:2).eq.'GC')) then
          write(output_unit_write,*) 'Radiation model: ',TRIM(params%GC_rad_model)
       end if
       write(output_unit_write,'("Collisions losses included: ",L1)') params%collisions
       if (params%collisions) then
          write(output_unit_write,*) 'Collision model: ',TRIM(params%collisions_model)
          write(output_unit_write,*) &
               'Bound electron model: ',TRIM(params%bound_electron_model)
       end if
       write(output_unit_write,'("Self-consistent E included: ",L1)') params%SC_E
       write(output_unit_write,'("* * * * * * * * * * * * * * * * * * * * *",/)')
    end if
  end subroutine load_korc_params

  subroutine initialize_korc_parameters(params)
    !! @note Interface for calling initialization subroutines @endnote
    TYPE(KORC_PARAMS), INTENT(INOUT) 	:: params
    !! Core KORC simulation parameters.
    INTEGER 							:: mpierr
    !! MPI error status.

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

    call read_namelist(params,params%path_to_inputs,.true.,params%path_to_outputs)

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

    call load_korc_params(params)

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  end subroutine initialize_korc_parameters

  subroutine define_time_step(params,F)
    !! @note Subroutine that defines or loads from restart file the time
    !! stepping parameters. @endnote
    TYPE(KORC_PARAMS), INTENT(INOUT) :: params
    TYPE(FIELDS), INTENT(INOUT) :: F
    !! Core KORC simulation parameters.

    if (params%restart) then
       call load_time_stepping_params(params)
       
    else if (params%proceed.or.params%reinit) then
       
       call load_prev_time(params)
       
       params%ito = 1_ip

       params%dt = params%dt*(2.0_rp*C_PI*params%cpp%time_r)

       params%t_steps = CEILING((params%simulation_time-params%init_time)/ &
            params%dt,ip)

       params%output_cadence = CEILING(params%snapshot_frequency/params%dt,ip)

       if (params%output_cadence.EQ.0_ip) params%output_cadence = 1_ip

       params%num_snapshots = params%t_steps/params%output_cadence

       if (params%t_steps.gt.params%output_cadence) then
          params%dt=params%snapshot_frequency/float(params%output_cadence)
       endif
       
       params%restart_output_cadence = CEILING(params%restart_overwrite_frequency/ &
            params%dt,ip)


       params%t_skip=min(params%t_steps,params%output_cadence)
       params%t_skip=max(1_ip,params%t_skip)


    else
       params%ito = 1_ip

       params%dt = params%dt*(2.0_rp*C_PI*params%cpp%time_r)

       params%t_steps = CEILING(params%simulation_time/params%dt,ip)
       
       params%output_cadence = CEILING(params%snapshot_frequency/params%dt,ip)

       if (params%output_cadence.EQ.0_ip) params%output_cadence = 1_ip

       params%num_snapshots = params%t_steps/params%output_cadence

       if (params%t_steps.gt.params%output_cadence) then
          params%dt=params%snapshot_frequency/float(params%output_cadence)
       endif

       params%restart_output_cadence = CEILING(params%restart_overwrite_frequency/ &
            params%dt,ip)

       params%t_skip=min(params%t_steps,params%output_cadence)
       params%t_skip=max(1_ip,params%t_skip)

    end if

    !    write(output_unit_write,*) 'dt',params%dt,'t_skip',params%t_skip

    if (params%mpi_params%rank .EQ. 0) then
       write(output_unit_write,'(/,"* * * * * TIME STEPPING PARAMETERS * * * * *")')
       write(output_unit_write,'("Simulation time: ",E17.10," s")') params%simulation_time
       write(output_unit_write,'("Initial time: ",E17.10," s")') params%init_time     
       write(output_unit_write,'("Output frequency: ",E17.10," s")') params%snapshot_frequency
       write(output_unit_write,'("Relativistic gyro-period: ",E17.10)') 2.0_rp*C_PI* &
            params%cpp%time_r
       write(output_unit_write,'("Time step: ",E17.10)') params%dt
       write(output_unit_write,'("Number of time steps: ",I16)') params%t_steps
       write(output_unit_write,'("Starting simulation at time step: ",I16)') params%ito
       write(output_unit_write,'("Output cadence: ",I16)') params%output_cadence
       write(output_unit_write,'("Restart cadence: ",I16)') params%restart_output_cadence
       write(output_unit_write,'("Number of outputs: ",I16)') params%num_snapshots
       write(output_unit_write,'("* * * * * * * * * * * * * * * * * * * * * * *",/)')
    end if
  end subroutine define_time_step

  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !
  ! * * * SUBROUTINES FOR INITIALIZING PARTICLES * * * !
  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !

  subroutine initialize_particles(params,F,P,spp)
    !! @note Subroutine that loads the information of the initial condition 
    !! of the different particle species. This subroutine calls
    !! the subroutine that generates the initial energy and pitch angle 
    !! distribution functions. @endnote
    TYPE(KORC_PARAMS), INTENT(IN) 				:: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN) 					:: F
    !! An instance of KORC's derived type FIELDS containing all the information 
    !! about the fields used in the simulation. See [[korc_types]]
    !!and [[korc_fields]].
    TYPE(PROFILES), INTENT(INOUT) 			:: P
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(OUT) 	:: spp
    !! An instance of KORC's derived type SPECIES containing all the information 
    !! of different electron species. See [[korc_types]].
    !REAL(rp), DIMENSION(:), ALLOCATABLE 				:: ppp
    !! Number of computational particles per MPI process 
    !! used to simulate each electron species.
    !REAL(rp), DIMENSION(:), ALLOCATABLE 				:: q
    !! Charge of each species.
    !REAL(rp), DIMENSION(:), ALLOCATABLE 				:: m
    !! Mass of each species.
    !REAL(rp), DIMENSION(:), ALLOCATABLE 				:: Eo
    !! Initial energy of each electron species in case of 
    !! using an initial mono-energetic distribution.
    !REAL(rp), DIMENSION(:), ALLOCATABLE 				:: etao
    !! Initial pitch-angle of each electron species in case of 
    !! using an initial mono-pitch-angle distribution.
    !REAL(rp), DIMENSION(:), ALLOCATABLE 				:: Eo_lims
    !! Minimum and maximum energy limits of a given initial
    !! non-mono-energetic distribution.
    !REAL(rp), DIMENSION(:), ALLOCATABLE 				:: etao_lims
    !! Minimum and maximum pitch-angle limits of a given initial
    !! non-mono-pitch-angle distribution.
    !LOGICAL, DIMENSION(:), ALLOCATABLE 				:: runaway
    !! Flag to decide whether a given electron is a runaway (runaway=T)
    !! or not (runaway=F).
    !CHARACTER(MAX_STRING_LENGTH),DIMENSION(:),ALLOCATABLE  :: spatial_distribution
    !! String describing the type of initial spatial distribution for
    !! each electron species.
    !CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: energy_distribution
    !! String describing the type of initial energy distribution for each
    !! electron species.
    !CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: pitch_distribution
    !! String describing the type of initial pitch-angle distribution
    !! for each electron species.
    !REAL(rp), DIMENSION(:), ALLOCATABLE 				:: Ro
    !! Radial position of the center of the electrons' initial
    !! spatial distribution.
    !REAL(rp), DIMENSION(:), ALLOCATABLE 				:: PHIo
    !! Azimuthal position of the electrons' initial spatial distribution, 
    !! in case of using a disk at a certain poloidal section.
    !REAL(rp), DIMENSION(:), ALLOCATABLE 				:: Zo
    !! Height of the center of the electrons' initial spatial distribution.
    !REAL(rp), DIMENSION(:), ALLOCATABLE 				:: r_inner
    !! Minimum minor radius of the electrons' initial spatial distribution.
    !REAL(rp), DIMENSION(:), ALLOCATABLE 				:: r_outter
    !! Maximum minor radius of the electrons' initial spatial distribution.
    !REAL(rp), DIMENSION(:), ALLOCATABLE 				:: falloff_rate
    !! Exponential falloff or standard deviation of a non-uniform radial
    !! distribution of electrons.
    !REAL(rp), DIMENSION(:), ALLOCATABLE 				:: shear_factor
    !! Shear factor used to generate an initial spatial 
    !! distribution with an elliptic poloidal cross section.
    !! See <em>Carbajal and del-Castillo-Negrete, Nuclear Fusion,
    !! submitted (2018)</em>.
    !REAL(rp), DIMENSION(:), ALLOCATABLE                           :: sigmaR
    !! Variance of the first dimension of a 2D Gaussian, spatial
    !! distribution function
    !REAL(rp), DIMENSION(:), ALLOCATABLE                           :: sigmaZ
    !! Variance of the second dimension of a 2D Gaussian, spatial
    !! distribution function
    !REAL(rp), DIMENSION(:), ALLOCATABLE                          :: theta_gauss
    !! Angle of counter-clockwise rotation (in degrees) of 2D Gaussian
    !! distribution relative to R,Z
    !REAL(rp), DIMENSION(:), ALLOCATABLE                          :: psi_max
    !! Maximum value of the argument of the 2D gaussian exponential, used for an
    !! indicator function that limits the region of MH sampling
    !REAL(rp), DIMENSION(:), ALLOCATABLE                          :: Spong_b
    !REAL(rp), DIMENSION(:), ALLOCATABLE                          :: Spong_w
    !REAL(rp), DIMENSION(:), ALLOCATABLE                          :: Spong_dlam
    !REAL(rp), DIMENSION(:), ALLOCATABLE                          :: dth,dgam
    !REAL(rp), DIMENSION(:), ALLOCATABLE                          :: dR
    !REAL(rp), DIMENSION(:), ALLOCATABLE                          :: dZ

    INTEGER 				                       	:: ii
    !! Iterator of spp structure.
    INTEGER 							:: mpierr
    !! MPI error status.
    !REAL(rp), DIMENSION(:), ALLOCATABLE :: Xtrace

    !NAMELIST /plasma_species/ ppp,q,m,Eo,etao,Eo_lims,etao_lims,runaway, &
    !     spatial_distribution,energy_distribution,pitch_distribution,Ro, &
    !     PHIo,Zo,r_inner,r_outter,falloff_rate,shear_factor,sigmaR,sigmaZ, &
    !     theta_gauss,psi_max,Xtrace,Spong_b,Spong_w,Spong_dlam,dth,dR,dZ,dgam

    ! Allocate array containing variables of particles for each species
    ALLOCATE(spp(params%num_species))

    !ALLOCATE(ppp(params%num_species))
    !ALLOCATE(q(params%num_species))
    !ALLOCATE(m(params%num_species))
    !ALLOCATE(Eo(params%num_species))
    !ALLOCATE(etao(params%num_species))
    !ALLOCATE(Eo_lims(2_idef*params%num_species))
    !ALLOCATE(etao_lims(2_idef*params%num_species))
    !ALLOCATE(runaway(params%num_species))
    !ALLOCATE(spatial_distribution(params%num_species))
    !ALLOCATE(energy_distribution(params%num_species))
    !ALLOCATE(pitch_distribution(params%num_species))
    !ALLOCATE(Ro(params%num_species))
    !ALLOCATE(PHIo(params%num_species))
    !ALLOCATE(Zo(params%num_species))
    !ALLOCATE(r_inner(params%num_species))
    !ALLOCATE(r_outter(params%num_species))
    !ALLOCATE(falloff_rate(params%num_species))
    !ALLOCATE(shear_factor(params%num_species))
    !ALLOCATE(sigmaR(params%num_species))
    !ALLOCATE(sigmaZ(params%num_species))
    !ALLOCATE(theta_gauss(params%num_species))
    !ALLOCATE(psi_max(params%num_species))
    !ALLOCATE(Spong_b(params%num_species))
    !ALLOCATE(Spong_w(params%num_species))
    !ALLOCATE(Spong_dlam(params%num_species))
    !ALLOCATE(dth(params%num_species))
    !ALLOCATE(dgam(params%num_species))
    !ALLOCATE(dR(params%num_species))
    !ALLOCATE(dZ(params%num_species))
    !ALLOCATE(Xtrace(3_idef*params%num_species))

    !open(unit=default_unit_open,file=TRIM(params%path_to_inputs), &
    !     status='OLD',form='formatted')
    !read(default_unit_open,nml=plasma_species)
    !close(default_unit_open)


    do ii=1_idef,params%num_species
       spp(ii)%runaway = runaway(ii)
       spp(ii)%spatial_distribution = TRIM(spatial_distribution(ii))
       spp(ii)%energy_distribution = TRIM(energy_distribution(ii))
       spp(ii)%pitch_distribution = TRIM(pitch_distribution(ii))
       spp(ii)%q = q(ii)*C_E
       spp(ii)%m = m(ii)*C_ME
       spp(ii)%ppp = ppp(ii)
       spp(ii)%pinit = pinit(ii)
       spp(ii)%pRE = pinit(ii)

       spp(ii)%Ro = Ro(ii)
       spp(ii)%PHIo = C_PI*PHIo(ii)/180.0_rp
       spp(ii)%Zo = Zo(ii)
       spp(ii)%r_inner = r_inner(ii)
       spp(ii)%r_outter = r_outter(ii)
       spp(ii)%falloff_rate = falloff_rate(ii)
       spp(ii)%shear_factor = shear_factor(ii)
       spp(ii)%sigmaR = sigmaR(ii)
       spp(ii)%sigmaZ = sigmaZ(ii)
       spp(ii)%theta_gauss = theta_gauss(ii)
       spp(ii)%psi_max = psi_max(ii)
       spp(ii)%Spong_w = Spong_w(ii)
       spp(ii)%Spong_b = Spong_b(ii)
       spp(ii)%Spong_dlam = Spong_dlam(ii)
       spp(ii)%dth = dth(ii)
       spp(ii)%dgam = dgam(ii)
       spp(ii)%dR = dR(ii)
       spp(ii)%dZ = dZ(ii)


       ! * * These values can change in initial_energy_pitch_dist * * !
       spp(ii)%Eo = Eno(ii)*C_E
       spp(ii)%Eo_lims = Eo_lims((ii-1_idef)*2_idef + 1_idef:2_idef*ii)*C_E
       spp(ii)%etao = etao(ii)
       spp(ii)%etao_lims = etao_lims((ii-1_idef)*2_idef + 1_idef:2_idef*ii)
       ! * * These values can change in initial_energy_pitch_dist * * !

       if (spp(ii)%spatial_distribution.eq.'TRACER') &
            spp(ii)%Xtrace = Xtrace((ii-1_idef)*3_idef + 1_idef:3_idef*ii)

       ALLOCATE( spp(ii)%vars%X(spp(ii)%ppp,3) )
       ALLOCATE( spp(ii)%vars%V(spp(ii)%ppp,3) )
       ALLOCATE( spp(ii)%vars%Rgc(spp(ii)%ppp,3) )
       ALLOCATE( spp(ii)%vars%Y(spp(ii)%ppp,3) )
       ALLOCATE( spp(ii)%vars%E(spp(ii)%ppp,3) )
       ALLOCATE( spp(ii)%vars%B(spp(ii)%ppp,3) )
       ALLOCATE( spp(ii)%vars%PSI_P(spp(ii)%ppp) )
       ALLOCATE( spp(ii)%vars%ne(spp(ii)%ppp) )
       ALLOCATE( spp(ii)%vars%ni(spp(ii)%ppp) )
       ALLOCATE( spp(ii)%vars%Te(spp(ii)%ppp) )
       ALLOCATE( spp(ii)%vars%Zeff(spp(ii)%ppp) )
       ALLOCATE( spp(ii)%vars%g(spp(ii)%ppp) )
       ALLOCATE( spp(ii)%vars%eta(spp(ii)%ppp) )
       ALLOCATE( spp(ii)%vars%mu(spp(ii)%ppp) )
       ALLOCATE( spp(ii)%vars%Prad(spp(ii)%ppp) )
       ALLOCATE( spp(ii)%vars%Pin(spp(ii)%ppp) )
       ALLOCATE( spp(ii)%vars%flagCon(spp(ii)%ppp) )       
       ALLOCATE( spp(ii)%vars%flagCol(spp(ii)%ppp) )
       ALLOCATE( spp(ii)%vars%initLCFS(spp(ii)%ppp) )
       ALLOCATE( spp(ii)%vars%flagRE(spp(ii)%ppp) )
       ALLOCATE( spp(ii)%vars%wt(spp(ii)%ppp) )
#ifdef FIO
       ALLOCATE( spp(ii)%vars%hint(spp(ii)%ppp))
#endif
       
       !     write(output_unit_write,'("0 size of PSI_P: ",I16)') size(spp(ii)%vars%PSI_P)

       spp(ii)%vars%X = 0.0_rp
       spp(ii)%vars%V = 0.0_rp
       spp(ii)%vars%Rgc = 0.0_rp
       spp(ii)%vars%Y = 0.0_rp
       spp(ii)%vars%E = 0.0_rp
       spp(ii)%vars%B = 0.0_rp
       spp(ii)%vars%PSI_P = 0.0_rp
       spp(ii)%vars%ne = 0.0_rp
       spp(ii)%vars%ni = 0.0_rp
       spp(ii)%vars%Te = 0.0_rp
       spp(ii)%vars%Zeff = 0.0_rp
       spp(ii)%vars%g = 0.0_rp
       spp(ii)%vars%eta = 0.0_rp
       spp(ii)%vars%mu = 0.0_rp
       spp(ii)%vars%Prad = 0.0_rp
       spp(ii)%vars%Pin = 0.0_rp
       spp(ii)%vars%flagCon(1:spp(ii)%pinit) = 1_is
       spp(ii)%vars%flagCol(1:spp(ii)%pinit) = 1_is
       spp(ii)%vars%flagRE(1:spp(ii)%pinit) = 1_is
       spp(ii)%vars%initLCFS(1:spp(ii)%pinit) = 1_is
       if (spp(ii)%pinit.lt.spp(ii)%ppp) then
          spp(ii)%vars%flagCon(spp(ii)%pinit+1:spp(ii)%ppp) = 0_is
          spp(ii)%vars%flagCol(spp(ii)%pinit+1:spp(ii)%ppp) = 0_is          
          spp(ii)%vars%flagRE(spp(ii)%pinit+1:spp(ii)%ppp) = 0_is
          spp(ii)%vars%initLCFS(spp(ii)%pinit+1:spp(ii)%ppp) = 0_is
       endif
       spp(ii)%vars%wt = 0.0_rp

       if (params%orbit_model(1:2).eq.'GC') then
          ALLOCATE( spp(ii)%vars%Y0(spp(ii)%ppp,3) )
          ALLOCATE( spp(ii)%vars%V0(spp(ii)%ppp) )
          ALLOCATE( spp(ii)%vars%k1(spp(ii)%ppp,4) )
          ALLOCATE( spp(ii)%vars%k2(spp(ii)%ppp,4) )
          ALLOCATE( spp(ii)%vars%k3(spp(ii)%ppp,4) )
          ALLOCATE( spp(ii)%vars%k4(spp(ii)%ppp,4) )
          ALLOCATE( spp(ii)%vars%k5(spp(ii)%ppp,4) )
          ALLOCATE( spp(ii)%vars%k6(spp(ii)%ppp,4) )
          if (params%orbit_model(3:5)=='pre'.or. &
               TRIM(params%field_model)=='M3D_C1'.or. &
               TRIM(params%field_model)=='NIMROD') then
             ALLOCATE( spp(ii)%vars%gradB(spp(ii)%ppp,3) )
             ALLOCATE( spp(ii)%vars%curlb(spp(ii)%ppp,3) )

             spp(ii)%vars%gradB = 0.0_rp
             spp(ii)%vars%curlb = 0.0_rp
          else if (params%orbit_model(3:6)=='grad') then
             ALLOCATE( spp(ii)%vars%BR(spp(ii)%ppp,3) )
             ALLOCATE( spp(ii)%vars%BPHI(spp(ii)%ppp,3) )
             ALLOCATE( spp(ii)%vars%BZ(spp(ii)%ppp,3) )

             spp(ii)%vars%BR = 0.0_rp
             spp(ii)%vars%BPHI = 0.0_rp
             spp(ii)%vars%BZ = 0.0_rp
          end if
          ALLOCATE( spp(ii)%vars%RHS(spp(ii)%ppp,5) )

          spp(ii)%vars%Y0 = 0.0_rp
          spp(ii)%vars%V0 = 0.0_rp
          spp(ii)%vars%k1 = 0.0_rp
          spp(ii)%vars%k2 = 0.0_rp
          spp(ii)%vars%k3 = 0.0_rp
          spp(ii)%vars%k4 = 0.0_rp
          spp(ii)%vars%k5 = 0.0_rp
          spp(ii)%vars%k6 = 0.0_rp
          spp(ii)%vars%RHS = 0.0_rp
       end if

    end do

    P%R0_RE=spp(1)%Ro
    P%Z0_RE=spp(1)%Zo
    P%n_REr0=max(sqrt(spp(1)%psi_max*2*spp(1)%sigmaR**2), &
         sqrt(spp(1)%psi_max*2*spp(1)%sigmaZ**2))  

    call initial_energy_pitch_dist(params,spp)


    DEALLOCATE(ppp)
    DEALLOCATE(pinit)
    DEALLOCATE(q)
    DEALLOCATE(m)
    DEALLOCATE(Eno)
    DEALLOCATE(etao)
    DEALLOCATE(Eo_lims)
    DEALLOCATE(etao_lims)
    DEALLOCATE(runaway)
    DEALLOCATE(spatial_distribution)
    DEALLOCATE(energy_distribution)
    DEALLOCATE(pitch_distribution)
    DEALLOCATE(Ro)
    DEALLOCATE(PHIo)
    DEALLOCATE(Zo)
    DEALLOCATE(r_inner)
    DEALLOCATE(r_outter)
    DEALLOCATE(falloff_rate)
    DEALLOCATE(shear_factor)
    DEALLOCATE(sigmaR)
    DEALLOCATE(sigmaZ)
    DEALLOCATE(theta_gauss)
    DEALLOCATE(psi_max)
    DEALLOCATE(Spong_b)
    DEALLOCATE(Spong_w)
    DEALLOCATE(Spong_dlam)
    DEALLOCATE(dth)
    DEALLOCATE(dgam)
    DEALLOCATE(dR)
    DEALLOCATE(dZ)
    DEALLOCATE(Xtrace)

  end subroutine initialize_particles


  subroutine set_up_particles_ic(params,F,spp,P)
    !! @note Subroutine with calls to subroutines to load particles' 
    !! information if it is a restarting simulation, or to initialize the
    !! spatial and velocity distribution of each species if it is a new  
    !! simulation. @endnote
    TYPE(KORC_PARAMS), INTENT(INOUT) 				:: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(INOUT) 					:: F
    !! An instance of KORC's derived type FIELDS containing all 
    !! the information about the fields used in the simulation. 
    !! See [[korc_types]] and [[korc_fields]].
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)       :: spp
    !! An instance of KORC's derived type SPECIES containing all 
    !! the information of different electron species. See [[korc_types]].
    TYPE(PROFILES), INTENT(IN)                                 :: P
    !! An instance of the KORC derived type PROFILES.
    INTEGER                                                    :: ii
    !! Species iterator.

    if (params%restart.OR.params%proceed.or.params%reinit) then
       call load_particles_ic(params,spp,F)

       if(params%LargeCollisions.and.(.not.params%load_balance)) then
          do ii=1_idef,params%num_species
             spp(ii)%pRE=int(sum(float(spp(ii)%vars%flagRE)))
          end do
       end if

       !write(6,*) 'flagRE',spp(1)%vars%flagRE
       !write(6,*) 'pRE',spp(1)%pRE
       
       call init_random_seed()              
    else

       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'("* * * * INITIALIZING SPATIAL DISTRIBUTION * * * *")')
          flush(output_unit_write)
       end if
       call intitial_spatial_distribution(params,spp,P,F)
       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'("* * * * * * * * * * * * * * * * * * * * * * * *",/)')
       end if

       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'("* * * * INITIALIZING VELOCITY COMPONENTS * * * *")')
       end if
       call initial_gyro_distribution(params,F,spp)
       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'("* * * * * * * * * * * * * * * * * * * * * * * *",/)')
       end if
    end if

  end subroutine set_up_particles_ic

end module korc_initialize
