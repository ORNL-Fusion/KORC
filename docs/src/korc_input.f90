module korc_input
  !! @note Module with subroutines to read in all namelists in supplied
  !! input file and broadcast to all mpi processes.
  USE korc_types
  USE korc_hpc
  
  IMPLICIT NONE


  !! Default values for all inputs
  !! -----------------------------------------------
  !! input_parameters
  !! -----------------------------------------------
  LOGICAL :: restart = .FALSE.
    ! Restart simulation that exited before simulation_time reached
  LOGICAL :: proceed = .FALSE.
  ! Append simulation results after previous simulation_time reached
  LOGICAL :: load_balance = .FALSE.
  LOGICAL :: reinit = .FALSE.
    ! Begin a new simulation, reinitializing from restart file state
  REAL(rp) :: simulation_time = 1.E-3
    ! Total aimed simulation time in seconds   
    ! Run 10 mu s If transients exist put 5 extra mu s.
  REAL(rp) :: snapshot_frequency = 1.E-5
    ! Time between snapshots in seconds
  REAL(rp) :: restart_overwrite_frequency = 1.E-1	
    ! Time between overwritting of restart file in seconds
  REAL(rp) :: dt = 1.E0
    ! Time step as fraction of relativistic gyro-period
  INTEGER :: num_species = 1
  REAL(rp) :: minimum_particle_energy = 1.0E5 
    ! Minimum energy of simulated particles in eV
  LOGICAL :: radiation = .FALSE.
  CHARACTER(30) :: GC_rad_model='SDE'
  LOGICAL :: collisions = .FALSE.
  LOGICAL :: LargeCollisions = .FALSE.
  CHARACTER(30) :: collisions_model = 'SINGLE_SPECIES' 
    ! Options are: 'NONE','SINGLE_SPECIES' and 'MULTIPLE_SPECIES'
  CHARACTER(30) :: bound_electron_model = 'HESSLOW' 
    ! Options are: 'NO_BOUND', 'HESSLOW', and 'ROSENBLUTH'
  CHARACTER(30) :: field_model = 'M3D_C1'
  CHARACTER(30) :: profile_model = 'M3D_C1'	
    ! The two options for this parameter are 'ANALYTICAL' or 'EXTERNAL'.
    ! For 'ANALYTICAL', the magnetic field is calculated based on
    ! the parameters given in the "analytic_mag_field_params" section.
    ! For 'EXTERNAL', the magnetic field is loaded from the file
    ! specified in "filename".
    ! 'UNIFORM' A uniform magnetic field used to advance only electrons' 
    ! velocity.
  CHARACTER(150) :: magnetic_field_filename = 'C1.h5'
  !  magnetic_field_filename = 'JFIT_D3D_164409_1405ms.h5'
  CHARACTER(150) ::  magnetic_field_directory
  CHARACTER(150) ::  magnetic_field_list
  INTEGER :: time_slice = 000
  REAL(rp) :: rmax =  1.60
  REAL(rp) :: rmin =  0.15
  REAL(rp) :: zmax =  1.65
  REAL(rp) :: zmin = -1.65
  CHARACTER(75) :: outputs_list = '{X,Y,V,B,E,g,eta,flagCon,flagCol,PSIp,ne}' 
    ! List of outputs
    !'{X,Y,V,E,B,g,mu,eta,Prad,Pin,flagCon,flagCol,gradB,curlb,ne,Te,Zeff,PSIp}'
  LOGICAL :: HDF5_error_handling = .TRUE.
  LOGICAL :: FO_GC_compare = .FALSE.
  CHARACTER(30) :: orbit_model = 'GC'
    ! 'FO' for full orbit, 'GCpre' for guiding center with pre-computed
    ! auxiliary fields, 'GCgrad' for guiding center with auxiliary
    ! fields computed with PSPLINE.
  CHARACTER(30) :: field_eval = 'interp'
    ! Set for plasma_model='ANALYTICAL'. Can be 'interp' or 'eqn',
    ! where 'eqn' evaluates particle fields at particle positions and
    ! 'interp' interpolates precomputed fields.
  LOGICAL :: FokPlan = .FALSE.
  LOGICAL :: SameRandSeed = .FALSE.
  LOGICAL :: SC_E = .FALSE.
  LOGICAL :: SC_E_add = .FALSE.
  INTEGER :: pchunk = 1

  !! -----------------------------------------------
  !! plasma_species
  !! As these inputs are vectors with dimension given by the number of species
  !! indicate default values for num_species=1 below, after the input_parameter
  !! namelist is read
  !! -----------------------------------------------
  LOGICAL, DIMENSION(:), ALLOCATABLE :: runaway 
    !! Flag to decide whether a given electron is a runaway (runaway=T)
    !! or not (runaway=F).
  INTEGER, DIMENSION(:), ALLOCATABLE :: ppp
  INTEGER, DIMENSION(:), ALLOCATABLE :: pinit    ! Number of particles per process (mpi)
  REAL(rp), DIMENSION(:), ALLOCATABLE :: q 
    ! Electric charge
  REAL(rp), DIMENSION(:), ALLOCATABLE :: m 
    ! In units of electron mass
  CHARACTER(30), DIMENSION(:), ALLOCATABLE :: spatial_distribution 
    !! String describing the type of initial spatial distribution for
    !! each electron species.
    ! Options are: 'UNIFORM', 'DISK', 'TORUS', 'EXPONENTIAL-TORUS',
    ! 'GAUSSIAN-TORUS', 'ELLIPTIC-TORUS', 'EXPONENTIAL-ELLIPTIC-TORUS',
    ! 'GAUSSIAN-ELLIPTICAL-TORUS', '2D-GAUSSIAN-ELLIPTIC-TORUS-MH',
    ! 'AVALANCHE-4D','TRACER','SPONG-3D','HOLLMANN-3D','HOLLMANN-3D-PSI',
    ! 'MH_psi'
  REAL(rp), DIMENSION(:), ALLOCATABLE :: Ro 
    !! Radial position of the center of the electrons' initial
    !! spatial distribution.
  REAL(rp), DIMENSION(:), ALLOCATABLE :: PHIo 
    !! Azimuthal position of the electrons' initial spatial distribution, 
    !! in case of using a disk at a certain poloidal section.
    ! In degrees
  REAL(rp), DIMENSION(:), ALLOCATABLE :: Zo 
    !! Height of the center of the electrons' initial spatial distribution.
  REAL(rp), DIMENSION(:), ALLOCATABLE :: r_inner 
  REAL(rp), DIMENSION(:), ALLOCATABLE :: r_outter 
  REAL(rp), DIMENSION(:), ALLOCATABLE :: shear_factor 
  REAL(rp), DIMENSION(:), ALLOCATABLE :: sigmaR 
  REAL(rp), DIMENSION(:), ALLOCATABLE :: sigmaZ 
  REAL(rp), DIMENSION(:), ALLOCATABLE :: theta_gauss 
  REAL(rp), DIMENSION(:), ALLOCATABLE :: psi_max
    !! Maximum value of the argument of the 2D gaussian exponential, used for an
    !! indicator function that limits the region of MH sampling
    ! goes as R^2 for HOLLMANN-3D, is psiN_max for HOLLMANN-3D-PSI
  REAL(rp), DIMENSION(:), ALLOCATABLE :: falloff_rate 
  CHARACTER(30), DIMENSION(:), ALLOCATABLE :: energy_distribution
    ! Options are: 'MONOENERGETIC', 'THERMAL', 'AVALANCHE', 
    ! 'EXPERIMENTAL', and 'UNIFORM'
  CHARACTER(30), DIMENSION(:), ALLOCATABLE :: pitch_distribution 
    ! Options are: 'MONOPITCH', 'AVALANCHE', 'EXPERIMENTAL', and 'UNIFORM'.
  REAL(rp), DIMENSION(:), ALLOCATABLE :: Eno 
    ! Initial energy in eV
  REAL(rp), DIMENSION(:), ALLOCATABLE :: etao 
    ! Initial pitch angle
  REAL(rp), DIMENSION(:), ALLOCATABLE :: Eo_lims 
    ! Lower and upper limit of simulated energy range, in eV.
  REAL(rp), DIMENSION(:), ALLOCATABLE :: etao_lims 
    ! Lower and upper limit of simulated pitch-angle range, in degrees.
  REAL(rp), DIMENSION(:), ALLOCATABLE  :: Xtrace 
    ! Initial position of tracer particle for debugging with
    ! spatial_distribution='TRACER'
  REAL(rp), DIMENSION(:), ALLOCATABLE :: Spong_b 
  REAL(rp), DIMENSION(:), ALLOCATABLE :: Spong_w 
  REAL(rp), DIMENSION(:), ALLOCATABLE :: Spong_dlam 
  REAL(rp), DIMENSION(:), ALLOCATABLE :: dth 
    ! Variance of sampling normal variate for pitch angle
  REAL(rp), DIMENSION(:), ALLOCATABLE :: dgam 
    ! Variance of sampling normal variate for pitch angle
  REAL(rp), DIMENSION(:), ALLOCATABLE :: dR 
    ! Variance of sampling normal variate for R location
  REAL(rp), DIMENSION(:), ALLOCATABLE :: dZ 
    ! Variance of sampling normal variate for Z location


  !! -----------------------------------------------
  !! analytical_fields_params
  !! -----------------------------------------------
  REAL(rp) :: Eo = 0. 
    ! In V/m
  CHARACTER(30) :: current_direction = 'ANTI-PARALLEL' 
    ! 'PARALLEL' or 'ANTI-PARALLEL'
  REAL(rp) :: Bo = 2.2 
    ! In Teslas. ITER: 5.42 DIII-D: 2.19
  REAL(rp) :: minor_radius = 0.7 
    ! Minor radius in meters. ITER: 1.5 DIII-D: 0.5
  REAL(rp) :: major_radius = 1.7 
    ! Major radius in meters. ITER: 6.5 DIII-D: 1.6
  REAL(rp) :: qa = 5 
    ! Safety factor at the separatrix (r=a)
  REAL(rp) :: qo = 1.5 
    ! Safety factor at the separatrix (r=a)
  REAL(rp) :: nR= 50
    ! Mesh points in R for analytical interpolation mesh
  REAL(rp) :: nZ= 50
  ! Mesh points in Z for analytical interpolation mesh
  REAL(rp) :: nPHI= 50
  ! Mesh points in PHI for analytical interpolation mesh
  CHARACTER(30) :: E_profile = 'none'
  REAL(rp) :: E_dyn = 0.	
  REAL(rp) :: E_pulse = 5.E-2
  REAL(rp) :: E_width = 2.5E-2
  REAL(rp) ::Ero=0
  ! amplitude of radial electric field in V/m
  REAL(rp) :: rmn =0.6
  !location of Er 
  REAL(rp) :: sigmamn=1.E-2
  ! half width of Er perturbation

  !! -----------------------------------------------
  !! externalPlasmaModel
  !! -----------------------------------------------
  LOGICAL :: Bfield = .FALSE.
  LOGICAL :: B1field = .FALSE.
  LOGICAL :: E1field = .FALSE.
  LOGICAL :: dBfield = .FALSE.
  LOGICAL :: axisymmetric_fields = .FALSE.
  LOGICAL :: Bflux = .FALSE.
  LOGICAL :: Bflux3D = .FALSE.
  LOGICAL :: Efield = .FALSE.
  LOGICAL :: Dim2x1t =.FALSE.
  LOGICAL :: E_2x1t = .FALSE.
  REAL(rp) :: t0_2x1t = 1.405
  INTEGER :: ind0_2x1t = 11
  LOGICAL :: ReInterp_2x1t = .FALSE.
  INTEGER :: res_double=0
  INTEGER :: dim_1D=50
  REAL(rp) :: dt_E_SC=1.E-7
  REAL(rp) :: Ip_exp=2E5
  REAL(rp) :: PSIp_lim=0.8446
  REAL(rp) :: PSIp_0=0.6
  REAL(rp) :: psip_conv=1.0
  REAL(rp)  :: MARS_AMP_Scale=1.0
  REAL(rp)  :: AORSA_AMP_Scale=1.0
  REAL(rp)  :: AORSA_freq=0.0
  REAL(rp)  :: AORSA_nmode=0.0
  CHARACTER(30) :: Analytic_IWL='NONE'
  INTEGER :: ntiles=42
  REAL(rp) :: circumradius=1.016
  LOGICAL :: useLCFS = .FALSE.
  
  !! -----------------------------------------------
  !! plasmaProfiles
  !! -----------------------------------------------
  LOGICAL :: axisymmetric = .TRUE.
  CHARACTER(30) :: filename = 'JFIT_D3D_157576_t1580_1.h5' !
  REAL(rp) :: radius_profile = 0.6
  CHARACTER(30) :: ne_profile = 'RE-EVO-PSIP-G' 
    ! Options are 'FLAT','POLYNOMIAL','RE-EVO','RE-EVO1','RE-EVO-PSI'
  REAL(rp) :: neo = 4.E20 
    ! In m^-3
  REAL(rp) :: n_ne = 2.5E19 !
  REAL(rp) :: n_shelf = 2.5E19
  !	
  REAL(rp), DIMENSION(4) :: a_ne = (/0.99713,0.047037,0.40023,-1.0466/)
  REAL(rp) :: n_REr0 = 0.4
  REAL(rp) :: n_tauion = 1.5e-2
  REAL(rp) :: n_tauin = 7.5e-3
  REAL(rp) :: n_tauout = 1.25e-2	
  REAL(rp) :: n_shelfdelay = 4.e-2
  REAL(rp) :: n_lamfront = 0.005
  REAL(rp) :: n_lamback = 0.005
  REAL(rp) :: n_lamshelf = 0.225
  REAL(rp) :: psiN_0 = 0.8
  
  
  CHARACTER(30) :: Te_profile = 'FLAT' 
    ! Options are 'FLAT' and 'POLYNOMIAL'
  REAL(rp) :: Teo = 1.5 
    ! In eV
  REAL(rp) :: n_Te = 0.1 !
  REAL(rp), DIMENSION(4) :: a_Te = (/1.0046,-0.076652,-2.6429,1.7415/)
  
  CHARACTER(30) :: Zeff_profile = 'FLAT' 
    ! Options are 'FLAT' and 'POLYNOMIAL'
  REAL(rp) :: Zeffo = 1. 
    ! In m^-3
  REAL(rp) :: n_Zeff = 3.0 !
  REAL(rp), DIMENSION(4) :: a_Zeff = (/1.0065,-0.12081,0.02834,-0.11796/)

  !! -----------------------------------------------
  !! CollisionParamsSingleSpecies
  !! -----------------------------------------------
  REAL(rp) :: Te_sing = 2.0 
    ! Background electron temperature in eV
  REAL(rp) :: Ti_sing = 2.0 
    ! Background ion temperature in eV
  REAL(rp) :: ne_sing = 4.E20 
    ! Background electron density in m^-3
  REAL(rp) :: Zeff_sing = 1. 
    ! Effective atomic number
  REAL(rp) :: dTau_sing = 5.E-2
  ! Subcycling time step in collisional time units (Tau)
  REAL(rp) :: p_therm = 1._rp
  LOGICAL :: ConserveLA = .TRUE.
  CHARACTER(30) :: Clog_model = 'HESSLOW'
  CHARACTER(30) :: min_secRE = 'CRIT'
  LOGICAL :: sample_test  = .FALSE.
  REAL(rp)  :: pmin_scale = 1._rp
  LOGICAL :: energy_diffusion = .TRUE.
  CHARACTER(30) :: LAC_gam_resolution = '2EXP'
  LOGICAL :: FP_bremsstrahlung = .TRUE.
  LOGICAL :: pitch_diffusion = .TRUE.
  INTEGER :: ngrid1 = 100
  REAL(rp)  :: Clog_const = 20._rp
  
  !! -----------------------------------------------
  !! CollisionParamsMultipleSpecies
  !! -----------------------------------------------
  INTEGER :: num_impurity_species = 1
  REAL(rp) :: Te_mult = 2.0 
    ! Background electron temperature in eV
  REAL(rp) :: ne_mult = 4.E20 
    ! Background electron density in 1/m^3
  REAL(rp), DIMENSION(10)  :: Zo_mult = 10.0 
    ! Full nuclear charge of each impurity: Z=1 for D, Z=10 for Ne, Z=18 for Ar
  REAL(rp), DIMENSION(10)  :: Zj_mult = 1.0  
    ! Average charge state of each impurity
  REAL(rp), DIMENSION(10)  :: nz_mult = 4.E20 
    ! Impurity densities
  REAL(rp), DIMENSION(10)  :: IZj_mult = 15.7596
    ! Ionization energy of impurity in eV
  CHARACTER(20) :: neut_prof = 'UNIFORM'

  !! -----------------------------------------------
  !! AvalancheGenerationPDF
  !! -----------------------------------------------
  REAL(rp) :: max_pitch_angle_aval = 10.0 
    ! Maximum pitch angle of sampled PDF in degrees
  REAL(rp) :: min_pitch_angle_aval = 0.0 !
  REAL(rp) :: dth_aval = 1.0
    ! Variance of sampling normal variate for pitch angle
  REAL(rp) :: dp_aval = 10.0
    ! Variance of sampling normal variate for momentum
  REAL(rp) :: dR_aval = 0.05
    ! Variance of sampling normal variate for R location
  REAL(rp) :: dZ_aval = 0.05
    ! Variance of sampling normal variate for Z location
  REAL(rp) :: max_energy_aval  = 60E6 
    ! Maximum energy of sampled PDF in MeV
  REAL(rp) :: min_energy_aval = 7.0E6 !
  REAL(rp) :: ne_aval = 8.5E18 
    ! Background electron density in m^-3
  REAL(rp) :: Zeff_aval = 1.0 
    ! Effective atomic number of ions
  REAL(rp) :: Epar_aval = 0.7427 
    ! Parallel electric field in V/m
  REAL(rp) :: Te_aval = 1.0 
  ! Background electron temperature in eV

  !! -----------------------------------------------
  !! ExperimentalPDF
  !! -----------------------------------------------
  REAL(rp) :: E_expt = 7.67042 
    ! Parallel electric field in units of Ec
  REAL(rp) :: Zeff_expt = 1.36632 
    ! Effective atomic number
  REAL(rp) :: max_pitch_angle_expt = 176.257 
    ! In degrees
  REAL(rp) :: min_pitch_angle_expt = 158.25 
    ! In degrees
  REAL(rp) :: min_energy_expt = 9.54997E6 
    ! In eV
  REAL(rp) :: max_energy_expt = 44.3241E6 
    ! In eV
  REAL(rp) :: k_expt = 11.89 
    ! Shape factor of Gamma distribution for energy
  REAL(rp) :: t_expt = 0.65 
    ! Scale factor of Gamma distribution for energy
  REAL(rp) :: Bo_expt = 2.2
    ! Characteristic magnetic field
  REAL(rp) :: lambda_expt = 4.0E-6 
    ! Characteristic wavelength
  REAL(rp) :: A_fact_expt=1.
  CHARACTER(30) :: filename_exp = 'Exp_PDF.h5' !
  
  !! -----------------------------------------------
  !! HollmannPDF
  !! -----------------------------------------------
  CHARACTER(30) :: filename_Hollmann = 'Hollmann_PDF_HR.h5' !
  INTEGER :: rho_ind = 7
  REAL(rp) :: Eo_Hollmann = 24.56
  ! Toroidal electric field from experimental diagnostics before SPI in
  ! physical units
  REAL(rp) :: E_Hollmann = 11.
    ! Parallel electric field in units of Ec
  REAL(rp) :: sigma_E_Hollmann=.2
  REAL(rp) :: Zeff_Hollmann = 5. 
    ! Effective atomic number
  REAL(rp) :: sigma_Z_Hollmann=10.	
  REAL(rp) :: max_pitch_angle_Hollmann = 40. 
    ! In degrees
  REAL(rp) :: min_pitch_angle_Hollmann = 0. 
    ! In degrees
  REAL(rp) :: min_energy_Hollmann = 1.E6 
    ! In eV
  REAL(rp) :: max_energy_Hollmann = 60.E6
    ! For Hollmann_PDF_HR.h5, needs to be less than 80 MeV
    ! In eV
  CHARACTER(30) :: current_direction_Hollmann = 'ANTICLOCKWISE'
  REAL(rp) :: Bo_Hollmann = 2.2 
    ! Characteristic magnetic field
  REAL(rp) :: lambda_Hollmann = 4.0E-6 
    ! Characteristic wavelengt
  REAL(rp) :: A_fact_Hollmann=1.
  LOGICAL :: gam_min_from_col=.FALSE.

  !! -----------------------------------------------
  !! EnergyGammaPDF
  !! -----------------------------------------------
  REAL(rp) :: min_energy_gamma = 1.0E6 
    ! In eV
  REAL(rp) :: max_energy_gamma = 30.0E6 
    ! In eV
  REAL(rp) :: k_gamma = 11.89 
    ! Shape factor of Gamma distribution for energy
  REAL(rp) :: t_gamma = 0.65 
    ! Scale factor of Gamma distribution for energy

  !! -----------------------------------------------
  !! SimpleEquilibriumPDF
  !! -----------------------------------------------
  REAL(rp) :: E_simple  = 4.5 
    ! Parallel electric field in units of Ec
  REAL(rp) :: Zeff_simple  = 4.0 
    ! Effective atomic number
  REAL(rp) :: max_pitch_angle_simple  = 60.0 
    ! In degrees
  REAL(rp) :: min_pitch_angle_simple  = 0.0 
    ! In degrees
  REAL(rp) :: Bo_simple  = 2.0 
    ! Characteristic magnetic field
  REAL(rp) :: lambda_simple = 890.0E-9 
    ! Characteristic wavelength

CONTAINS

  subroutine read_namelist(params,infile,echo_in,outdir)

    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    CHARACTER(*), INTENT(IN) :: infile,outdir
    LOGICAL, INTENT(IN) :: echo_in

    INTEGER :: read_stat,nc
    INTEGER :: number_of_namelists=0,il,inst
    INTEGER, DIMENSION(20) :: namel_order=0
    CHARACTER(20) :: tempfile
    CHARACTER(128) :: ctmp
    CHARACTER(128) :: outfile
    LOGICAL :: reading
    INTEGER :: mpierr
    INTEGER :: tmp
    
    !! Namelist declarations
    NAMELIST /input_parameters/ restart,field_model,magnetic_field_filename, &
         simulation_time,snapshot_frequency,dt,num_species,radiation, &
         collisions,collisions_model,outputs_list,minimum_particle_energy, &
         HDF5_error_handling,orbit_model,field_eval,proceed,profile_model, &
         restart_overwrite_frequency,FokPlan,GC_rad_model,bound_electron_model,&
         FO_GC_compare,SameRandSeed,SC_E,reinit,SC_E_add,time_slice,rmax, &
         rmin,zmax,zmin,pchunk,magnetic_field_directory,magnetic_field_list,&
         LargeCollisions,load_balance
    NAMELIST /plasma_species/ ppp,q,m,Eno,etao,Eo_lims,etao_lims,runaway, &
         spatial_distribution,energy_distribution,pitch_distribution,Ro, &
         PHIo,Zo,r_inner,r_outter,falloff_rate,shear_factor,sigmaR,sigmaZ, &
         theta_gauss,psi_max,Xtrace,Spong_b,Spong_w,Spong_dlam,dth,dR,dZ,dgam,&
         pinit
    NAMELIST /analytical_fields_params/ Bo,minor_radius,major_radius,&
         qa,qo,Eo,current_direction,nR,nZ,nPHI,dim_1D,dt_E_SC,Ip_exp, &
         E_dyn,E_pulse,E_width,E_profile,Ero,rmn,sigmamn
    NAMELIST /externalPlasmaModel/ Efield, Bfield, Bflux,Bflux3D,dBfield, &
         axisymmetric_fields, Eo,E_dyn,E_pulse,E_width,res_double, &
         dim_1D,dt_E_SC,Ip_exp,PSIp_lim,Dim2x1t,t0_2x1t,E_2x1t,ReInterp_2x1t, &
         ind0_2x1t,PSIp_0,B1field,psip_conv,MARS_AMP_Scale,Analytic_IWL, &
         ntiles,circumradius,AORSA_AMP_Scale,AORSA_freq,AORSA_nmode,E1field, &
         useLCFS
    NAMELIST /plasmaProfiles/ radius_profile,ne_profile,neo,n_ne,a_ne, &
         Te_profile,Teo,n_Te,a_Te,n_REr0,n_tauion,n_lamfront,n_lamback, &
         Zeff_profile,Zeffo,n_Zeff,a_Zeff,filename,axisymmetric, &
         n_lamshelf,n_shelfdelay,n_tauin,n_tauout,n_shelf,psiN_0
    NAMELIST /CollisionParamsSingleSpecies/ Te_sing,Ti_sing,ne_sing, &
         Zeff_sing,dTau_sing,p_therm,ConserveLA,Clog_model,sample_test,&
         min_secRE,pmin_scale,energy_diffusion,LAC_gam_resolution, &
         FP_bremsstrahlung,pitch_diffusion,ngrid1,Clog_const
    NAMELIST /CollisionParamsMultipleSpecies/ num_impurity_species,Te_mult, &
         ne_mult,Zo_mult,Zj_mult,nz_mult,IZj_mult,neut_prof
    NAMELIST /AvalancheGenerationPDF/ max_pitch_angle_aval, &
         min_pitch_angle_aval,max_energy_aval,min_energy_aval,ne_aval, &
         Zeff_aval,Epar_aval,Te_aval,dth_aval,dp_aval,dR_aval,dZ_aval
    NAMELIST /ExperimentalPDF/ max_pitch_angle_expt,min_pitch_angle_expt, &
         max_energy_expt,min_energy_expt,Zeff_expt,E_expt,k_expt,t_expt, &
         Bo_expt,lambda_expt,A_fact_expt,filename_exp
    NAMELIST /HollmannPDF/ E_Hollmann,Zeff_Hollmann,max_pitch_angle_Hollmann, &
         min_pitch_angle_Hollmann,max_energy_Hollmann, &
         min_energy_Hollmann,filename_Hollmann,Bo_Hollmann,lambda_Hollmann, &
         current_direction_Hollmann,A_fact_Hollmann,sigma_E_Hollmann, &
         sigma_Z_Hollmann,Eo_Hollmann,rho_ind,gam_min_from_col
    NAMELIST /SimpleEquilibriumPDF/ max_pitch_angle_simple, &
         min_pitch_angle_simple,Zeff_simple,E_simple, &
         Bo_simple,lambda_simple
    NAMELIST /EnergyGammaPDF/ max_energy_gamma,min_energy_gamma,k_gamma,t_gamma

!!-----------------------------------------------------------------------
!!     open input file.
!!     Remove comments from input file and put into temporary file.
!!-----------------------------------------------------------------------
    tempfile='tempinput.korc'
    if (params%mpi_params%rank.eq.0) then
       CALL rmcoment(infile,tempfile)
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)   
    OPEN(UNIT=default_unit_open,FILE=tempfile,STATUS='OLD',POSITION='REWIND')
!!-----------------------------------------------------------------------
!!    check namelist file for namelist order and number.
!!-----------------------------------------------------------------------
    DO
       READ(UNIT=default_unit_open,FMT='(a)',IOSTAT=read_stat) ctmp 
       IF (read_stat/=0) EXIT
       nc=LEN_TRIM(ctmp)
       IF (nc<1) CYCLE
       ctmp=ADJUSTL(ctmp)
       reading=.false.
       IF (ctmp(1:1)=='&') THEN
          number_of_namelists=number_of_namelists+1
!!-----------------------------------------------------------------------
!!         trim all but the namelist name.
!!-----------------------------------------------------------------------
          DO il=2,nc+1
             IF (ctmp(il:il)/=' ') THEN
                IF (.NOT.reading) inst=il
                reading=.true.
                CYCLE
             ENDIF
             IF (ctmp(il:il)==' '.AND.reading) THEN
                ctmp=ctmp(inst:il-1)
                EXIT
             ENDIF
          ENDDO
          BACKSPACE(default_unit_open)
!!-----------------------------------------------------------------------
!!         select and read namelist.
!!-----------------------------------------------------------------------
          SELECT CASE(TRIM(ctmp))
          CASE('input_parameters')
             READ(UNIT=default_unit_open,NML=input_parameters,IOSTAT=read_stat)
             
          CASE('plasma_species')

             !write(6,*) 'reading plasma_species namelist'
             ALLOCATE(runaway(num_species))
             ALLOCATE(ppp(num_species))
             ALLOCATE(pinit(num_species))
             ALLOCATE(q(num_species))
             ALLOCATE(m(num_species))
             ALLOCATE(spatial_distribution(num_species))
             ALLOCATE(Ro(num_species))
             ALLOCATE(PHIo(num_species))
             ALLOCATE(Zo(num_species))
             ALLOCATE(r_inner(num_species))
             ALLOCATE(r_outter(num_species))
             ALLOCATE(shear_factor(num_species))
             ALLOCATE(sigmaR(num_species))
             ALLOCATE(sigmaZ(num_species))
             ALLOCATE(theta_gauss(num_species))
             ALLOCATE(psi_max(num_species))
             ALLOCATE(falloff_rate(num_species))
             ALLOCATE(energy_distribution(num_species))
             ALLOCATE(pitch_distribution(num_species))
             ALLOCATE(Eno(num_species))
             ALLOCATE(etao(num_species))
             ALLOCATE(Eo_lims(2_idef*num_species))
             ALLOCATE(etao_lims(2_idef*num_species))
             ALLOCATE(Xtrace(3_idef*num_species))
             ALLOCATE(Spong_b(num_species))
             ALLOCATE(Spong_w(num_species))
             ALLOCATE(Spong_dlam(num_species))
             ALLOCATE(dth(num_species))
             ALLOCATE(dgam(num_species))
             ALLOCATE(dR(num_species))
             ALLOCATE(dZ(num_species))

             if (num_species.eq.1) then
                runaway = .FALSE.
                ppp = 1E0
                pinit = 0
                q = -1.0 
                m = 1.0 
                spatial_distribution = 'TRACER'
                Ro = 1.1
                PHIo = 0.0
                Zo = -0.05
                r_inner = 0.0 
                r_outter = 0.6 
                shear_factor = 0.35 
                sigmaR = 1.e6
                sigmaZ = 0.2
                theta_gauss = 0.0
                psi_max=.8446
                falloff_rate = 0.0 
                energy_distribution = 'MONOENERGETIC'   
                pitch_distribution = 'MONOPITCH'	
                Eno = 10.0E6 
                etao = 1.0  
                Eo_lims = (/1.0E6,50.0E6/) 
                etao_lims = (/0.0,20.0 /)
                Xtrace =(/1.1,0.0,0.0/)
                Spong_b = 0.2
                Spong_w = 0.1
                Spong_dlam = 0.1
                dth = 3.
                dgam = 3.	  
                dR = 0.1
                dZ = 0.1
             else
                if (params%mpi_params%rank .EQ. 0) then
                   write(output_unit_write,'("Need to supply all inputs for num_species .gt. 1")')
                end if
             end if
             
             READ(UNIT=default_unit_open,NML=plasma_species,IOSTAT=read_stat)

             if (pinit(1).eq.0) pinit(:)=ppp(:)
             ! set pinit equal to ppp if no pinit input
             
          CASE('analytical_fields_params')
             READ(UNIT=default_unit_open,NML=analytical_fields_params,IOSTAT=read_stat)
          CASE('externalPlasmaModel')
             READ(UNIT=default_unit_open,NML=externalPlasmaModel,IOSTAT=read_stat)
          CASE('plasmaProfiles')
             READ(UNIT=default_unit_open,NML=plasmaProfiles,IOSTAT=read_stat)
          CASE('CollisionParamsSingleSpecies')
             READ(UNIT=default_unit_open,NML=CollisionParamsSingleSpecies,IOSTAT=read_stat)
          CASE('CollisionParamsMultipleSpecies')
             READ(UNIT=default_unit_open,NML=CollisionParamsMultipleSpecies,IOSTAT=read_stat)
          CASE('AvalancheGenerationPDF')
             READ(UNIT=default_unit_open,NML=AvalancheGenerationPDF,IOSTAT=read_stat)
          CASE('ExperimentalPDF')
             READ(UNIT=default_unit_open,NML=ExperimentalPDF,IOSTAT=read_stat)
          CASE('HollmannPDF')
             READ(UNIT=default_unit_open,NML=HollmannPDF,IOSTAT=read_stat)
          CASE('EnergyGammaPDF')
             READ(UNIT=default_unit_open,NML=EnergyGammaPDF,IOSTAT=read_stat)
          CASE('SimpleEquilibriumPDF')
             READ(UNIT=default_unit_open,NML=SimpleEquilibriumPDF,IOSTAT=read_stat)
          CASE DEFAULT
             write(output_unit_write,*) (TRIM(ctmp)//' is an unrecognized &
                  &namelist.')
             call korc_abort(13)
          END SELECT
          IF (read_stat/=0) then
             write(output_unit_write,*) ('Error reading namelist '//TRIM(ctmp)//'.')             
             call korc_abort(13)
          end if
       ENDIF
    ENDDO

!!-----------------------------------------------------------------------
!!     close input file.
!!       Delete it since it is the temporary file
!!-----------------------------------------------------------------------
    if (params%mpi_params%rank.ne.0) then
       CLOSE(default_unit_open)
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
    if (params%mpi_params%rank.eq.0) then
       CLOSE(default_unit_open,STATUS='DELETE')
    end if
!!-----------------------------------------------------------------------
!!     echo the input parameters to the output file.
!!-----------------------------------------------------------------------


    IF (echo_in) THEN
       if (params%mpi_params%rank .EQ. 0) then

          WRITE(output_unit_write,'(a,/)') 'VALUE OF ALL INPUTS:'
          WRITE(UNIT=output_unit_write,NML=input_parameters)
          WRITE(output_unit_write,'(/)') 
          WRITE(UNIT=output_unit_write,NML=plasma_species)
          WRITE(output_unit_write,'(/)') 
          WRITE(UNIT=output_unit_write,NML=analytical_fields_params)
          WRITE(output_unit_write,'(/)') 
          WRITE(UNIT=output_unit_write,NML=externalPlasmaModel)
          WRITE(output_unit_write,'(/)') 
          WRITE(UNIT=output_unit_write,NML=plasmaProfiles)
          WRITE(output_unit_write,'(/)') 
          WRITE(UNIT=output_unit_write,NML=CollisionParamsSingleSpecies)
          WRITE(output_unit_write,'(/)') 
          WRITE(UNIT=output_unit_write,NML=CollisionParamsMultipleSpecies)
          WRITE(output_unit_write,'(/)') 
          WRITE(UNIT=output_unit_write,NML=AvalancheGenerationPDF)
          WRITE(output_unit_write,'(/)') 
          WRITE(UNIT=output_unit_write,NML=ExperimentalPDF)
          WRITE(output_unit_write,'(/)') 
          WRITE(UNIT=output_unit_write,NML=HollmannPDF)
          WRITE(output_unit_write,'(/)') 
          WRITE(UNIT=output_unit_write,NML=EnergyGammaPDF)
          WRITE(output_unit_write,'(/)') 
          WRITE(UNIT=output_unit_write,NML=SimpleEquilibriumPDF)
          WRITE(output_unit_write,'(/)')
             
       end if
    end if

!!---------------------------------------------------------
!!     some tests
!!---------------------------------------------------------

    !write(6,*) TRIM(outputs_list),len(TRIM(outputs_list))
    
    tmp=len(TRIM(outputs_list))
    if (outputs_list(tmp:tmp).ne.'}') then
       if(params%mpi_params%rank.eq.0) then
          write(6,*) &
               'Check that enough characters are allocated for&
               & outputs list!'
       end if
       call korc_abort(13)
    end if

    !write(6,*) TRIM(magnetic_field_filename),len(TRIM(magnetic_field_filename))

    tmp=len(TRIM(magnetic_field_filename))
    if (magnetic_field_filename(tmp-2:tmp).ne.'.h5'.and. &
         magnetic_field_filename(tmp-5:tmp-5).ne.'.') then
       if(params%mpi_params%rank.eq.0) then
          write(6,*) &
               'Check that enough characters are allocated for&
               & magnetic field filename!'
       end if
       call korc_abort(13)
    end if 
      
    end subroutine read_namelist

    SUBROUTINE rmcoment(fileold,filenew)

    CHARACTER(*), INTENT(IN) :: fileold,filenew
    CHARACTER(128) :: line
    INTEGER, PARAMETER :: nold=55,nnew=56
    INTEGER cmax, ios
    LOGICAL :: file_exist
!!-----------------------------------------------------------------------
!!     Open files, but make sure the old one exists first.
!!-----------------------------------------------------------------------
    INQUIRE(FILE=fileold,EXIST=file_exist)
    IF(.NOT. file_exist) THEN
       PRINT *,'The file "',fileold,'" could not be found.'
       STOP
    ENDIF

    OPEN(UNIT=default_unit_open,FILE=fileold,status="OLD",form='formatted')
    OPEN(UNIT=default_unit_write,FILE=filenew,status='REPLACE')

!!-----------------------------------------------------------------------
!!     Strip comments.     Note: line lengths limited to 127 characters
!!-----------------------------------------------------------------------
    DO
       READ(UNIT=default_unit_open,FMT='(a)',IOSTAT=ios) line
       IF (ios /= 0) EXIT
       cmax=1
       DO WHILE(line(cmax:cmax).NE.'!' .AND. cmax .LE. 127)
          cmax=cmax+1
       ENDDO
       IF(cmax .GT. 1) WRITE(default_unit_write,'(a)') line(1:cmax-1)
    ENDDO

!!-----------------------------------------------------------------------
!!     Close files and exit
!!-----------------------------------------------------------------------
    CLOSE(default_unit_open)
    CLOSE(default_unit_write)

  END SUBROUTINE rmcoment
  
end module korc_input
