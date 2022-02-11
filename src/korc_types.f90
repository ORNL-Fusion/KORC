module korc_types
  !! @note Module containing the definition of KORC derived types and
  !! KORC variables, the building blocks of the code. @endnote
#ifdef FIO
  USE, INTRINSIC :: iso_c_binding
#endif
  implicit none

  ! * * * * * * * * * * * * * * * * * * * * !
  ! * * * Real and integer precisions * * * !
  ! * * * * * * * * * * * * * * * * * * * * !

  INTEGER, PUBLIC, PARAMETER 	:: is = KIND(INT(1,1)) 
  !! Definition of 1 Byte (8 bits) Fortran KORC integer type.
  INTEGER, PUBLIC, PARAMETER 	:: ip = KIND(INT(1,8)) 
  !! Definition of 8 Bytes (64 bits) Fortran KORC integer type.
  INTEGER, PUBLIC, PARAMETER 	:: idef = KIND(1) 
  !! Definition of the default KORC integer type on the system where
  !! KORC is compiled.
  INTEGER, PUBLIC, PARAMETER 	:: rdef = KIND(1.0) 
  !! Definition of the default KORC real type on the system where
  !! KORC is compiled.
#ifdef DOUBLE_PRECISION
  INTEGER, PUBLIC, PARAMETER 	:: rp = KIND(0.d0) 
  !! Definition of the KORC double precision real type.
#elif SINGLE_PRECISION
  INTEGER, PUBLIC, PARAMETER 	:: rp = KIND(1.0) 
  !! Definition of the KORC single precision real type.
#endif
  REAL(rp), PUBLIC, PARAMETER :: korc_zero = 1.0E-15 
  !! Definition of the zero in KORC.

  ! * * * * * * * * * * * * * * * * * * * * !
  ! * * * Real and integer precisions * * * !
  ! * * * * * * * * * * * * * * * * * * * * !

  INTEGER, PUBLIC, PARAMETER 	:: MAX_STRING_LENGTH = 1000 
  !! Default length of a KORC_STRING variable.
  INTEGER, PUBLIC, PARAMETER 	:: default_unit_open = 101 
  !! Default file unit for opening and reading from an external text file.
  INTEGER, PUBLIC, PARAMETER 	:: default_unit_write = 201 
  !! Default file unit for opening and writing to external an external text file.
  INTEGER, PUBLIC, PARAMETER 	:: output_unit_write = 202 
  !! Default file unit for opening and writing to external an external text file.

  !> @note KORC string type of length MAX_STRING_LENGTH. @endnote
  TYPE, PUBLIC :: KORC_STRING 
     CHARACTER(MAX_STRING_LENGTH) :: str
  END TYPE KORC_STRING

  ! * * * * * * * * * * * * * * * * * * * * !
  ! * * * Basic korc array structures * * * !
  ! * * * * * * * * * * * * * * * * * * * * !

  TYPE, PUBLIC :: V_FIELD_3D
     !! @note KORC 3-D vector field type @endnote
     !! This KORC type represents a 3-D vector field varible in
     !! cylindrical coordinates. For example, this could be the 3-D magnetic
     !! field, which can be written as $$\mathbf{B}(R,\phi,Z) = B_R(R,\phi,Z)
     !! \hat{R} + B_\phi(R,\phi,Z) \hat{\phi} + B_Z(R,\phi,Z) \hat{Z}.$$
     !! All the members (components) of the V_FIELD_3D type follow the
     !! following index convention:
     !! (\(R\) index,\(\phi\) index, \(Z\) index)

     REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: R 
     !! \(R\) component of the vector field variable.
     REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: PHI 
     !! \(\phi\) component of the vector field variable.
     REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: Z 
     !! \(Z\) component of the vector field variable.
  END TYPE V_FIELD_3D

  TYPE, PUBLIC :: V_FIELD_2DX
     !! @note KORC 2-D vector field type @endnote
     !! This KORC type represents a 2-D vector field varible in
     !! cartesian coordinates. For example, this could be the magnetic
     !! field in an axisymmetirc plasma, which can be written as
     !! $$\mathbf{B}(R,Z) = B_X(R,Z)
     !! \hat{X} + B_Y(R,Z) \hat{Y} + B_Z(R,Z) \hat{Z}.$$
     !! All the members (components) of the V_FIELD_3D type follow the
     !! following index convention:
     !! (\(X\) index,\(Y\) index, \(Z\) index)

     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: X 
     !! \(R\) component of the vector field variable.
     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: Y 
     !! \(\phi\) component of the vector field variable.
     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: Z 
     !! \(Z\) component of the vector field variable.
  END TYPE V_FIELD_2DX

  TYPE, PUBLIC :: V_FIELD_2D
     !! @note KORC 2-D vector field type @endnote
     !! This KORC type represents a 2-D vector field varible in cylindrical
     !! coordinates. For example, this could be the magnetic
     !! field in an axisymmetric plasma, which can be written as
     !! $$\mathbf{B}(R,Z) = B_R(R,Z) \hat{R} + B_\phi(R,Z) \hat{\phi} + B_Z(R,Z)
     !! \hat{Z}.$$
     !! All the members (components) of the V_FIELD_2D type follow the
     !! following index convention:
     !! (\(R\) index,\(Z\) index).
     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: R 
     !! \(R \) component of the vector field variable.
     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: PHI 
     !! \(\phi \) component of the vector field variable.
     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: Z 
     !! \(Z \) component of the vector field variable.
  END TYPE V_FIELD_2D

  TYPE, PUBLIC :: V_FIELD_1D
     !! @note KORC 1-D vector field type @endnote
     !! This KORC type represents a 1-D vector field varible in cylindrical
     !! coordinates. For example, this could be the magnetic
     !! field in an axisymmetric plasma, which can be written as
     !! $$\mathbf{B}(r) = B_R(r) \hat{R} + B_\phi(r) \hat{\phi} + B_Z(r)
     !! \hat{Z}.$$
     !! All the members (components) of the V_FIELD_1D type follow the
     !! following index convention:
     !! (\(r\) index).
     REAL(rp), DIMENSION(:), ALLOCATABLE :: R 
     !! \(R \) component of the vector field variable.
     REAL(rp), DIMENSION(:), ALLOCATABLE :: PHI 
     !! \(\phi \) component of the vector field variable.
     REAL(rp), DIMENSION(:), ALLOCATABLE :: Z 
     !! \(Z \) component of the vector field variable.
  END TYPE V_FIELD_1D
  
  ! * * * * * * * * * * * * * * * * * * * * !
  ! * * * Basic korc array structures * * * !
  ! * * * * * * * * * * * * * * * * * * * * !

  TYPE, PRIVATE :: KORC_MPI
     !! @note KORC derived type to keep relevant MPI parameters. @endnote
     INTEGER :: nmpi 
     !! Number of MPI processes.
     INTEGER :: nmpi_prev 
     INTEGER :: rank 
     !! Rank in WORLD COMMON communicator.
     INTEGER :: rank_topo 
     !! Rank in mpi_topo communicator
     INTEGER :: mpi_topo 
     !! MPI communicator for a certain topology.
  END TYPE KORC_MPI


  TYPE, PUBLIC :: CHARCS_PARAMS
     !! @note KORC derived type containing characteristic scales used in the normalization of KORC variables. @endnote
     !! These characteristic scales are problem-dependent quantities. They are calculated in [[korc_units(module)]] using the input
     !! parameters of a given KORC simulation.

     REAL(rp) :: time
     !! @note Characteristic non-relativistic time scale given by \(1/\omega_{ce}\), where \(\omega_{ce}=e B_0/m_e\) is the
     !! largest electron cyclotron frequency in the simulation. @endnote
     REAL(rp) :: time_r
     !! @note Characteristic relativistic time scale given by \(1/\omega_{ce}\), where \(\omega_{ce}=e B_0/\gamma m_e\) is the
     !! largest relativistic electron cyclotron frequency in the simulation. @endnote
     REAL(rp) :: velocity 
     !! Characteristic velocity. This is fixed to the speed of \(c\).
     REAL(rp) :: length 
     !! Characteristic length scale calculated as \(c\) times the relativistic time scale.
     REAL(rp) :: mass 
     !! Characteristic particle mass. This is equal to the electron mass \(m_e\).
     REAL(rp) :: charge 
     !! Characteristic particle charge. This is equal to the electron charge \(q_e\).
     REAL(rp) :: density 
     !! Characteristic particle density. This is equal to \(1/l^3\), with \(l\) the characteristic length.
     REAL(rp) :: Eo 
     !! Characteristic electric field \(E_0\). Usually \(E_0\) at the magnetic axis.
     REAL(rp) :: Bo 
     !! Characteristic magnetic field \(B_0\). Usually \(B_0\) at the magnetic axis.
     REAL(rp) :: energy 
     !! Characteristic energy. This is equal to \(m_e c^2\).
     REAL(rp) :: pressure 
     !! Characteristic pressure. @todo This needs to be defined.
     REAL(rp) :: temperature 
     !! Characteristic plasma temperature (Joules). This is equal to \(m_e c^2\).
  END TYPE CHARCS_PARAMS


  TYPE, PUBLIC :: KORC_PARAMS
     !! @note Core KORC parameters. @endnote
     !!  This KORC derived type contains the variables that control KORC's
     !! core behavior.

     CHARACTER(MAX_STRING_LENGTH) :: path_to_inputs 
     !! Absolute path to KORC's input file.
     CHARACTER(MAX_STRING_LENGTH) :: path_to_outputs 
     !! Absolute path to the outputs' folder.
     INTEGER 			:: num_omp_threads 
     !! Number of open MP threads per MPI process used in the simulation.
     LOGICAL 			:: restart 
     !! Flag to indicate if the simulations proceeds (restart=T) or not
     !! (restart=F). Restart simulation that exited before simulation_time
     !! reached.
     LOGICAL 			:: proceed 
     !! Flag to indicate if the simulations continues (proceed=T) or not
     !! (proceed=F). Append simulation results after previous simulation_time
     !! reached.
     LOGICAL 			:: load_balance
     LOGICAL  :: reinit
     !! Flag to begin a new simulation, reinitializing from restart file state
     REAL(rp) 			:: simulation_time 
     !! Total simulation time in seconds.
     REAL(rp) 			:: snapshot_frequency 
     !! Time between snapshots in time of the simulation.
     REAL(rp) 			:: restart_overwrite_frequency 
     !! Time between overwrites of restart file in time of the simulation.
     REAL(rp) 			:: dt 
     !! Time step in the simulation as a fraction of the relativistic electron
     !! gyro-period \(\tau_e = 2\pi\gamma m_e/eB_0\).
     REAL(rp) 			:: time = 0.0_rp 
     !! Current physical time in the simulation.
     INTEGER(ip) 			:: ito = 0_ip 
     !! Initial time iteration in the simulation, this is different from zero
     !! in case is a restarting simulation.
     INTEGER(ip) 			:: it = 0_ip 
     !! Current time iteration in the simulation, this is different from zero
     !! in case is a restarting simulation.
     REAL(rp)  :: init_time = 0.0_rp
     !! Time at the beginning of a run with proceed=T     
     INTEGER(ip) 			:: t_steps
     INTEGER(ip) 			:: prev_iter_2x1t 
     !! Number of time steps needed for evolving the electrons up to
     !! "simulation_time".
     INTEGER(ip) 			:: t_skip
     INTEGER(ip) 			:: t_it_SC=1_ip 
     INTEGER(ip) 			:: output_cadence
     INTEGER(ip) 			:: coll_per_dump
     INTEGER(ip) 			:: orbits_per_coll
     REAL(rp) 			:: coll_per_dump_dt
     !! Time iteration offset used to decide when the outputs are generated.
     INTEGER(ip) 			:: restart_output_cadence 
     !! Time iteration offset used to decide when the restart files are
     !! generated.
     INTEGER(ip) 			:: coll_cadence 
     INTEGER(ip) 			:: num_snapshots 
     !! Number of snapshots in time for generating the output files.
     INTEGER 			:: num_species 
     !! Number of different populations of simulated relativistic electrons
     !! in KORC.
     REAL(rp) 			:: minimum_particle_energy 
     !! Minimum allowed energy of simulated electrons. 
     !! @todo To improve the criterium to decide when an electron will not
     !! be followed anymore in the simulation.
     REAL(rp) 			:: minimum_particle_g 
     !! Minimum allowed relativistic factor \(\gamma\) of simulated electrons.
     LOGICAL 			:: radiation 
     !! Flag to indicate if synchrotron radiation losses are included
     !! (radiation=T) or not (radiation=F).
     LOGICAL 			:: collisions 
     !! Flag to indicate if collisionsare included (collisions=T) or not
     !! (collisions=F).
     LOGICAL 			:: LargeCollisions 
     CHARACTER(MAX_STRING_LENGTH) :: GC_rad_model
     CHARACTER(MAX_STRING_LENGTH) :: collisions_model 
     !! String with the name of the collisions model to be used in the
     !! simulation.
     CHARACTER(MAX_STRING_LENGTH) :: bound_electron_model 
     CHARACTER(MAX_STRING_LENGTH) :: field_model
     CHARACTER(MAX_STRING_LENGTH) :: profile_model 
     !! String with the name of the model for the fields and plasma profiles.
     CHARACTER(MAX_STRING_LENGTH) :: magnetic_field_filename
     CHARACTER(MAX_STRING_LENGTH) :: magnetic_field_directory
     CHARACTER(MAX_STRING_LENGTH) :: magnetic_field_list
     CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: magnetic_field_filenames
     REAL(rp), DIMENSION(:), ALLOCATABLE :: time_of_filenames
     !! String with the name of the model for the fields and plasma profiles.
     CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: outputs_list 
     !! List of electron variables to include in the outputs.
     INTEGER 			:: HDF5_error_handling 
     !! Flag to indicate whether we allow HDF5 to show warnings 
     !! during runtime (HDF5_error_handling=1) or not (HDF5_error_handling=0)
     TYPE(KORC_MPI) 		:: mpi_params 
     !! An instance of the KORC_MPI derived type.
     TYPE(CHARCS_PARAMS) 		:: cpp 
     !! An instance of the CHARCS_PARAMS derived type.
     LOGICAL :: FO_GC_compare
     CHARACTER(MAX_STRING_LENGTH) :: orbit_model
     !! String with the name of the orbit model ('FO' or 'GC').
     CHARACTER(MAX_STRING_LENGTH) :: field_eval
     !! String with the name of the field evaluation method for
     !! analytical fields ('interp' or 'eqn')
     LOGICAL :: GC_coords
     !! Flag to [[get_fields]] to control whether cartesian to cylindrical
     !! coordinate transformation needs to be performed
     LOGICAL :: FokPlan
     !! Flag to decouple spatial-dependence of evolution
     LOGICAL :: SameRandSeed
     LOGICAL :: SC_E
     LOGICAL :: SC_E_add
     INTEGER                      :: time_slice !< M3D-C1 time slice to use.
     REAL(rp)                     :: rmax !< Maximum r for M3D-C1 fields.
     REAL(rp)                     :: rmin !< Minimum r for M3D-C1 fields.
     REAL(rp)                     :: zmax !< Maximum z for M3D-C1 fields.
     REAL(rp)                     :: zmin !< Minimum z for M3D-C1 fields.
     INTEGER  :: pchunk
     !! number of particles per vectorized chunk
     INTEGER  :: num_impurity_species
     REAL(rp), DIMENSION(:), ALLOCATABLE        :: Zj
     REAL(rp) :: gam_min
  END TYPE KORC_PARAMS


  TYPE, PUBLIC :: PARTICLES
     !! @note Derived type containing all the electrons' variables
     !!in the simulation. @endnote


     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: X 
     !! Cartesian coordinates of the electrons' position.
     !! dim(X) = (3,SPECIES::ppp).
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: V 
     !! Cartesian components of the electrons' velocity. dim(V) = dim(X).
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: Rgc 
     !! Cartesian coordinates of the electrons' guiding-center position.
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: Y 
     !! Coordinates of the electrons' position in cylindrical or toroidal
     !! coordinates.
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: Y0 
     !! Placeholder coordinates of the electrons' position in cylindrical
     !! coordinates for GC orbit model.
     REAL(rp), DIMENSION(:), ALLOCATABLE 	:: V0 
     !! Placeholder of the electrons' parallel momentum for the GC orbit model
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: E 
     !! Cartesian components of the electric field felt by each electron.
     !! dim(E) = dim(X).
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: B 
     !! Cartesian components of the magnetic field felt by each electron.
     !! dim(B) = dim(X).
     REAL(rp), DIMENSION(:), ALLOCATABLE 	:: PSI_P 
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: BR 
     !! Cartesian components of the gradient of the R-component of the
     !! magnetic field felt by each electron. dim(B) = dim(X).
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: BPHI 
     !! Cartesian components of the gradient of the PHI-component of the
     !! magnetic field felt by each electron. dim(B) = dim(X).
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: BZ 
     !! Cartesian components of the gradient of the Z-component of the
     !! magnetic field felt by each electron. dim(B) = dim(X).
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: gradB 
     !! Cylindrical components of the gradient of magnitude of magnetic
     !! field felt by each electron. dim(B) = dim(X).
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: curlb 
     !! Cylindrical components of the curl of the magnetic field unit
     !! vector felt by each electron. dim(B) = dim(X).
     REAL(rp), DIMENSION(:,:), ALLOCATABLE      :: RHS
     !! RHS of equations of motion for GC orbit model
     REAL(rp), DIMENSION(:,:), ALLOCATABLE      :: k1
     !! Cash-Karp Runge-Kutta coefficient for GC orbit model
     REAL(rp), DIMENSION(:,:), ALLOCATABLE      :: k2
     !! Cash-Karp Runge-Kutta coefficient for GC orbit model
     REAL(rp), DIMENSION(:,:), ALLOCATABLE      :: k3
     !! Cash-Karp Runge-Kutta coefficient for GC orbit model
     REAL(rp), DIMENSION(:,:), ALLOCATABLE      :: k4
     !! Cash-Karp Runge-Kutta coefficient for GC orbit model
     REAL(rp), DIMENSION(:,:), ALLOCATABLE      :: k5
     !! Cash-Karp Runge-Kutta coefficient for GC orbit model
     REAL(rp), DIMENSION(:,:), ALLOCATABLE      :: k6
     !! Cash-Karp Runge-Kutta coefficient for GC orbit model
     REAL(rp), DIMENSION(:), ALLOCATABLE 	:: ne 
     !! Electron density seen by each electron. dim(ne) = (1,SPECIES::ppp).
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: nimp
     REAL(rp), DIMENSION(:), ALLOCATABLE 	:: ni 
     REAL(rp), DIMENSION(:), ALLOCATABLE 	:: Te 
     !! Electron temperature seen by each electron. dim(Te) = (1,SPECIES::ppp).
     REAL(rp), DIMENSION(:), ALLOCATABLE 	:: Zeff 
     !! Zeff seen by each electron. dim(Zeff) = (1,SPECIES::ppp).
     REAL(rp), DIMENSION(:), ALLOCATABLE 	:: g 
     !! Instantaneous relativistic \(\gamma = 1/\sqrt{1 - v^2/c^2}\) factor
     !! of each electron in the simulation.
     REAL(rp), DIMENSION(:), ALLOCATABLE 	:: eta 
     !! Instantaneous pitch angle of each electron in the simulation.
     REAL(rp), DIMENSION(:), ALLOCATABLE 	:: mu 
     !! Magnetic moment of each electron in the simulation.
     REAL(rp), DIMENSION(:), ALLOCATABLE 	:: Prad 
     !! Instantaneous radiated power by each electron in the simulation.
     REAL(rp), DIMENSION(:), ALLOCATABLE 	:: Pin 
     !! Instantaneous input power of each electron due to the electric
     !! field acceleration.
     INTEGER(is), DIMENSION(:), ALLOCATABLE 	:: flagCon
     !! Flag for each particle to decide whether it is contained
     INTEGER(is), DIMENSION(:), ALLOCATABLE 	:: flagCol
     !! Flag for each particle to decide whether it is part of the high
     !! energy population
     INTEGER(is), DIMENSION(:), ALLOCATABLE 	:: flagRE 
     !! Flag for each particle to decide whether it is initialized
     REAL(rp), DIMENSION(:), ALLOCATABLE 	:: AUX 
     !! An auxiliary scalar variable for each electron.
     REAL(rp), DIMENSION(:), ALLOCATABLE 	:: wt 
     !! Weight of each electron. This is used when sampling weighted
     !! PDFs and in the synthetic camera diagnostic.
     INTEGER(is), DIMENSION(:), ALLOCATABLE 	:: initLCFS
#ifdef FIO
     TYPE(C_PTR), DIMENSION(:), ALLOCATABLE :: hint
     !! Hint for M3D_C1 interpolation.
#endif
     LOGICAL                                :: cart
  END TYPE PARTICLES


  TYPE, PUBLIC :: SPECIES
     !! @note Derived type containing the initial parameters of each electron
     !! ensemble followed in a KORC simulation. @endnote

     TYPE(PARTICLES) 			:: vars 
     !! An instance of the PARTICLES derived type.
     LOGICAL 				:: runaway 
     !! Flag to decide whether a given electron is a runaway (runaway=T)
     !! or not (runaway=F).
     CHARACTER(MAX_STRING_LENGTH) 	:: spatial_distribution 
     !! String describing the type of initial spatial distribution for
     !! each electron species.
     CHARACTER(MAX_STRING_LENGTH) 	:: energy_distribution 
     !! String describing the type of initial energy distribution for
     !! each electron species.
     CHARACTER(MAX_STRING_LENGTH) 	:: pitch_distribution 
     !! String describing the type of initial pitch-angle distribution for
     !! each electron species.
     REAL(rp) 				:: Eo 
     !! Initial energy of each electron species in case of using an initial
     !! mono-energetic distribution.
     REAL(rp) 				:: go 
     !! Corresponding relativisitc factor of each electron species in case
     !! of using an initial mono-energetic distribution.
     REAL(rp) 				:: etao 
     !! Initial pitch-angle of each electron species in case of using an
     !! initial mono-pitch-angle distribution.
     REAL(rp), DIMENSION(2) 		:: Eo_lims 
     !! Minimum and maximum energy limits of a given initial
     !! non-mono-energetic distribution.
     REAL(rp), DIMENSION(2) 		:: etao_lims 
     !! Minimum and maximum pitch-angle limits of a given initial
     !! non-mono-pitch-angle distribution.
     REAL(rp) 				:: wc 
     !! The mean electron cyclotron frequency of each electron species.
     REAL(rp) 				:: wc_r 
     !! The mean relativistic electron cyclotron frequency of each electron
     !! species.
     REAL(rp) 				:: q 
     !! Charge of each species. @note This was left explicit to allow 
     !! KORC to follow electrons and ions in the future. @endnote
     REAL(rp) 				:: m 
     !! Mass of each species. @note This was left explicit to allow KORC 
     !! to follow electrons and ions in the future. @endnote
     INTEGER 				:: ppp 
     !! Number of computational particles used to simulate each electron
     !! species.
     INTEGER 				:: pinit 
     !! Number of computational particles initialized for each electron species
     !! to give room for sources of additional electrons
     INTEGER 				:: pRE 
     REAL(rp) 				:: Ro 
     !! Radial position of the center of the electrons' initial spatial
     !! distribution.
     REAL(rp) 				:: PHIo 
     !! Azimuthal position of the electrons' initial spatial 
     !! distribution, in case of using a disk at a certain poloidal section.
     REAL(rp) 				:: Zo 
     !! Height of the center of the electrons' initial spatial distribution.
     REAL(rp) 				:: r_inner 
     !! Minimum minor radius of the electrons' initial spatial distribution.
     REAL(rp) 				:: r_outter 
     !! Maximum minor radius of the electrons' initial spatial distribution.
     REAL(rp) 				:: falloff_rate 
     !! Exponential falloff or standard deviation of a non-uniform radial
     !! distribution of electrons.
     REAL(rp) 				:: shear_factor 
     !! Shear factor used to generate an initial spatial distribution with an
     !! elliptic poloidal cross section. 
     !! @note See <em>Carbajal and del-Castillo-Negrete, Nuclear Fusion,
     !! submitted (2018)</em>. @endnote
     REAL(rp)                            :: sigmaR
     !! Variance of the first dimension of a 2D Gaussian, spatial
     !! distribution function
     REAL(rp)                            :: sigmaZ
     !! Variance of the second dimension of a 2D Gaussian, spatial
     !! distribution function
     REAL(rp)                            :: theta_gauss
     !! Angle of counter-clockwise rotation (in degrees) of 2D Gaussian
     !! distribution relative to R,Z
     REAL(rp)                            :: psi_max
     !! Maximum value of the argument of the 2D gaussian exponential, used
     !! for an indicator function that limits the region of MH sampling
     REAL(rp)                            :: Spong_b
     REAL(rp)                            :: Spong_w
     REAL(rp)                            :: Spong_dlam
     REAL(rp)                            :: dth
     REAL(rp)                            :: dgam
     REAL(rp)                            :: dR
     REAL(rp)                            :: dZ
     REAL(rp),DIMENSION(3) :: Xtrace
     !! Initial position in Cartesian coordinates for tracer particle
  END TYPE SPECIES


  TYPE, PRIVATE :: A_FIELD
     !! @note Derived type having all the parameters of the analytical
     !! magnetic field included in KORC. @endnote
     !! The analytical magnetic field is given by:
     !! $$\mathbf{B}(r,\vartheta) = \frac{1}{1 + \eta \cos{\vartheta}}
     !! \left[ B_0 \hat{e}_\zeta  + B_\vartheta(r) \hat{e}_\vartheta \right],$$
     !! where \(\eta = r/R_0\) is the aspect ratio, the constant \(B_0\)
     !! denotes the magnitude of the toroidal magnetic field, 
     !! and \(B_\vartheta(r) = \eta B_0/q(r)\) is the poloidal magnetic
     !! field with safety factor 
     !! \(q(r) = q_0\left( 1 + \frac{r^2}{\lambda^2} \right)\). The
     !! constant \(q_0\) is the safety factor at 
     !! the magnetic axis and the constant \(\lambda\) is obtained from
     !! the values of \(q_0\) and \(q(r)\) 
     !! at the plasma edge \(r=r_{edge}\).

     REAL(rp) 			:: Bo 
     !! Magnitude of the toroidal magnetic field \(B_0\).
     REAL(rp) 			:: a 
     !! Plasma edge \(r_{edge}\) as measured from the magnetic axis.
     REAL(rp) 			:: Ro 
     !! Radial position of the magnetic axis \(R_0\)
     REAL(rp) 			:: qa 
     !! Safety factor at the plasma edge.
     REAL(rp) 			:: qo 
     !! Safety factor at the magnetic axis \(q_0\).
     REAL(rp) 			:: lambda 
     !! \(\lambda\) parameter of \(q(r)\).
     REAL(rp) 			:: Bpo 
     !! @deprecated Parameter not used anymore. @todo Delete parameter.
     REAL(rp) 			:: Bp_sign 
     !! Sign of \(B_\vartheta(r)\). This depends on current_direction,
     !! Bp_sign=1 for 
     !! current_direction='PARALLEL', and Bp_sign=-1 for
     !! current_direction='ANTI-PARALLEL'.
     CHARACTER(MAX_STRING_LENGTH) :: current_direction 
     !! Direction of plasma current: PARALLEL or ANTI-PARALLEL to the
     !! toroidal magnetic field.
     REAL(rp) 			:: Ero 
     !! radial electric field amplitude.
     REAL(rp) 			:: rmn
     !! mode location rmn
     REAL(rp)			:: sigmamn
     !! mode width sigmamn
  END TYPE A_FIELD

  TYPE, PRIVATE :: MESH
     !! Derived type with the cylindrical coordinates of the grid nodes 
     !! at which the pre-computed plasma profiles and fields are known.

     REAL(rp), DIMENSION(:), ALLOCATABLE :: R 
     !! Radial grid.
     REAL(rp), DIMENSION(:), ALLOCATABLE :: PHI 
     !! Azimuthal grid.
     REAL(rp), DIMENSION(:), ALLOCATABLE :: Z 
     !! Z grid.
  END TYPE MESH


  TYPE, PUBLIC :: FIELDS
     !! @note Derived type with all the variables and data of analytical 
     !! and pre-computed electric and magnetic fields. @endnote

     TYPE(A_FIELD) 	 			:: AB 
     !! An instance of the KORC derived data type A_FIELD.
     TYPE(V_FIELD_3D) 				:: E_3D 
     !! KORC 3-D vector field of the pre-computed electric field.
     TYPE(V_FIELD_3D) 				:: B_3D
     TYPE(V_FIELD_3D) 				:: dBdR_3D
     TYPE(V_FIELD_3D) 				:: dBdPHI_3D
     TYPE(V_FIELD_3D) 				:: dBdZ_3D 
     !! KORC 3-D vector field of the pre-computed magnetic field.
     TYPE(V_FIELD_2D) 				:: E_2D 
     !! KORC 2-D vector field of the pre-computed electric field.
     TYPE(V_FIELD_2D) 				:: B_2D
     TYPE(V_FIELD_2D) 				:: B1Re_2D
     TYPE(V_FIELD_2D) 				:: B1Im_2D
     TYPE(V_FIELD_2DX) 				:: B1Re_2DX
     TYPE(V_FIELD_2DX) 				:: E1Im_2DX
     TYPE(V_FIELD_2DX) 				:: E1Re_2DX
     TYPE(V_FIELD_2DX) 				:: B1Im_2DX
     TYPE(V_FIELD_2D) 				:: dBdR_2D
     TYPE(V_FIELD_2D) 				:: dBdPHI_2D
     TYPE(V_FIELD_2D) 				:: dBdZ_2D 
     !! KORC 3-D vector field of the pre-computed magnetic field.
     TYPE(V_FIELD_2D) 				:: gradB_2D
     TYPE(V_FIELD_3D) 				:: gradB_3D 
     !! KORC 3-D vector field of the gradient of the magnitude of the
     !! pre-computed magnetic field.
     TYPE(V_FIELD_2D) 				:: curlb_2D
     TYPE(V_FIELD_3D) 				:: curlb_3D 
     !! KORC 3-D vector field of the curl of the unit vector in the
     !! direction of the pre-computed magnetic field.
     TYPE(V_FIELD_1D) 				:: E_SC_1D
     TYPE(V_FIELD_1D) 				:: J0_SC_1D
     TYPE(V_FIELD_1D) 				:: J1_SC_1D
     TYPE(V_FIELD_1D) 				:: J2_SC_1D
     TYPE(V_FIELD_1D) 				:: J3_SC_1D
     TYPE(V_FIELD_1D) 				:: A1_SC_1D
     TYPE(V_FIELD_1D) 				:: A2_SC_1D
     TYPE(V_FIELD_1D) 				:: A3_SC_1D
     
     REAL(rp), DIMENSION(:), ALLOCATABLE :: r_1D
     REAL(rp), DIMENSION(:), ALLOCATABLE :: PSIP_1D
     REAL(rp), DIMENSION(:), ALLOCATABLE :: dMagPsiSqdPsiP
     REAL(rp), DIMENSION(:), ALLOCATABLE :: ddMagPsiSqdPsiPSq
     TYPE(MESH) 		 		:: X 
     !! An instance of the KORC derived type MESH.
     CHARACTER(MAX_STRING_LENGTH) :: E_model
     !! Name for dynamical, analytic, electric field model to be added to
     REAL(rp)  :: E_dyn
     REAL(rp)  :: E_pulse
     REAL(rp)  :: E_width
     CHARACTER(30) :: E_profile
     REAL(rp)  :: PSIP_min
     REAL(rp)  :: PSIp_lim,PSIp_0
     REAL(rp)  :: AMP
     REAL(rp)  :: MARS_AMP_Scale
     REAL(rp)  :: AORSA_AMP_Scale
     REAL(rp)  :: AORSA_freq
     REAL(rp)  :: AORSA_nmode
     !! interpolated E field
     CHARACTER(30) :: Analytic_IWL
     INTEGER :: ntiles
     REAL(rp) :: circumradius
     INTEGER 			:: res_double
     INTEGER, DIMENSION(3) 			:: dims 
     !! Dimensions of the KORC vector field. dims=(number of grid 
     !! nodes along \(R\), number of grid nodes along \(\phi\), 
     !! number of grid nodes along \(Z\)).
     INTEGER 			:: dim_1D
     INTEGER 			:: subcycle_E_SC
     REAL(rp)  :: dt_E_SC,Ip_exp,Ip0
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: PSIp
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: PSIp_FS
     REAL(rp), DIMENSION(:,:,:), ALLOCATABLE 	:: PSIp3D 
     !! 2-D array for storing the data of the poloidal magnetic flux.
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: FLAG2D 
     !! 2-D array defining the simulation domain where pre-computed data exist.
     REAL(rp), DIMENSION(:,:,:), ALLOCATABLE      :: FLAG3D 
     !! 3-D array defining the simulation domain where pre-computed data exist.
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: LCFS2D
     REAL(rp), DIMENSION(:,:,:), ALLOCATABLE    :: LCFS3D 
     REAL(rp) 					:: Eo 
     !! Characteristic electric field.
     REAL(rp) 					:: Bo 
     !! Characteristic magnetic field.
     REAL(rp) 					:: Ro 
     !! Radial position of the magnetic axis.
     REAL(rp) 					:: Zo 
     !! \(Z\) position of the magnetic axis.
     LOGICAL 					:: Bfield 
     !! Flag to indicate whether a pre-computed magnetic field will be
     !! used (Bfield=T) or not (Bfield=F).
     LOGICAL 					:: dBfield 
     !! Flag to indicate whether a pre-computed magnetic field will be
     !! used (Bfield=T) or not (Bfield=F).
     LOGICAL 					:: Bflux
     LOGICAL 					:: Bflux3D
     LOGICAL 					:: B1field
     LOGICAL 					:: E1field
     !! Flag to indicate whether a pre-computed poloidal magnetic flux will
     !! be used (Bflux=T) or not (Bflux=F).
     LOGICAL 					:: Efield 
     !! Flag to indicate whether a pre-computed electric field will be used
     !! (Efield=T) or not (Efield=F).
     LOGICAL 					:: Bfield_in_file 
     !! Flag to indicate if a pre-computed magnetic field is in the input file.
     LOGICAL 					:: dBfield_in_file
     LOGICAL 					:: B1field_in_file
     LOGICAL 					:: E1field_in_file 
     !! Flag to indicate if a pre-computed magnetic field is in the input file.
     LOGICAL 					:: Bflux_in_file 
     !! Flag to indicate if a pre-computed poloidal magnetic flux is in the
     !! input file.
     LOGICAL 					:: Efield_in_file 
     !! Flag to indicate if a pre-computed electric field is in the input file.
     LOGICAL 					:: axisymmetric_fields 
     !! Flag to indicate if the pre-computed fields are axisymmetric.
     LOGICAL 					:: Dim2x1t
     LOGICAL 					:: E_2x1t,ReInterp_2x1t
     REAL(rp)  :: t0_2x1t
     INTEGER  :: ind0_2x1t,ind_2x1t
#ifdef FIO
     INTEGER  :: isrc
     INTEGER (C_INT)                         :: FIO_B
     !! An M3D-C1 magnetic field.
     INTEGER (C_INT)                         :: FIO_E
     !! An M3D-C1 Electric field.
     INTEGER (C_INT)                         :: FIO_A
     !! An M3D-C1 vector potential.
     real(c_double)  ::  time0,time1
#endif
     REAL(rp)  :: psip_conv
     LOGICAL :: useLCFS
     
  END TYPE FIELDS


  TYPE, PUBLIC :: PROFILES
     !! @note KORC derived data type having information about the plasma
     !! profiles. 
     !! See [[korc_profiles.f90("file")]] for more information. @endnote
     !! KORC can run using either analytical and pre-computed plasma profiles. 
     !! Pre-computed plasma profiles, 
     !! as in the case of pre-computed electric or magnetic fields, are
     !! interpolated 
     !! to electrons' position in [[korc_profiles]].
     !!
     !! There are two types of analytical plsama profiles that can be used
     !! in KORC: 
     !! 3rd degree polynomial radial plasma profiles,
     !! $$f(r) = a_3r^3 + a_2r^2 +a_1r + a_0,$$
     !! and radial plasma profiles with a \(\tanh(r)\) dependency:
     !! $$f(r) = f_0\left[1 - \tanh^n\left(\frac{2r}{a}\right)\right],$$
     !! where \(f_0\) is a given plasma parameter at the magnetic axis,
     !! and \(a\) is 
     !! the plasma radius as measured 
     !! from the magnetic axis to the last closed flux surface. Notice that the 
     !! larger \(n\) is, the more uniform the radial profiles are.

     TYPE(MESH) 				        :: X 
     !! An instance of the KORC derived data type MESH.
     REAL(rp) 				        :: a 
     !! Plasma radius as measured from the magnetic axis
     REAL(rp) 				        :: R0
     REAL(rp) 				        :: Z0

     REAL(rp) 				        :: R0_RE
     REAL(rp) 				        :: Z0_RE
     
     INTEGER, DIMENSION(3) 			:: dims 
     !! Dimensions of the arrays containing the pre-computed profiles data. dims=(number of grid nodes along \(R\), 
     !! number of grid nodes along \(\phi\), number of grid nodes along \(Z\)).
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: FLAG2D 
     !! 2-D array defining the simulation domain where pre-computed data exist.
     REAL(rp), DIMENSION(:,:,:), ALLOCATABLE      :: FLAG3D 
     !! 3-D array defining the simulation domain where pre-computed data exist.

     REAL(rp) 					:: n_ne 
     !! \(n\) used in \(\tanh^n(r)\) of the electron density profile.
     REAL(rp) 					:: n_Te 
     !! \(n\) used in \(\tanh^n(r)\) of the electron temperature profile.
     REAL(rp) 					:: n_Zeff 
     !! \(n\) used in \(\tanh^n(r)\) of the \(Z_{eff}\) profile.

     REAL(rp)  ::  n_REr0=0._rp
     REAL(rp)  ::  n_tauion=0._rp
     REAL(rp)  ::  n_tauin=0._rp
     REAL(rp)  ::  n_tauout=0._rp
     REAL(rp)  ::  n_shelfdelay=0._rp
     REAL(rp)  ::  n_lamfront=0._rp
     REAL(rp)  ::  n_lamback=0._rp
     REAL(rp)  ::  n_lamshelf=0._rp
     REAL(rp)  ::  n_shelf=0._rp
     REAL(rp)  ::  psiN_0=1._rp

     
     REAL(rp), DIMENSION(4) 			:: a_ne 
     !! Coefficients of the polynomial electron density profile. 
     !! See detailed description above, a_ne=(\(a_{0}\),\(a_{2}\),\(a_{3}\),\(a_{4}\)).
     REAL(rp), DIMENSION(4) 			:: a_Te 
     !! Coefficients of the polynomial electron temperature profile. 
     !! See detailed description above, a_ne=(\(a_{0}\),\(a_{2}\),\(a_{3}\),\(a_{4}\)).
     REAL(rp), DIMENSION(4) 			:: a_Zeff 
     !! Coefficients of the \(Z_{eff}\) profile. 
     !! See detailed description above, a_ne=(\(a_{0}\),\(a_{2}\),\(a_{3}\),\(a_{4}\)).

     ! Zeff
     CHARACTER(MAX_STRING_LENGTH) 		:: Zeff_profile 
     !! String containing the type of \(Z_{eff}\) profile to be used in the simulation.
     REAL(rp) 					:: Zeffo 
     !! \(Z_{eff}\) at the magnetic axis.
     REAL(rp), DIMENSION(:,:,:), ALLOCATABLE      :: Zeff_3D 
     !! 3-D array for keeping the pre-computed data of the \(Z_{eff}\) profile.
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: Zeff_2D 
     !! 2-D array for keeping the pre-computed data of the \(Z_{eff}\) profile.

     ! Density
     CHARACTER(MAX_STRING_LENGTH) 		:: ne_profile 
     !! String containing the type of electron density profile to be used in the simulation.
     REAL(rp) 					:: neo 
     !! Electron density at the magnetic axis
     REAL(rp), DIMENSION(:,:,:), ALLOCATABLE      :: ne_3D 
     !! 3-D array for keeping the pre-computed data of the electron density profile.
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: ne_2D 
     !! 2-D array for keeping the pre-computed data of the electron density profile.

     !Temperature
     CHARACTER(MAX_STRING_LENGTH) 		:: Te_profile 
     !! String containing the type of electron temperature profile to be used in the simulation.
     REAL(rp) 					:: Teo 
     !! Electron temperature at the magnetic axis
     REAL(rp), DIMENSION(:,:,:), ALLOCATABLE      :: Te_3D 
     !! 3-D array for keeping the pre-computed data of the electron density profile.
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: Te_2D 
     !! 2-D array for keeping the pre-computed data of the electron density profile.

     CHARACTER(MAX_STRING_LENGTH) 		:: filename 
     !! Full path to the HDF5 file containing the pre-computed plasma profiles.
     LOGICAL 					:: axisymmetric 
     !! Flag to indicate if the plasma profiles are axisymmetric.
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: RHON
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: nRE_2D
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: nAr0_2D
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: nAr1_2D
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: nAr2_2D
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: nAr3_2D
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: nD_2D
     REAL(rp), DIMENSION(:,:), ALLOCATABLE 	:: nD1_2D 
     
#ifdef FIO
     INTEGER (C_INT)                         :: FIO_ne
     INTEGER (C_INT)                         :: FIO_ni
     INTEGER (C_INT)                         :: FIO_te
     INTEGER (C_INT), DIMENSION(:), ALLOCATABLE  :: FIO_nimp
     INTEGER (C_INT)                         :: FIO_zeff
#endif
  END TYPE PROFILES

end module korc_types
