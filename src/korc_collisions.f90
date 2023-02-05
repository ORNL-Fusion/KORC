module korc_collisions
  use korc_types
  use korc_constants
  use korc_HDF5
  use korc_interp
  use korc_profiles
  use korc_fields
  use korc_input

#ifdef FIO
  use korc_fio
#endif

#ifdef PARALLEL_RANDOM
  use korc_random
#endif

  IMPLICIT NONE

  CHARACTER(LEN=*), PRIVATE, PARAMETER 	:: MODEL1 = 'SINGLE_SPECIES'
  CHARACTER(LEN=*), PRIVATE, PARAMETER 	:: MODEL2 = 'MULTIPLE_SPECIES'
  REAL(rp), PRIVATE, PARAMETER 			:: infinity = HUGE(1.0_rp)

  TYPE, PRIVATE :: PARAMS_MS
     INTEGER 					:: num_impurity_species
     REAL(rp) 					:: Te
     ! Background electron temperature in eV
     REAL(rp) 					:: ne
     ! Background electron density in 1/m^3
     REAL(rp) 					:: nH
     ! Background proton density in 1/m^3
     REAL(rp) 					:: nef
     ! Free electron density in 1/m^3
     REAL(rp), DIMENSION(:), ALLOCATABLE        :: neb
     ! Bound electron density in 1/m^3
     REAL(rp), DIMENSION(:), ALLOCATABLE        :: Zi
     ! Atomic number of (majority) background ions
     REAL(rp), DIMENSION(:), ALLOCATABLE        :: Zo
     ! Full nuclear charge of each impurity: Z=1 for D, Z=10 for Ne
     REAL(rp), DIMENSION(:), ALLOCATABLE        :: Zj
     ! Atomic number of each impurity: Z=1 for D, Z=10 for Ne
     REAL(rp), DIMENSION(:), ALLOCATABLE        :: nz
     ! Impurity densities
     REAL(rp), DIMENSION(:), ALLOCATABLE        :: IZj,aZj
     ! Ionization energy of impurity in eV
     REAL(rp), DIMENSION(:), ALLOCATABLE        :: Ee_IZj
     ! me*c^2/IZj dimensionless parameter

     REAL(rp) 					:: rD
     ! Debye length
     REAL(rp) 					:: re
     ! Classical electron radius
     REAL(rp) 					:: Gammac_min


     REAL(rp), DIMENSION(2) :: aH=(/274._rp,0._rp/)
     REAL(rp), DIMENSION(7) :: aC=(/144._rp,118._rp,95._rp,70._rp, &
          42._rp,39._rp,0._rp/)
     REAL(rp), DIMENSION(11) :: aNe=(/111._rp,100._rp,90._rp,80._rp, &
          71._rp,62._rp,52._rp,40._rp,24._rp,23._rp,0._rp/)
     REAL(rp), DIMENSION(19) :: aAr=(/96._rp,90._rp,84._rp,78._rp,72._rp, &
          65._rp,59._rp,53._rp,47._rp,44._rp,41._rp,38._rp,25._rp,32._rp, &
          27._rp,21._rp,13._rp,13._rp,0._rp/)


     REAL(rp), DIMENSION(2) :: IH=(/14.9916_rp,huge(1._rp)/)
     REAL(rp), DIMENSION(7) :: IC=(/65.9_rp,92.6_rp,134.8_rp,214.2_rp, &
          486.2_rp,539.5_rp,huge(1._rp)/)
     REAL(rp), DIMENSION(11) :: INe=(/137.2_rp,165.2_rp,196.9_rp,235.2_rp, &
          282.8_rp,352.6_rp,475.0_rp,696.8_rp,1409.2_rp,1498.4_rp,huge(1._rp)/)
     REAL(rp), DIMENSION(19) :: IAr=(/188.5_rp,219.4_rp,253.8_rp,293.4_rp, &
          339.1_rp,394.5_rp,463.4_rp,568.0_rp,728.0_rp,795.9_rp,879.8_rp, &
          989.9_rp,1138.1_rp,1369.5_rp,1791.2_rp,2497.0_rp,4677.2_rp, &
          4838.2_rp,huge(1._rp)/)

     CHARACTER(30) :: neut_prof
     REAL(rp)  :: neut_edge_fac
     REAL(rp) 			:: Ec,Ec_min
     ! Critical electric field
     LOGICAL  :: LargeCollisions
     LOGICAL :: lowKE_REs
     REAL(rp)  :: lowKE_LAC_not_ionized

  END TYPE PARAMS_MS

  TYPE, PRIVATE :: PARAMS_SS
     REAL(rp) 			:: Te
     ! Electron temperature
     REAL(rp) 			:: Ti
     ! Ion temperature
     REAL(rp) 			:: ne
     ! Background electron density
     REAL(rp) 			:: Zeff
     ! Effective atomic number of ions
     REAL(rp) 			:: rD
     ! Debye radius
     REAL(rp) 			:: re
     ! Classical electron radius
     REAL(rp) 			:: CoulombLogee,CoulombLogei
     ! Coulomb logarithm
     REAL(rp) 			:: CLog1, CLog2,CLog0_1, CLog0_2
     REAL(rp) 			:: VTe
     ! Thermal velocity of background electrons
     REAL(rp) 			:: VTeo
     REAL(rp) 			:: delta
     ! delta parameter
     REAL(rp) 			:: deltao
     REAL(rp) 			:: Gammac
     ! Collisional Gamma factor
     REAL(rp) 			:: Gammaco
     ! Collisional gamma factor normalized for SDE for dp
     REAL(rp) 			:: Tau
     ! Collisional time of relativistic particles
     REAL(rp) 			:: Tauc
     ! Collisional time of thermal particles
     REAL(rp) 			:: taur
     ! radiation timescale
     REAL(rp) 			:: Ec
     ! Critical electric field
     REAL(rp) 			:: ED
     ! Dreicer electric field
     REAL(rp) 			:: dTau
     ! Subcycling time step in collisional time units (Tau)
     INTEGER(ip)		:: subcycling_iterations,ngrid1
     REAL(rp) :: coll_per_dump_dt,Clog_const
     REAL(rp) :: p_min,p_crit,p_therm,gam_min,gam_crit,gam_therm,pmin_scale
     LOGICAL :: ConserveLA,sample_test,avalanche,energy_diffusion,FP_bremsstrahlung,pitch_diffusion,always_aval
     CHARACTER(30) :: Clog_model,min_secRE,LAC_gam_resolution

     REAL(rp), DIMENSION(3) 	:: x = (/1.0_rp,0.0_rp,0.0_rp/)
     REAL(rp), DIMENSION(3) 	:: y = (/0.0_rp,1.0_rp,0.0_rp/)
     REAL(rp), DIMENSION(3) 	:: z = (/0.0_rp,0.0_rp,1.0_rp/)

     TYPE(PROFILES) 			   :: P

     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: rnd_num
     INTEGER 				   :: rnd_num_count
     INTEGER 				   :: rnd_dim = 40000000_idef

     LOGICAL :: slowing_down

  END TYPE PARAMS_SS

  TYPE(PARAMS_MS), PRIVATE :: cparams_ms
  TYPE(PARAMS_SS), PRIVATE :: cparams_ss

  PUBLIC :: initialize_collision_params,&
       normalize_collisions_params,&
       collision_force,&
       deallocate_collisions_params,&
       save_collision_params,&
       include_CoulombCollisions_GC_p,&
       include_CoulombCollisionsLA_GC_p,&
       include_CoulombCollisions_FO_p,&
       check_collisions_params,&
       define_collisions_time_step,&
       large_angle_source
  PRIVATE :: load_params_ms,&
       load_params_ss,&
       normalize_params_ms,&
       normalize_params_ss,&
       save_params_ms,&
       save_params_ss,&
       deallocate_params_ms,&
       cross,&
       CA,&
       CB_ee,&
       CB_ei,&
       CB_ei_FIO,&
       CF,&
       CF_FIO,&
       fun,&
       nu_S,&
       nu_S_FIO,&
       nu_par,&
       nu_D,&
       nu_D_FIO,&
       Gammac_wu,&
       CLog_wu,&
       VTe_wu,&
       Gammacee,&
       CLog,&
       VTe,&
       CA_SD,&
       CB_ee_SD,&
       CB_ei_SD,&
       CF_SD,&
       delta,&
       unitVectorsC,&
       unitVectors_p

contains

  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !
  ! * SUBROUTINES FOR INITIALIZING COLLISIONS PARAMS * !
  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !


  subroutine load_params_ms(params)
    TYPE(KORC_PARAMS), INTENT(INOUT) 	:: params
    !REAL(rp) 				:: Te
    ! Background electron temperature in eV
    !REAL(rp) 				:: ne
    ! Background electron density in 1/m^3
    !INTEGER 				:: num_impurity_species
    !REAL(rp), DIMENSION(10) 		:: Zo
    ! Full nuclear charge of each impurity: Z=1 for D, Z=10 for Ne
    !REAL(rp), DIMENSION(10) 		:: Zj
    ! Atomic number of each impurity: Z=1 for D, Z=10 for Ne
    !REAL(rp), DIMENSION(10) 		:: nz
    ! Impurity densities
    !REAL(rp), DIMENSION(10) 		:: IZj
    ! Ionization energy of impurity in eV
    !REAL(rp), DIMENSION(10) 		:: aZj
    INTEGER :: i

    !NAMELIST /CollisionParamsMultipleSpecies/ num_impurity_species,Te,ne, &
    !     Zo,Zj,nz,IZj


    !open(unit=output_unit_write,file=TRIM(params%path_to_inputs), &
    !     status='OLD',form='formatted')
    !read(output_unit_write,nml=CollisionParamsMultipleSpecies)
    !close(output_unit_write)

    if (params%profile_model(10:10).eq.'H') then
       cparams_ms%num_impurity_species = 6
       params%num_impurity_species = 6
    else
       cparams_ms%num_impurity_species = num_impurity_species
       params%num_impurity_species = num_impurity_species
    endif


    ALLOCATE(cparams_ms%Zj(cparams_ms%num_impurity_species))
    ALLOCATE(cparams_ms%Zo(cparams_ms%num_impurity_species))
    ALLOCATE(cparams_ms%nz(cparams_ms%num_impurity_species))
    ALLOCATE(cparams_ms%neb(cparams_ms%num_impurity_species))
    ALLOCATE(cparams_ms%IZj(cparams_ms%num_impurity_species))
    ALLOCATE(cparams_ms%aZj(cparams_ms%num_impurity_species))
    ALLOCATE(cparams_ms%Ee_IZj(cparams_ms%num_impurity_species))

    cparams_ms%Te = Te_mult*C_E
    cparams_ms%ne = ne_mult
    cparams_ms%nH = ne_mult

    if  (params%profile_model.eq.'M3D_C1') then
       do i=1,cparams_ms%num_impurity_species
          cparams_ms%Zj(i) = real(i)-1._rp
          cparams_ms%Zo(i) = Zo_mult(1)
       end do
       cparams_ms%nz(1) = nz_mult(1)

       params%Zj=cparams_ms%Zj
    elseif (params%profile_model(10:10).eq.'H') then
       cparams_ms%Zj(1)=0
       cparams_ms%Zo(1)=18
       cparams_ms%Zj(2)=1
       cparams_ms%Zo(2)=18
       cparams_ms%Zj(3)=2
       cparams_ms%Zo(3)=18
       cparams_ms%Zj(4)=3
       cparams_ms%Zo(4)=18
       cparams_ms%Zj(5)=0
       cparams_ms%Zo(5)=1
       cparams_ms%Zj(6)=1
       cparams_ms%Zo(6)=1

       cparams_ms%nz(1) = nz_mult(1)

       params%Zj=cparams_ms%Zj
    else
       cparams_ms%Zj = Zj_mult(1:cparams_ms%num_impurity_species)
       cparams_ms%Zo = Zo_mult(1:cparams_ms%num_impurity_species)
       cparams_ms%nz = nz_mult(1:cparams_ms%num_impurity_species)
    endif

    if (.not.(params%profile_model(10:10).eq.'H')) then
       do i=1,cparams_ms%num_impurity_species
          if (int(cparams_ms%Zo(i)).eq.1) then
             cparams_ms%IZj(i) = C_E*cparams_ms%IH(int(cparams_ms%Zj(i)+1))
             cparams_ms%aZj(i) = cparams_ms%aH(int(cparams_ms%Zj(i)+1))
          else if (int(cparams_ms%Zo(i)).eq.6) then
             cparams_ms%IZj(i) = C_E*cparams_ms%IC(int(cparams_ms%Zj(i)+1))
             cparams_ms%aZj(i) = cparams_ms%aC(int(cparams_ms%Zj(i)+1))
          else if (int(cparams_ms%Zo(i)).eq.10) then
             cparams_ms%IZj(i) = C_E*cparams_ms%INe(int(cparams_ms%Zj(i)+1))
             cparams_ms%aZj(i) = cparams_ms%aNe(int(cparams_ms%Zj(i)+1))
          else if (int(cparams_ms%Zo(i)).eq.18) then
             cparams_ms%IZj(i) = C_E*cparams_ms%IAr(int(cparams_ms%Zj(i)+1))
             cparams_ms%aZj(i) = cparams_ms%aAr(int(cparams_ms%Zj(i)+1))
          else
             if (params%mpi_params%rank .EQ. 0) then
                write(output_unit_write,'("Atomic number not defined!")')
             end if
             exit
          end if
       end do
    else
       cparams_ms%IZj(1) = C_E*cparams_ms%IAr(1)
       cparams_ms%aZj(1) = cparams_ms%aAr(1)
       cparams_ms%IZj(2) = C_E*cparams_ms%IAr(2)
       cparams_ms%aZj(2) = cparams_ms%aAr(2)
       cparams_ms%IZj(3) = C_E*cparams_ms%IAr(3)
       cparams_ms%aZj(3) = cparams_ms%aAr(3)
       cparams_ms%IZj(4) = C_E*cparams_ms%IAr(4)
       cparams_ms%aZj(4) = cparams_ms%aAr(4)
       cparams_ms%IZj(5) = C_E*cparams_ms%IH(1)
       cparams_ms%aZj(5) = cparams_ms%aH(1)
       cparams_ms%IZj(6) = C_E*cparams_ms%IH(2)
       cparams_ms%aZj(6) = cparams_ms%aH(2)
    endif

    cparams_ms%nef = ne_mult + sum(cparams_ms%Zj*cparams_ms%nz)
    cparams_ms%neb = (cparams_ms%Zo-cparams_ms%Zj)*cparams_ms%nz

    cparams_ms%rD = SQRT( C_E0*cparams_ms%Te/(cparams_ms%ne*C_E**2) )
    cparams_ms%re = C_RE
    cparams_ms%Ee_IZj = C_ME*C_C**2/cparams_ms%IZj

    cparams_ms%neut_prof=neut_prof
    cparams_ms%neut_edge_fac=neut_edge_fac
    cparams_ms%lowKE_REs=lowKE_REs
    cparams_ms%lowKE_LAC_not_ionized=lowKE_LAC_not_ionized

    cparams_ms%Gammac_min = Gammac_wu(params,cparams_ss%P%n_ne,cparams_ss%Te)
    cparams_ms%LargeCollisions = LargeCollisions

    if (params%mpi_params%rank .EQ. 0) then
       write(output_unit_write,'("Number of impurity species: ",I16)')&
            cparams_ms%num_impurity_species
       do i=1,cparams_ms%num_impurity_species
          if (cparams_ms%Zo(i).eq.1) then
             write(output_unit_write,'("H with charge state: ",I16)') int(cparams_ms%Zj(i))
             write(output_unit_write,'("Mean excitation energy I (eV)",E17.10)') &
                  cparams_ms%IZj(i)/C_E
             write(output_unit_write,'("Effective ion length scale a (a_0)",E17.10)') &
                  cparams_ms%aZj(i)
          else if (cparams_ms%Zo(i).eq.6) then
             write(output_unit_write,'("C with charge state: ",I16)') int(cparams_ms%Zj(i))
             write(output_unit_write,'("Mean excitation energy I (eV)",E17.10)') &
                  cparams_ms%IZj(i)/C_E
             write(output_unit_write,'("Effective ion length scale a (a_0)",E17.10)') &
                  cparams_ms%aZj(i)
          else if (cparams_ms%Zo(i).eq.10) then
             write(output_unit_write,'("Ne with charge state: ",I16)') int(cparams_ms%Zj(i))
             write(output_unit_write,'("Mean excitation energy I (eV)",E17.10)') &
                  cparams_ms%IZj(i)/C_E
             write(output_unit_write,'("Effective ion length scale a (a_0)",E17.10)') &
                  cparams_ms%aZj(i)
          else if (cparams_ms%Zo(i).eq.18) then
             write(output_unit_write,'("Ar with charge state: ",I16)') int(cparams_ms%Zj(i))
             write(output_unit_write,'("Mean excitation energy I (eV)",E17.10)') &
                  cparams_ms%IZj(i)/C_E
             write(output_unit_write,'("Effective ion length scale a (a_0)",E17.10)') &
                  cparams_ms%aZj(i)
          end if
       end do
    end if

  end subroutine load_params_ms


  subroutine load_params_ss(params)
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    !REAL(rp) 				:: Te
    ! Electron temperature
    !REAL(rp) 				:: Ti
    ! Ion temperature
    !REAL(rp) 				:: ne
    ! Background electron density
    !REAL(rp) 				:: Zeff
    ! Effective atomic number of ions
    !REAL(rp) 				:: dTau
    ! Subcycling time step in collisional time units (Tau)
    !CHARACTER(MAX_STRING_LENGTH) 	:: ne_profile
    !CHARACTER(MAX_STRING_LENGTH) 	:: Te_profile
    !CHARACTER(MAX_STRING_LENGTH) 	:: Zeff_profile
    !CHARACTER(MAX_STRING_LENGTH) 	:: filename
    !REAL(rp) 				:: radius_profile
    !REAL(rp) 				:: neo
    !REAL(rp) 				:: Teo
    !REAL(rp) 				:: Zeffo
    !REAL(rp) 				:: n_ne
    !REAL(rp) 				:: n_Te
    !REAL(rp) 				:: n_Zeff
    !REAL(rp), DIMENSION(4) 		:: a_ne
    !REAL(rp), DIMENSION(4) 		:: a_Te
    !REAL(rp), DIMENSION(4) 		:: a_Zeff
    !LOGICAL 				:: axisymmetric
    !REAL(rp)  ::  n_REr0
    !REAL(rp)  ::  n_tauion
    !REAL(rp)  ::  n_lamfront,psiN_0
    !REAL(rp)  ::  n_lamback,n_lamshelf,n_shelfdelay,n_tauin,n_tauout,n_shelf

    !NAMELIST /CollisionParamsSingleSpecies/ Te, Ti, ne, Zeff, dTau

    !NAMELIST /plasmaProfiles/ radius_profile,ne_profile,neo,n_ne,a_ne,&
    !     Te_profile,Teo,n_Te,a_Te,n_REr0,n_tauion,n_lamfront,n_lamback, &
    !     Zeff_profile,Zeffo,n_Zeff,a_Zeff,filename,axisymmetric, &
    !     n_lamshelf,n_shelfdelay,n_tauin,n_tauout,n_shelf,psiN_0


    !open(unit=output_unit_write,file=TRIM(params%path_to_inputs), &
    !     status='OLD',form='formatted')
    !read(output_unit_write,nml=CollisionParamsSingleSpecies)
    !close(output_unit_write)

    cparams_ss%Te = Te_sing*C_E
    cparams_ss%Ti = Ti_sing*C_E
    cparams_ss%ne = ne_sing
    cparams_ss%Zeff = Zeff_sing
    cparams_ss%dTau = dTau_sing
    cparams_ss%p_therm = p_therm
    cparams_ss%ConserveLA = ConserveLA
    cparams_ss%sample_test = sample_test
    cparams_ss%always_aval = always_aval
    cparams_ss%Clog_model = Clog_model
    cparams_ss%min_secRE = min_secRE
    cparams_ss%pmin_scale = pmin_scale
    cparams_ss%energy_diffusion = energy_diffusion
    cparams_ss%pitch_diffusion = pitch_diffusion
    cparams_ss%LAC_gam_resolution = LAC_gam_resolution
    cparams_ss%FP_bremsstrahlung = FP_bremsstrahlung
    cparams_ss%ngrid1 = ngrid1
    cparams_ss%Clog_const = Clog_const

    cparams_ss%gam_therm = sqrt(1+p_therm*p_therm)
    cparams_ss%gam_min = cparams_ss%gam_therm
    cparams_ss%p_min = cparams_ss%p_therm

    cparams_ss%rD = SQRT(C_E0*cparams_ss%Te/(cparams_ss%ne*C_E**2*(1.0_rp + &
         cparams_ss%Te/cparams_ss%Ti)))

    cparams_ss%re = C_E**2/(4.0_rp*C_PI*C_E0*C_ME*C_C**2)
    cparams_ss%CoulombLogee = CLogee_wu(params,cparams_ss%ne,cparams_ss%Te)
    cparams_ss%CoulombLogei = CLogei_wu(params,cparams_ss%ne,cparams_ss%Te)

    cparams_ss%VTe = VTe_wu(cparams_ss%Te)
    cparams_ss%delta = cparams_ss%VTe/C_C
    cparams_ss%Gammaco = C_E**4/(4.0_rp*C_PI*C_E0**2)
    cparams_ss%Gammac = Gammac_wu(params,cparams_ss%ne,cparams_ss%Te)

    cparams_ss%Tauc = C_ME**2*cparams_ss%VTe**3/cparams_ss%Gammac
    cparams_ss%Tau = C_ME**2*C_C**3/cparams_ss%Gammac

    cparams_ss%Ec = C_ME*C_C/(C_E*cparams_ss%Tau)
    cparams_ss%ED = cparams_ss%ne*C_E**3*cparams_ss%CoulombLogee/ &
         (4.0_rp*C_PI*C_E0**2*cparams_ss%Te)

    cparams_ss%taur=6*C_PI*C_E0*(C_ME*C_C)**3/(C_E**4*params%cpp%Bo**2)

    !	ALLOCATE(cparams_ss%rnd_num(3,cparams_ss%rnd_dim))
    !	call RANDOM_NUMBER(cparams_ss%rnd_num)
    cparams_ss%rnd_num_count = 1_idef

    !open(unit=output_unit_write,file=TRIM(params%path_to_inputs), &
    !     status='OLD',form='formatted')
    !read(output_unit_write,nml=plasmaProfiles)
    !close(output_unit_write)

    cparams_ss%P%a = radius_profile
    cparams_ss%P%ne_profile = TRIM(ne_profile)
    cparams_ss%P%neo = neo
    cparams_ss%P%n_ne = n_ne
    cparams_ss%P%a_ne = a_ne

    cparams_ss%P%Te_profile = TRIM(Te_profile)
    cparams_ss%P%Teo = Teo*C_E
    cparams_ss%P%n_Te = n_Te
    cparams_ss%P%a_Te = a_Te

    cparams_ss%P%Zeff_profile = TRIM(Zeff_profile)
    cparams_ss%P%Zeffo = Zeffo
    cparams_ss%P%n_Zeff = n_Zeff
    cparams_ss%P%a_Zeff = a_Zeff

    cparams_ss%slowing_down = slowing_down

  end subroutine load_params_ss


  subroutine initialize_collision_params(params,spp,P,F,init)
    TYPE(KORC_PARAMS), INTENT(INOUT) :: params
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)       :: spp
    TYPE(PROFILES), INTENT(INOUT)  :: P
    TYPE(FIELDS), INTENT(IN)                :: F
    LOGICAL, INTENT(IN) :: init
    INTEGER 				                       	:: ii
    REAL(rp) :: p_crit,gam_crit,maxEinterp,minEinterp

    if (params%collisions.or.((TRIM(params%field_model).eq.'M3D_C1'.or. &
         TRIM(params%field_model).eq.'NIMROD').and. &
         params%radiation)) then

       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'(/,"* * * * * * * INITIALIZING COLLISIONS * * * * * * *")')
       end if

       if (init) then

          SELECT CASE (TRIM(params%collisions_model))
          CASE (MODEL1)
             call load_params_ss(params)

             SELECT CASE(TRIM(params%bound_electron_model))
             CASE ('NO_BOUND')
                call load_params_ms(params)

                cparams_ms%Ec=cparams_ss%Ec
                if (.not.(P%ne_profile(1:6).eq.'RE-EVO')) then
                   cparams_ms%Ec_min=cparams_ms%Ec
                else
                   cparams_ms%Ec_min=cparams_ms%Ec* &
                        cparams_ms%Gammac_min/cparams_ss%Gammac
                endif

             CASE('HESSLOW')
                call load_params_ms(params)

                cparams_ms%Ec=cparams_ss%Ec

                if (.not.(cparams_ms%lowKE_REs)) then
                  cparams_ms%Ec=cparams_ms%Ec* &
                        (1._rp+sum((cparams_ms%Zo-cparams_ms%Zj)* &
                        cparams_ms%nz)/cparams_ss%ne)
                end if

                if (.not.(P%ne_profile(1:6).eq.'RE-EVO')) then
                   cparams_ms%Ec_min=cparams_ms%Ec
                else
                   cparams_ms%Ec_min=cparams_ms%Ec* &
                        cparams_ms%Gammac_min/cparams_ss%Gammac
                end if

             CASE('ROSENBLUTH')
                call load_params_ms(params)

                cparams_ms%Ec=cparams_ss%Ec* &
                     (1._rp+sum((cparams_ms%Zo-cparams_ms%Zj)* &
                     cparams_ms%nz)/cparams_ss%ne)
                if (.not.(P%ne_profile(1:6).eq.'RE-EVO')) then
                   cparams_ms%Ec_min=cparams_ms%Ec
                else
                   cparams_ms%Ec_min=cparams_ms%Ec* &
                        cparams_ms%Gammac_min/cparams_ss%Gammac
                end if

             CASE DEFAULT
                write(output_unit_write,'("Default case")')
             END SELECT

             do ii=1_idef,params%num_species
                ALLOCATE( spp(ii)%vars%nimp(spp(ii)%ppp, &
                     cparams_ms%num_impurity_species) )
                spp(ii)%vars%nimp = 0.0_rp
             end do

#ifdef FIO
             if (TRIM(params%field_model) .eq. 'M3D_C1') then
                call initialize_m3d_c1_imp(params,F,P, &
                     cparams_ms%num_impurity_species,.true.)
             endif
#endif

          CASE (MODEL2)
             call load_params_ms(params)
          CASE DEFAULT
             write(output_unit_write,'("Default case")')
          END SELECT

       end if

       if (params%LargeCollisions) then

          if (params%mpi_params%rank .EQ. 0) then
             write(output_unit_write,'(/,"* * * * * * * LARGE ANGLE COLLISIONS * * * * * * *")')
          end if

          if (TRIM(params%field_model) .eq. 'ANALYTICAL') then

             !write(6,*) 'Eo',F%Eo
             !write(6,*) 'Ec',cparams_ss%Ec
             !write(6,*) 'Ec_min',cparams_ms%Ec_min

             cparams_ss%avalanche=.TRUE.
             if (TRIM(params%collisions_model).eq.'NO_BOUND') then
                if (abs(F%Eo).lt.cparams_ss%Ec) then
                   cparams_ss%avalanche=.FALSE.
                end if
             else
                if (abs(F%Eo).lt.cparams_ms%Ec_min) then
                   cparams_ss%avalanche=.FALSE.
                end if
             end if

             if (cparams_ss%avalanche) then

                if (TRIM(params%collisions_model).eq.'NO_BOUND') then
                   p_crit=1/sqrt(abs(F%Eo)/cparams_ss%Ec-1._rp)
                else
                   p_crit=1/sqrt(abs(F%Eo)/cparams_ms%Ec_min-1._rp)
                end if
             end if

          else if ((TRIM(params%field_model) .eq. 'EXTERNAL-PSI')) then
             if (F%ReInterp_2x1t) then
                maxEinterp=maxval(F%E_3D%PHI(:,F%ind_2x1t,:)* &
                     F%FLAG3D(:,F%ind_2x1t,:))

                minEinterp=minval(F%E_3D%PHI(:,F%ind_2x1t,:)* &
                     F%FLAG3D(:,F%ind_2x1t,:))
             else if (F%E_profile.eq.'MST_FSA') then
                maxEinterp=F%E_dyn
                minEinterp=0._rp
             end if

             cparams_ss%avalanche=.TRUE.
             if (TRIM(params%collisions_model).eq.'NO_BOUND') then
                if ((abs(maxEinterp).lt.cparams_ss%Ec).and. &
                     (abs(minEinterp).lt.cparams_ss%Ec)) &
                     cparams_ss%avalanche=.FALSE.
             else
                if ((abs(maxEinterp).lt.cparams_ms%Ec_min).and. &
                     (abs(minEinterp).lt.cparams_ms%Ec_min)) &
                     cparams_ss%avalanche=.FALSE.
             end if

             !write(6,*) 'maxEinterp',maxEinterp,'minEinterp',minEinterp, &
             !      'E_c',cparams_ms%Ec,'E_c,min',cparams_ms%Ec_min, &
             !      cparams_ss%avalanche

             if (cparams_ss%avalanche) then

                if (abs(maxEinterp).gt.abs(minEinterp)) then
                   if (TRIM(params%collisions_model).eq.'NO_BOUND') then
                      p_crit=1/sqrt(abs(maxEinterp)/cparams_ss%Ec-1._rp)
                   else

                      p_crit=1/sqrt(abs(maxEinterp)/cparams_ms%Ec_min-1._rp)
                   end if
                else
                   if (TRIM(params%collisions_model).eq.'NO_BOUND') then
                      p_crit=1/sqrt(abs(minEinterp)/cparams_ss%Ec-1._rp)
                   else
                      p_crit=1/sqrt(abs(minEinterp)/cparams_ms%Ec_min-1._rp)
                   end if
                end if

             end if

          else
             write(6,*) 'Need to set p_crit!'
             call korc_abort(25)
          end if

          !if (cparams_ss%always_aval) then
          !   cparams_ss%avalanche=.TRUE.
          !   p_crit = 1.53073
          !endif

          if (cparams_ss%avalanche) then

             cparams_ss%p_crit=p_crit

             gam_crit=sqrt(1+p_crit*p_crit)

             cparams_ss%gam_crit=gam_crit

             cparams_ss%gam_therm=(gam_crit+1._rp)/2._rp
             cparams_ss%p_therm=sqrt(cparams_ss%gam_therm*cparams_ss%gam_therm-1)

             if(cparams_ss%min_secRE.eq.'THERM') then
                cparams_ss%p_min=min(cparams_ss%p_therm,cparams_ss%p_min)
                cparams_ss%gam_min=sqrt(1+cparams_ss%p_min*cparams_ss%p_min)
             else
                cparams_ss%p_min=p_crit
                cparams_ss%gam_min=gam_crit
             end if

             !write(6,*) p_crit,gam_crit,cparams_ss%p_therm,cparams_ss%gam_therm,cparams_ss%p_min,cparams_ss%gam_min



             if (params%mpi_params%rank .EQ. 0) then
                write(output_unit_write,*) 'Minimum energy of secondary RE is ',&
                     cparams_ss%min_secRE
                write(output_unit_write,*) 'p_crit/(me*c) and gam_crit are: ',p_crit,gam_crit
                write(output_unit_write,*) 'p_min/(me*c) and gam_min are: ', &
                     cparams_ss%p_min,cparams_ss%gam_min
                if(.not.init) then
                   if (TRIM(params%field_model) .eq. 'ANALYTICAL') then
                         write(output_unit_write,*) 'Maximum E_PHI : ',F%Eo*params%cpp%Eo,'V/m'
                   else if ((TRIM(params%field_model) .eq. 'EXTERNAL-PSI') &
                        .AND.(F%ReInterp_2x1t)) then
                      if (abs(maxEinterp).gt.abs(minEinterp)) then
                         write(output_unit_write,*) 'Maximum E_PHI : ',maxEinterp*params%cpp%Eo,'V/m'
                      else
                         write(output_unit_write,*) 'Maximum E_PHI : ',minEinterp*params%cpp%Eo,'V/m'
                      end if
                   endif

                   if (TRIM(params%collisions_model).eq.'NO_BOUND') then
                      write(output_unit_write,*) 'E_CH is: ',cparams_ss%Ec*params%cpp%Eo,'V/m'
                   else
                      write(output_unit_write,*) 'E_CH is: ',cparams_ms%Ec_min*params%cpp%Eo,'V/m'
                   end if
                   write(output_unit_write,*) 'tau_c,rel is: ',cparams_ss%Tau*params%cpp%time,'s'
                else
                   if (TRIM(params%field_model) .eq. 'ANALYTICAL') then
                      write(output_unit_write,*) 'Maximum E_PHI : ',F%Eo,'V/m'
                   else if ((TRIM(params%field_model) .eq. 'EXTERNAL-PSI') &
                        .AND.(F%ReInterp_2x1t)) then
                      if (abs(maxEinterp).gt.abs(minEinterp)) then
                         write(output_unit_write,*) 'Maximum E_PHI : ',maxEinterp,'V/m'
                      else
                         write(output_unit_write,*) 'Maximum E_PHI : ',minEinterp,'V/m'
                      end if
                   end if

                   if (TRIM(params%collisions_model).eq.'NO_BOUND') then
                      write(output_unit_write,*) 'E_CH is: ',cparams_ss%Ec,'V/m'
                   else
                      write(output_unit_write,*) 'E_CH is: ',cparams_ms%Ec_min,'V/m'
                   end if
                   write(output_unit_write,*) 'tau_c,rel is: ',cparams_ss%Tau,'s'
                end if
                write(output_unit_write,'("* * * * * * * * * * * * * * * * * * * * * * * * * *",/)')
             end if

          else

             if (params%mpi_params%rank .EQ. 0) then
                if(.not.init) then
                   if (TRIM(params%field_model) .eq. 'ANALYTICAL') then
                      write(output_unit_write,*) 'Maximum E_PHI : ',F%Eo*params%cpp%Eo,'V/m'
                   else if ((TRIM(params%field_model) .eq. 'EXTERNAL-PSI') &
                        .AND.(F%ReInterp_2x1t)) then
                      if (abs(maxEinterp).gt.abs(minEinterp)) then
                         write(output_unit_write,*) 'Maximum E_PHI : ',maxEinterp*params%cpp%Eo,'V/m'
                      else
                         write(output_unit_write,*) 'Maximum E_PHI : ',minEinterp*params%cpp%Eo,'V/m'
                      end if
                   end if

                   if (TRIM(params%collisions_model).eq.'NO_BOUND') then
                      write(output_unit_write,*) 'E_CH is: ',cparams_ss%Ec*params%cpp%Eo,'V/m'
                   else
                      write(output_unit_write,*) 'E_CH is: ',cparams_ms%Ec_min*params%cpp%Eo,'V/m'
                   end if

                else
                   if (TRIM(params%field_model) .eq. 'ANALYTICAL') then
                      write(output_unit_write,*) 'Maximum E_PHI : ',F%Eo,'V/m'
                   else if ((TRIM(params%field_model) .eq. 'EXTERNAL-PSI') &
                        .AND.(F%ReInterp_2x1t)) then
                      if (abs(maxEinterp).gt.abs(minEinterp)) then
                         write(output_unit_write,*) 'Maximum E_PHI : ',maxEinterp,'V/m'
                      else
                         write(output_unit_write,*) 'Maximum E_PHI : ',minEinterp,'V/m'
                      end if
                   end if

                   if (TRIM(params%collisions_model).eq.'NO_BOUND') then
                         write(output_unit_write,*) 'E_CH is: ',cparams_ss%Ec,'V/m'
                   else
                      write(output_unit_write,*) 'E_CH is: ',cparams_ms%Ec_min,'V/m'
                   end if
                end if
                write(output_unit_write,*) 'No secondary REs will be calculated in this interval'
                write(output_unit_write,*) 'p_min from initial or last time interval used'
                write(output_unit_write,*) 'to calculate collision time scales'

                write(output_unit_write,'("* * * * * * * * * * * * * * * * * * * * * * * * * *",/)')
             end if

          end if

       end if

       ! for passing minimum gamma to Hollmann sampling routines
       params%gam_min=sqrt(1+cparams_ss%p_min*cparams_ss%p_min* &
            cparams_ss%pmin_scale*cparams_ss%pmin_scale)

       !write(6,*) params%gam_min

       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'("* * * * * * * * * * * * * * * * * * * * * * * * * *",/)')
       end if

    end if


  end subroutine initialize_collision_params


  subroutine normalize_params_ms(params)
    TYPE(KORC_PARAMS), INTENT(IN) :: params

    cparams_ms%Te = cparams_ms%Te/params%cpp%temperature
    cparams_ms%ne = cparams_ms%ne/params%cpp%density
    cparams_ms%nH = cparams_ms%nH/params%cpp%density
    cparams_ms%nef = cparams_ms%nef/params%cpp%density
    cparams_ms%neb = cparams_ms%neb/params%cpp%density
    if (ALLOCATED(cparams_ms%nz)) cparams_ms%nz = cparams_ms%nz/ &
         params%cpp%density
    if (ALLOCATED(cparams_ms%IZj)) cparams_ms%IZj = cparams_ms%IZj/ &
         params%cpp%energy
    cparams_ms%rD = cparams_ms%rD/params%cpp%length
    cparams_ms%re = cparams_ms%re/params%cpp%length
    cparams_ms%Ec = cparams_ms%Ec/params%cpp%Eo
    cparams_ms%Ec_min = cparams_ms%Ec_min/params%cpp%Eo
    cparams_ms%Gammac_min = cparams_ms%Gammac_min*params%cpp%time/ &
         (params%cpp%mass**2*params%cpp%velocity**3)
  end subroutine normalize_params_ms


  subroutine normalize_params_ss(params)
    !! Calculate constant quantities used in various functions within
    !! this module
    TYPE(KORC_PARAMS), INTENT(IN) :: params

    cparams_ss%Clog1 = -1.15_rp*LOG10(1.0E-6_rp*params%cpp%density)
    cparams_ss%Clog2 = 2.3_rp*LOG10(params%cpp%temperature/C_E)
    cparams_ss%Clog0_1 = -LOG(1.0E-20_rp*params%cpp%density)/2._rp
    cparams_ss%Clog0_2 = LOG(1.0E-3 *params%cpp%temperature/C_E)
    cparams_ss%Gammaco = cparams_ss%Gammaco*params%cpp%density* &
         params%cpp%time/(params%cpp%mass**2*params%cpp%velocity**3)
    cparams_ss%VTeo = SQRT(params%cpp%temperature/C_ME)/params%cpp%velocity
    cparams_ss%deltao = params%cpp%velocity/C_C

    cparams_ss%Te = cparams_ss%Te/params%cpp%temperature
    cparams_ss%Ti = cparams_ss%Ti/params%cpp%temperature
    cparams_ss%ne = cparams_ss%ne/params%cpp%density
    cparams_ss%rD = cparams_ss%rD/params%cpp%length
    cparams_ss%re = cparams_ss%re/params%cpp%length
    cparams_ss%VTe = cparams_ss%VTe/params%cpp%velocity
    cparams_ss%Gammac = cparams_ss%Gammac*params%cpp%time/ &
         (params%cpp%mass**2*params%cpp%velocity**3)
    cparams_ss%Tau = cparams_ss%Tau/params%cpp%time
    cparams_ss%Tauc = cparams_ss%Tauc/params%cpp%time
    cparams_ss%Ec = cparams_ss%Ec/params%cpp%Eo
    cparams_ss%ED = cparams_ss%ED/params%cpp%Eo

    cparams_ss%taur=cparams_ss%taur/params%cpp%time

    cparams_ss%P%a = cparams_ss%P%a/params%cpp%length
    cparams_ss%P%neo = cparams_ss%P%neo/params%cpp%density
    cparams_ss%P%n_ne = cparams_ss%P%n_ne/params%cpp%density
    cparams_ss%P%Teo = cparams_ss%P%Teo/params%cpp%temperature
  end subroutine normalize_params_ss


  subroutine normalize_collisions_params(params)
    TYPE(KORC_PARAMS), INTENT(IN) :: params

    if (params%collisions) then
       SELECT CASE (TRIM(params%collisions_model))
       CASE (MODEL1)
          call normalize_params_ss(params)

          SELECT CASE(TRIM(params%bound_electron_model))
          CASE ('NO_BOUND')
             call normalize_params_ms(params)
          CASE('HESSLOW')
             call normalize_params_ms(params)
          CASE('ROSENBLUTH')
             call normalize_params_ms(params)
          CASE DEFAULT
             write(output_unit_write,'("Default case")')
          END SELECT

       CASE (MODEL2)
          call normalize_params_ms(params)
       CASE DEFAULT
          write(output_unit_write,'("Default case")')
       END SELECT
    end if
  end subroutine normalize_collisions_params


  subroutine collision_force(spp,U,Fcoll)
    !! For multiple-species collisions
    !! J. R. Martin-Solis et al. PoP 22, 092512 (2015)
    !! if (params%collisions .AND. (TRIM(params%collisions_model) .EQ.
    !! 'MULTIPLE_SPECIES')) then call collision_force(spp(ii),U_os,Fcoll)
    !!	U_RC = U_RC + a*Fcoll/spp(ii)%q end if

    TYPE(SPECIES), INTENT(IN) 		:: spp
    REAL(rp), DIMENSION(3), INTENT(IN) 	:: U
    REAL(rp), DIMENSION(3), INTENT(OUT) :: Fcoll
    REAL(rp), DIMENSION(3) 		:: V
    REAL(rp), DIMENSION(3) 		:: Fcolle
    REAL(rp), DIMENSION(3) 		:: Fcolli
    REAL(rp) 				:: gamma
    REAL(rp) 				:: tmp
    REAL(rp) 				:: ae
    REAL(rp) 				:: ai
    REAL(rp) 				:: Clog_ef
    REAL(rp) 				:: Clog_eb
    REAL(rp) 				:: Clog_eH
    REAL(rp) 				:: Clog_eZj
    REAL(rp) 				:: Clog_eZo
    INTEGER 				:: ppi

    gamma = SQRT(1.0_rp + DOT_PRODUCT(U,U))
    V = U/gamma

    tmp = (gamma - 1.0_rp)*SQRT(gamma + 1.0_rp)
    Clog_ef = log(0.5_rp*tmp*(cparams_ms%rD/cparams_ms%re)/gamma)
    ae = cparams_ms%nef*Clog_ef
    do ppi=1_idef,cparams_ms%num_impurity_species
       Clog_eb = log(tmp*cparams_ms%Ee_IZj(ppi))
       ae = ae + cparams_ms%neb(ppi)*Clog_eb
    end do

    tmp = (gamma**2 - 1.0_rp)/gamma
    Clog_eH = log( tmp*(cparams_ms%rD/cparams_ms%re) )
    ai = cparams_ms%nH*Clog_eH
    do ppi=1_idef,cparams_ms%num_impurity_species
       Clog_eZj = log( cparams_ms%rD/(cparams_ms%Zj(ppi)* &
            cparams_ms%re*cparams_ms%Ee_IZj(ppi)) )
       Clog_eZo = log(tmp*cparams_ms%Ee_IZj(ppi))
       ai = ai + &
            cparams_ms%nz(ppi)*(Clog_eZj*cparams_ms%Zj(ppi)**2 + &
            Clog_eZo*cparams_ms%Zo(ppi)**2)
    end do

    tmp = gamma*(gamma + 1.0_rp)/(SQRT(DOT_PRODUCT(U,U))**3)
    Fcolle = -4.0_rp*C_PI*ae*spp%m*(cparams_ms%re**2)*tmp*U

    tmp = gamma/(SQRT(DOT_PRODUCT(U,U))**3)
    Fcolli = -4.0_rp*C_PI*ai*spp%m*(cparams_ms%re**2)*tmp*U

    Fcoll = Fcolle + Fcolli
  end subroutine collision_force


  subroutine define_collisions_time_step(params,F,init)
    TYPE(KORC_PARAMS), INTENT(INOUT) 	:: params
    TYPE(FIELDS), INTENT(IN) :: F
    LOGICAL, INTENT(IN)  :: init
    INTEGER(ip) 			:: iterations
    REAL(rp) 				:: E,E_min
    REAL(rp) 				:: v
    REAL(rp) 				:: Tau
    REAL(rp), DIMENSION(3) 		:: nu
    REAL(rp) 				:: num_collisions_in_simulation


    if (params%collisions) then
       E = C_ME*C_C**2 + params%minimum_particle_energy*params%cpp%energy


       E_min=sqrt((cparams_ss%p_min*cparams_ss%pmin_scale* &
            params%cpp%mass*params%cpp%velocity* &
            C_C)**2+(C_ME*C_C**2)**2)

       !write(6,'("E_min (MeV)",E17.10)') E/(10**6*C_E)
       !write(6,'("E_min (MeV)",E17.10)') E_min/(10**6*C_E)

       !if (.not.params%LargeCollisions) then
       !   v = SQRT(1.0_rp - (C_ME*C_C**2/E)**2)
          !write(6,*) 'v_min',v
       !else
       v = SQRT(1.0_rp - (C_ME*C_C**2/E_min)**2)
          !write(6,*) 'v_therm',v
       !end if



       if ((params%profile_model.eq.'M3D_C1').or. &
            (params%profile_model(10:10).eq.'H')) then
          nu = (/nu_S_FIO(params,v),nu_D_FIO(params,v),nu_par(v)/)
       else
          nu = (/nu_S(params,v),nu_D(params,v),nu_par(v)/)
       endif

       if (.not.cparams_ss%slowing_down) nu(1)=tiny(0._rp)
       if (.not.cparams_ss%pitch_diffusion) nu(2)=tiny(0._rp)
       if (.not.cparams_ss%energy_diffusion) nu(3)=tiny(0._rp)

       Tau = MINVAL( 1.0_rp/nu )


       !write(output_unit_write,'("collision freqencies ",F25.12)') nu(3)
       !write(6,*) 'collision times',1/nu*params%cpp%time
       !write(6,*) 'p_min',cparams_ss%p_min

       cparams_ss%subcycling_iterations = ceiling(cparams_ss%dTau*Tau/ &
            params%dt,ip)
       params%coll_cadence=cparams_ss%subcycling_iterations

       if (params%LargeCollisions.and.params%snapshot_frequency.gt.0._rp) then

          !write(6,*) 'params%snapshot_frequency',params%snapshot_frequency*params%cpp%time
          !write(6,*) 'cparams_ss%dTau*Tau',cparams_ss%dTau*Tau*params%cpp%time
          !write(6,*) 'FLOOR(params%snapshot_frequency/cparams_ss%dTau*Tau)', &
          !     FLOOR(params%snapshot_frequency/ &
          !     (cparams_ss%dTau*Tau),ip)

          params%coll_per_dump=ceiling(params%snapshot_frequency/ &
               (cparams_ss%dTau*Tau))

          cparams_ss%coll_per_dump_dt=params%snapshot_frequency/params%coll_per_dump

          params%coll_per_dump_dt=cparams_ss%coll_per_dump_dt

          if (params%coll_per_dump.gt.params%t_skip) then
             write(6,*) 'more collisional iterations than orbit iterations, decrease orbit timestep!'
             call korc_abort(26)
          endif

          params%orbits_per_coll=ceiling(cparams_ss%coll_per_dump_dt/ &
               params%dt)

          params%dt=cparams_ss%coll_per_dump_dt/float(params%orbits_per_coll)

       end if

       if (init) num_collisions_in_simulation = params%simulation_time/Tau

       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'("* * * * * * * * * * * SUBCYCLING FOR  &
               COLLISIONS * * * * * * * * * * *")')

         write(output_unit_write,'("Slowing down freqency (CF): ",E17.10)') &
               nu(1)/params%cpp%time
          write(output_unit_write,'("Pitch angle scattering freqency (CB): ",E17.10)') &
               nu(2)/params%cpp%time
          write(output_unit_write,'("Speed diffusion freqency (CA): ",E17.10)') &
               nu(3)/params%cpp%time

!          write(6,*) Tau

          write(output_unit_write,'("The shorter collisional time in the simulations  &
               is: ",E17.10," s")') Tau*params%cpp%time
          write(output_unit_write,'("Number of KORC iterations per collision: ",I16)')  &
               cparams_ss%subcycling_iterations
          if (init) then
             write(output_unit_write,'("Number of collisions in simulated time: ",E17.10)')  &
                  num_collisions_in_simulation
          end if

          if (params%LargeCollisions) then

             write(output_unit_write,'("Number of collision steps per dump step: ",I16)') params%coll_per_dump

             write(output_unit_write,'("Collision time step: ",E17.10)') params%coll_per_dump_dt*params%cpp%time

             write(output_unit_write,'("Number of orbit steps per collision step: ",I16)') params%orbits_per_coll

             write(output_unit_write,'("Orbit time step: ",E17.10)') params%dt*params%cpp%time

          end if


          write(output_unit_write,'("* * * * * * * * * * * * * * * * * * * * &
               * * * * * * * * * * * * * * *",/)')
       end if
    end if
  end subroutine define_collisions_time_step


  ! * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * !
  ! * FUNCTIONS OF COLLISION OPERATOR FOR SINGLE-SPECIES PLASMAS * !
  ! * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * !

  ! *_wu functions have physical units!

  function VTe_wu(Te)
    REAL(rp), INTENT(IN) 	:: Te
    !! In Joules
    REAL(rp) 			:: VTe_wu

    VTe_wu = SQRT(2.0_rp*Te/C_ME)
  end function VTe_wu


  function VTe(Te)
    !! Dimensionless temperature
    REAL(rp), INTENT(IN) 	:: Te
    REAL(rp) 				:: VTe

    VTe = SQRT(2.0_rp*Te)*cparams_ss%VTeo
  end function VTe


  function Gammac_wu(params,ne,Te)
    !! With units
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    REAL(rp), INTENT(IN) 	:: ne
    REAL(rp), INTENT(IN) 	:: Te ! In Joules
    REAL(rp) 				:: Gammac_wu

    Gammac_wu = ne*CLogee_wu(params,ne,Te)*cparams_ss%Gammaco
  end function Gammac_wu


  function Gammacee(v,ne,Te)
    !! Dimensionless ne and Te
    REAL(rp), INTENT(IN) 	:: v
    REAL(rp), INTENT(IN) 	:: ne
    REAL(rp), INTENT(IN) 	:: Te
    REAL(rp) 				:: Gammacee

    Gammacee = ne*CLogee(v,ne,Te)*cparams_ss%Gammaco
  end function Gammacee

  function CLog_wu(ne,Te)
    !! With units
    REAL(rp), INTENT(IN) 	:: ne
    !! ne is in m^-3 and below is converted to cm^-3
    REAL(rp), INTENT(IN) 	:: Te ! In Joules
    REAL(rp) 				:: CLog_wu

    CLog_wu = 25.3_rp - 1.15_rp*LOG10(1E-6_rp*ne) + 2.3_rp*LOG10(Te/C_E)

  end function CLog_wu

  function CLog0_wu(ne,Te)
    !! With units
    REAL(rp), INTENT(IN) 	:: ne
    !! ne is in m^-3 and below is converted to cm^-3
    REAL(rp), INTENT(IN) 	:: Te ! In Joules
    REAL(rp) 				:: CLog0_wu

    CLog0_wu = 14.9_rp - LOG(1E-20_rp*ne)/2._rp + LOG(1E-3_rp*Te/C_E)

  end function CLog0_wu

  function CLogee_wu(params,ne,Te)

    !! With units
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    REAL(rp), INTENT(IN) 	:: ne
    !! ne is in m^-3 and below is converted to cm^-3
    REAL(rp), INTENT(IN) 	:: Te ! In Joules
    REAL(rp) 				:: CLogee_wu
    REAL(rp)  :: k=5._rp

    if (cparams_ss%Clog_model.eq.'HESSLOW') then
       CLogee_wu = CLog0_wu(ne,Te)+ &
            log(1+(2*(params%minimum_particle_g-1)/ &
            (VTe_wu(Te)/C_C)**2)**(k/2._rp))/k

    else if (cparams_ss%Clog_model.eq.'CONSTANT') then
       CLogee_wu = cparams_ss%Clog_const

    else if (cparams_ss%Clog_model.eq.'MCDEVITT') then
       CLogee_wu = CLog0_wu(ne,Te)+ &
            log(1+(2*(params%minimum_particle_g-1)/ &
            (VTe_wu(Te)/C_C)**2)**(k/2._rp))/k
    end if

  end function CLogee_wu

  function CLogei_wu(params,ne,Te)

    !! With units
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    REAL(rp), INTENT(IN) 	:: ne
    !! ne is in m^-3 and below is converted to cm^-3
    REAL(rp), INTENT(IN) 	:: Te ! In Joules
    REAL(rp) 				:: CLogei_wu
    REAL(rp)  :: k=5._rp
    REAL(rp)  :: p

    p=sqrt(params%minimum_particle_g**2-1)

    if (cparams_ss%Clog_model.eq.'CONSTANT') then
       CLogei_wu = cparams_ss%Clog_const
    else
       CLogei_wu = CLog0_wu(ne,Te)+ &
            log(1+(2*p/(VTe_wu(Te)/C_C))**k)/k
    end if
  end function CLogei_wu

  function CLog(ne,Te) ! Dimensionless ne and Te
    REAL(rp), INTENT(IN) 	:: ne
    REAL(rp), INTENT(IN) 	:: Te
    REAL(rp) 				:: CLog

    CLog = 25.3_rp - 1.15_rp*LOG10(ne) + 2.3_rp*LOG10(Te) + &
         cparams_ss%CLog1 + cparams_ss%CLog2
  end function CLog

  function CLog0(ne,Te) ! Dimensionless ne and Te
    REAL(rp), INTENT(IN) 	:: ne
    REAL(rp), INTENT(IN) 	:: Te
    REAL(rp) 				:: CLog0

    CLog0 = 14.9_rp - LOG(ne)/2._rp + LOG(Te) + &
         cparams_ss%CLog0_1 + cparams_ss%CLog0_2
  end function CLog0

  function CLogee(v,ne,Te)

    REAL(rp), INTENT(IN) 	:: v
    REAL(rp), INTENT(IN) 	:: ne
    !! ne is in m^-3 and below is converted to cm^-3
    REAL(rp), INTENT(IN) 	:: Te ! In Joules
    REAL(rp) 				:: CLogee
    REAL(rp)  :: k=5._rp
    REAL(rp)  :: gam
    REAL(rp) :: gam_min

    gam=1/sqrt(1-v**2)
    gam_min=cparams_ss%gam_min

    if (cparams_ss%Clog_model.eq.'HESSLOW') then
       CLogee = CLog0(ne,Te)+ &
            log(1+(2*(gam-1)/VTe(Te)**2)**(k/2._rp))/k

    else if (cparams_ss%Clog_model.eq.'CONSTANT') then
       CLogee = cparams_ss%Clog_const

    else if (cparams_ss%Clog_model.eq.'MCDEVITT') then
       CLogee = CLog0(ne,Te)+ &
            log(1+(2*(gam-1)/VTe(Te)**2)**(k/2._rp))/k+ &
            log(sqrt(2*(gam_min-1._rp)/(gam-1._rp)))
    end if

!    write(output_unit_write,*) gam,CLogee
  end function CLogee

  function CLogei(v,ne,Te)

    REAL(rp), INTENT(IN) 	:: v
    REAL(rp), INTENT(IN) 	:: ne
    !! ne is in m^-3 and below is converted to cm^-3
    REAL(rp), INTENT(IN) 	:: Te ! In Joules
    REAL(rp) 				:: CLogei
    REAL(rp)  :: k=5._rp
    REAL(rp)  :: gam,p

    gam=1/sqrt(1-v**2)
    p=gam*v

    if (cparams_ss%Clog_model.eq.'CONSTANT') then
       CLogei = cparams_ss%Clog_const
    else
       CLogei = CLog0(ne,Te)+log(1+(2*p/VTe(Te))**k)/k
    end if

  end function CLogei

  function delta(Te)
    REAL(rp), INTENT(IN) 	:: Te
    REAL(rp) 				:: delta

    delta = VTe(Te)*cparams_ss%deltao
  end function delta


  function psi(x)
    REAL(rp), INTENT(IN) 	:: x
    REAL(rp) 				:: psi

    psi = 0.5_rp*(ERF(x) - 2.0_rp*x*EXP(-x**2)/SQRT(C_PI))/x**2
  end function psi


  function CA(v)
    REAL(rp), INTENT(IN) 	:: v
    REAL(rp) 				:: CA
    REAL(rp) 				:: x

    x = v/cparams_ss%VTe

    if ((.not.cparams_ms%LargeCollisions).or. &
         (.not.(cparams_ss%P%ne_profile(1:6).eq.'RE-EVO'))) then
       CA  = cparams_ss%Gammac*psi(x)/v
    else
       CA  = cparams_ms%Gammac_min*psi(x)/v
    endif
  end function CA

  function CA_SD(v,ne,Te)
    REAL(rp), INTENT(IN) 	:: v
    REAL(rp), INTENT(IN) 	:: ne
    REAL(rp), INTENT(IN) 	:: Te
    REAL(rp) 				:: CA_SD
    REAL(rp) 				:: x

 !   write(6,*) ne,Te

    x = v/VTe(Te)
    CA_SD  = Gammacee(v,ne,Te)*psi(x)/v

!    write(output_unit_write,'("ne, "E17.10)') ne
!    write(output_unit_write,'("Te, "E17.10)') Te

!    write(output_unit_write,'("x, "E17.10)') x
!    write(output_unit_write,'("psi, "E17.10)') psi(x)
!    write(output_unit_write,'("Gammac, "E17.10)') Gammac(ne,Te)

  end function CA_SD

  function dCA_SD(v,me,ne,Te)
    REAL(rp), INTENT(IN) 	:: v
    REAL(rp), INTENT(IN) 	:: me
    REAL(rp), INTENT(IN) 	:: ne
    REAL(rp), INTENT(IN) 	:: Te
    REAL(rp) 				:: dCA_SD
    REAL(rp) 				:: x
    real(rp)  :: gam

    gam=1/sqrt(1-v**2)
    x = v/VTe(Te)
    dCA_SD  = Gammacee(v,ne,Te)*((2*(gam*v)**2-1)*psi(x)+ &
         2.0_rp*x*EXP(-x**2)/SQRT(C_PI))/(gam**3*me*v**2)
  end function dCA_SD

  function CF(params,v)
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    REAL(rp), INTENT(IN) 	:: v
    REAL(rp) 				:: CF
    REAL(rp) 				:: CF_temp
    REAL(rp) 				:: x
    INTEGER :: i
    REAL(rp)  :: k=5._rp

    x = v/cparams_ss%VTe

    if ((.not.cparams_ms%LargeCollisions).or. &
         (.not.(cparams_ss%P%ne_profile(1:6).eq.'RE-EVO'))) then
       CF  = cparams_ss%Gammac*psi(x)/cparams_ss%Te
    else
       CF  = cparams_ms%Gammac_min*psi(x)/cparams_ss%Te
    endif

    if (params%bound_electron_model.eq.'HESSLOW') then
       CF_temp=CF
       do i=1,cparams_ms%num_impurity_species
          if ((.not.cparams_ms%LargeCollisions).or. &
               (.not.(cparams_ss%P%ne_profile(1:6).eq.'RE-EVO'))) then
             CF_temp=CF_temp+CF*cparams_ms%nz(i)/cparams_ms%ne* &
                  (cparams_ms%Zo(i)-cparams_ms%Zj(i))/ &
                  CLogee(v,cparams_ss%ne,cparams_ss%Te)* &
                  (log(1+h_j(i,v)**k)/k-v**2)
          else
             CF_temp=CF_temp+CF*cparams_ms%nz(i)/cparams_ss%P%n_ne* &
                  (cparams_ms%Zo(i)-cparams_ms%Zj(i))/ &
                  CLogee(v,cparams_ss%P%n_ne,cparams_ss%Te)* &
                  (log(1+h_j(i,v)**k)/k-v**2)
          end if
       end do
       CF=CF_temp

    else if (params%bound_electron_model.eq.'ROSENBLUTH') then
       CF_temp=CF
       do i=1,cparams_ms%num_impurity_species
          CF_temp=CF_temp+CF*cparams_ms%nz(i)/cparams_ms%ne* &
               (cparams_ms%Zo(i)-cparams_ms%Zj(i))/2._rp
       end do
       CF=CF_temp

    end if

  end function CF

  function CF_FIO(params,v)
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    REAL(rp), INTENT(IN) 	:: v
    REAL(rp) 				:: CF_FIO
    REAL(rp) 				:: CF_temp
    REAL(rp) 				:: x
    INTEGER :: i
    REAL(rp)  :: k=5._rp

    x = v/cparams_ss%VTe
    CF_FIO  = cparams_ss%Gammac*psi(x)/cparams_ss%Te

    if (params%bound_electron_model.eq.'HESSLOW') then
       CF_temp=CF_FIO
       do i=1,1
          CF_temp=CF_temp+CF_FIO*cparams_ms%nz(i)/cparams_ms%ne* &
               (cparams_ms%Zo(i)-cparams_ms%Zj(i))/ &
               CLogee(v,cparams_ss%ne,cparams_ss%Te)* &
               (log(1+h_j(i,v)**k)/k-v**2)
       end do
       CF_FIO=CF_temp

    end if

  end function CF_FIO

  function CF_SD(params,v,ne,Te,P,Y_R,Y_Z)
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    REAL(rp), INTENT(IN) 	:: v
    REAL(rp), INTENT(IN) 	:: ne
    REAL(rp), INTENT(IN) 	:: Te
    REAL(rp) 				:: CF_SD
    REAL(rp) 				:: CF_temp
    REAL(rp) 				:: x
    INTEGER :: i
    REAL(rp)  :: k=5._rp,ra
    TYPE(PROFILES), INTENT(IN)  :: P
    REAL(rp), INTENT(IN) 			:: Y_R,Y_Z

    x = v/VTe(Te)
    CF_SD  = Gammacee(v,ne,Te)*psi(x)/Te

    if (params%bound_electron_model.eq.'HESSLOW') then
       CF_temp=CF_SD
       if ((cparams_ms%Zj(1).eq.0.0).and. &
            (neut_prof.eq.'UNIFORM')) then
          CF_temp=CF_temp+CF_SD*cparams_ms%nz(1)/ne* &
               (cparams_ms%Zo(1)-cparams_ms%Zj(1))/ &
               CLogee(v,ne,Te)*(log(1+h_j(1,v)**k)/k-v**2)
       else if ((cparams_ms%Zj(1).eq.0.0).and. &
            (neut_prof.eq.'HOLLOW')) then
          CF_temp=CF_temp+CF_SD*max(cparams_ms%nz(1)-ne,0._rp)/ne* &
               (cparams_ms%Zo(1)-cparams_ms%Zj(1))/ &
               CLogee(v,ne,Te)*(log(1+h_j(1,v)**k)/k-v**2)
       else if ((cparams_ms%Zj(1).eq.0.0).and. &
            (neut_prof.eq.'EDGE')) then
          ra=sqrt((Y_R-P%R0)**2+(Y_Z-P%Z0)**2)/P%a
          CF_temp=CF_temp+CF_SD*cparams_ms%nz(1)*ra**cparams_ms%neut_edge_fac/ne* &
               (cparams_ms%Zo(1)-cparams_ms%Zj(1))/ &
               CLogee(v,ne,Te)*(log(1+h_j(1,v)**k)/k-v**2)

          !write(6,*) 'ra',ra,'nimp',cparams_ms%nz(1)*ra**cparams_ms%neut_edge_fac* &
         !     params%cpp%density

       else
          CF_temp=CF_temp+CF_SD*cparams_ms%nz(1)/cparams_ms%ne* &
               (cparams_ms%Zo(1)-cparams_ms%Zj(1))/ &
               CLogee(v,ne,Te)*(log(1+h_j(1,v)**k)/k-v**2)
       endif

       do i=2,cparams_ms%num_impurity_species
          CF_temp=CF_temp+CF_SD*cparams_ms%nz(i)/cparams_ms%ne* &
               (cparams_ms%Zo(i)-cparams_ms%Zj(i))/ &
               CLogee(v,ne,Te)*(log(1+h_j(i,v)**k)/k-v**2)
       end do
       CF_SD=CF_temp

    else if (params%bound_electron_model.eq.'ROSENBLUTH') then
       CF_temp=CF_SD
       do i=1,cparams_ms%num_impurity_species
          CF_temp=CF_temp+CF_SD*cparams_ms%nz(i)/cparams_ms%ne* &
               (cparams_ms%Zo(i)-cparams_ms%Zj(i))/2._rp
       end do
       CF_SD=CF_temp

    end if

  end function CF_SD

  function CF_SD_FIO(params,v,ne,Te,nimp)
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    REAL(rp), INTENT(IN) 	:: v
    REAL(rp), INTENT(IN) 	:: ne
    REAL(rp), INTENT(IN) 	:: Te
    REAL(rp),DIMENSION(cparams_ms%num_impurity_species), INTENT(IN) 	:: nimp
    REAL(rp) 				:: CF_SD_FIO
    REAL(rp) 				:: CF_temp
    REAL(rp) 				:: x
    INTEGER :: i
    REAL(rp)  :: k=5._rp

    x = v/VTe(Te)
    CF_SD_FIO  = Gammacee(v,ne,Te)*psi(x)/Te

    if (params%bound_electron_model.eq.'HESSLOW') then
       CF_temp=CF_SD_FIO
       do i=1,cparams_ms%num_impurity_species
          CF_temp=CF_temp+CF_SD_FIO*nimp(i)/ne* &
               (cparams_ms%Zo(i)-cparams_ms%Zj(i))/ &
               CLogee(v,ne,Te)*(log(1+h_j(i,v)**k)/k-v**2)
       end do
       CF_SD_FIO=CF_temp
    end if

  end function CF_SD_FIO

  function CB_ee(v)
    REAL(rp), INTENT(IN) 	:: v
    REAL(rp) 				:: CB_ee
    REAL(rp) 				:: x

    x = v/cparams_ss%VTe

    if ((.not.cparams_ms%LargeCollisions).or. &
         (.not.(cparams_ss%P%ne_profile(1:6).eq.'RE-EVO'))) then
       CB_ee  = (0.5_rp*cparams_ss%Gammac/v)*(ERF(x) - &
            psi(x) + 0.5_rp*cparams_ss%delta**4*x**2 )
    else
       CB_ee  = (0.5_rp*cparams_ms%Gammac_min/v)*(ERF(x) - &
            psi(x) + 0.5_rp*cparams_ss%delta**4*x**2 )
    endif


  end function CB_ee

  function CB_ei(params,v)
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    REAL(rp), INTENT(IN) 	:: v
    REAL(rp) 				:: CB_ei
    REAL(rp) 				:: CB_ei_temp
    REAL(rp) 				:: x
    INTEGER :: i

    x = v/cparams_ss%VTe
    if ((.not.cparams_ms%LargeCollisions).or. &
         (.not.(cparams_ss%P%ne_profile(1:6).eq.'RE-EVO'))) then
       CB_ei  = (0.5_rp*cparams_ss%Gammac/v)*(cparams_ss%Zeff* &
            CLogei(v,cparams_ss%ne,cparams_ss%Te)/ &
            CLogee(v,cparams_ss%ne,cparams_ss%Te))
    else
       CB_ei  = (0.5_rp*cparams_ms%Gammac_min/v)*(cparams_ss%Zeff* &
            CLogei(v,cparams_ss%P%n_ne,cparams_ss%Te)/ &
            CLogee(v,cparams_ss%P%n_ne,cparams_ss%Te))
    end if


    if (params%bound_electron_model.eq.'HESSLOW') then
       CB_ei_temp=CB_ei
       do i=1,cparams_ms%num_impurity_species
          if ((.not.cparams_ms%LargeCollisions).or. &
               (.not.(cparams_ss%P%ne_profile(1:6).eq.'RE-EVO'))) then
             CB_ei_temp=CB_ei_temp+CB_ei*cparams_ms%nz(i)/(cparams_ms%ne* &
                  cparams_ss%Zeff*CLogei(v,cparams_ss%ne,cparams_ss%Te))* &
                  g_j(i,v)
          else
             CB_ei_temp=CB_ei_temp+CB_ei*cparams_ms%nz(i)/(cparams_ms%ne* &
                  cparams_ss%Zeff*CLogei(v,cparams_ss%P%n_ne,cparams_ss%Te))* &
                  g_j(i,v)
          end if
       end do
       CB_ei=CB_ei_temp

    else if (params%bound_electron_model.eq.'ROSENBLUTH') then
       CB_ei_temp=CB_ei
       do i=1,cparams_ms%num_impurity_species
          CB_ei_temp=CB_ei_temp+CB_ei*cparams_ms%nz(i)/cparams_ms%ne* &
               (cparams_ms%Zo(i)-cparams_ms%Zj(i))/2._rp
       end do
       CB_ei=CB_ei_temp

    end if

  end function CB_ei

  function CB_ei_FIO(params,v)
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    REAL(rp), INTENT(IN) 	:: v
    REAL(rp) 				:: CB_ei_FIO
    REAL(rp) 				:: CB_ei_temp
    REAL(rp) 				:: x
    INTEGER :: i

    x = v/cparams_ss%VTe
    CB_ei_FIO  = (0.5_rp*cparams_ss%Gammac/v)*(cparams_ss%Zeff* &
         CLogei(v,cparams_ss%ne,cparams_ss%Te)/ &
         CLogee(v,cparams_ss%ne,cparams_ss%Te))


    if (params%bound_electron_model.eq.'HESSLOW') then
       CB_ei_temp=CB_ei_FIO
       do i=1,1
          CB_ei_temp=CB_ei_temp+CB_ei_FIO*cparams_ms%nz(i)/(cparams_ms%ne* &
               cparams_ss%Zeff*CLogei(v,cparams_ss%ne,cparams_ss%Te))* &
               g_j(i,v)
       end do
       CB_ei_FIO=CB_ei_temp
    end if

  end function CB_ei_FIO

  function CB_ee_SD(v,ne,Te,Zeff)
    REAL(rp), INTENT(IN) 	:: v
    REAL(rp), INTENT(IN) 	:: ne
    REAL(rp), INTENT(IN) 	:: Te
    REAL(rp), INTENT(IN) 	:: Zeff
    REAL(rp) 				:: CB_ee_SD
    REAL(rp) 				:: x

    x = v/VTe(Te)
    CB_ee_SD  = (0.5_rp*Gammacee(v,ne,Te)/v)* &
         (ERF(x) - psi(x) + &
         0.5_rp*delta(Te)**4*x**2 )
  end function CB_ee_SD

  function CB_ei_SD(params,v,ne,Te,Zeff,P,Y_R,Y_Z)
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    REAL(rp), INTENT(IN) 	:: v
    REAL(rp), INTENT(IN) 	:: ne
    REAL(rp), INTENT(IN) 	:: Te
    REAL(rp), INTENT(IN) 	:: Zeff
    REAL(rp) 				:: CB_ei_SD
    REAL(rp) 				:: CB_ei_temp
    REAL(rp) 				:: x,ra
    INTEGER :: i
    TYPE(PROFILES), INTENT(IN)  :: P
    REAL(rp), INTENT(IN) 			:: Y_R,Y_Z

    x = v/VTe(Te)
    CB_ei_SD  = (0.5_rp*Gammacee(v,ne,Te)/v)* &
         (Zeff*CLogei(v,ne,Te)/CLogee(v,ne,Te))

    if (params%bound_electron_model.eq.'HESSLOW') then
       CB_ei_temp=CB_ei_SD
       if ((cparams_ms%Zj(1).eq.0.0).and. &
            (neut_prof.eq.'UNIFORM')) then
          CB_ei_temp=CB_ei_temp+CB_ei_SD*cparams_ms%nz(1)/(ne* &
               Zeff*CLogei(v,ne,Te))*g_j(1,v)
       else if ((cparams_ms%Zj(1).eq.0.0).and. &
            (neut_prof.eq.'HOLLOW')) then
          CB_ei_temp=CB_ei_temp+CB_ei_SD*max(cparams_ms%nz(1)-ne,0._rp)/(ne* &
               Zeff*CLogei(v,ne,Te))*g_j(1,v)
      else if ((cparams_ms%Zj(1).eq.0.0).and. &
            (neut_prof.eq.'EDGE')) then
          ra=sqrt((Y_R-P%R0)**2+(Y_Z-P%Z0)**2)/P%a
          CB_ei_temp=CB_ei_temp+CB_ei_SD*cparams_ms%nz(1)*ra**cparams_ms%neut_edge_fac/(ne* &
               Zeff*CLogei(v,ne,Te))*g_j(1,v)
       else
          CB_ei_temp=CB_ei_temp+CB_ei_SD*cparams_ms%nz(1)/(cparams_ms%ne* &
               Zeff*CLogei(v,ne,Te))*g_j(1,v)
       endif

       do i=2,cparams_ms%num_impurity_species
          CB_ei_temp=CB_ei_temp+CB_ei_SD*cparams_ms%nz(i)/(cparams_ms%ne* &
               Zeff*CLogei(v,ne,Te))*g_j(i,v)
       end do
       CB_ei_SD=CB_ei_temp

    else if (params%bound_electron_model.eq.'ROSENBLUTH') then
       CB_ei_temp=CB_ei_SD
       do i=1,cparams_ms%num_impurity_species
          CB_ei_temp=CB_ei_temp+CB_ei_SD*cparams_ms%nz(i)/cparams_ms%ne* &
               (cparams_ms%Zo(i)-cparams_ms%Zj(i))/2._rp
       end do
       CB_ei_SD=CB_ei_temp

    end if

  end function CB_ei_SD

  function CB_ei_SD_FIO(params,v,ne,Te,nimp,Zeff)
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    REAL(rp), INTENT(IN) 	:: v
    REAL(rp), INTENT(IN) 	:: ne
    REAL(rp), INTENT(IN) 	:: Te
    REAL(rp), INTENT(IN) 	:: Zeff
    REAL(rp), DIMENSION(cparams_ms%num_impurity_species), INTENT(IN) 	:: nimp
    REAL(rp) 				:: CB_ei_SD_FIO
    REAL(rp) 				:: CB_ei_temp
    REAL(rp) 				:: x
    INTEGER :: i

    x = v/VTe(Te)
    CB_ei_SD_FIO  = (0.5_rp*Gammacee(v,ne,Te)/v)* &
         (Zeff*CLogei(v,ne,Te)/CLogee(v,ne,Te))

    if (params%bound_electron_model.eq.'HESSLOW') then
       CB_ei_temp=CB_ei_SD_FIO
       do i=1,cparams_ms%num_impurity_species
          CB_ei_temp=CB_ei_temp+CB_ei_SD_FIO*nimp(i)/(ne* &
               Zeff*CLogei(v,ne,Te))*g_j(i,v)
       end do
       CB_ei_SD_FIO=CB_ei_temp

    end if

  end function CB_ei_SD_FIO

  function nu_S(params,v)
    ! Slowing down collision frequency
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    REAL(rp), INTENT(IN) 	:: v
      ! Normalised particle speed
    REAL(rp) 				:: nu_S
    REAL(rp) 				:: nu_S_temp
    REAL(rp) 				:: p

    p = v/SQRT(1.0_rp - v**2)
    nu_S = CF(params,v)/p

  end function nu_S

  function nu_S_FIO(params,v)
    ! Slowing down collision frequency
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    REAL(rp), INTENT(IN) 	:: v
      ! Normalised particle speed
    REAL(rp) 				:: nu_S_FIO
    REAL(rp) 				:: nu_S_temp
    REAL(rp) 				:: p

    p = v/SQRT(1.0_rp - v**2)
    nu_S_FIO = CF_FIO(params,v)/p

  end function nu_S_FIO

  function h_j(i,v)
    INTEGER, INTENT(IN) 	:: i
    REAL(rp), INTENT(IN) 	:: v
    REAL(rp)  :: gam
    REAL(rp)  :: p
    REAL(rp)  :: h_j

    gam=1/sqrt(1-v**2)
    p=v*gam

    h_j=p*sqrt(gam-1)/cparams_ms%IZj(i)

  end function h_j

  function g_j(i,v)
    INTEGER, INTENT(IN) 	:: i
    REAL(rp), INTENT(IN) 	:: v
    REAL(rp)  :: gam
    REAL(rp)  :: p
    REAL(rp)  :: g_j

    gam=1/sqrt(1-v**2)
    p=v*gam

    g_j=2._rp/3._rp*((cparams_ms%Zo(i)**2-cparams_ms%Zj(i)**2)* &
         log((p*cparams_ms%aZj(i))**(3._rp/2._rp)+1)- &
         (cparams_ms%Zo(i)-cparams_ms%Zj(i))**2* &
         (p*cparams_ms%aZj(i))**(3._rp/2._rp)/ &
         ((p*cparams_ms%aZj(i))**(3._rp/2._rp)+1))

!    write(output_unit_write,'("g_j: ",E17.10)') g_j

  end function g_j

  function nu_D(params,v)
    ! perpendicular diffusion (pitch angle scattering) collision frequency
    REAL(rp), INTENT(IN) 	:: v
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
      ! Normalised particle speed
    REAL(rp) 				:: nu_D
    REAL(rp) 				:: p

    p = v/SQRT(1.0_rp - v**2)
    nu_D = 2.0_rp*(CB_ee(v)+CB_ei(params,v))/p**2
  end function nu_D

    function nu_D_FIO(params,v)
    ! perpendicular diffusion (pitch angle scattering) collision frequency
    REAL(rp), INTENT(IN) 	:: v
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
      ! Normalised particle speed
    REAL(rp) 				:: nu_D_FIO
    REAL(rp) 				:: p

    p = v/SQRT(1.0_rp - v**2)
    nu_D_FIO = 2.0_rp*(CB_ee(v)+CB_ei_FIO(params,v))/p**2
  end function nu_D_FIO


  function nu_par(v)
    ! parallel (speed diffusion) collision frequency
    REAL(rp), INTENT(IN) 	:: v
      ! Normalised particle speed
    REAL(rp) 				:: nu_par
    REAL(rp) 				:: p

    p = v/SQRT(1.0_rp - v**2)
    nu_par = 2.0_rp*CA(v)/p**2
  end function nu_par


  function fun(v)
    REAL(rp), INTENT(IN) 	:: v
    REAL(rp) 				:: fun
    REAL(rp) 				:: x

    x = v/cparams_ss%VTe
    fun = 2.0_rp*( 1.0_rp/x + x )*EXP(-x**2)/SQRT(C_PI) - ERF(x)/x**2 - psi(v)
  end function fun


  function cross(a,b)
    REAL(rp), DIMENSION(3), INTENT(IN) 	:: a
    REAL(rp), DIMENSION(3), INTENT(IN) 	:: b
    REAL(rp), DIMENSION(3) 				:: cross

    cross(1) = a(2)*b(3) - a(3)*b(2)
    cross(2) = a(3)*b(1) - a(1)*b(3)
    cross(3) = a(1)*b(2) - a(2)*b(1)
  end function cross

  subroutine unitVectorsC(B,b1,b2,b3)
    REAL(rp), DIMENSION(3), INTENT(IN) 	:: B
    REAL(rp), DIMENSION(3), INTENT(OUT) :: b1
    REAL(rp), DIMENSION(3), INTENT(OUT) :: b2
    REAL(rp), DIMENSION(3), INTENT(OUT) :: b3

    b1 = B/SQRT(DOT_PRODUCT(B,B))

    b2 = cross(b1,(/0.0_rp,0.0_rp,1.0_rp/))
    b2 = b2/SQRT(DOT_PRODUCT(b2,b2))

    b3 = cross(b1,b2)
    b3 = b3/SQRT(DOT_PRODUCT(b3,b3))
  end subroutine unitVectorsC

  subroutine unitVectors_p(pchunk,b_unit_X,b_unit_Y,b_unit_Z,b1_X,b1_Y,b1_Z, &
       b2_X,b2_Y,b2_Z,b3_X,b3_Y,b3_Z)
    INTEGER, INTENT(IN)  :: pchunk
    REAL(rp), DIMENSION(pchunk), INTENT(IN) 	:: b_unit_X,b_unit_Y,b_unit_Z
    REAL(rp), DIMENSION(pchunk), INTENT(OUT) :: b1_X,b1_Y,b1_Z
    REAL(rp), DIMENSION(pchunk), INTENT(OUT) :: b2_X,b2_Y,b2_Z
    REAL(rp), DIMENSION(pchunk), INTENT(OUT) :: b3_X,b3_Y,b3_Z
    REAL(rp), DIMENSION(pchunk) :: b2mag,b3mag
    integer(ip) :: cc

    !$OMP SIMD
!    !$OMP& aligned(b1_X,b1_Y,b1_Z,b_unit_X,b_unit_Y,b_unit_Z, &
!    !$OMP& b2_X,b2_Y,b2_Z,b2mag,b3_X,b3_Y,b3_Z,b3mag)
    do cc=1_idef,pchunk
       b1_X(cc) = b_unit_X(cc)
       b1_Y(cc) = b_unit_Y(cc)
       b1_Z(cc) = b_unit_Z(cc)

       b2_X(cc) = b1_Y(cc)
       b2_Y(cc) = -b1_X(cc)
       b2_Z(cc) = 0._rp

       b2mag(cc)=sqrt(b2_X(cc)*b2_X(cc)+b2_Y(cc)*b2_Y(cc)+b2_Z(cc)*b2_Z(cc))

       b2_X(cc) = b2_X(cc)/b2mag(cc)
       b2_Y(cc) = b2_Y(cc)/b2mag(cc)
       b2_Z(cc) = b2_Z(cc)/b2mag(cc)

       b3_X(cc)=b1_Y(cc)*b2_Z(cc)-b1_Z(cc)*b2_Y(cc)
       b3_Y(cc)=b1_Z(cc)*b2_X(cc)-b1_X(cc)*b2_Z(cc)
       b3_Z(cc)=b1_X(cc)*b2_Y(cc)-b1_Y(cc)*b2_X(cc)

       b3mag(cc)=sqrt(b3_X(cc)*b3_X(cc)+b3_Y(cc)*b3_Y(cc)+b3_Z(cc)*b3_Z(cc))

       b3_X(cc) = b3_X(cc)/b3mag(cc)
       b3_Y(cc) = b3_Y(cc)/b3mag(cc)
       b3_Z(cc) = b3_Z(cc)/b3mag(cc)
    end do
    !$OMP END SIMD

  end subroutine unitVectors_p

  subroutine check_collisions_params(spp)
#ifdef PARALLEL_RANDOM
    USE omp_lib
#endif
    TYPE(SPECIES), INTENT(IN) :: spp
    INTEGER aux

    aux = cparams_ss%rnd_num_count + 2_idef*INT(spp%ppp,idef)

    if (aux.GE.cparams_ss%rnd_dim) then
#ifdef PARALLEL_RANDOM
       cparams_ss%rnd_num = get_random()
#else
       call RANDOM_NUMBER(cparams_ss%rnd_num)
#endif
       cparams_ss%rnd_num_count = 1_idef
    end if
  end subroutine check_collisions_params

  ! * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * !
  ! * FUNCTIONS OF COLLISION OPERATOR FOR SINGLE-SPECIES PLASMAS * !
  ! * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * !




  subroutine include_CoulombCollisions_FO_p(tt,params,X_X,X_Y,X_Z, &
       U_X,U_Y,U_Z,B_X,B_Y,B_Z,me,P,F,flagCon,flagCol,PSIp)
    !! This subroutine performs a Stochastic collision process consistent
    !! with the Fokker-Planck model for relativitic electron colliding with
    !! a thermal (Maxwellian) plasma. The collision operator is in spherical
    !! coordinates of the form found in Papp et al., NF (2011). CA
    !! corresponds to the parallel (speed diffusion) process, CF corresponds
    !! to a slowing down (momentum loss) process, and CB corresponds to a
    !! perpendicular diffusion process. Ordering of the processes are
    !! $$ \sqrt{CB}\gg CB \gg CF \sim \sqrt{CA} \gg CA,$$
    !! and only the dominant terms are kept.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)      :: F
    TYPE(KORC_PARAMS), INTENT(IN) 		:: params
    REAL(rp), DIMENSION(params%pchunk), INTENT(IN) 	:: X_X,X_Y,X_Z,PSIp
    REAL(rp), DIMENSION(params%pchunk)  	:: Y_R,Y_PHI,Y_Z
    REAL(rp), DIMENSION(params%pchunk), INTENT(INOUT) 	:: U_X,U_Y,U_Z

    REAL(rp), DIMENSION(params%pchunk) 			:: ne,Te,Zeff
    INTEGER(is), DIMENSION(params%pchunk), INTENT(INOUT) 			:: flagCol
    INTEGER(is), DIMENSION(params%pchunk), INTENT(INOUT) 			:: flagCon
    REAL(rp), INTENT(IN)  :: me

    INTEGER(ip), INTENT(IN) 			:: tt

    REAL(rp), DIMENSION(params%pchunk), INTENT(IN) 		:: B_X,B_Y,B_Z

    REAL(rp), DIMENSION(params%pchunk) 		:: b_unit_X,b_unit_Y,b_unit_Z
    REAL(rp), DIMENSION(params%pchunk) 		:: b1_X,b1_Y,b1_Z
    REAL(rp), DIMENSION(params%pchunk) 		:: b2_X,b2_Y,b2_Z
    REAL(rp), DIMENSION(params%pchunk) 		:: b3_X,b3_Y,b3_Z
    REAL(rp), DIMENSION(params%pchunk) 		:: Bmag


    REAL(rp), DIMENSION(params%pchunk,3) 			:: dW
    !! 3D Weiner process
    REAL(rp), DIMENSION(params%pchunk,3) 			:: rnd1

    REAL(rp) 					:: dt,time
    REAL(rp), DIMENSION(params%pchunk) 					:: um
    REAL(rp), DIMENSION(params%pchunk) 					:: dpm
    REAL(rp), DIMENSION(params%pchunk) 					:: vm
    REAL(rp), DIMENSION(params%pchunk) 					:: pm

    REAL(rp),DIMENSION(params%pchunk) 			:: Ub_X,Ub_Y,Ub_Z
    REAL(rp), DIMENSION(params%pchunk) 			:: xi
    REAL(rp), DIMENSION(params%pchunk) 			:: dxi
    REAL(rp), DIMENSION(params%pchunk)  			:: phi
    REAL(rp), DIMENSION(params%pchunk)  			:: dphi
    !! speed of particle
    REAL(rp),DIMENSION(params%pchunk) 					:: CAL
    REAL(rp),DIMENSION(params%pchunk) 					:: dCAL
    REAL(rp),DIMENSION(params%pchunk) 					:: CFL
    REAL(rp),DIMENSION(params%pchunk) 					:: CBL

    integer :: cc,pchunk

    pchunk=params%pchunk

    if (MODULO(params%it+tt,cparams_ss%subcycling_iterations) .EQ. 0_ip) then
       dt = REAL(cparams_ss%subcycling_iterations,rp)*params%dt
       time=params%init_time+(params%it-1+tt)*params%dt
       ! subcylcling iterations a fraction of fastest collision frequency,
       ! where fraction set by dTau in namelist &CollisionParamsSingleSpecies

       call cart_to_cyl_p(pchunk,X_X,X_Y,X_Z,Y_R,Y_PHI,Y_Z)

       if (params%profile_model(1:10).eq.'ANALYTICAL') then
          call analytical_profiles_p(pchunk,time,params,Y_R,Y_Z,P,F,ne,Te,Zeff,PSIp)
       else  if (params%profile_model(1:8).eq.'EXTERNAL') then
#ifdef PSPLINE
          call interp_FOcollision_p(pchunk,Y_R,Y_PHI,Y_Z,ne,Te,Zeff,flagCon)
#endif
       end if

       !$OMP SIMD
!       !$OMP& aligned(um,pm,vm,U_X,U_Y,U_Z,Bmag,B_X,B_Y,B_Z, &
!       !$OMP& b_unit_X,b_unit_Y,b_unit_Z,xi)
       do cc=1_idef,pchunk

          um(cc) = SQRT(U_X(cc)*U_X(cc)+U_Y(cc)*U_Y(cc)+U_Z(cc)*U_Z(cc))
          pm(cc)=me*um(cc)
          vm(cc) = um(cc)/SQRT(1.0_rp + um(cc)*um(cc))
          ! um is gamma times v, this solves for v

          Bmag(cc)= SQRT(B_X(cc)*B_X(cc)+B_Y(cc)*B_Y(cc)+B_Z(cc)*B_Z(cc))

          b_unit_X(cc)=B_X(cc)/Bmag(cc)
          b_unit_Y(cc)=B_Y(cc)/Bmag(cc)
          b_unit_Z(cc)=B_Z(cc)/Bmag(cc)

          xi(cc)=(U_X(cc)*b_unit_X(cc)+U_Y(cc)*b_unit_Y(cc)+ &
               U_Z(cc)*b_unit_Z(cc))/um(cc)

          ! pitch angle in b_unit reference frame
       end do
       !$OMP END SIMD

!       write(output_unit_write,'("vm: ",E17.10)') vm
!       write(output_unit_write,'("xi: ",E17.10)') xi

       call unitVectors_p(pchunk,b_unit_X,b_unit_Y,b_unit_Z,b1_X,b1_Y,b1_Z, &
            b2_X,b2_Y,b2_Z,b3_X,b3_Y,b3_Z)
          ! b1=b_unit, (b1,b2,b3) is right-handed

       !$OMP SIMD
!       !$OMP& aligned(phi,U_X,U_Y,U_Z,b3_X,b3_Y,b3_Z,b2_X,b2_Y,b2_Z)
       do cc=1_idef,pchunk
          phi(cc) = atan2((U_X(cc)*b3_X(cc)+U_Y(cc)*b3_Y(cc)+ &
               U_Z(cc)*b3_Z(cc)), &
               (U_X(cc)*b2_X(cc)+U_Y(cc)*b2_Y(cc)+U_Z(cc)*b2_Z(cc)))
          ! azimuthal angle in b_unit refernce frame
       end do
       !$OMP END SIMD

!       write(output_unit_write,'("phi: ",E17.10)') phi

       !$OMP SIMD
!       !$OMP& aligned(rnd1,dW,CAL,dCAL,CFL,CBL,vm,ne,Te,Zeff,dpm, &
!       !$OMP& flagCon,flagCol,dxi,xi,pm,dphi,um,Ub_X,Ub_Y,Ub_Z,U_X,U_Y,U_Z, &
!       !$OMP& b1_X,b1_Y,b1_Z,b2_X,b2_Y,b2_Z,b3_X,b3_Y,b3_Z)
       do cc=1_idef,pchunk

#ifdef PARALLEL_RANDOM
          ! uses C library to generate normal_distribution random variables,
          ! preserving parallelization where Fortran random number generator
          ! does not
          rnd1(cc,1) = get_random()
          rnd1(cc,2) = get_random()
          rnd1(cc,3) = get_random()
#else
          call RANDOM_NUMBER(rnd1)
#endif

          dW(cc,1) = SQRT(3*dt)*(-1+2*rnd1(cc,1))
          dW(cc,2) = SQRT(3*dt)*(-1+2*rnd1(cc,2))
          dW(cc,3) = SQRT(3*dt)*(-1+2*rnd1(cc,3))
          ! 3D Weiner process

          CAL(cc) = CA_SD(vm(cc),ne(cc),Te(cc))
          dCAL(cc)= dCA_SD(vm(cc),me,ne(cc),Te(cc))
          CFL(cc) = CF_SD(params,vm(cc),ne(cc),Te(cc),P,Y_R(cc),Y_Z(cc))
          CBL(cc) = (CB_ee_SD(vm(cc),ne(cc),Te(cc),Zeff(cc))+ &
               CB_ei_SD(params,vm(cc),ne(cc),Te(cc),Zeff(cc),P,Y_R(cc),Y_Z(cc)))

         if (.not.cparams_ss%slowing_down) CFL(cc)=0._rp
         if (.not.cparams_ss%pitch_diffusion) CBL(cc)=0._rp
         if (.not.cparams_ss%energy_diffusion) THEN
            CAL(cc)=0._rp
            dCAL(cc)=0._rp
         ENDIF

          dpm(cc)=REAL(flagCol(cc))*REAL(flagCon(cc))* &
               ((-CFL(cc)+dCAL(cc))*dt+ &
               sqrt(2.0_rp*CAL(cc))*dW(cc,1))
          dxi(cc)=REAL(flagCol(cc))*REAL(flagCon(cc))* &
               (-2*xi(cc)*CBL(cc)/(pm(cc)*pm(cc))*dt- &
               sqrt(2.0_rp*CBL(cc)*(1-xi(cc)*xi(cc)))/pm(cc)*dW(cc,2))
          dphi(cc)=REAL(flagCol(cc))*REAL(flagCon(cc))* &
               (sqrt(2*CBL(cc))/(pm(cc)* &
               sqrt(1-xi(cc)*xi(cc)))*dW(cc,3))

          pm(cc)=pm(cc)+dpm(cc)
          xi(cc)=xi(cc)+dxi(cc)
          phi(cc)=phi(cc)+dphi(cc)

!          if (pm(cc)<0) pm(cc)=-pm(cc)

          ! Keep xi between [-1,1]
          if (xi(cc)>1) then
             xi(cc)=1-mod(xi(cc),1._rp)
          else if (xi(cc)<-1) then
             xi(cc)=-1-mod(xi(cc),-1._rp)
          endif

          ! Keep phi between [0,pi]
!          if (phi(cc)>C_PI) then
!             phi(cc)=C_PI-mod(phi(cc),C_PI)
!          else if (phi(cc)<0) then
!             phi(cc)=mod(-phi(cc),C_PI)
!          endif

          um(cc)=pm(cc)/me

          Ub_X(cc)=um(cc)*xi(cc)
          Ub_Y(cc)=um(cc)*sqrt(1-xi(cc)*xi(cc))*cos(phi(cc))
          Ub_Z(cc)=um(cc)*sqrt(1-xi(cc)*xi(cc))*sin(phi(cc))

          U_X(cc) = Ub_X(cc)*b1_X(cc)+Ub_Y(cc)*b2_X(cc)+Ub_Z(cc)*b3_X(cc)
          U_Y(cc) = Ub_X(cc)*b1_Y(cc)+Ub_Y(cc)*b2_Y(cc)+Ub_Z(cc)*b3_Y(cc)
          U_Z(cc) = Ub_X(cc)*b1_Z(cc)+Ub_Y(cc)*b2_Z(cc)+Ub_Z(cc)*b3_Z(cc)

       end do
       !$OMP END SIMD

!       if (tt .EQ. 1_ip) then
!          write(output_unit_write,'("CA: ",E17.10)') CAL(1)
!          write(output_unit_write,'("dCA: ",E17.10)') dCAL(1)
!          write(output_unit_write,'("CF ",E17.10)') CFL(1)
!          write(output_unit_write,'("CB: ",E17.10)') CBL(1)
!       end if


       do cc=1_idef,pchunk
          if (pm(cc).lt.0) then
             write(output_unit_write,'("Momentum less than zero")')
             stop
          end if
       end do

    end if
  end subroutine include_CoulombCollisions_FO_p

#ifdef FIO
  subroutine include_CoulombCollisions_FOfio_p(tt,params,X_X,X_Y,X_Z, &
       U_X,U_Y,U_Z,B_X,B_Y,B_Z,me,P,F,flagCon,flagCol,PSIp,hint)
    !! This subroutine performs a Stochastic collision process consistent
    !! with the Fokker-Planck model for relativitic electron colliding with
    !! a thermal (Maxwellian) plasma. The collision operator is in spherical
    !! coordinates of the form found in Papp et al., NF (2011). CA
    !! corresponds to the parallel (speed diffusion) process, CF corresponds
    !! to a slowing down (momentum loss) process, and CB corresponds to a
    !! perpendicular diffusion process. Ordering of the processes are
    !! $$ \sqrt{CB}\gg CB \gg CF \sim \sqrt{CA} \gg CA,$$
    !! and only the dominant terms are kept.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)      :: F
    TYPE(KORC_PARAMS), INTENT(IN) 		:: params
    REAL(rp), DIMENSION(params%pchunk), INTENT(IN) 	:: X_X,X_Y,X_Z,PSIp
    REAL(rp), DIMENSION(params%pchunk)  	:: Y_R,Y_PHI,Y_Z
    REAL(rp), DIMENSION(params%pchunk), INTENT(INOUT) 	:: U_X,U_Y,U_Z

    REAL(rp), DIMENSION(params%pchunk) 			:: ne,Te,Zeff,ni
    REAL(rp), DIMENSION(params%pchunk,params%num_impurity_species) 	:: nimp
    INTEGER(is), DIMENSION(params%pchunk), INTENT(INOUT) 			:: flagCol
    INTEGER(is), DIMENSION(params%pchunk), INTENT(INOUT) 			:: flagCon
    REAL(rp), INTENT(IN)  :: me

    INTEGER(ip), INTENT(IN) 			:: tt

    REAL(rp), DIMENSION(params%pchunk), INTENT(IN) 		:: B_X,B_Y,B_Z
    TYPE(C_PTR),DIMENSION(params%pchunk) :: hint

    REAL(rp), DIMENSION(params%pchunk) 		:: b_unit_X,b_unit_Y,b_unit_Z
    REAL(rp), DIMENSION(params%pchunk) 		:: b1_X,b1_Y,b1_Z
    REAL(rp), DIMENSION(params%pchunk) 		:: b2_X,b2_Y,b2_Z
    REAL(rp), DIMENSION(params%pchunk) 		:: b3_X,b3_Y,b3_Z
    REAL(rp), DIMENSION(params%pchunk) 		:: Bmag


    REAL(rp), DIMENSION(params%pchunk,3) 			:: dW
    !! 3D Weiner process
    REAL(rp), DIMENSION(params%pchunk,3) 			:: rnd1

    REAL(rp) 					:: dt,time
    REAL(rp), DIMENSION(params%pchunk) 					:: um
    REAL(rp), DIMENSION(params%pchunk) 					:: dpm
    REAL(rp), DIMENSION(params%pchunk) 					:: vm
    REAL(rp), DIMENSION(params%pchunk) 					:: pm

    REAL(rp),DIMENSION(params%pchunk) 			:: Ub_X,Ub_Y,Ub_Z
    REAL(rp), DIMENSION(params%pchunk) 			:: xi
    REAL(rp), DIMENSION(params%pchunk) 			:: dxi
    REAL(rp), DIMENSION(params%pchunk)  			:: phi
    REAL(rp), DIMENSION(params%pchunk)  			:: dphi
    !! speed of particle
    REAL(rp),DIMENSION(params%pchunk) 					:: CAL
    REAL(rp),DIMENSION(params%pchunk) 					:: dCAL
    REAL(rp),DIMENSION(params%pchunk) 					:: CFL
    REAL(rp),DIMENSION(params%pchunk) 					:: CBL

    integer :: cc,pchunk

    pchunk=params%pchunk

    if (MODULO(params%it+tt,cparams_ss%subcycling_iterations) .EQ. 0_ip) then
       dt = REAL(cparams_ss%subcycling_iterations,rp)*params%dt
       time=params%init_time+(params%it-1+tt)*params%dt
       ! subcylcling iterations a fraction of fastest collision frequency,
       ! where fraction set by dTau in namelist &CollisionParamsSingleSpecies

       call cart_to_cyl_p(pchunk,X_X,X_Y,X_Z,Y_R,Y_PHI,Y_Z)

       call get_fio_profile_p(params,P,Y_R,Y_PHI,Y_Z, &
            ne,Te,flagCon,hint)

       call get_fio_ion_p(params,P,Y_R,Y_PHI,Y_Z, &
            ne,ni,nimp,Zeff,flagCon,hint)

       !write(6,*) ne,Te,nimp,Zeff

       !$OMP SIMD
!       !$OMP& aligned(um,pm,vm,U_X,U_Y,U_Z,Bmag,B_X,B_Y,B_Z, &
!       !$OMP& b_unit_X,b_unit_Y,b_unit_Z,xi)
       do cc=1_idef,pchunk

          um(cc) = SQRT(U_X(cc)*U_X(cc)+U_Y(cc)*U_Y(cc)+U_Z(cc)*U_Z(cc))
          pm(cc)=me*um(cc)
          vm(cc) = um(cc)/SQRT(1.0_rp + um(cc)*um(cc))
          ! um is gamma times v, this solves for v

          Bmag(cc)= SQRT(B_X(cc)*B_X(cc)+B_Y(cc)*B_Y(cc)+B_Z(cc)*B_Z(cc))

          b_unit_X(cc)=B_X(cc)/Bmag(cc)
          b_unit_Y(cc)=B_Y(cc)/Bmag(cc)
          b_unit_Z(cc)=B_Z(cc)/Bmag(cc)

          xi(cc)=(U_X(cc)*b_unit_X(cc)+U_Y(cc)*b_unit_Y(cc)+ &
               U_Z(cc)*b_unit_Z(cc))/um(cc)

          ! pitch angle in b_unit reference frame
       end do
       !$OMP END SIMD

!       write(output_unit_write,'("vm: ",E17.10)') vm
!       write(output_unit_write,'("xi: ",E17.10)') xi

       call unitVectors_p(pchunk,b_unit_X,b_unit_Y,b_unit_Z,b1_X,b1_Y,b1_Z, &
            b2_X,b2_Y,b2_Z,b3_X,b3_Y,b3_Z)
          ! b1=b_unit, (b1,b2,b3) is right-handed

       !$OMP SIMD
!       !$OMP& aligned(phi,U_X,U_Y,U_Z,b3_X,b3_Y,b3_Z,b2_X,b2_Y,b2_Z)
       do cc=1_idef,pchunk
          phi(cc) = atan2((U_X(cc)*b3_X(cc)+U_Y(cc)*b3_Y(cc)+ &
               U_Z(cc)*b3_Z(cc)), &
               (U_X(cc)*b2_X(cc)+U_Y(cc)*b2_Y(cc)+U_Z(cc)*b2_Z(cc)))
          ! azimuthal angle in b_unit refernce frame
       end do
       !$OMP END SIMD

!       write(output_unit_write,'("phi: ",E17.10)') phi

       !$OMP SIMD
!       !$OMP& aligned(rnd1,dW,CAL,dCAL,CFL,CBL,vm,ne,Te,Zeff,nimp,dpm, &
!       !$OMP& flagCon,flagCol,dxi,xi,pm,dphi,um,Ub_X,Ub_Y,Ub_Z,U_X,U_Y,U_Z, &
!       !$OMP& b1_X,b1_Y,b1_Z,b2_X,b2_Y,b2_Z,b3_X,b3_Y,b3_Z)
       do cc=1_idef,pchunk

#ifdef PARALLEL_RANDOM
          ! uses C library to generate normal_distribution random variables,
          ! preserving parallelization where Fortran random number generator
          ! does not
          rnd1(cc,1) = get_random()
          rnd1(cc,2) = get_random()
          rnd1(cc,3) = get_random()
#else
          call RANDOM_NUMBER(rnd1)
#endif

          dW(cc,1) = SQRT(3*dt)*(-1+2*rnd1(cc,1))
          dW(cc,2) = SQRT(3*dt)*(-1+2*rnd1(cc,2))
          dW(cc,3) = SQRT(3*dt)*(-1+2*rnd1(cc,3))
          ! 3D Weiner process

          CAL(cc) = CA_SD(vm(cc),ne(cc),Te(cc))
          dCAL(cc)= dCA_SD(vm(cc),me,ne(cc),Te(cc))
          CFL(cc) = CF_SD(params,vm(cc),ne(cc),Te(cc),P,Y_R(cc),Y_Z(cc))
          CBL(cc) = (CB_ee_SD(vm(cc),ne(cc),Te(cc),Zeff(cc))+ &
               CB_ei_SD(params,vm(cc),ne(cc),Te(cc),Zeff(cc),P,Y_R(cc),Y_Z(cc)))


          dpm(cc)=REAL(flagCol(cc))*REAL(flagCon(cc))* &
               ((-CFL(cc)+dCAL(cc))*dt+ &
               sqrt(2.0_rp*CAL(cc))*dW(cc,1))
          dxi(cc)=REAL(flagCol(cc))*REAL(flagCon(cc))* &
               (-2*xi(cc)*CBL(cc)/(pm(cc)*pm(cc))*dt- &
               sqrt(2.0_rp*CBL(cc)*(1-xi(cc)*xi(cc)))/pm(cc)*dW(cc,2))
          dphi(cc)=REAL(flagCol(cc))*REAL(flagCon(cc))* &
               (sqrt(2*CBL(cc))/(pm(cc)* &
               sqrt(1-xi(cc)*xi(cc)))*dW(cc,3))

          pm(cc)=pm(cc)+dpm(cc)
          xi(cc)=xi(cc)+dxi(cc)
          phi(cc)=phi(cc)+dphi(cc)

!          if (pm(cc)<0) pm(cc)=-pm(cc)

          ! Keep xi between [-1,1]
          if (xi(cc)>1) then
             xi(cc)=1-mod(xi(cc),1._rp)
          else if (xi(cc)<-1) then
             xi(cc)=-1-mod(xi(cc),-1._rp)
          endif

          ! Keep phi between [0,pi]
!          if (phi(cc)>C_PI) then
!             phi(cc)=C_PI-mod(phi(cc),C_PI)
!          else if (phi(cc)<0) then
!             phi(cc)=mod(-phi(cc),C_PI)
!          endif

          um(cc)=pm(cc)/me

          Ub_X(cc)=um(cc)*xi(cc)
          Ub_Y(cc)=um(cc)*sqrt(1-xi(cc)*xi(cc))*cos(phi(cc))
          Ub_Z(cc)=um(cc)*sqrt(1-xi(cc)*xi(cc))*sin(phi(cc))

          U_X(cc) = Ub_X(cc)*b1_X(cc)+Ub_Y(cc)*b2_X(cc)+Ub_Z(cc)*b3_X(cc)
          U_Y(cc) = Ub_X(cc)*b1_Y(cc)+Ub_Y(cc)*b2_Y(cc)+Ub_Z(cc)*b3_Y(cc)
          U_Z(cc) = Ub_X(cc)*b1_Z(cc)+Ub_Y(cc)*b2_Z(cc)+Ub_Z(cc)*b3_Z(cc)

       end do
       !$OMP END SIMD

!       if (tt .EQ. 1_ip) then
!          write(output_unit_write,'("CA: ",E17.10)') CAL(1)
!          write(output_unit_write,'("dCA: ",E17.10)') dCAL(1)
!          write(output_unit_write,'("CF ",E17.10)') CFL(1)
!          write(output_unit_write,'("CB: ",E17.10)') CBL(1)
!       end if


       do cc=1_idef,pchunk
          if (pm(cc).lt.0) then
             write(output_unit_write,'("Momentum less than zero")')
             stop
          end if
       end do

    end if
  end subroutine include_CoulombCollisions_FOfio_p
#endif

  subroutine include_CoulombCollisions_GC_p(tt,params,Y_R,Y_PHI,Y_Z, &
       Ppll,Pmu,me,flagCon,flagCol,F,P,E_PHI,ne,PSIp)

    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    TYPE(KORC_PARAMS), INTENT(INOUT) 		:: params
    REAL(rp), DIMENSION(params%pchunk), INTENT(INOUT) 	:: Ppll
    REAL(rp), DIMENSION(params%pchunk), INTENT(INOUT) 	:: Pmu
    REAL(rp), DIMENSION(params%pchunk) 			:: Bmag
    REAL(rp), DIMENSION(params%pchunk) 	:: B_R,B_PHI,B_Z
    REAL(rp), DIMENSION(params%pchunk) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp), DIMENSION(params%pchunk) :: gradB_R,gradB_PHI,gradB_Z
    REAL(rp), DIMENSION(params%pchunk) 	:: E_R,E_Z
    REAL(rp), DIMENSION(params%pchunk), INTENT(OUT) 	:: E_PHI,ne,PSIp
    REAL(rp), DIMENSION(params%pchunk), INTENT(IN) 			:: Y_R,Y_PHI,Y_Z
    INTEGER(is), DIMENSION(params%pchunk), INTENT(INOUT) 	:: flagCol
    INTEGER(is), DIMENSION(params%pchunk), INTENT(INOUT) 	:: flagCon
    REAL(rp), INTENT(IN) 			:: me
    REAL(rp), DIMENSION(params%pchunk) 			:: Te,Zeff
    REAL(rp), DIMENSION(params%pchunk) 			:: nAr0,nAr1,nAr2,nAr3
    REAL(rp), DIMENSION(params%pchunk) 			:: nD,nD1
    REAL(rp), DIMENSION(params%pchunk,2) 			:: dW
    REAL(rp), DIMENSION(params%pchunk,2) 			:: rnd1
    REAL(rp) 					:: dt,time
    REAL(rp), DIMENSION(params%pchunk) 					:: pm
    REAL(rp), DIMENSION(params%pchunk)  					:: dp
    REAL(rp), DIMENSION(params%pchunk)  					:: xi
    REAL(rp), DIMENSION(params%pchunk)  					:: dxi
    REAL(rp), DIMENSION(params%pchunk)  					:: v,gam
    !! speed of particle
    REAL(rp), DIMENSION(params%pchunk) 					:: CAL
    REAL(rp) , DIMENSION(params%pchunk)					:: dCAL
    REAL(rp), DIMENSION(params%pchunk) 					:: CFL
    REAL(rp), DIMENSION(params%pchunk) 					:: CBL
    REAL(rp), DIMENSION(params%pchunk) 	:: SC_p,SC_mu,BREM_p
    REAL(rp) 					:: kappa
    integer :: cc,pchunk
    integer(ip),INTENT(IN) :: tt
    REAL(rp), DIMENSION(params%pchunk,params%num_impurity_species) 	:: nimp
    REAL(rp), DIMENSION(params%pchunk) 	:: E_PHI_tmp

    pchunk=params%pchunk

    if (MODULO(params%it+tt,cparams_ss%subcycling_iterations) .EQ. 0_ip) then
       dt = REAL(cparams_ss%subcycling_iterations,rp)*params%dt
       time=params%init_time+(params%it-1+tt)*params%dt


       if (params%field_eval.eq.'eqn') then
          call analytical_fields_GC_p(pchunk,F,Y_R,Y_PHI, &
               Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
               gradB_R,gradB_PHI,gradB_Z,PSIp)
       else if (params%field_eval.eq.'interp') then
#ifdef PSPLINE
          if (F%axisymmetric_fields) then
             if (F%Bflux) then
                if (.not.params%SC_E) then
                   call calculate_GCfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
                        E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                        gradB_R,gradB_PHI,gradB_Z,flagCon,PSIp)
                else
                   call calculate_GCfields_p_FS(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
                        E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                        gradB_R,gradB_PHI,gradB_Z,flagCon,PSIp)
                end if

             else if (F%ReInterp_2x1t) then
                call calculate_GCfieldswE_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
                     E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                     gradB_R,gradB_PHI,gradB_Z,flagCon,PSIp)
             else if (F%Bflux3D) then
                call calculate_GCfields_2x1t_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
                     E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                     gradB_R,gradB_PHI,gradB_Z,flagCon,PSIp,time)
             else if (F%dBfield) then
                call calculate_2DBdBfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
                     E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                     gradB_R,gradB_PHI,gradB_Z,flagCon,PSIp)
             else
                call interp_fields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
                     E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                     gradB_R,gradB_PHI,gradB_Z,flagCon)
             end if
          else
             if (F%dBfield) then
                call calculate_2DBdBfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
                     E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                     gradB_R,gradB_PHI,gradB_Z,flagCon,PSIp)
             else
                call interp_fields_3D_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
                     E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                     gradB_R,gradB_PHI,gradB_Z,flagCon)
             end if
          endif
          call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)
#endif
       end if


       if (params%profile_model(1:10).eq.'ANALYTICAL') then
          call analytical_profiles_p(pchunk,time,params,Y_R,Y_Z,P,F,ne,Te,Zeff,PSIp)
       else if (params%profile_model(1:8).eq.'EXTERNAL') then
#ifdef PSPLINE
          if (params%profile_model(10:10).eq.'H') then
             call interp_Hcollision_p(pchunk,Y_R,Y_PHI,Y_Z,ne,Te,Zeff, &
                  nAr0,nAr1,nAr2,nAr3,nD,nD1,flagCon)
             do cc=1_idef,pchunk
                nimp(cc,1)=nAr0(cc)
                nimp(cc,2)=nAr1(cc)
                nimp(cc,3)=nAr2(cc)
                nimp(cc,4)=nAr3(cc)
                nimp(cc,5)=nD(cc)
                nimp(cc,6)=nD1(cc)
             end do
          else
             call interp_FOcollision_p(pchunk,Y_R,Y_PHI,Y_Z,ne,Te,Zeff,flagCon)
          endif
#endif
       end if

       E_PHI_tmp=E_PHI
       if (.not.params%FokPlan) E_PHI=0._rp

       !$OMP SIMD
!       !$OMP& aligned (pm,xi,v,Ppll,Bmag,Pmu)
       do cc=1_idef,pchunk
          Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))
          ! Transform p_pll,mu to P,eta
          pm(cc) = SQRT(Ppll(cc)*Ppll(cc)+2*me*Bmag(cc)*Pmu(cc))
          xi(cc) = Ppll(cc)/pm(cc)

          gam(cc) = sqrt(1+pm(cc)*pm(cc))

          v(cc) = pm(cc)/gam(cc)
          ! normalized speed (v_K=v_P/c)
       end do
       !$OMP END SIMD

!       write(output_unit_write,'("ne: "E17.10)') ne
!       write(output_unit_write,'("Te: "E17.10)') Te
!       write(output_unit_write,'("Bmag: "E17.10)') Bmag
!       write(output_unit_write,'("v: ",E17.10)') v
!       write(output_unit_write,'("xi: ",E17.10)') xi
       !       write(output_unit_write,'("size(E_PHI_GC): ",I16)') size(E_PHI)


       !$OMP SIMD
!       !$OMP& aligned(rnd1,dW,CAL,dCAL,CFL,CBL,v,ne,Te,Zeff,dp, &
!       !$OMP& flagCon,flagCol,dxi,xi,pm,Ppll,Pmu,Bmag)
       do cc=1_idef,pchunk

#ifdef PARALLEL_RANDOM
          rnd1(cc,1) = get_random()
          rnd1(cc,2) = get_random()
          !       rnd1(:,1) = get_random_mkl()
          !       rnd1(:,2) = get_random_mkl()
#else
          call RANDOM_NUMBER(rnd1)
#endif

          dW(cc,1) = SQRT(3*dt)*(-1+2*rnd1(cc,1))
          dW(cc,2) = SQRT(3*dt)*(-1+2*rnd1(cc,2))

!          write(output_unit_write,'("dW1: ",E17.10)') dW(cc,1)
!          write(output_unit_write,'("dW2: ",E17.10)') dW(cc,2)

          if (params%profile_model(10:10).eq.'H') then
             CAL(cc) = CA_SD(v(cc),ne(cc),Te(cc))
             dCAL(cc)= dCA_SD(v(cc),me,ne(cc),Te(cc))
             CFL(cc) = CF_SD_FIO(params,v(cc),ne(cc),Te(cc),nimp(cc,:))
             CBL(cc) = (CB_ee_SD(v(cc),ne(cc),Te(cc),Zeff(cc))+ &
                  CB_ei_SD_FIO(params,v(cc),ne(cc),Te(cc),nimp(cc,:),Zeff(cc)))
          else
             CAL(cc) = CA_SD(v(cc),ne(cc),Te(cc))
             dCAL(cc)= dCA_SD(v(cc),me,ne(cc),Te(cc))
             CFL(cc) = CF_SD(params,v(cc),ne(cc),Te(cc),P,Y_R(cc),Y_Z(cc))
             CBL(cc) = (CB_ee_SD(v(cc),ne(cc),Te(cc),Zeff(cc))+ &
                  CB_ei_SD(params,v(cc),ne(cc),Te(cc),Zeff(cc),P,Y_R(cc),Y_Z(cc)))
          endif

          if (.not.cparams_ss%slowing_down) CFL(cc)=0._rp
          if (.not.cparams_ss%pitch_diffusion) CBL(cc)=0._rp
          if (.not.cparams_ss%energy_diffusion) THEN
             CAL(cc)=0._rp
             dCAL(cc)=0._rp
          ENDIF

          dp(cc)=REAL(flagCol(cc))*REAL(flagCon(cc))* &
               ((-CFL(cc)+dCAL(cc)+E_PHI(cc)*xi(cc))*dt+ &
               sqrt(2.0_rp*CAL(cc))*dW(cc,1))

          dxi(cc)=REAL(flagCol(cc))*REAL(flagCon(cc))* &
               ((-2*xi(cc)*CBL(cc)/(pm(cc)*pm(cc))+ &
               E_PHI(cc)*(1-xi(cc)*xi(cc))/pm(cc))*dt- &
               sqrt(2.0_rp*CBL(cc)*(1-xi(cc)*xi(cc)))/pm(cc)*dW(cc,2))

!          write(output_unit_write,'("dp: ",E17.10)') dp(cc)
!          write(output_unit_write,'("dxi: ",E17.10)') dxi(cc)

       end do
       !$OMP END SIMD

       if (params%FokPlan.and.params%radiation) then
          if(params%GC_rad_model.eq.'SDE') then

             !$OMP SIMD
             do cc=1_idef,pchunk

                SC_p(cc)=-gam(cc)*pm(cc)*(1-xi(cc)*xi(cc))/ &
                     (cparams_ss%taur/Bmag(cc)**2)
                SC_mu(cc)=xi(cc)*(1-xi(cc)*xi(cc))/ &
                     ((cparams_ss%taur/Bmag(cc)**2)*gam(cc))

                kappa=2._rp*C_PI*C_RE**2._rp*C_ME*C_C**2._rp
                BREM_p(cc)=-2._rp*ne(cc)*kappa*Zeff(cc)*(Zeff(cc)+1._rp)* &
                     C_a/C_PI*(gam(cc)-1._rp)*(log(2._rp*gam(cc))-1._rp/3._rp)


                dp(cc)=dp(cc)+(SC_p(cc)+BREM_p(cc))*dt* &
                     REAL(flagCol(cc))*REAL(flagCon(cc))
                dxi(cc)=dxi(cc)+(SC_mu(cc))*dt* &
                     REAL(flagCol(cc))*REAL(flagCon(cc))

             end do
             !$OMP END SIMD

          end if
       end if

       !$OMP SIMD
       do cc=1_idef,pchunk

          pm(cc)=pm(cc)+dp(cc)
          xi(cc)=xi(cc)+dxi(cc)
       end do
       !$OMP END SIMD

       do cc=1_idef,pchunk
!          if (pm(cc)<0) pm(cc)=-pm(cc)

          ! Keep xi between [-1,1]
          if (xi(cc)>1) then
!             write(output_unit_write,'("High xi at: ",E17.10," with dxi: ",E17.10)') &
!                  time*params%cpp%time, dxi(cc)
             xi(cc)=1-mod(xi(cc),1._rp)
          else if (xi(cc)<-1) then
             xi(cc)=-1-mod(xi(cc),-1._rp)
!             write(output_unit_write,'("Low xi at: ",E17.10," with dxi: ",E17.10)') &
!                  time*params%cpp%time, dxi(cc)
          endif

       end do

       !$OMP SIMD
       do cc=1_idef,pchunk
          ! Transform P,xi to p_pll,mu
          Ppll(cc)=pm(cc)*xi(cc)
          Pmu(cc)=(pm(cc)*pm(cc)-Ppll(cc)*Ppll(cc))/(2*me*Bmag(cc))
       end do
       !$OMP END SIMD


!       write(output_unit_write,'("rnd1: ",E17.10)') rnd1
!       write(output_unit_write,'("flag: ",I16)') flag
!       write(output_unit_write,'("CA: ",E17.10)') CAL
!       write(output_unit_write,'("dCA: ",E17.10)') dCAL
!       write(output_unit_write,'("CF ",E17.10)') CFL
!       write(output_unit_write,'("CB: ",E17.10)') CBL
!       write(output_unit_write,'("dp: ",E17.10)') dp
!       write(output_unit_write,'("dxi: ",E17.10)') dxi
!       write(output_unit_write,'("Ppll: ",E17.10)') Ppll
!      write(output_unit_write,'("Pmu: ",E17.10)') Pmu
!       write(output_unit_write,'("E_PHI_COL: ",E17.10)') E_PHI

       do cc=1_idef,pchunk
          if ((pm(cc).lt.min(cparams_ss%p_min*cparams_ss%pmin_scale, &
               p_therm)).and.flagCol(cc).eq.1_ip) then
!             write(output_unit_write,'("Momentum less than zero")')
             !             stop
!             write(output_unit_write,'("Particle not tracked at: ",E17.10," &
!                  & with xi: ",E17.10)') time*params%cpp%time, xi(cc)
             flagCol(cc)=0_ip
          end if
       end do

!       if (tt .EQ. 1_ip) then
!          write(output_unit_write,'("dp_rad: ",E17.10)') &
!               -gam(1)*pm(1)*(1-xi(1)*xi(1))/ &
!               (cparams_ss%taur/Bmag(1)**2)*dt
!          write(output_unit_write,'("dxi_rad: ",E17.10)') &
!               xi(1)*(1-xi(1)*xi(1))/ &
!               ((cparams_ss%taur/Bmag(1)**2)*gam(1))*dt
!       end if

!       if (tt .EQ. 1_ip) then
!          write(output_unit_write,'("CA: ",E17.10)') CAL(1)
!          write(output_unit_write,'("dCA: ",E17.10)') dCAL(1)
!          write(output_unit_write,'("CF ",E17.10)') CFL(1)
!          write(output_unit_write,'("CB: ",E17.10)') CBL(1)
!       end if

       if (.not.params%FokPlan) E_PHI=E_PHI_tmp

    end if

  end subroutine include_CoulombCollisions_GC_p

  subroutine include_CoulombCollisionsLA_GC_p(spp,achunk,tt,params, &
       Y_R,Y_PHI,Y_Z,Ppll,Pmu,me,flagCon,flagCol,F,P,E_PHI,ne,Te,PSIp)

    TYPE(SPECIES), INTENT(INOUT)    :: spp
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    TYPE(KORC_PARAMS), INTENT(INOUT) 		:: params
    INTEGER, INTENT(IN) :: achunk
    REAL(rp), DIMENSION(achunk), INTENT(INOUT) 	:: Ppll
    REAL(rp), DIMENSION(achunk), INTENT(INOUT) 	:: Pmu
    REAL(rp), DIMENSION(achunk) 			:: Bmag
    REAL(rp), DIMENSION(achunk) 	:: B_R,B_PHI,B_Z
    REAL(rp), DIMENSION(achunk) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp), DIMENSION(achunk) :: gradB_R,gradB_PHI,gradB_Z,ntot
    REAL(rp), DIMENSION(achunk) 	:: E_R,E_Z,E_PHI_LAC
    REAL(rp), DIMENSION(achunk), INTENT(OUT) 	:: E_PHI,ne,Te,PSIp
    REAL(rp), DIMENSION(achunk), INTENT(IN) 			:: Y_R,Y_PHI,Y_Z
    INTEGER(is), DIMENSION(achunk), INTENT(INOUT) 	:: flagCol
    INTEGER(is), DIMENSION(achunk), INTENT(INOUT) 	:: flagCon
    REAL(rp), INTENT(IN) 			:: me
    REAL(rp), DIMENSION(achunk) 			:: Zeff
    REAL(rp), DIMENSION(achunk) 			:: nAr0,nAr1,nAr2,nAr3
    REAL(rp), DIMENSION(achunk) 			:: nD,nD1
    REAL(rp), DIMENSION(achunk,2) 			:: dW
    REAL(rp), DIMENSION(achunk,2) 			:: rnd1
    REAL(rp) 					:: dt,time
    REAL(rp), DIMENSION(achunk) 	:: pm,pm0
    REAL(rp), DIMENSION(achunk)  	:: dp
    REAL(rp), DIMENSION(achunk)  	:: xi,xi0
    REAL(rp), DIMENSION(achunk)  	:: dxi
    REAL(rp), DIMENSION(achunk)  					:: v,gam
    !! speed of particle
    REAL(rp), DIMENSION(achunk) 					:: CAL
    REAL(rp) , DIMENSION(achunk)					:: dCAL
    REAL(rp), DIMENSION(achunk) 					:: CFL
    REAL(rp), DIMENSION(achunk) 					:: CBL
    REAL(rp), DIMENSION(achunk) 	:: SC_p,SC_xi,BREM_p
    REAL(rp) 					:: kappa,ra
    integer :: cc,ii
    integer(ip),INTENT(IN) :: tt
    REAL(rp), DIMENSION(achunk,params%num_impurity_species) 	:: nimp


    dt=cparams_ss%coll_per_dump_dt
    time=params%init_time+(params%it-1)*params%dt+ &
         tt*cparams_ss%coll_per_dump_dt


    if (params%field_eval.eq.'eqn') then
       call analytical_fields_GC_p(achunk,F,Y_R,Y_PHI, &
            Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
            gradB_R,gradB_PHI,gradB_Z,PSIp)
    else if (params%field_eval.eq.'interp') then
#ifdef PSPLINE
       if (F%axisymmetric_fields) then
          if (F%Bflux) then
             if (.not.params%SC_E) then
                call calculate_GCfields_p(achunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
                     E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                     gradB_R,gradB_PHI,gradB_Z,flagCon,PSIp)
             else
                call calculate_GCfields_p_FS(achunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
                     E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                     gradB_R,gradB_PHI,gradB_Z,flagCon,PSIp)
             end if

          else if (F%ReInterp_2x1t) then
             call calculate_GCfieldswE_p(achunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
                  E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                  gradB_R,gradB_PHI,gradB_Z,flagCon,PSIp)
          else if (F%Bflux3D) then
             call calculate_GCfields_2x1t_p(achunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
                  E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                  gradB_R,gradB_PHI,gradB_Z,flagCon,PSIp,time)
          else if (F%dBfield) then
             call calculate_2DBdBfields_p(achunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
                  E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                  gradB_R,gradB_PHI,gradB_Z,flagCon,PSIp)
          else
             call interp_fields_p(achunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
                  E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                  gradB_R,gradB_PHI,gradB_Z,flagCon)
          end if
       else
          if (F%dBfield) then
             call calculate_2DBdBfields_p(achunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
                  E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                  gradB_R,gradB_PHI,gradB_Z,flagCon,PSIp)
          else
             call interp_fields_3D_p(achunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
                  E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                  gradB_R,gradB_PHI,gradB_Z,flagCon)
          end if
       endif
#endif
       if(.not.F%ReInterp_2x1t) call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)
    end if


    if (params%profile_model(1:10).eq.'ANALYTICAL') then
       call analytical_profiles_p(achunk,time,params,Y_R,Y_Z,P,F, &
            ne,Te,Zeff,PSIp)
    else if (params%profile_model(1:8).eq.'EXTERNAL') then
#ifdef PSPLINE
       if (params%profile_model(10:10).eq.'H') then
          call interp_Hcollision_p(achunk,Y_R,Y_PHI,Y_Z,ne,Te,Zeff, &
               nAr0,nAr1,nAr2,nAr3,nD,nD1,flagCon)
          do cc=1_idef,achunk
             nimp(cc,1)=nAr0(cc)
             nimp(cc,2)=nAr1(cc)
             nimp(cc,3)=nAr2(cc)
             nimp(cc,4)=nAr3(cc)
             nimp(cc,5)=nD(cc)
             nimp(cc,6)=nD1(cc)
          end do
       else
          call interp_FOcollision_p(achunk,Y_R,Y_PHI,Y_Z,ne,Te,Zeff,flagCon)
       endif
#endif
    end if

    E_PHI_LAC=E_PHI
    if (.not.params%FokPlan) E_PHI=0._rp

    !$OMP SIMD
    !       !$OMP& aligned (pm,xi,v,Ppll,Bmag,Pmu)
    do cc=1_idef,achunk
       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))
       ! Transform p_pll,mu to P,eta
       pm(cc) = SQRT(Ppll(cc)*Ppll(cc)+2*me*Bmag(cc)*Pmu(cc))
       pm0(cc)=pm(cc)
       xi(cc) = Ppll(cc)/pm(cc)
       xi0(cc)=xi(cc)

       gam(cc) = sqrt(1+pm(cc)*pm(cc))

       v(cc) = pm(cc)/gam(cc)
       ! normalized speed (v_K=v_P/c)
    end do
    !$OMP END SIMD

    !       write(output_unit_write,'("ne: "E17.10)') ne
    !       write(output_unit_write,'("Te: "E17.10)') Te
    !       write(output_unit_write,'("Bmag: "E17.10)') Bmag
    !       write(output_unit_write,'("v: ",E17.10)') v
    !       write(output_unit_write,'("xi: ",E17.10)') xi
    !       write(output_unit_write,'("size(E_PHI_GC): ",I16)') size(E_PHI)


    !$OMP SIMD
    !       !$OMP& aligned(rnd1,dW,CAL,dCAL,CFL,CBL,v,ne,Te,Zeff,dp, &
    !       !$OMP& flagCon,flagCol,dxi,xi,pm,Ppll,Pmu,Bmag)
    do cc=1_idef,achunk

#ifdef PARALLEL_RANDOM
       rnd1(cc,1) = get_random()
       rnd1(cc,2) = get_random()
       !       rnd1(:,1) = get_random_mkl()
       !       rnd1(:,2) = get_random_mkl()
#else
       call RANDOM_NUMBER(rnd1)
#endif

       dW(cc,1) = SQRT(3*dt)*(-1+2*rnd1(cc,1))
       dW(cc,2) = SQRT(3*dt)*(-1+2*rnd1(cc,2))

       !          write(output_unit_write,'("dW1: ",E17.10)') dW(cc,1)
       !          write(output_unit_write,'("dW2: ",E17.10)') dW(cc,2)

       if (params%profile_model(10:10).eq.'H') then
          CAL(cc) = CA_SD(v(cc),ne(cc),Te(cc))
          dCAL(cc)= dCA_SD(v(cc),me,ne(cc),Te(cc))
          CFL(cc) = CF_SD_FIO(params,v(cc),ne(cc),Te(cc),nimp(cc,:))
          CBL(cc) = (CB_ee_SD(v(cc),ne(cc),Te(cc),Zeff(cc))+ &
               CB_ei_SD_FIO(params,v(cc),ne(cc),Te(cc),nimp(cc,:),Zeff(cc)))
          ! using *_FIO routine to include additions due to Hollmann species
          ! profile
       else
          CAL(cc) = CA_SD(v(cc),ne(cc),Te(cc))
          dCAL(cc)= dCA_SD(v(cc),me,ne(cc),Te(cc))
          CFL(cc) = CF_SD(params,v(cc),ne(cc),Te(cc),P,Y_R(cc),Y_Z(cc))
          CBL(cc) = (CB_ee_SD(v(cc),ne(cc),Te(cc),Zeff(cc))+ &
               CB_ei_SD(params,v(cc),ne(cc),Te(cc),Zeff(cc),P,Y_R(cc),Y_Z(cc)))
       endif

       if (.not.cparams_ss%slowing_down) CFL(cc)=0._rp
       if (.not.cparams_ss%pitch_diffusion) CBL(cc)=0._rp
       if (.not.cparams_ss%energy_diffusion) THEN
          CAL(cc)=0._rp
          dCAL(cc)=0._rp
       ENDIF

       dp(cc)=REAL(flagCol(cc))*REAL(flagCon(cc))* &
            ((-CFL(cc)+dCAL(cc)+E_PHI(cc)*xi(cc))*dt+ &
            sqrt(2.0_rp*CAL(cc))*dW(cc,1))

       dxi(cc)=REAL(flagCol(cc))*REAL(flagCon(cc))* &
            ((-2*xi(cc)*CBL(cc)/(pm(cc)*pm(cc))+ &
            E_PHI(cc)*(1-xi(cc)*xi(cc))/pm(cc))*dt- &
            sqrt(2.0_rp*CBL(cc)*(1-xi(cc)*xi(cc)))/pm(cc)*dW(cc,2))

       if(cparams_ss%sample_test) then
          dp(cc)=0._rp
          dxi(cc)=0._rp
       end if

       !          write(output_unit_write,'("dp: ",E17.10)') dp(cc)
       !          write(output_unit_write,'("dxi: ",E17.10)') dxi(cc)

       !write(6,*) 'gam,xi',gam(cc),xi(cc)

       !write(6,*) 'dpE',E_PHI(cc)*xi(cc)*dt
       !write(6,*) 'dpCF',CFL(cc)*dt
       !write(6,*) 'dpCA',sqrt(2*CAL(cc)*dt)

       !write(6,*) 'dxiE',E_PHI(cc)*(1-xi(cc)*xi(cc))/pm(cc)*dt
       !write(6,*) 'dxiCB',2*xi(cc)*CBL(cc)/(pm(cc)*pm(cc))*dt

    end do
    !$OMP END SIMD

    if (params%FokPlan.and.params%radiation) then
       if(params%GC_rad_model.eq.'SDE') then

          !$OMP SIMD
          do cc=1_idef,achunk

             SC_p(cc)=-gam(cc)*pm(cc)*(1-xi(cc)*xi(cc))/ &
                  (cparams_ss%taur/Bmag(cc)**2)
             SC_xi(cc)=xi(cc)*(1-xi(cc)*xi(cc))/ &
                  ((cparams_ss%taur/Bmag(cc)**2)*gam(cc))

             kappa=2._rp*C_PI*C_RE**2._rp*C_ME*C_C**2._rp/ &
                  (params%cpp%length**2._rp*params%cpp%energy)
             BREM_p(cc)=-2._rp*ne(cc)*kappa*Zeff(cc)*(Zeff(cc)+1._rp)* &
                  C_a/C_PI*(gam(cc)-1._rp)*(log(2._rp*gam(cc))-1._rp/3._rp)

             !write(6,*) 'dpR',SC_p(cc)*dt
             !write(6,*) 'dpB',BREM_p(cc)*dt

             !write(6,*) 'dxiR',SC_xi(cc)*dt

             if (.not.FP_bremsstrahlung) BREM_p(cc)=0._rp

             dp(cc)=dp(cc)+(SC_p(cc)+BREM_p(cc))*dt* &
                  REAL(flagCol(cc))*REAL(flagCon(cc))
             dxi(cc)=dxi(cc)+(SC_xi(cc))*dt* &
                  REAL(flagCol(cc))*REAL(flagCon(cc))

          end do
          !$OMP END SIMD

       end if
    end if

#if DBG_CHECK
    do cc=1_idef,achunk
       if (dp(cc).gt.pm(cc)) then
          write(6,*) 'small angle collision'
          write(6,*) 'p0,xi0',pm0(cc),xi0(cc)
          write(6,*) 'p,xi',pm(cc),xi(cc)
          write(6,*) 'dp,dxi',dp(cc),dxi(cc)
          write(6,*) 'CBL',CBL(cc)
          write(6,*) 'v,ne,Te,Zeff',v(cc),ne(cc)*params%cpp%density,Te(cc)*params%cpp%temperature,Zeff(cc)
          write(6,*) 'ppll,pmu,Bmag',Ppll(cc),Pmu(cc),Bmag(cc)
          call korc_abort(24)
       endif
    end do
#endif

    !$OMP SIMD
    do cc=1_idef,achunk

       pm(cc)=pm(cc)+dp(cc)
       xi(cc)=xi(cc)+dxi(cc)
    end do
    !$OMP END SIMD

    do cc=1_idef,achunk
       !          if (pm(cc)<0) pm(cc)=-pm(cc)

       ! Keep xi between [-1,1]
       if (xi(cc)>1) then
          !             write(output_unit_write,'("High xi at: ",E17.10," with dxi: ",E17.10)') &
          !                  time*params%cpp%time, dxi(cc)
          xi(cc)=1-mod(xi(cc),1._rp)
       else if (xi(cc)<-1) then
          xi(cc)=-1-mod(xi(cc),-1._rp)
          !             write(output_unit_write,'("Low xi at: ",E17.10," with dxi: ",E17.10)') &
          !                  time*params%cpp%time, dxi(cc)
       endif

       if ((pm(cc).lt.min(cparams_ss%p_min*cparams_ss%pmin_scale, &
            p_therm)).and.flagCol(cc).eq.1_ip) then
          !             write(output_unit_write,'("Momentum less than zero")')
          !             stop
          !             write(output_unit_write,'("Particle not tracked at: ",E17.10," &
          !                  & with xi: ",E17.10)') time*params%cpp%time, xi(cc)
          flagCol(cc)=0_ip
       end if

    end do


    !       write(output_unit_write,'("rnd1: ",E17.10)') rnd1
    !       write(output_unit_write,'("flag: ",I16)') flag
    !       write(output_unit_write,'("CA: ",E17.10)') CAL
    !       write(output_unit_write,'("dCA: ",E17.10)') dCAL
    !       write(output_unit_write,'("CF ",E17.10)') CFL
    !       write(output_unit_write,'("CB: ",E17.10)') CBL
    !       write(output_unit_write,'("dp: ",E17.10)') dp
    !       write(output_unit_write,'("dxi: ",E17.10)') dxi
    !       write(output_unit_write,'("Ppll: ",E17.10)') Ppll
    !      write(output_unit_write,'("Pmu: ",E17.10)') Pmu
    !       write(output_unit_write,'("E_PHI_COL: ",E17.10)') E_PHI


    !       if (tt .EQ. 1_ip) then
    !          write(output_unit_write,'("dp_rad: ",E17.10)') &
    !               -gam(1)*pm(1)*(1-xi(1)*xi(1))/ &
    !               (cparams_ss%taur/Bmag(1)**2)*dt
    !          write(output_unit_write,'("dxi_rad: ",E17.10)') &
    !               xi(1)*(1-xi(1)*xi(1))/ &
    !               ((cparams_ss%taur/Bmag(1)**2)*gam(1))*dt
    !       end if

    !       if (tt .EQ. 1_ip) then
    !          write(output_unit_write,'("CA: ",E17.10)') CAL(1)
    !          write(output_unit_write,'("dCA: ",E17.10)') dCAL(1)
    !          write(output_unit_write,'("CF ",E17.10)') CFL(1)
    !          write(output_unit_write,'("CB: ",E17.10)') CBL(1)
    !       end if

#if DBG_CHECK
    do cc=1_idef,achunk
       if (ISNAN(xi(cc)).or.(abs(xi(cc)).gt.1._rp)) then
          write(6,*) 'xi is NaN or >1 before LAC'
          write(6,*) 'p0,xi0',pm0(cc),xi0(cc)
          write(6,*) 'p,xi',pm(cc),xi(cc)
          write(6,*) 'dp,dxi',dp(cc),dxi(cc)
          write(6,*) 'CBL',CBL(cc)
          write(6,*) 'v,ne,Te,Zeff',v(cc),ne(cc)*params%cpp%density,Te(cc)*params%cpp%temperature,Zeff(cc)
          write(6,*) 'ppll,pmu,Bmag',Ppll(cc),Pmu(cc),Bmag(cc)
          call korc_abort(24)
       end if
    end do
# endif

    if (cparams_ss%avalanche) then

       !$OMP SIMD
       do cc=1_idef,achunk
          ntot(cc)=ne(cc)
          if (params%bound_electron_model.eq.'HESSLOW') then

             if (.not.cparams_ms%lowKE_REs) then

                if ((cparams_ms%Zj(1).eq.0.0).and. &
                     (neut_prof.eq.'UNIFORM')) then
                   ntot(cc)=ntot(cc)+cparams_ms%nz(1)* &
                        (cparams_ms%Zo(1)-cparams_ms%Zj(1))
                else if ((cparams_ms%Zj(1).eq.0.0).and. &
                     (neut_prof.eq.'HOLLOW')) then
                   ntot(cc)=ntot(cc)+max(cparams_ms%nz(1)-ne(cc),0._rp)* &
                        (cparams_ms%Zo(1)-cparams_ms%Zj(1))
                else if ((cparams_ms%Zj(1).eq.0.0).and. &
                     (neut_prof.eq.'EDGE')) then
                   ra=sqrt((Y_R(cc)-P%R0)**2+(Y_Z(cc)-P%Z0)**2)/P%a
                   ntot(cc)=ntot(cc)+cparams_ms%nz(1)*ra**cparams_ms%neut_edge_fac* &
                        (cparams_ms%Zo(1)-cparams_ms%Zj(1))
                else
                   ntot(cc)=ntot(cc)+ne(cc)*cparams_ms%nz(1)/cparams_ms%ne* &
                        (cparams_ms%Zo(1)-cparams_ms%Zj(1))
                endif

                do ii=2,cparams_ms%num_impurity_species
                   ntot(cc)=ntot(cc)+ne(cc)*cparams_ms%nz(ii)/cparams_ms%ne* &
                        (cparams_ms%Zo(ii)-cparams_ms%Zj(ii))
                end do
             else
                if ((cparams_ms%Zj(1).eq.0.0).and. &
                     (neut_prof.eq.'UNIFORM')) then
                   ntot(cc)=ntot(cc)+cparams_ms%nz(1)* &
                        (cparams_ms%Zo(1)-cparams_ms%Zj(1) &
                              -cparams_ms%lowKE_LAC_not_ionized)
                else if ((cparams_ms%Zj(1).eq.0.0).and. &
                     (neut_prof.eq.'HOLLOW')) then
                   ntot(cc)=ntot(cc)+max(cparams_ms%nz(1)-ne(cc),0._rp)* &
                        (cparams_ms%Zo(1)-cparams_ms%Zj(1) &
                              -cparams_ms%lowKE_LAC_not_ionized)
                else if ((cparams_ms%Zj(1).eq.0.0).and. &
                     (neut_prof.eq.'EDGE')) then
                   ntot(cc)=ntot(cc)+cparams_ms%nz(1)*ra**cparams_ms%neut_edge_fac* &
                        (cparams_ms%Zo(1)-cparams_ms%Zj(1) &
                        -cparams_ms%lowKE_LAC_not_ionized)
                else
                   ntot(cc)=ntot(cc)+ne(cc)*cparams_ms%nz(1)/cparams_ms%ne* &
                        (cparams_ms%Zo(1)-cparams_ms%Zj(1) &
                              -cparams_ms%lowKE_LAC_not_ionized)
                endif

                do ii=2,cparams_ms%num_impurity_species
                   ntot(cc)=ntot(cc)+ne(cc)*cparams_ms%nz(ii)/cparams_ms%ne* &
                        (cparams_ms%Zo(ii)-cparams_ms%Zj(ii) &
                              -cparams_ms%lowKE_LAC_not_ionized)
                end do
             endif

          end if
       end do
       !$OMP END SIMD

       !write(6,*) 'ntot',ntot*params%cpp%density

       call large_angle_source(spp,params,achunk,F,Y_R,Y_PHI,Y_Z, &
            pm,xi,ne,ntot,Te,Bmag,E_PHI_LAC,me,flagCol,flagCon,B_R,B_PHI,B_Z)

#if DBG_CHECK
       do cc=1_idef,achunk
          if (abs(pm(cc)-pm0(cc)).gt.pm0(cc)) then
             write(6,*) 'large angle collision'
             write(6,*) 'p0,xi0',pm0(cc),xi0(cc)
             write(6,*) 'p,xi',pm(cc),xi(cc)
             write(6,*) 'dp,dxi',dp(cc),dxi(cc)
             write(6,*) 'CBL',CBL(cc)
             write(6,*) 'v,ne,Te,Zeff',v(cc),ne(cc)*params%cpp%density,Te(cc)*params%cpp%temperature/params%cpp%charge,Zeff(cc)
             write(6,*) 'ppll,pmu,Bmag',Ppll(cc),Pmu(cc),Bmag(cc)
             call korc_abort(24)
          endif
       end do
#endif

    end if

#if DBG_CHECK
    do cc=1_idef,achunk
       if (ISNAN(xi(cc)).or.(abs(xi(cc)).gt.1._rp)) then
          write(6,*) 'xi is NaN or >1 after LAC'
          write(6,*) 'p0,xi0',pm0(cc),xi0(cc)
          write(6,*) 'p,xi',pm(cc),xi(cc)
          write(6,*) 'dp,dxi',dp(cc),dxi(cc)
          write(6,*) 'CBL',CBL(cc)
          write(6,*) 'v,ne,Te,Zeff',v(cc),ne(cc)*params%cpp%density,Te(cc)*params%cpp%temperature,Zeff(cc)
          call korc_abort(24)
       end if
    end do
# endif

    !$OMP SIMD
    do cc=1_idef,achunk
       ! Transform P,xi to p_pll,mu
       Ppll(cc)=pm(cc)*xi(cc)
       Pmu(cc)=(pm(cc)*pm(cc)-Ppll(cc)*Ppll(cc))/(2*me*Bmag(cc))
    end do
    !$OMP END SIMD

#if DBG_CHECK
    do cc=1_idef,achunk
       if (Pmu(cc).lt.0._rp) then
          write(6,*) 'mu is negative'
          write(6,*) 'p,xi',pm(cc),xi(cc)
          write(6,*) 'V_PLL,V_MU',Ppll(cc),Pmu(cc)
          call korc_abort(24)
       endif
    end do
#endif

    E_PHI=E_PHI_LAC

  end subroutine include_CoulombCollisionsLA_GC_p

#ifdef FIO
  subroutine include_CoulombCollisions_GCfio_p(tt,params,Y_R,Y_PHI,Y_Z, &
       Ppll,Pmu,me,flagCon,flagCol,F,P,E_PHI,ne,ni,Te,Zeff,nimp,PSIp,hint)

    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    TYPE(KORC_PARAMS), INTENT(INOUT) 		:: params
    REAL(rp), DIMENSION(params%pchunk), INTENT(INOUT) 	:: Ppll
    REAL(rp), DIMENSION(params%pchunk), INTENT(INOUT) 	:: Pmu
    REAL(rp), DIMENSION(params%pchunk) 			:: Bmag
    REAL(rp), DIMENSION(params%pchunk) 	:: B_R,B_PHI,B_Z
    REAL(rp), DIMENSION(params%pchunk) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp), DIMENSION(params%pchunk) :: gradB_R,gradB_PHI,gradB_Z
    REAL(rp), DIMENSION(params%pchunk) 	:: E_R,E_Z
    REAL(rp), DIMENSION(params%pchunk) 	:: E_PHI,PSIp,E_PHI0
    REAL(rp), DIMENSION(params%pchunk), INTENT(INOUT) 	:: ne,Te,ni
    REAL(rp), DIMENSION(params%pchunk,params%num_impurity_species),&
         & INTENT(INOUT) 	:: nimp
    REAL(rp), DIMENSION(params%pchunk), INTENT(IN) 	:: Y_R,Y_PHI&
         &,Y_Z
    INTEGER(is), DIMENSION(params%pchunk), INTENT(INOUT)  :: flagCol
    INTEGER(is), DIMENSION(params%pchunk), INTENT(INOUT)  :: flagCon
    TYPE(C_PTR), DIMENSION(params%pchunk), INTENT(INOUT)  :: hint
    REAL(rp), INTENT(IN) 			:: me
    REAL(rp), DIMENSION(params%pchunk), INTENT(OUT) 	:: Zeff
    REAL(rp), DIMENSION(params%pchunk,2) 			:: dW
    REAL(rp), DIMENSION(params%pchunk,2) 			:: rnd1
    REAL(rp) 					:: dt,time
    REAL(rp), DIMENSION(params%pchunk) 	:: pm
    REAL(rp), DIMENSION(params%pchunk)  :: dp
    REAL(rp), DIMENSION(params%pchunk)  :: xi
    REAL(rp), DIMENSION(params%pchunk)  :: dxi
    REAL(rp), DIMENSION(params%pchunk)  :: v,gam
    !! speed of particle
    REAL(rp), DIMENSION(params%pchunk) 	:: CAL
    REAL(rp) , DIMENSION(params%pchunk)	:: dCAL
    REAL(rp), DIMENSION(params%pchunk) 	:: CFL
    REAL(rp), DIMENSION(params%pchunk) 	:: CBL
    REAL(rp), DIMENSION(params%pchunk) 	:: SC_p,SC_mu,BREM_p
    REAL(rp) 					:: kappa
    integer :: cc,pchunk,ii
    integer(ip),INTENT(IN) :: tt

    pchunk=params%pchunk

    if (MODULO(params%it+tt,cparams_ss%subcycling_iterations) .EQ. 0_ip) then
       dt = REAL(cparams_ss%subcycling_iterations,rp)*params%dt

       call get_fio_GCmagnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
            B_R,B_PHI,B_Z,gradB_R,gradB_PHI,gradB_Z, &
            curlb_R,curlb_PHI,curlb_Z,flagCon,hint)

       if (F%FIO_E .ge. 0) then
          call get_fio_GCelectric_fields_p(params,F, &
               Y_R,Y_PHI,Y_Z,E_R,E_PHI,E_Z,flagCon,hint)
       end if
       call get_fio_vector_potential_p(params,F,Y_R,Y_PHI,Y_Z, &
            PSIp,flagCon,hint)

       call get_fio_profile_p(params,P,Y_R,Y_PHI,Y_Z, &
            ne,Te,flagCon,hint)

       call get_fio_ion_p(params,P,Y_R,Y_PHI,Y_Z, &
            ne,ni,nimp,Zeff,flagCon,hint)


       !write(6,*) ne,Te,nimp,Zeff

       !$OMP SIMD
!       !$OMP& aligned (pm,xi,v,Ppll,Bmag,Pmu)
       do cc=1_idef,pchunk
          Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))
          ! Transform p_pll,mu to P,eta
          pm(cc) = SQRT(Ppll(cc)*Ppll(cc)+2*me*Bmag(cc)*Pmu(cc))
          xi(cc) = Ppll(cc)/pm(cc)

          gam(cc) = sqrt(1+pm(cc)*pm(cc))

          v(cc) = pm(cc)/gam(cc)
          ! normalized speed (v_K=v_P/c)

          E_PHI0(cc)=E_PHI(cc)
       end do
       !$OMP END SIMD

       !write(6,*) 'R',Y_R*params%cpp%length
       !write(6,*) 'PHI',Y_PHI
       !write(6,*) 'Z',Y_Z*params%cpp%length
       !write(6,*) 'ne',ne*params%cpp%density
       !write(6,*) 'ni',ni*params%cpp%density
       !write(6,*) 'nimp',nimp(:,1:2)*params%cpp%density
       !write(6,*) 'Zeff',Zeff

       if (.not.params%FokPlan) E_PHI=0._rp

!       write(output_unit_write,'("ne: "E17.10)') ne
!       write(output_unit_write,'("Te: "E17.10)') Te
!       write(output_unit_write,'("Bmag: "E17.10)') Bmag
!       write(output_unit_write,'("v: ",E17.10)') v
!       write(output_unit_write,'("xi: ",E17.10)') xi
       !       write(output_unit_write,'("size(E_PHI_GC): ",I16)') size(E_PHI)


       !$OMP SIMD
!       !$OMP& aligned(rnd1,dW,CAL,dCAL,CFL,CBL,v,ne,Te,Zeff,dp, &
!       !$OMP& flagCon,flagCol,dxi,xi,pm,Ppll,Pmu,Bmag)
       do cc=1_idef,pchunk

#ifdef PARALLEL_RANDOM
          rnd1(cc,1) = get_random()
          rnd1(cc,2) = get_random()
          !       rnd1(:,1) = get_random_mkl()
          !       rnd1(:,2) = get_random_mkl()
#else
          call RANDOM_NUMBER(rnd1)
#endif

          dW(cc,1) = SQRT(3*dt)*(-1+2*rnd1(cc,1))
          dW(cc,2) = SQRT(3*dt)*(-1+2*rnd1(cc,2))

!          write(output_unit_write,'("dW1: ",E17.10)') dW(cc,1)
!          write(output_unit_write,'("dW2: ",E17.10)') dW(cc,2)

          CAL(cc) = CA_SD(v(cc),ne(cc),Te(cc))
          dCAL(cc)= dCA_SD(v(cc),me,ne(cc),Te(cc))
          CFL(cc) = CF_SD_FIO(params,v(cc),ne(cc),Te(cc),nimp(cc,:))
          CBL(cc) = (CB_ee_SD(v(cc),ne(cc),Te(cc),Zeff(cc))+ &
               CB_ei_SD_FIO(params,v(cc),ne(cc),Te(cc),nimp(cc,:),Zeff(cc)))


          dp(cc)=REAL(flagCol(cc))*REAL(flagCon(cc))* &
               ((-CFL(cc)+dCAL(cc)+E_PHI(cc)*xi(cc))*dt+ &
               sqrt(2.0_rp*CAL(cc))*dW(cc,1))

          dxi(cc)=REAL(flagCol(cc))*REAL(flagCon(cc))* &
               ((-2*xi(cc)*CBL(cc)/(pm(cc)*pm(cc))+ &
               E_PHI(cc)*(1-xi(cc)*xi(cc))/pm(cc))*dt- &
               sqrt(2.0_rp*CBL(cc)*(1-xi(cc)*xi(cc)))/pm(cc)*dW(cc,2))

!          write(output_unit_write,'("dp: ",E17.10)') dp(cc)
!          write(output_unit_write,'("dxi: ",E17.10)') dxi(cc)

       end do
       !$OMP END SIMD

       if (params%FokPlan.and.params%radiation) then
          if(params%GC_rad_model.eq.'SDE') then

             !$OMP SIMD
             do cc=1_idef,pchunk

                SC_p(cc)=-gam(cc)*pm(cc)*(1-xi(cc)*xi(cc))/ &
                     (cparams_ss%taur/Bmag(cc)**2)
                SC_mu(cc)=xi(cc)*(1-xi(cc)*xi(cc))/ &
                     ((cparams_ss%taur/Bmag(cc)**2)*gam(cc))

                kappa=2._rp*C_PI*C_RE**2._rp*C_ME*C_C**2._rp
                BREM_p(cc)=-2._rp*ne(cc)*kappa*Zeff(cc)*(Zeff(cc)+1._rp)* &
                     C_a/C_PI*(gam(cc)-1._rp)*(log(2._rp*gam(cc))-1._rp/3._rp)


                dp(cc)=dp(cc)+(SC_p(cc)+BREM_p(cc))*dt* &
                     REAL(flagCol(cc))*REAL(flagCon(cc))
                dxi(cc)=dxi(cc)+(SC_mu(cc))*dt* &
                     REAL(flagCol(cc))*REAL(flagCon(cc))

             end do
             !$OMP END SIMD

          end if
       end if

       !$OMP SIMD
       do cc=1_idef,pchunk

          pm(cc)=pm(cc)+dp(cc)
          xi(cc)=xi(cc)+dxi(cc)

!          if (pm(cc)<0) pm(cc)=-pm(cc)

          ! Keep xi between [-1,1]
          if (xi(cc)>1) then
!             write(output_unit_write,'("High xi at: ",E17.10," with dxi: ",E17.10)') &
!                  time*params%cpp%time, dxi(cc)
             xi(cc)=1-mod(xi(cc),1._rp)
          else if (xi(cc)<-1) then
             xi(cc)=-1-mod(xi(cc),-1._rp)
!             write(output_unit_write,'("Low xi at: ",E17.10," with dxi: ",E17.10)') &
!                  time*params%cpp%time, dxi(cc)
          endif

          ! Transform P,xi to p_pll,mu
          Ppll(cc)=pm(cc)*xi(cc)
          Pmu(cc)=(pm(cc)*pm(cc)-Ppll(cc)*Ppll(cc))/(2*me*Bmag(cc))
       end do
       !$OMP END SIMD

#if DBG_CHECK
       do cc=1_idef,pchunk
          if((isnan(Ppll(cc)).or.isnan(Pmu(cc))).and. &
             ((flagCol(cc).eq.1_is).and.(flagCon(cc).eq.1_is))) then
             write(6,*) 'End collision'
             write(6,*) 'Ppll',Ppll(cc)
             write(6,*) 'Pmu',Pmu(cc)
             write(6,*) 'Bmag',Bmag(cc)
             write(6,*) 'pm',pm(cc)
             write(6,*) 'xi',xi(cc)
             write(6,*) 'dp',dp(cc)
             write(6,*) 'dxi',dxi(cc)
             write(6,*) 'CFL',CFL(cc)
             write(6,*) 'CBL',CBL(cc)
             write(6,*) 'dCAL',dCAL(cc)
             write(6,*) 'CAL',CAL(cc)
             write(6,*) 'v',v(cc)
             write(6,*) 'ne',ne(cc)
             write(6,*) 'Te',Te(cc)
             write(6,*) 'nimp',nimp(cc,:)
             write(6,*) 'Zeff',Zeff(cc)
             write(6,*) 'Y_R',Y_R(cc)
             write(6,*) 'Y_PHI',Y_PHI(cc)
             write(6,*) 'Y_Z',Y_Z(cc)

             stop 'Ppll or Pmu is NaN'
          endif
       end do
#endif

!       write(output_unit_write,'("rnd1: ",E17.10)') rnd1
!       write(output_unit_write,'("flag: ",I16)') flag
!       write(output_unit_write,'("CA: ",E17.10)') CAL
!       write(output_unit_write,'("dCA: ",E17.10)') dCAL
!       write(output_unit_write,'("CF ",E17.10)') CFL
!       write(output_unit_write,'("CB: ",E17.10)') CBL
!       write(output_unit_write,'("dp: ",E17.10)') dp
!       write(output_unit_write,'("dxi: ",E17.10)') dxi
!       write(output_unit_write,'("Ppll: ",E17.10)') Ppll
!      write(output_unit_write,'("Pmu: ",E17.10)') Pmu
!       write(output_unit_write,'("E_PHI_COL: ",E17.10)') E_PHI

       do cc=1_idef,pchunk
          if ((pm(cc).lt.min(cparams_ss%p_min*cparams_ss%pmin_scale, &
               p_therm)).and.flagCol(cc).eq.1_ip) then
!             write(output_unit_write,'("Momentum less than zero")')
             !             stop
!             write(output_unit_write,'("Particle not tracked at: ",E17.10," &
!                  & with xi: ",E17.10)') time*params%cpp%time, xi(cc)
             flagCol(cc)=0_ip
          end if
       end do

!       if (tt .EQ. 1_ip) then
!          write(output_unit_write,'("dp_rad: ",E17.10)') &
!               -gam(1)*pm(1)*(1-xi(1)*xi(1))/ &
!               (cparams_ss%taur/Bmag(1)**2)*dt
!          write(output_unit_write,'("dxi_rad: ",E17.10)') &
!               xi(1)*(1-xi(1)*xi(1))/ &
!               ((cparams_ss%taur/Bmag(1)**2)*gam(1))*dt
!       end if

!       if (tt .EQ. 1_ip) then
!          write(output_unit_write,'("CA: ",E17.10)') CAL(1)
!          write(output_unit_write,'("dCA: ",E17.10)') dCAL(1)
!          write(output_unit_write,'("CF ",E17.10)') CFL(1)
!          write(output_unit_write,'("CB: ",E17.10)') CBL(1)
!       end if

       if (.not.params%FokPlan) E_PHI=E_PHI0

    end if

  end subroutine include_CoulombCollisions_GCfio_p
#endif

  subroutine large_angle_source(spp,params,achunk,F,Y_R,Y_PHI,Y_Z, &
       pm,xi,ne,netot,Te,Bmag,E_PHI,me,flagCol,flagCon,B_R,B_PHI,B_Z)
    TYPE(SPECIES), INTENT(INOUT)    :: spp
    TYPE(KORC_PARAMS), INTENT(IN) 			:: params
    TYPE(FIELDS), INTENT(IN)                                   :: F
    INTEGER, INTENT(IN) :: achunk
    REAL(rp), INTENT(INOUT), DIMENSION(achunk)  :: pm,xi
    REAL(rp), INTENT(IN), DIMENSION(achunk)  :: Y_R,Y_PHI,Y_Z
    REAL(rp), INTENT(IN), DIMENSION(achunk)  :: B_R,B_PHI,B_Z
    REAL(rp), INTENT(IN), DIMENSION(achunk)  :: ne,netot,Te
    REAL(rp), INTENT(IN), DIMENSION(achunk)  :: Bmag,E_PHI
    INTEGER(is), INTENT(IN), DIMENSION(achunk)  :: flagCol,flagCon
    REAL(rp), INTENT(IN)  :: me
    REAL(rp), DIMENSION(achunk)  :: gam,prob0,prob1,pm0,xi0,gam0
    REAL(rp) :: gam_min,p_min,gammax,dt,gamsecmax,psecmax,ptrial
    REAL(rp) :: gamtrial,cosgam1,xirad,xip,xim,xitrial,sinsq1,cossq1,pitchprob1
    REAL(rp) :: dsigdgam1,S_LAmax,S_LA1,tmppm,gamvth,vmin,E_C,p_c,gam_c,pRE
    INTEGER :: ngam1,neta1
    INTEGER :: ii,jj,cc,seciter
    REAL(rp), DIMENSION(cparams_ss%ngrid1) :: gam1,pm1,tmpgam1,tmpcosgam,tmpdsigdgam,tmpsecthreshgam,probtmp,intpitchprob
    REAL(rp), DIMENSION(cparams_ss%ngrid1-1) :: dpm1
    REAL(rp), DIMENSION(cparams_ss%ngrid1) :: eta1,tmpsinsq,tmpcossq
    REAL(rp), DIMENSION(cparams_ss%ngrid1-1) :: deta1
    REAL(rp), DIMENSION(cparams_ss%ngrid1,cparams_ss%ngrid1) :: cosgam,sinsq,cossq,tmpcossq1,pitchprob,dsigdgam
    REAL(rp), DIMENSION(cparams_ss%ngrid1,cparams_ss%ngrid1) :: secthreshgam,pm11,gam11,eta11,S_LA,pitchrad
    LOGICAL :: accepted


    dt=cparams_ss%coll_per_dump_dt*params%cpp%time
    ngam1=cparams_ss%ngrid1
    neta1=cparams_ss%ngrid1


    !$OMP SIMD
    do cc=1_idef,achunk
       pm0(cc)=pm(cc)
       xi0(cc)=xi(cc)

       gam(cc) = sqrt(1+pm(cc)*pm(cc))
       gam0(cc)=gam(cc)

#ifdef PARALLEL_RANDOM
       prob0(cc) = get_random()
#else
       call RANDOM_NUMBER(prob0)
#endif

    end do
    !$OMP END SIMD

    vmin=1/sqrt(1+1/(cparams_ss%p_min*cparams_ss%pmin_scale)**2)

    !! For each primary RE, calculating probability to generate a secondary RE

    do cc=1_idef,achunk

       if ((flagCol(cc).lt.1).or.(flagCon(cc).lt.1)) cycle

       if (.not.(cparams_ms%lowKE_REs)) then
          E_C=netot(cc)/ne(cc)*Gammacee(vmin,ne(cc),Te(cc))
       else
          E_C=Gammacee(vmin,ne(cc),Te(cc))
       end if



       !write(6,*) 'E',E_PHI*params%cpp%Eo
       !write(6,*) 'E_C',E_C*params%cpp%Eo
       !write(6,*) 'E_c,min',cparams_ms%Ec_min*params%cpp%Eo
       !write(6,*) 'ne',ne(cc)*params%cpp%density
       !write(6,*) 'netot',netot(cc)*params%cpp%density
       !write(6,*) 'Te',Te(cc)*params%cpp%temperature
       !write(6,*) 'Clog',CLogee_wu(params,ne(cc)*params%cpp%density,Te(cc)*params%cpp%temperature)

       if (.not.(cparams_ss%always_aval)) then
          if (E_C.gt.abs(E_PHI(cc))) cycle

          p_c=cparams_ss%pmin_scale/sqrt(abs(E_PHI(cc))/E_C-1)
          gam_c=sqrt(1+p_c**2)

          if(cparams_ss%min_secRE.eq.'THERM') then
             gam_min=(gam_c+1)/2
             p_min=sqrt(gam_min**2-1)
          else
             gam_min=gam_c
             p_min=p_c
          end if
       else
          p_min=cparams_ss%p_min*cparams_ss%pmin_scale
          gam_min=sqrt(1+p_min**2)
       endif

       gammax=(gam(cc)+1._rp)/2._rp

       if (gam_min.eq.1._rp) then
          write(6,*) 'R',Y_R(cc)*params%cpp%length,'Z',Y_Z(cc)*params%cpp%length
          write(6,*) 'vmin',vmin,'netot',netot(cc)*params%cpp%density,'Te',Te(cc)*params%cpp%temperature/C_E
          write(6,*) 'E',E_PHI*params%cpp%Eo,'E_c',E_C*params%cpp%Eo
          write(6,*) 'p_c',p_c,'gam_c',gam_c
          write(6,*) 'p_min',p_min,'gam_min',gam_min
          !write(6,*) 'LAC_gam_resolution: ',TRIM(cparams_ss%LAC_gam_resolution)
          !write(6,*) 'gam_min,gammax',gam_min,gammax
       end if

       !! Generating 1D and 2D ranges for secondary RE distribution

       if (TRIM(cparams_ss%LAC_gam_resolution).eq.'LIN') then
          do ii=1,ngam1
             gam1(ii)=gam_min+(gammax-gam_min)* &
                  REAL(ii-1)/REAL(ngam1-1)
          end do
       elseif (TRIM(cparams_ss%LAC_gam_resolution).eq.'EXP') then
          do ii=1,ngam1
             tmpgam1(ii)=log10(gam_min)+ &
                  (log10(gammax)-log10(gam_min))* &
                  REAL(ii-1)/REAL(ngam1-1)
          end do
          gam1=10**(tmpgam1)
       elseif (TRIM(cparams_ss%LAC_gam_resolution).eq.'2EXP') then
          do ii=1,ngam1
             tmpgam1(ii)=log10(log10(gam_min))+ &
                  (log10(log10(gammax))-log10(log10(gam_min)))* &
                  REAL(ii-1)/REAL(ngam1-1)
          end do
          gam1=10**(10**(tmpgam1))
       end if

       !write(6,*) 'tmpgam1',tmpgam1
       !write(6,*) 'gam1',gam1


       pm1=sqrt(gam1**2-1)

       do ii=1,ngam1-1
          dpm1(ii)=pm1(ii+1)-pm1(ii)
       end do

       do ii=1,neta1
          eta1(ii)=C_PI*(ii-1)/(neta1-1)
       end do

       !write(6,*) 'eta1',eta1

       do ii=1,neta1-1
          deta1(ii)=eta1(ii+1)-eta1(ii)
       end do

       do ii=1,neta1
          pm11(:,ii)=pm1
          gam11(:,ii)=gam1
       end do

       do ii=1,ngam1
          eta11(ii,:)=eta1
       end do


       tmpcosgam=sqrt(((gam(cc)+1)*(gam1-1))/((gam(cc)-1)*(gam1+1)))
       tmpdsigdgam=2*C_PI*C_RE**2/(gam(cc)**2-1)* &
            (((gam(cc)-1)**2*gam(cc)**2)/((gam1-1)**2*(gam(cc)-gam1)**2)- &
            (2*gam(cc)**2+2*gam(cc)-1)/((gam1-1)*(gam(cc)-gam1))+1)
       tmpsecthreshgam=1._rp
       where(gam1.gt.(gam(cc)+1)/2._rp) tmpsecthreshgam=0._rp

       !write(6,*) 'tmpcosgam',tmpcosgam
       !write(6,*) 'tmpdsigdgam',tmpdsigdgam
       !write(6,*) 'tmpsecthreshgam',tmpsecthreshgam


       do ii=1,neta1
          cosgam(:,ii)=tmpcosgam
          dsigdgam(:,ii)=tmpdsigdgam
          secthreshgam(:,ii)=tmpsecthreshgam
       end do

       !if (cc.eq.1) then
          !write(6,*) cosgam
       !end if

       tmpsinsq=(1-xi(cc)**2)*sin(eta1)**2
       tmpcossq=xi(cc)*cos(eta1)

       do ii=1,ngam1
          sinsq(ii,:)=tmpsinsq
          tmpcossq1(ii,:)=tmpcossq
       end do

       !if (cc.eq.1) then
          !write(6,*) sinsq
       !end if

       cossq=(cosgam-tmpcossq1)**2

       pitchrad=sinsq-cossq
       where(pitchrad.lt.0) pitchrad=tiny(0._rp)

       !if (cc.eq.1) then
          !write(6,*) cossq
          !write(6,*) 'sinsq-cossq',sqrt(sinsq-cossq)
          !write(6,*) 'pitchrad',pitchrad
       !end if

       pitchprob=1/(C_PI*sqrt(pitchrad))

       where(pitchprob.eq.1/(C_PI*sqrt(tiny(0._rp)))) pitchprob=0._rp

       !if (cc.eq.1) then
          !write(6,*) 'pitchprob',pitchprob
       !end if

       S_LA=netot(cc)*params%cpp%density*C_C/(2*C_PI)* &
            (pm11/gam11)*(pm(cc)/gam(cc))* &
            pitchprob*dsigdgam*secthreshgam

       !! Saving maximum secondary RE source for use in rejection-acceptance
       !! sampling algorithm
       S_LAmax=maxval(S_LA)

       S_LA=S_LA*sin(eta11)

       !! Trapezoidal integration of secondary RE source to find probabilty

       do ii=1,ngam1
          probtmp(ii)=S_LA(ii,1)*deta1(1)/2+S_LA(ii,neta1)*deta1(neta1-1)/2
          !intpitchprob(ii)=pitchprob(ii,1)*sin(eta1(1))*deta1(1)/2+ &
          !     pitchprob(ii,neta1)*sin(eta1(neta1))*deta1(neta1-1)/2
          do jj=2,neta1-1
             probtmp(ii)=probtmp(ii)+S_LA(ii,jj)*(deta1(jj)+deta1(jj-1))/2
          !   intpitchprob(ii)=intpitchprob(ii)+ &
          !        pitchprob(ii,jj)*sin(eta1(jj))*(deta1(jj)+deta1(jj-1))/2
          end do
       end do

       prob1(cc)=probtmp(1)*dpm1(1)/2+probtmp(ngam1)*dpm1(ngam1-1)/2
       do jj=2,ngam1-1
          prob1(cc)=prob1(cc)+probtmp(jj)*(dpm1(jj)+dpm1(jj-1))/2
       end do

       !write(6,*) 'prob1pre',prob1(cc),'flagCol',flagCol(cc),'flagCon',flagCon(cc),'dt',dt

       !write(6,*) 'intpitchprob',intpitchprob

       prob1(cc)=prob1(cc)*dt*2*C_PI

       if (ISNAN(prob1(cc))) then
          write(6,*) 'NaN probability from secondary RE source'
          write(6,*) 'p,xi',pm(cc),xi(cc)
          write(6,*) 'gam_min,gammax',gam_min,gammax
          write(6,*) 'E',E_PHI(cc)*params%cpp%Eo
          write(6,*) 'E_C',E_C*params%cpp%Eo
          !write(6,*) 'pitchprob',pitchprob
          !write(6,*) 'S_LA',S_LA
          call korc_abort(24)
       end if

       if (prob1(cc).gt.1._rp) then
          write(6,*) 'Multiple secondary REs generated in a collisional time step'
          write(6,*) 'p,xi',pm(cc),xi(cc)
          write(6,*) 'gam_min,gammax',gam_min,gammax
          write(6,*) 'E',E_PHI(cc)*params%cpp%Eo
          write(6,*) 'E_C',E_C*params%cpp%Eo
          call korc_abort(24)
       end if

       !write(6,*) 'gam',gam(cc),'xi',xi(cc)
       !write(6,*) 'prob1',prob1(cc),'prob0',prob0(cc)

       if (prob1(cc).gt.prob0(cc)) then

          !! If secondary RE generated, begin acceptance-rejection sampling
          !! algorithm

          !write(6,*) 'secondary RE from ',prob1,prob0

          accepted=.false.

          gamsecmax=(gam(cc)+1)/2
          psecmax=sqrt(gamsecmax**2-1)

          seciter=0
          do while (.not.accepted)

             ptrial=p_min+(psecmax-p_min)*get_random()
             gamtrial=sqrt(1+ptrial*ptrial)

             cosgam1=sqrt(((gam(cc)+1)*(gamtrial-1))/((gam(cc)-1)*(gamtrial+1)))

             !xirad=sqrt((cosgam1*xi(cc))**2-(xi(cc)**2+cosgam1**2-1))

             !if (isnan(xirad)) then
             !   write(6,*) 'Sample not in allowable region of phase space'
             !   call korc_abort(24)
             !end if

             !xip=cosgam1*xi(cc)+xirad
             !xim=cosgam1*xi(cc)-xirad

             !xitrial=xim+(xip-xim)*get_random()
             xitrial=-1+2*get_random()

             sinsq1=(1-xi(cc)*xi(cc))*(1-xitrial*xitrial)
             cossq1=(cosgam1-xi(cc)*xitrial)**2

             if ((sinsq1-cossq1).le.0._rp) cycle

             pitchprob1=1/(C_PI*sqrt(sinsq1-cossq1))

             !if (isnan(pitchprob1)) cycle

             dsigdgam1=2*C_PI*C_RE**2/(gam(cc)**2-1)* &
                  (((gam(cc)-1)**2*gam(cc)**2)/ &
                  ((gamtrial-1)**2*(gam(cc)-gamtrial)**2)- &
                  (2*gam(cc)**2+2*gam(cc)-1)/ &
                  ((gamtrial-1)*(gam(cc)-gamtrial))+1)

             S_LA1=netot(cc)*params%cpp%density*C_C/(2*C_PI)* &
                  (ptrial/gamtrial)*(pm(cc)/gam(cc))* &
                  pitchprob1*dsigdgam1

             if (S_LA1/S_LAmax.gt.get_random()) accepted=.true.

             seciter=seciter+1

             !if (mod(seciter,100).eq.0) then
             !   write(6,*) 'iteration',seciter
             !end if

          end do

          !! Write secondary RE degrees of freedom to particle derived type


          !$OMP ATOMIC UPDATE
          spp%pRE=spp%pRE+1
          !$OMP ATOMIC WRITE
          spp%vars%flagRE(spp%pRE)=1
          !$OMP ATOMIC WRITE
          spp%vars%flagCon(spp%pRE)=1
          !$OMP ATOMIC WRITE
          spp%vars%flagCol(spp%pRE)=1
          !$OMP ATOMIC WRITE
          spp%vars%V(spp%pRE,1)=ptrial*xitrial
          !$OMP ATOMIC WRITE
          spp%vars%V(spp%pRE,2)=ptrial*ptrial*(1-xitrial*xitrial)/ &
               (2*me*Bmag(cc))
          !$OMP ATOMIC WRITE
          spp%vars%Y(spp%pRE,1)=Y_R(cc)
          !$OMP ATOMIC WRITE
          spp%vars%Y(spp%pRE,2)=Y_PHI(cc)
          !$OMP ATOMIC WRITE
          spp%vars%Y(spp%pRE,3)=Y_Z(cc)
          !$OMP ATOMIC WRITE
          spp%vars%Yborn(spp%pRE,1)=Y_R(cc)
          !$OMP ATOMIC WRITE
          spp%vars%Yborn(spp%pRE,2)=Y_PHI(cc)
          !$OMP ATOMIC WRITE
          spp%vars%Yborn(spp%pRE,3)=Y_Z(cc)
          !$OMP ATOMIC WRITE
          spp%vars%B(spp%pRE,1)=B_R(cc)
          !$OMP ATOMIC WRITE
          spp%vars%B(spp%pRE,2)=B_PHI(cc)
          !$OMP ATOMIC WRITE
          spp%vars%B(spp%pRE,3)=B_Z(cc)


          !! Write changes to primary RE degrees of freedom to temporary
          !! arrays for passing back out to particle derived type

#if DBG_CHECK
          if (spp%vars%V(spp%pRE,2).lt.0._rp) then
             write(6,*) 'mu is negative for secondary'
             write(6,*) 'ppll,mu',spp%vars%V(spp%pRE,1),spp%vars%V(spp%pRE,2)
             write(6,*) 'pm,xi',ptrial,xitrial
             write(6,*) 'xim,xip',xim,xip
             write(6,*) 'pmin,psecmax',p_min,psecmax
             call korc_abort(24)
          end if
# endif

          if (cparams_ss%ConserveLA) then
             tmppm=pm(cc)
             gamvth=1/sqrt(1-2*Te(cc))
             gam(cc)=gam(cc)-gamtrial+gamvth
             pm(cc)=sqrt(gam(cc)*gam(cc)-1)
             xi(cc)=(tmppm*xi(cc)-ptrial*xitrial)/pm(cc)
          end if

#if DBG_CHECK
          if (abs(pm(cc)-pm0(cc)).gt.pm0(cc)) then
             write(6,*) 'Conservative update'
             write(6,*) 'p0,gam0,xi0',pm0(cc),gam0(cc),xi0(cc)
             write(6,*) 'ptrial,gamtrial,xitrial',ptrial,gamtrial,xitrial
             write(6,*) 'gamvth,Te',gamvth,Te(cc)*params%cpp%temperature/params%cpp%charge
             write(6,*) 'p,xi',pm(cc),xi(cc)
             call korc_abort(24)
          endif

#endif

          !$OMP ATOMIC READ
          pRE=spp%pRE

          if (pRE.eq.spp%ppp) then
             write(6,*) 'All REs allocated on proc',params%mpi_params%rank
             call korc_abort(24)
          end if

       end if

    end do



  end subroutine large_angle_source


  subroutine save_params_ms(params)
    TYPE(KORC_PARAMS), INTENT(IN) 			:: params
    CHARACTER(MAX_STRING_LENGTH) 			:: filename
    CHARACTER(MAX_STRING_LENGTH) 			:: gname
    CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
    CHARACTER(MAX_STRING_LENGTH) 			:: dset
    CHARACTER(MAX_STRING_LENGTH) 			:: attr
    INTEGER(HID_T) 					:: h5file_id
    INTEGER(HID_T) 					:: group_id
    INTEGER 						:: h5error
    REAL(rp) 						:: units

    if (params%mpi_params%rank .EQ. 0) then
       filename = TRIM(params%path_to_outputs) // "simulation_parameters.h5"
       call h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, h5file_id, h5error)

       gname = "collisions_ms"
       call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

       ALLOCATE(attr_array(cparams_ms%num_impurity_species))

       dset = TRIM(gname) // "/model"
       call save_string_parameter(h5file_id,dset,(/params%collisions_model/))

       dset = TRIM(gname) // "/num_impurity_species"
       attr = "Number of impurity species"
       call save_to_hdf5(h5file_id,dset,cparams_ms%num_impurity_species,attr)

       dset = TRIM(gname) // "/Te"
       attr = "Background electron temperature in eV"
       units = params%cpp%temperature/C_E
       call save_to_hdf5(h5file_id,dset,units*cparams_ms%Te,attr)

       dset = TRIM(gname) // "/ne"
       attr = "Background electron density in m^-3"
       units = params%cpp%density
       call save_to_hdf5(h5file_id,dset,units*cparams_ms%ne,attr)

       dset = TRIM(gname) // "/nH"
       attr = "Background proton density in m^-3"
       units = params%cpp%density
       call save_to_hdf5(h5file_id,dset,units*cparams_ms%nH,attr)

       dset = TRIM(gname) // "/nef"
       attr = "Free electron density in m^-3"
       units = params%cpp%density
       call save_to_hdf5(h5file_id,dset,units*cparams_ms%nef,attr)

       dset = TRIM(gname) // "/neb"
       attr_array(1) = "Bound electron density per impurity in m^-3"
       units = params%cpp%density
       call save_1d_array_to_hdf5(h5file_id,dset,units*cparams_ms%neb, &
            attr_array)

       dset = TRIM(gname) // "/Zo"
       attr_array(1) = "Full nuclear charge of impurities"
       call save_1d_array_to_hdf5(h5file_id,dset,cparams_ms%Zo,attr_array)

       dset = TRIM(gname) // "/Zj"
       attr_array(1) = "Average charge state of impurities"
       call save_1d_array_to_hdf5(h5file_id,dset,cparams_ms%Zj,attr_array)

       dset = TRIM(gname) // "/nz"
       attr_array(1) = "Density of impurities in m^-3"
       units = params%cpp%density
       call save_1d_array_to_hdf5(h5file_id,dset,units*cparams_ms%nz,attr_array)

       dset = TRIM(gname) // "/IZj"
       attr_array(1) = " Ionization energy of impurities in eV"
       units = params%cpp%energy/C_E
       call save_1d_array_to_hdf5(h5file_id,dset,units*cparams_ms%IZj, &
            attr_array)

       dset = TRIM(gname) // "/rD"
       attr = "Debye length in m"
       units = params%cpp%length
       call save_to_hdf5(h5file_id,dset,units*cparams_ms%rD,attr)

       dset = TRIM(gname) // "/re"
       attr = "Classical electron radius in m"
       units = params%cpp%length
       call save_to_hdf5(h5file_id,dset,units*cparams_ms%re,attr)

       dset = TRIM(gname) // "/Ec"
       attr = "Critical electric field with impurities"
       units = params%cpp%Eo
       call save_to_hdf5(h5file_id,dset,units*cparams_ms%Ec_min,attr)

       DEALLOCATE(attr_array)

       call h5gclose_f(group_id, h5error)

       call h5fclose_f(h5file_id, h5error)
    end if
  end subroutine save_params_ms


  subroutine save_params_ss(params)
    TYPE(KORC_PARAMS), INTENT(IN) 				:: params
    CHARACTER(MAX_STRING_LENGTH) 				:: filename
    CHARACTER(MAX_STRING_LENGTH) 				:: gname
    CHARACTER(MAX_STRING_LENGTH) 				:: subgname
    CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
    CHARACTER(MAX_STRING_LENGTH) 				:: dset
    CHARACTER(MAX_STRING_LENGTH) 				:: attr
    INTEGER(HID_T) 						:: h5file_id
    INTEGER(HID_T) 						:: group_id
    INTEGER(HID_T) 						:: subgroup_id
    INTEGER 							:: h5error
    REAL(rp) 							:: units


    if (params%mpi_params%rank .EQ. 0) then
       filename = TRIM(params%path_to_outputs) // "simulation_parameters.h5"
       call h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, h5file_id, h5error)

       gname = "collisions_ss"
       call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

       ALLOCATE(attr_array(cparams_ms%num_impurity_species))

       dset = TRIM(gname) // "/collisions_model"
       call save_string_parameter(h5file_id,dset,(/params%collisions_model/))

       dset = TRIM(gname) // "/Te"
       attr = "Background electron temperature in eV"
       units = params%cpp%temperature/C_E
       call save_to_hdf5(h5file_id,dset,units*cparams_ss%Te,attr)

       dset = TRIM(gname) // "/Ti"
       attr = "Background ion temperature in eV"
       units = params%cpp%temperature/C_E
       call save_to_hdf5(h5file_id,dset,units*cparams_ss%Ti,attr)

       dset = TRIM(gname) // "/ne"
       attr = "Background electron density in m^-3"
       units = params%cpp%density
       call save_to_hdf5(h5file_id,dset,units*cparams_ss%ne,attr)

       dset = TRIM(gname) // "/Zeff"
       attr = "Effective nuclear charge of impurities"
       call save_to_hdf5(h5file_id,dset,cparams_ss%Zeff,attr)

       dset = TRIM(gname) // "/rD"
       attr = "Debye length in m"
       units = params%cpp%length
       call save_to_hdf5(h5file_id,dset,units*cparams_ss%rD,attr)

       dset = TRIM(gname) // "/re"
       attr = "Classical electron radius in m"
       units = params%cpp%length
       call save_to_hdf5(h5file_id,dset,units*cparams_ss%re,attr)

       dset = TRIM(gname) // "/Clogee"
       attr = "Coulomb logarithm"
       call save_to_hdf5(h5file_id,dset,cparams_ss%CoulombLogee,attr)

       dset = TRIM(gname) // "/Clogei"
       attr = "Coulomb logarithm"
       call save_to_hdf5(h5file_id,dset,cparams_ss%CoulombLogei,attr)

       dset = TRIM(gname) // "/VTe"
       attr = "Background electron temperature"
       units = params%cpp%velocity
       call save_to_hdf5(h5file_id,dset,units*cparams_ss%VTe,attr)

       dset = TRIM(gname) // "/delta"
       attr = "Delta parameter VTe/C"
       call save_to_hdf5(h5file_id,dset,cparams_ss%delta,attr)

       dset = TRIM(gname) // "/Gamma"
       attr = "Gamma coefficient"
       units = (params%cpp%mass**2*params%cpp%velocity**3)/params%cpp%time
       call save_to_hdf5(h5file_id,dset,units*cparams_ss%Gammac,attr)

       dset = TRIM(gname) // "/Tau_rad"
       attr = "Synchroton damping time in s"
       units = params%cpp%time
       call save_to_hdf5(h5file_id,dset,units*cparams_ss%taur,attr)

       dset = TRIM(gname) // "/Tau"
       attr = "Relativistic collisional time in s"
       units = params%cpp%time
       call save_to_hdf5(h5file_id,dset,units*cparams_ss%Tau,attr)

       dset = TRIM(gname) // "/Tauc"
       attr = "Thermal collisional time in s"
       units = params%cpp%time
       call save_to_hdf5(h5file_id,dset,units*cparams_ss%Tauc,attr)

       dset = TRIM(gname) // "/dTau"
       attr = "Subcycling time step in s"
       units = params%cpp%time
       call save_to_hdf5(h5file_id,dset,units*cparams_ss%dTau* &
            cparams_ss%Tau,attr)

       dset = TRIM(gname) // "/subcycling_iterations"
       attr = "KORC iterations per collision"
       call save_to_hdf5(h5file_id,dset,cparams_ss%subcycling_iterations,attr)

       dset = TRIM(gname) // "/Ec"
       attr = "Critical electric field"
       units = params%cpp%Eo
       call save_to_hdf5(h5file_id,dset,units*cparams_ss%Ec,attr)

       dset = TRIM(gname) // "/ED"
       attr = "Dreicer electric field"
       units = params%cpp%Eo
       call save_to_hdf5(h5file_id,dset,units*cparams_ss%ED,attr)

       call h5gclose_f(group_id, h5error)

       call h5fclose_f(h5file_id, h5error)
    end if
  end subroutine save_params_ss


  subroutine save_collision_params(params)
    TYPE(KORC_PARAMS), INTENT(IN) :: params

    if (.NOT.(params%restart)) then

       if (params%collisions) then
          SELECT CASE (TRIM(params%collisions_model))
          CASE (MODEL1)
             call save_params_ss(params)

             SELECT CASE(TRIM(params%bound_electron_model))
             CASE ('NO_BOUND')
                call save_params_ms(params)
             CASE('HESSLOW')
                call save_params_ms(params)
             CASE('ROSENBLUTH')
                call save_params_ms(params)
             CASE DEFAULT
                write(output_unit_write,'("Default case")')
             END SELECT

          CASE (MODEL2)
             call save_params_ms(params)
          CASE DEFAULT
             write(output_unit_write,'("Default case")')
          END SELECT
       end if

    end if
  end subroutine save_collision_params


  subroutine deallocate_params_ms()
    if (ALLOCATED(cparams_ms%Zj)) DEALLOCATE(cparams_ms%Zj)
    if (ALLOCATED(cparams_ms%Zo)) DEALLOCATE(cparams_ms%Zo)
    if (ALLOCATED(cparams_ms%nz)) DEALLOCATE(cparams_ms%nz)
    if (ALLOCATED(cparams_ms%neb)) DEALLOCATE(cparams_ms%neb)
    if (ALLOCATED(cparams_ms%IZj)) DEALLOCATE(cparams_ms%IZj)
    if (ALLOCATED(cparams_ms%Zj)) DEALLOCATE(cparams_ms%Ee_IZj)
  end subroutine deallocate_params_ms


  subroutine deallocate_collisions_params(params)
    TYPE(KORC_PARAMS), INTENT(IN) :: params

    if (params%collisions) then
       SELECT CASE (TRIM(params%collisions_model))
       CASE (MODEL1)
          !				write(output_unit_write,'("Something to be done")')

          SELECT CASE(TRIM(params%bound_electron_model))
          CASE ('NO_BOUND')
             call deallocate_params_ms()
          CASE('HESSLOW')
             call deallocate_params_ms()
          CASE('ROSENBLUTH')
             call deallocate_params_ms()
          CASE DEFAULT
             write(output_unit_write,'("Default case")')
          END SELECT

       CASE (MODEL2)
          call deallocate_params_ms()
       CASE DEFAULT
          write(output_unit_write,'("Default case")')
       END SELECT
    end if
  end subroutine deallocate_collisions_params

end module korc_collisions
