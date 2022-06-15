MODULE korc_avalanche
  USE korc_types
  USE korc_constants
  USE korc_HDF5
  USE korc_hpc
  USE korc_fields
  USE korc_profiles
  USE korc_coords
  USE korc_input
  
  IMPLICIT NONE

  TYPE, PRIVATE :: AVALANCHE_PDF_PARAMS
     REAL(rp) :: max_pitch_angle
     !! Maximum pitch angle of sampled PDF in degrees
     REAL(rp) :: min_pitch_angle
     !! Minimum pitch angle of sampled PDF in degrees
     REAL(rp) :: min_energy
     !! Minimum energy of sampled PDF in MeV
     REAL(rp) :: max_energy
     !! Maximum energy of sampled PDF in MeV
     REAL(rp) :: min_p
     !! Minimum momentum of sampled PDF
     REAL(rp) :: max_p
     !! Maximum momentum of sampled PDF
     REAL(rp) :: ne
     !! Background electron density in m^-3
     REAL(rp) :: Zeff
     !! Effective atomic number of ions
     REAL(rp) :: Ec
     !! Critical electric field in V/m
     REAL(rp) :: Epar
     !! Parallel electric field in V/m
     REAL(rp) :: Ebar
     !! Epar/Ec
     REAL(rp) :: Te
     !! Background electron temperature in eV
     REAL(rp) :: lD
     !! Debye length
     REAL(rp) :: bmin
     !! Maximum approach radius
     REAL(rp) :: CoulombLog
     !! Coulomb Logarithm
     REAL(rp) :: Tau
     !! Collisional time

     REAL(rp) :: dth
     !! Variance of sampling normal variate for pitch angle
     REAL(rp) :: dp
     !! Variance of sampling normal variate for momentum
     REAL(rp) :: dR
     !! Variance of sampling normal variate for R location
     REAL(rp) :: dZ
     !! Variance of sampling normal variate for Z location
     
     REAL(rp) :: fo
     REAL(rp) :: alpha
     REAL(rp) :: cz
     REAL(rp) :: C1
     REAL(rp) :: C2
  END TYPE AVALANCHE_PDF_PARAMS

  TYPE(AVALANCHE_PDF_PARAMS), PRIVATE :: aval_params
  REAL(rp), PRIVATE, PARAMETER :: minmax_buffer_size = 10.0_rp

  PUBLIC :: get_avalanche_distribution,&
       get_Avalanche_4D
  PRIVATE :: initialize_avalanche_params,&
       save_avalanche_params,&
       fRE,&
       sample_distribution,&
       indicator,&
       PSI_ROT,&
       random_norm,&
       update_avalanche_params,&
       Avalanche_4D

CONTAINS

  SUBROUTINE get_avalanche_distribution(params,g,eta,go,etao)
    TYPE(KORC_PARAMS), INTENT(IN) :: params
    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: g
    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: eta
    REAL(rp), INTENT(OUT) :: go
    REAL(rp), INTENT(OUT) :: etao

    call initialize_avalanche_params(params)

    call save_avalanche_params(params)

    call sample_distribution(params,g,eta,go,etao)
  END SUBROUTINE get_avalanche_distribution

  SUBROUTINE get_Avalanche_4D(params,spp,P,F)
    TYPE(KORC_PARAMS), INTENT(IN) :: params
    TYPE(SPECIES), INTENT(INOUT)  :: spp
    TYPE(PROFILES), INTENT(IN)    :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    !! An instance of the KORC derived type FIELDS.

    call initialize_avalanche_params(params)

    call save_avalanche_params(params)

    call Avalanche_4D(params,spp,P,F)
  END SUBROUTINE get_Avalanche_4D


  SUBROUTINE initialize_avalanche_params(params)
    TYPE(KORC_PARAMS), INTENT(IN) :: params
    !REAL(rp) :: max_pitch_angle
    !REAL(rp) :: min_pitch_angle
    !REAL(rp) :: max_energy
    !REAL(rp) :: min_energy
    !REAL(rp) :: ne
    !REAL(rp) :: Zeff
    !REAL(rp) :: Epar
    !REAL(rp) :: Te
    !REAL(rp) :: dth
    !REAL(rp) :: dp
    !REAL(rp) :: dR
    !REAL(rp) :: dZ
    !NAMELIST /AvalancheGenerationPDF/ max_pitch_angle,min_pitch_angle, &
    !     max_energy,min_energy,ne,Zeff,Epar,Te,dth,dp,dR,dZ

    !open(unit=default_unit_open,file=TRIM(params%path_to_inputs), &
    !     status='OLD',form='formatted')
    !read(default_unit_open,nml=AvalancheGenerationPDF)
    !close(default_unit_open)

    aval_params%dth = dth_aval 
    aval_params%dp  = dp_aval 
    aval_params%dR  = dR_aval/params%cpp%length
    aval_params%dZ  = dZ_aval/params%cpp%length
    
    aval_params%max_pitch_angle = max_pitch_angle_aval
    aval_params%min_pitch_angle = min_pitch_angle_aval
    aval_params%max_energy = max_energy_aval*C_E ! In Joules
    aval_params%ne = ne_aval
    aval_params%Zeff = Zeff_aval
    aval_params%Te = Te_aval*C_E ! In Joules

    aval_params%lD = SQRT(C_E0*aval_params%Te/(aval_params%ne*C_E**2))
    aval_params%bmin = aval_params%Zeff/(12.0_rp*C_PI*aval_params%ne* &
         aval_params%lD**2)
    aval_params%CoulombLog = LOG(aval_params%lD/aval_params%bmin)
    aval_params%Tau = 1.0_rp/(4.0_rp*C_PI*C_C*C_RE**2*aval_params%ne* &
         aval_params%CoulombLog)

    aval_params%Ec = C_ME*C_C/(C_E*aval_params%Tau)
    aval_params%Epar = Epar_aval
    aval_params%Ebar = aval_params%Epar/aval_params%Ec

    if (min_energy_aval .EQ. 0.0_rp) then
       aval_params%max_p = SQRT((aval_params%max_energy/(C_ME*C_C**2))**2 &
            - 1.0_rp)
       ! In units of mec^2
       aval_params%min_p = SQRT(aval_params%Ebar - 1.0_rp)
       ! In units of mec^2

       aval_params%min_energy = SQRT(1.0_rp + aval_params%min_p**2)*C_ME*C_C**2
    else
       aval_params%min_energy = min_energy_aval*C_E
       ! In Joules

       aval_params%max_p = SQRT((aval_params%max_energy/(C_ME*C_C**2))**2 - 1.0_rp)
       ! In units of mec^2
       aval_params%min_p = SQRT((aval_params%min_energy/(C_ME*C_C**2))**2 - 1.0_rp)
       ! In units of mec^2
    end if

    aval_params%alpha = (aval_params%Ebar - 1.0_rp)/(1.0_rp + aval_params%Zeff)
    aval_params%cz = SQRT(3.0_rp*(aval_params%Zeff + 5.0_rp)/C_PI)* &
         aval_params%CoulombLog
    aval_params%fo = aval_params%alpha/aval_params%cz
    aval_params%C1 = 0.5_rp*aval_params%alpha
    aval_params%C2 = 1.0_rp/aval_params%cz - aval_params%C1
  END SUBROUTINE initialize_avalanche_params


  FUNCTION deg2rad(x)
    REAL(rp), INTENT(IN) :: x
    REAL(rp) :: deg2rad

    deg2rad = C_PI*x/180.0_rp
  END FUNCTION deg2rad


  FUNCTION fRE(x,p)
    REAL(rp), INTENT(IN) :: x ! x = cos(pitch)
    REAL(rp), INTENT(IN) :: p ! momentum
    REAL(rp) :: fRE

    fRE = aval_params%fo*p*EXP(-p*(aval_params%C2*x + aval_params%C1/x))/x
  END FUNCTION fRE

  FUNCTION log10fRE(x,p)
    REAL(rp), INTENT(IN) :: x ! x = cos(pitch)
    REAL(rp), INTENT(IN) :: p ! momentum
    REAL(rp) :: log10fRE

    log10fRE = LOG(fRE(x,p))
  END FUNCTION log10fRE


  SUBROUTINE sample_distribution(params,g,eta,go,etao)
    TYPE(KORC_PARAMS), INTENT(IN) :: params
    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: g
    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: eta
    REAL(rp), INTENT(OUT) :: go
    REAL(rp), INTENT(OUT) :: etao
    REAL(rp), DIMENSION(:), ALLOCATABLE :: p
    REAL(rp) :: chi, chi_test
    REAL(rp) :: p_buffer, p_test
    REAL(rp) :: eta_buffer, eta_test
    REAL(rp) :: ratio, rand_unif
    REAL(rp), DIMENSION(:), ALLOCATABLE :: p_samples
    REAL(rp), DIMENSION(:), ALLOCATABLE :: eta_samples
    REAL(rp), DIMENSION(:), ALLOCATABLE :: p_tmp
    REAL(rp), DIMENSION(:), ALLOCATABLE :: eta_tmp
    REAL(rp) :: minmax,min_p, max_p, min_pitch_angle, max_pitch_angle
    REAL(4), DIMENSION(2) :: tarray
    REAL(4) :: time_elapsed
    REAL(rp) :: deta
    REAL(rp) :: dp
    LOGICAL :: lp,leta
    INTEGER :: num_accepted
    INTEGER :: ii,jj,ppp,nsamples
    INTEGER :: mpierr


    ppp = SIZE(g)
    nsamples = ppp*params%mpi_params%nmpi
    ALLOCATE(p(ppp))

    deta = (aval_params%max_pitch_angle - aval_params%min_pitch_angle)/100.0_rp
    dp = (aval_params%max_p - aval_params%min_p)/100.0_rp

    do jj=1_idef,INT(minmax_buffer_size,idef)
       minmax = aval_params%min_p - REAL(jj,rp)*dp
       if (minmax.GT.0.0_rp) then
          min_p = minmax
       end if
    end do

    max_p = aval_params%max_p + minmax_buffer_size*dp

    if (aval_params%min_pitch_angle.GE.korc_zero) then
       do jj=1_idef,INT(minmax_buffer_size,idef)
          minmax = aval_params%min_pitch_angle -  REAL(jj,rp)*deta
          if (minmax.GT.0.0_rp) then
             min_pitch_angle = minmax
          end if
       end do
    else
       min_pitch_angle = aval_params%min_pitch_angle
    end if

    do jj=1_idef,INT(minmax_buffer_size,idef)
       minmax = aval_params%max_pitch_angle + REAL(jj,rp)*deta
       if (minmax.LE.90.0_rp) then
          max_pitch_angle = minmax
       else
          max_pitch_angle = aval_params%max_pitch_angle
          EXIT
       end if
    end do

    if (params%mpi_params%rank.EQ.0_idef) then
       ALLOCATE(p_samples(nsamples))
       ! Number of samples to distribute among all MPI processes
       ALLOCATE(eta_samples(nsamples))
       ! Number of samples to distribute among all MPI processes
       ALLOCATE(p_tmp(nsamples))
       ! Number of samples to distribute among all MPI processes
       ALLOCATE(eta_tmp(nsamples))
       ! Number of samples to distribute among all MPI processes

       ! Transient !
       call RANDOM_SEED()

       call RANDOM_NUMBER(rand_unif)
       eta_buffer = aval_params%min_pitch_angle + (aval_params%max_pitch_angle  &
            - aval_params%min_pitch_angle)*rand_unif

       call RANDOM_NUMBER(rand_unif)
       p_buffer = aval_params%min_p + (aval_params%max_p - aval_params%min_p)* &
            rand_unif

       ii=2_idef
       do while (ii .LE. 1000_idef)
          eta_test = eta_buffer + random_norm(0.0_rp,deta)
          do while ((ABS(eta_test) .GT. aval_params%max_pitch_angle).OR. &
               (ABS(eta_test) .LT. aval_params%min_pitch_angle))
             eta_test = eta_buffer + random_norm(0.0_rp,deta)
          end do
          chi_test = COS(deg2rad(eta_test))
          chi = COS(deg2rad(eta_buffer))

          p_test = p_buffer + random_norm(0.0_rp,dp)
          do while ((p_test.LT.aval_params%min_p).OR.(p_test.GT.aval_params%max_p))
             p_test = p_buffer + random_norm(0.0_rp,dp)
          end do

          ratio = fRE(chi_test,p_test)/fRE(chi,p_buffer)

          if (ratio .GE. 1.0_rp) then
             p_buffer = p_test
             eta_buffer = eta_test
             ii = ii + 1_idef
          else 
             call RANDOM_NUMBER(rand_unif)
             if (rand_unif .LT. ratio) then
                p_buffer = p_test
                eta_buffer = eta_test
                ii = ii + 1_idef
             end if
          end if
       end do
       ! Transient !


       call RANDOM_SEED()
       call RANDOM_NUMBER(rand_unif)

       eta_tmp(1) = eta_buffer
       p_tmp(1) = p_buffer

       num_accepted = 0_idef
       do while(num_accepted.LT.nsamples)
          ii=2_idef
          do while (ii .LE. nsamples)
             eta_test = eta_tmp(ii-1) + random_norm(0.0_rp,deta)
             do while ((ABS(eta_test) .GT. max_pitch_angle).OR.(ABS(eta_test) &
                  .LT. min_pitch_angle))
                eta_test = eta_tmp(ii-1) + random_norm(0.0_rp,deta)
             end do
             chi_test = COS(deg2rad(eta_test))
             chi = COS(deg2rad(eta_tmp(ii-1)))

             p_test = p_tmp(ii-1) + random_norm(0.0_rp,dp)
             do while ((p_test.LT.min_p).OR.(p_test.GT.max_p))
                p_test = p_tmp(ii-1) + random_norm(0.0_rp,dp)
             end do

             ratio = fRE(chi_test,p_test)/fRE(chi,p_tmp(ii-1))

             if (ratio .GE. 1.0_rp) then
                p_tmp(ii) = p_test
                eta_tmp(ii) = eta_test
                ii = ii + 1_idef
             else 
                call RANDOM_NUMBER(rand_unif)
                if (rand_unif .LT. ratio) then
                   p_tmp(ii) = p_test
                   eta_tmp(ii) = eta_test
                   ii = ii + 1_idef
                end if
             end if
          end do

          eta_tmp = ABS(eta_tmp)

          ii = 1_idef
          do while ( (ii.LT.nsamples).AND.(num_accepted.LT.nsamples) )
             lp = (p_tmp(ii).LE.aval_params%max_p).AND.(p_tmp(ii).GE. &
                  aval_params%min_p)
             leta = (eta_tmp(ii).LE.aval_params%max_pitch_angle).AND. &
                  (eta_tmp(ii).GE.aval_params%min_pitch_angle)
             if (lp.AND.leta) then
                num_accepted = num_accepted + 1_idef
                p_samples(num_accepted) = p_tmp(ii)
                eta_samples(num_accepted) = eta_tmp(ii)
             end if
             ii = ii + 1_idef
          end do

          eta_tmp(1) = eta_tmp(ii)
          p_tmp(1) = p_tmp(ii)
       end do

       go = SUM(SQRT(1.0_rp + p_samples**2))/nsamples
       etao = SUM(eta_samples)/nsamples

    end if

    CALL MPI_SCATTER(p_samples,ppp,MPI_REAL8,p,ppp,MPI_REAL8,0, &
         MPI_COMM_WORLD,mpierr)

    CALL MPI_SCATTER(eta_samples,ppp,MPI_REAL8,eta,ppp,MPI_REAL8, &
         0,MPI_COMM_WORLD,mpierr)

    CALL MPI_BCAST(go,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

    CALL MPI_BCAST(etao,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

    g = SQRT(1.0_rp + p**2)

    DEALLOCATE(p)
    if (params%mpi_params%rank.EQ.0_idef) then
       DEALLOCATE(p_samples)
       DEALLOCATE(eta_samples)
    end if
  END SUBROUTINE sample_distribution

  FUNCTION PSI_ROT(R,R0,sigR,Z,Z0,sigZ,theta)
    !! Calculates value of argument of 2D Gaussian spatial distribution with
    !! with counter-clockwise rotation.
    REAL(rp), INTENT(IN) 	:: R
    !! R-coordinate of MH sampled location
    REAL(rp), INTENT(IN) 	:: R0
    !! R-coordinate of center of 2D Gaussian
    REAL(rp), INTENT(IN) 	:: sigR
    !! Variance of first dimension of 2D Gaussian
    REAL(rp), INTENT(IN) 	:: Z
    !! Z-coordinate of MH sampled location
    REAL(rp), INTENT(IN) 	:: Z0
    !! Z-coordinate of center of 2D Gaussian
    REAL(rp), INTENT(IN) 	:: sigZ
    !! Variance of second dimension of 2D Gaussian
    REAL(rp), INTENT(IN) 	:: theta
    !! Angle of counter-clockwise rotation (in radians), of 2D Gaussian
    !! distribution relative to R,Z
    REAL(rp) 		:: PSI_ROT
    !! Argument of exponential comprising 2D Gaussian distribution 

    PSI_ROT=(R-R0)**2*((cos(theta))**2/(2*sigR**2)+(sin(theta))**2/ &
         (2*sigZ**2))+2*(R-R0)*(Z-Z0)*cos(theta)*sin(theta)*(1/ &
         (2*sigR**2)-1/(2*sigZ**2))+(Z-Z0)**2*((sin(theta))**2/ &
         (2*sigR**2)+(cos(theta))**2/(2*sigZ**2))

  END FUNCTION PSI_ROT


  FUNCTION indicator(psi,psi_max)
    !! Compares argument psi to chosen psi_max, returning step function.
    REAL(rp), INTENT(IN)  :: psi
    REAL(rp), INTENT(IN)  :: psi_max
    REAL(rp)              :: indicator

    IF (psi.LT.psi_max) THEN
       indicator=1
    ELSE
       indicator=0
    END IF

  END FUNCTION indicator


  FUNCTION random_norm(mean,sigma)
    REAL(rp), INTENT(IN) :: mean
    REAL(rp), INTENT(IN) :: sigma
    REAL(rp)             :: random_norm
    REAL(rp)             :: rand1, rand2

    call RANDOM_NUMBER(rand1)
    call RANDOM_NUMBER(rand2)

    random_norm = mean+sigma*SQRT(-2.0_rp*LOG(rand1))*COS(2.0_rp*C_PI*rand2);
  END FUNCTION random_norm


SUBROUTINE update_avalanche_params(params,prtcls)
  !! Updates the avalanche parameters aval_params% at each step
  !! in the MCMC after the profiles are interpolated at the sampled
  !! R,Z location.
  TYPE(KORC_PARAMS), INTENT(IN) 	:: params
  !! Core KORC simulation parameters.
  TYPE(PARTICLES), INTENT(IN) :: prtcls
  !! An instance of PARTICLES containing the variables of a
  !! given species. Call to this subroutine generally passes spp%vars.

  aval_params%ne = prtcls%ne(1)*params%cpp%density
  aval_params%Zeff = prtcls%Zeff(1)
  aval_params%Te = prtcls%Te(1)*params%cpp%temperature ! In Joules
  
  aval_params%lD = SQRT(C_E0*aval_params%Te/(aval_params%ne*C_E**2))
  aval_params%bmin = aval_params%Zeff/(12.0_rp*C_PI*aval_params%ne* &
       aval_params%lD**2)
  aval_params%CoulombLog = LOG(aval_params%lD/aval_params%bmin)
  aval_params%Tau = 1.0_rp/(4.0_rp*C_PI*C_C*C_RE**2*aval_params%ne* &
       aval_params%CoulombLog)

  aval_params%Ec = C_ME*C_C/(C_E*aval_params%Tau)
  aval_params%Ebar = aval_params%Epar/aval_params%Ec

  aval_params%alpha = (aval_params%Ebar - 1.0_rp)/(1.0_rp +&
       aval_params%Zeff)
  aval_params%cz = SQRT(3.0_rp*(aval_params%Zeff + 5.0_rp)/C_PI)* &
       aval_params%CoulombLog
  aval_params%fo = aval_params%alpha/aval_params%cz
  aval_params%C1 = 0.5_rp*aval_params%alpha
  aval_params%C2 = 1.0_rp/aval_params%cz - aval_params%C1

  
END SUBROUTINE update_avalanche_params

subroutine Avalanche_4D(params,spp,P,F)
  !! @note Subroutine that generates a 2D Gaussian distribution in an 
  !! elliptic torus as the initial spatial condition of a given particle 
  !! species in the simulation. @endnote
  TYPE(KORC_PARAMS), INTENT(IN) 	:: params
  !! Core KORC simulation parameters.
  TYPE(SPECIES), INTENT(INOUT) 		:: spp
  !! An instance of the derived type SPECIES containing all the parameters
  !! and simulation variables of the different species in the simulation.
  TYPE(PROFILES), INTENT(IN)            :: P
  !! An instance of the KORC derived type PROFILES.
  TYPE(FIELDS), INTENT(IN)                                   :: F
  !! An instance of the KORC derived type FIELDS.
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: R_samples
  !! Major radial location of all samples
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: PHI_samples
  !! Azimuithal angle of all samples
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: Z_samples
  !! Vertical location of all samples
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: P_samples
  !! Magnitude of momentum of all samples
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: T_samples
  !! Pitch angle of all samples
  REAL(rp), DIMENSION(:), ALLOCATABLE   :: mom
  !! Magnitude of momentum of samples on individual mpi process
  REAL(rp) 				:: psi_max_buff
  !! Value of buffer above desired maximum argument of 2D Gaussian spatial
  !! profile
  REAL(rp) 				:: minmax
  !! Temporary variable used for setting buffers
  REAL(rp) 				:: max_p
  !! Maximum domain for momentum sampling including buffer
  REAL(rp) 				:: min_p
  !! Minimum domain for momentum sampling including buffer
  REAL(rp) 				:: max_pitch_angle
  !! Maximum domain for pitch angle sampling including buffer
  REAL(rp) 				:: min_pitch_angle
  !! Minimum domain for pitch angle sampling including buffer
  REAL(rp) 				:: theta_rad
  !! Angle of rotation of 2D Gaussian spatial distribution in radians
  REAL(rp) 				:: R_buffer
  !! Previous sample of R location
  REAL(rp) 				:: Z_buffer
  !! Previous sample of Z location
  REAL(rp) 				:: P_buffer
  !! Previous sample of momentum
  REAL(rp) 				:: T_buffer
  !! Previous sample of pitch angle
  REAL(rp) 				:: R_test
  !! Present sample of R location
  REAL(rp) 				:: Z_test
  !! Present sample of Z location
  REAL(rp) 				:: P_test
  !! Present sample of momentum
  REAL(rp) 				:: T_test
  !! Present sample of pitch angle
  REAL(rp) 				:: psi0
  !! Previous value of 2D Gaussian argument based on R_buffer, Z_buffer
  REAL(rp) 				:: psi1
  !! Present value of 2D Gaussian argument based on R_test, Z_test
  REAL(rp) 				:: fRE0
  !! Evaluation of Avalanche distribution with previous sample
  REAL(rp) 				:: fRE1
  !! Evaluation of Avalanche distribution with present sample
  REAL(rp) 				:: rand_unif
  !! Uniform random variable [0,1]
  REAL(rp) 				:: ratio
  !! MH selection criteria
  INTEGER				:: nsamples
  !! Total number of samples to be distributed over all mpi processes
  INTEGER 				:: ii
  !! Sample iterator.
  INTEGER 				:: mpierr
  !! mpi error indicator
  
  ALLOCATE(mom(spp%ppp))
  
  nsamples = spp%ppp*params%mpi_params%nmpi

  psi_max_buff = spp%psi_max*1.1_rp

  theta_rad=C_PI*spp%theta_gauss/180.0_rp

  ! buffer at minimum p boundary
  do ii=1_idef,INT(minmax_buffer_size,idef)
     minmax = aval_params%min_p - REAL(ii,rp)*(aval_params%max_p- &
       aval_params%min_p)/100_rp
     if (minmax.GT.0.0_rp) then
        min_p = minmax
     end if
  end do

  ! buffer at maximum p boundary
  max_p = aval_params%max_p + minmax_buffer_size*(aval_params%max_p- &
       aval_params%min_p)/100_rp

  ! buffer at minimum pitch angle boundary  
  if (aval_params%min_pitch_angle.GE.korc_zero) then
     do ii=1_idef,INT(minmax_buffer_size,idef)
        minmax = aval_params%min_pitch_angle - REAL(ii,rp)* &
             (aval_params%max_pitch_angle-aval_params%min_pitch_angle)/100_rp
        if (minmax.GT.0.0_rp) then
           min_pitch_angle = minmax
        end if
     end do
  else
     min_pitch_angle = aval_params%min_pitch_angle
  end if

  ! buffer at maximum pitch angle boundary
  do ii=1_idef,INT(minmax_buffer_size,idef)
     minmax = aval_params%max_pitch_angle + REAL(ii,rp)* &
          (aval_params%max_pitch_angle-aval_params%min_pitch_angle)/100_rp
     if (minmax.LE.180.0_rp) then
        max_pitch_angle = minmax
     end if
  end do

  
  if (params%mpi_params%rank.EQ.0_idef) then
     ALLOCATE(R_samples(nsamples))
     ! Number of samples to distribute among all MPI processes
     ALLOCATE(PHI_samples(nsamples))
     ! Number of samples to distribute among all MPI processes
     ALLOCATE(Z_samples(nsamples))
     ! Number of samples to distribute among all MPI processes
     ALLOCATE(P_samples(nsamples))
     ! Number of samples to distribute among all MPI processes
     ALLOCATE(T_samples(nsamples))
     ! Number of samples to distribute among all MPI processes
     
     ! Transient !

     R_buffer = spp%Ro
     Z_buffer = spp%Zo

     
     call RANDOM_NUMBER(rand_unif)
     T_buffer = min_pitch_angle + (max_pitch_angle  &
          - min_pitch_angle)*rand_unif

     call RANDOM_NUMBER(rand_unif)
     P_buffer = min_p + (max_p - min_p)*rand_unif

!     write(output_unit_write,'("length norm: ",E17.10)') params%cpp%length
     
     ii=1_idef
     do while (ii .LE. 1000_idef)

!        write(output_unit_write,'("burn:",I15)') ii
        
        R_test = R_buffer + random_norm(0.0_rp,aval_params%dR)
        Z_test = Z_buffer + random_norm(0.0_rp,aval_params%dZ)
        P_test = P_buffer + random_norm(0.0_rp,aval_params%dp)
        T_test = T_buffer + random_norm(0.0_rp,aval_params%dth)


        ! Test that pitch angle and momentum are within chosen boundary
        do while ((T_test .GT. aval_params%max_pitch_angle).OR. &
             (T_test .LT. aval_params%min_pitch_angle))
           T_test = T_buffer + random_norm(0.0_rp,aval_params%dth)
        end do

        do while ((P_test.LT.aval_params%min_p).OR. &
             (P_test.GT.aval_params%max_p))
           P_test = P_buffer + random_norm(0.0_rp,aval_params%dp)
        end do

        ! initialize 2D gaussian argument and distribution function, or
        ! copy from previous sample
        if (ii==1) then
           psi0=PSI_ROT(R_buffer,spp%Ro,spp%sigmaR,Z_buffer,spp%Zo, &
                spp%sigmaZ,theta_rad)

           ! Repetedly put location of sample into first entry of spp%vars%Y,
           ! beacause spp%vas%X isn't filled until the end of this subrouting.
           ! This allows us to use the machinery for calling PSPLINES
           ! interpolation routines here also.
           spp%vars%Y(1,1)=R_buffer
           spp%vars%Y(1,2)=0
           spp%vars%Y(1,3)=Z_buffer

           call get_profiles(params,spp%vars,P,F)          

!           write(output_unit_write,'("ne",E17.10)') spp%vars%ne(1)
!           write(output_unit_write,'("Te",E17.10)') spp%vars%Te(1)
!           write(output_unit_write,'("Zeff",E17.10)') spp%vars%Zeff(1)
           
           ! Update avalanche parameters with interpolated fields to be used
           ! in call to avalanche distribution function
           call update_avalanche_params(params,spp%vars)
           
           fRE0=fRE(cos(deg2rad(T_buffer)),P_buffer)          
        else
           psi0=psi1
           fRE0=fRE1
        end if
        
        psi1=PSI_ROT(R_test,spp%Ro,spp%sigmaR,Z_test,spp%Zo, &
             spp%sigmaZ,theta_rad)       
        
        spp%vars%Y(1,1)=R_test
        spp%vars%Y(1,2)=0
        spp%vars%Y(1,3)=Z_test
        
        call get_profiles(params,spp%vars,P,F)

        call update_avalanche_params(params,spp%vars)
        
        fRE1=fRE(COS(deg2rad(T_test)),P_test)

!        write(output_unit_write,'("psi0: ",E17.10)') psi0
!        write(output_unit_write,'("psi1: ",E17.10)') psi1
        
        ! Calculate acceptance ratio for MH algorithm. fRE function
        ! incorporates p^2 factor of spherical coordinate Jacobian
        ! for velocity phase space, factors of sin(pitch angle) for velocity
        ! phase space and cylindrical coordinate Jacobian R for spatial
        ! phase space incorporated here.
        ratio = indicator(psi1,spp%psi_max)*R_test*EXP(-psi1)* &
             sin(deg2rad(T_test))*fRE1/(R_buffer*EXP(-psi0)* &
             sin(deg2rad(T_buffer))*fRE0)
        
        if (ratio .GE. 1.0_rp) then
           R_buffer = R_test
           Z_buffer = Z_test
           P_buffer = P_test
           T_buffer = T_test
           ii = ii + 1_idef
        else
           call RANDOM_NUMBER(rand_unif)
           if (rand_unif .LT. ratio) then
              R_buffer = R_test
              Z_buffer = Z_test
              P_buffer = P_test
              T_buffer = T_test
              ii = ii + 1_idef
           end if
        end if
     end do
     ! Transient !

     ii=1_idef
     do while (ii .LE. nsamples)

!        write(output_unit_write,'("sample:",I15)') ii
        
        if (modulo(ii,10000).eq.0) then
           write(output_unit_write,'("Sample: ",I10)') ii
        end if
        
        R_test = R_buffer + random_norm(0.0_rp,aval_params%dR)
        Z_test = Z_buffer + random_norm(0.0_rp,aval_params%dZ)
        P_test = P_buffer + random_norm(0.0_rp,aval_params%dp)
        T_test = T_buffer + random_norm(0.0_rp,aval_params%dth)

        ! Selection boundary is set with buffer region
        do while ((T_test .GT. max_pitch_angle).OR. &
             (T_test .LT. min_pitch_angle))
           if (T_test.lt.0) then
              T_test=abs(T_test)
              exit
           end if
           T_test = T_buffer + random_norm(0.0_rp,aval_params%dth)
        end do

        do while ((P_test.LT.min_p).OR.(P_test.GT.max_p))
           if (P_test.lt.0) then
              P_test=abs(P_test)
              exit
           end if
           P_test = P_buffer + random_norm(0.0_rp,aval_params%dp)
        end do

        psi0=psi1
        fRE0=fRE1
        
        psi1=PSI_ROT(R_test,spp%Ro,spp%sigmaR,Z_test,spp%Zo, &
             spp%sigmaZ,theta_rad)

        spp%vars%Y(1,1)=R_test
        spp%vars%Y(1,2)=0
        spp%vars%Y(1,3)=Z_test

        call get_profiles(params,spp%vars,P,F)

        call update_avalanche_params(params,spp%vars)
        
        fRE1=fRE(COS(deg2rad(T_test)),P_test)
        
        ratio = indicator(psi1,psi_max_buff)*R_test*EXP(-psi1)* &
             sin(deg2rad(T_test))*fRE1/(R_buffer*EXP(-psi0)* &
             sin(deg2rad(T_buffer))*fRE0)

        if (ratio .GE. 1.0_rp) then
           R_buffer = R_test
           Z_buffer = Z_test
           P_buffer = P_test
           T_buffer = T_test
        else
           call RANDOM_NUMBER(rand_unif)
           if (rand_unif .LT. ratio) then
              R_buffer = R_test
              Z_buffer = Z_test
              P_buffer = P_test
              T_buffer = T_test
           end if
        end if
        
        ! Only accept sample if it is within desired boundary, but
        ! add to MC above if within buffer. This helps make the boundary
        ! more defined.
        IF ((INT(indicator(psi1,spp%psi_max)).EQ.1).AND. &
             (P_buffer.LE.aval_params%max_p).AND. &
             (P_buffer.GE.aval_params%min_p).AND. &
             (T_buffer.LE.aval_params%max_pitch_angle).AND. &
             (T_buffer.GE.aval_params%min_pitch_angle)) THEN
           R_samples(ii) = R_buffer
           Z_samples(ii) = Z_buffer
           P_samples(ii) = P_buffer
           T_samples(ii) = T_buffer
           ! Sample phi location uniformly
           call RANDOM_NUMBER(rand_unif)
           PHI_samples(ii) = 2.0_rp*C_PI*rand_unif
           ii = ii + 1_idef 
        END IF

        
     end do

!  if (minval(R_samples(:)).lt.1._rp/params%cpp%length) stop 'error with sample'
!  write(output_unit_write,'("R_sample: ",E17.10)') R_samples(:)*params%cpp%length
  
  end if

  CALL MPI_SCATTER(R_samples*cos(PHI_samples),spp%ppp,MPI_REAL8, &
       spp%vars%X(:,1),spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  CALL MPI_SCATTER(R_samples*sin(PHI_samples),spp%ppp,MPI_REAL8, &
       spp%vars%X(:,2),spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  CALL MPI_SCATTER(Z_samples,spp%ppp,MPI_REAL8, &
       spp%vars%X(:,3),spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  CALL MPI_SCATTER(P_samples,spp%ppp,MPI_REAL8, &
       mom,spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  CALL MPI_SCATTER(T_samples,spp%ppp,MPI_REAL8, &
       spp%vars%eta,spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  
  
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

!  write(output_unit_write,'("X_X: ",E17.10)') spp%vars%X(:,1)*params%cpp%length
  
  ! gamma is kept for each particle, not the momentum
  spp%vars%g = SQRT(1.0_rp + mom**2)

  if (params%orbit_model(1:2).eq.'GC') call cart_to_cyl(spp%vars%X,spp%vars%Y)

!  write(output_unit_write,'("Y_R: ",E17.10)') spp%vars%Y(:,1)*params%cpp%length
  
!  if (minval(spp%vars%Y(:,1)).lt.1._rp/params%cpp%length) stop 'error with avalanche'
  
  DEALLOCATE(mom)
  if (params%mpi_params%rank.EQ.0_idef) then
     DEALLOCATE(R_samples)
     DEALLOCATE(Z_samples)
     DEALLOCATE(PHI_samples)
     DEALLOCATE(P_samples)
     DEALLOCATE(T_samples)
  end if
  
  
end subroutine Avalanche_4D

  
SUBROUTINE save_avalanche_params(params)
  TYPE(KORC_PARAMS), INTENT(IN) :: params
  CHARACTER(MAX_STRING_LENGTH) :: filename
  CHARACTER(MAX_STRING_LENGTH) :: gname
  CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
  CHARACTER(MAX_STRING_LENGTH) :: dset
  CHARACTER(MAX_STRING_LENGTH) :: attr
  INTEGER(HID_T) :: h5file_id
  INTEGER(HID_T) :: group_id
  INTEGER :: h5error
  REAL(rp) :: units

  if (params%mpi_params%rank .EQ. 0) then
     filename = TRIM(params%path_to_outputs) // "avalanche_parameters.h5"
     call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)

     gname = "pdf_params"
     call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

     dset = TRIM(gname) // "/max_pitch_angle"
     attr = "Maximum pitch angle in avalanche PDF (degrees)"
     call save_to_hdf5(h5file_id,dset,aval_params%max_pitch_angle,attr)

     dset = TRIM(gname) // "/min_pitch_angle"
     attr = "Minimum pitch angle in avalanche PDF (degrees)"
     call save_to_hdf5(h5file_id,dset,aval_params%min_pitch_angle,attr)

     dset = TRIM(gname) // "/min_energy"
     attr = "Minimum energy in avalanche PDF (eV)"
     units = 1.0_rp/C_E
     call save_to_hdf5(h5file_id,dset,units*aval_params%min_energy,attr)

     dset = TRIM(gname) // "/max_energy"
     attr = "Maximum energy in avalanche PDF (eV)"
     units = 1.0_rp/C_E
     call save_to_hdf5(h5file_id,dset,units*aval_params%max_energy,attr)

     dset = TRIM(gname) // "/max_p"
     attr = "Maximum momentum in avalanche PDF (me*c^2)"
     call save_to_hdf5(h5file_id,dset,aval_params%max_p,attr)

     dset = TRIM(gname) // "/min_p"
     attr = "Maximum momentum in avalanche PDF (me*c^2)"
     call save_to_hdf5(h5file_id,dset,aval_params%min_p,attr)

     dset = TRIM(gname) // "/ne"
     attr = "Background electron density (m^-3)"
     call save_to_hdf5(h5file_id,dset,aval_params%ne,attr)

     dset = TRIM(gname) // "/Zeff"
     attr = "Effective atomic number of ions."
     call save_to_hdf5(h5file_id,dset,aval_params%Zeff,attr)

     dset = TRIM(gname) // "/Ec"
     attr = "Critical electric field in (V/m)"
     call save_to_hdf5(h5file_id,dset,aval_params%Ec,attr)

     dset = TRIM(gname) // "/Epar"
     attr = "Parallel electric field in (V/m)"
     call save_to_hdf5(h5file_id,dset,aval_params%Epar,attr)

     dset = TRIM(gname) // "/Te"
     attr = "Background electron temperature (eV)"
     units = 1.0_rp/C_E
     call save_to_hdf5(h5file_id,dset,units*aval_params%Te,attr)

     dset = TRIM(gname) // "/lambda_D"
     attr = "Debye length (m)"
     call save_to_hdf5(h5file_id,dset,aval_params%lD,attr)

     dset = TRIM(gname) // "/bmin"
     attr = "Maximum approach radius (m)"
     call save_to_hdf5(h5file_id,dset,aval_params%bmin,attr)

     dset = TRIM(gname) // "/Clog"
     attr = "Coulomb logarithm"
     call save_to_hdf5(h5file_id,dset,aval_params%CoulombLog,attr)

     dset = TRIM(gname) // "/Tau"
     attr = "Collision time (s)"
     call save_to_hdf5(h5file_id,dset,aval_params%Tau,attr)

     call h5gclose_f(group_id, h5error)

     call h5fclose_f(h5file_id, h5error)
  end if
END SUBROUTINE save_avalanche_params

END MODULE korc_avalanche
