MODULE korc_velocity_distribution
  !! @note Module containing subroutines to initialize the velocity 
  !! distribution of the simulated particles. @endnote
  USE korc_types
  USE korc_constants
  USE korc_HDF5
  USE korc_hpc
  use korc_fields
  use korc_rnd_numbers
  use korc_hammersley_generator

  use korc_avalanche
  use korc_experimental_pdf
  use korc_energy_pdfs
  use korc_simple_equilibrium_pdf

  IMPLICIT NONE

  PUBLIC :: initial_gyro_distribution,&
       thermal_distribution,&
       initial_energy_pitch_dist
  PRIVATE :: fth_3V,&
       random_norm,&
       gyro_distribution

CONTAINS


  FUNCTION fth_3V(Vth,V)
    !! @note Function used to sample the probability density function of a 
    !! thermal plasma in the 3-dimensional velocity space. @endnote
    !! This function returns \(f_{T_e}(v) = \exp{\left( v^2/2v_{T_e}^2 \right)}\), 
    !! where \(v_{T_e} = \sqrt{T_e/m_e}\) is
    !! the temperature of the thermal electrons, and \(v = |\mathbf{v}|\) 
    !! is the speed of the sampled electron.
    REAL(rp), DIMENSION(3), INTENT(IN) 	:: V
    !! Velocity of the sampled electron \(\mathbf{v}\).
    REAL(rp), INTENT(IN) 			:: Vth
    !! Thermal velocity of the background electrons \(v_{T_e}\).
    REAL(rp) 				:: fth_3V
    !! Value of \(f_{T_e}(v)\).

    fth_3V = EXP(-0.5_rp*DOT_PRODUCT(V,V)/Vth**2.0_rp)
  END FUNCTION fth_3V


  FUNCTION random_norm(mu,sigma)
    !! @note Gaussian random number generator. @endnote
    !! This function returns a deviate of a Gaussian distribution
    !! $$f_G(x;\mu,\sigma) = 
    !! \frac{1}{\sigma\sqrt{2\pi}} \exp{\left( -(x-\mu)^2/2\sigma^2 \right)},$$
    !!
    !! with mean \(\mu\), and standard deviation \(\sigma\).
    !!
    !! We use the Inverse Transform Sampling Method for sampling \(x\). 
    !! With this method we get \(x = \sqrt{-2\log{(1-y)}}\cos(2\pi z)\),
    !! where \(y\) and \(z\) are uniform random numbers in the interval \([0,1]\).
    REAL(rp), INTENT(IN) 	:: mu
    !! Mean value \(\mu\) of the Gaussian distribution.
    REAL(rp), INTENT(IN) 	:: sigma
    !! Standard deviation \(\sigma\) of the Gaussian distribution.
    REAL(rp) 				:: random_norm
    !! Sampled number \(x\) from the Gaussian distribution \(f_G(x;\mu,\sigma)\).
    REAL(rp) 				:: rand1
    !! Uniform random number in the interval \([0,1]\).
    REAL(rp) 				:: rand2
    !! Uniform random number in the interval \([0,1]\).

    call RANDOM_NUMBER(rand1)
    call RANDOM_NUMBER(rand2)

    random_norm = SQRT(-2.0_rp*LOG(1.0_rp-rand1))*COS(2.0_rp*C_PI*rand2);
  END FUNCTION random_norm


  subroutine thermal_distribution(params,spp)
    !! @note Subroutine that samples a thermal distribution function
    !! of electrons for generating the initial condition of a set of
    !! simulated particles. @endnote
    !! This subroutine uses the Inverse Transform Sampling Method along
    !! with the  Metropolis-Hastings algorithm to generate an
    !! initial condition of the velocity distribution that follows a
    !! 3-dimensional (in velocity space) thermal distribution.
    !! @todo Check that the gyro-distribution is initialized right in
    !! this function.
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    !! Core KORC simulation parameters.
    TYPE(SPECIES), INTENT(INOUT) 	:: spp
    !! An instance of the derived type SPECIES containing all the
    !! parameters and simulation variables of the different species
    !! in the simulation.
    REAL(rp) 				:: Vmax
    !! Velocity cutoff where we stop sampling the tail of the thermal
    !! distribution.
    REAL(rp) 				:: vsq
    !! 
    REAL(rp) 				:: Vth
    !! Thermal velocity of the sampled distribution \(v_{T_e} =
    !! \sqrt{T_e/m_e}\).
    REAL(rp) 				:: sv
    !! Step to sample the velocity space using the Metropolis-Hastings
    !! algorithm.
    REAL(rp) 				:: ratio
    !! Ratio used to accept or reject a sampling in the Metropolis-Hastings
    !! algorithm.
    REAL(rp) 				:: rand_unif
    !! Uniform random deviate in the interval  \([0,1]\).
    REAL(rp), DIMENSION(3) 		:: V
    !! Sampled velocity.
    REAL(rp), DIMENSION(3) 		:: U
    !! Sampled velocity.
    REAL(rp), DIMENSION(3) 		:: b = (/1.0_rp,0.0_rp,0.0_rp/)
    !! Temporary variable representing a unit vector along the \(x\)-axis.
    INTEGER 				:: ii
    !! Iterator.
    INTEGER 				:: ppp
    !! Number of particles per species.
    
    Vmax = 0.9_rp
    Vth = SQRT(spp%Eo*ABS(spp%q)/spp%m)
    ppp = spp%ppp

    V = (/0.0_rp,0.0_rp,0.0_rp/)
    sv = Vth/10.0_rp

    ii=2_idef
    do while (ii .LE. 1000_idef)
       U(1) = V(1) + random_norm(0.0_rp,sv)
       do while (ABS(U(1)) .GT. Vmax)
          U(1) = V(1) + random_norm(0.0_rp,sv)
       end do

       U(2) = V(2) + random_norm(0.0_rp,sv)
       do while (ABS(U(2)) .GT. Vmax)
          U(2) = V(2) + random_norm(0.0_rp,sv)
       end do

       U(3) = V(3) + random_norm(0.0_rp,sv)
       do while (ABS(U(3)) .GT. Vmax)
          U(3) = V(3) + random_norm(0.0_rp,sv)
       end do

       ratio = fth_3V(Vth,U)/fth_3V(Vth,V)

       if (ratio .GE. 1.0_rp) then
          V = U
          ii = ii + 1_idef
       else
          call RANDOM_NUMBER(rand_unif)
          if (ratio .GT. rand_unif) then
             V = U
             ii = ii + 1_idef
          end if
       end if
    end do

    spp%vars%V(1,1) = V(1)
    spp%vars%V(1,2) = V(2)
    spp%vars%V(1,3) = V(3)
    ii=2_idef
    do while (ii .LE. ppp)
       U(1) = spp%vars%V(ii-1,1) + random_norm(0.0_rp,sv)
       do while (ABS(U(1)) .GT. Vmax)
          U(1) = spp%vars%V(ii-1,1) + random_norm(0.0_rp,sv)
       end do
       U(2) = spp%vars%V(ii-1,2) + random_norm(0.0_rp,sv)
       do while (ABS(U(2)) .GT. Vmax)
          U(2) = spp%vars%V(ii-1,2) + random_norm(0.0_rp,sv)
       end do
       U(3) = spp%vars%V(ii-1,3) + random_norm(0.0_rp,sv)
       do while (ABS(U(3)) .GT. Vmax)
          U(3) = spp%vars%V(ii-1,3) + random_norm(0.0_rp,sv)
       end do

       ratio = fth_3V(Vth,U)/fth_3V(Vth,spp%vars%V(ii-1,:))

       if (ratio .GE. 1.0_rp) then
          spp%vars%V(ii,1) = U(1)
          spp%vars%V(ii,2) = U(2)
          spp%vars%V(ii,3) = U(3)
          ii = ii + 1_idef
       else
          call RANDOM_NUMBER(rand_unif)
          if (ratio .GT. rand_unif) then
             spp%vars%V(ii,1) = U(1)
             spp%vars%V(ii,2) = U(2)
             spp%vars%V(ii,3) = U(3)
             ii = ii + 1_idef
          end if
       end if
    end do

    do ii=1_idef,ppp
       vsq = spp%vars%V(ii,1)*spp%vars%V(ii,1)             &
            + spp%vars%V(ii,2)*spp%vars%V(ii,2)             &
            + spp%vars%V(ii,3)*spp%vars%V(ii,3)
       spp%vars%g(ii) = 1.0_rp/SQRT(1.0_rp - vsq)
       spp%vars%eta(ii) = ACOS(spp%vars%V(ii,1)/SQRT(vsq))
    end do

    spp%go = spp%Eo/(spp%m*C_C**2)
    spp%etao = 90.0_rp
  end subroutine thermal_distribution


  subroutine initial_energy_pitch_dist(params,spp)
    !! @note Subroutine that calls subroutines of different modules to 
    !! initialize the energy and pitch-angle distribution in various ways. @endnote
    TYPE(KORC_PARAMS), INTENT(IN) 				:: params
    !! Core KORC simulation parameters.
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)       :: spp
    !! An instance of the derived type SPECIES containing all the parameters and 
    !! simulation variables of the different species in the simulation.
    INTEGER 							:: ii
    !! Species iterator.
    INTEGER 							:: mpierr
    !! MPI error status.

    do ii=1_idef,params%num_species

       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'(/,"* * * * * SPECIES: ",I2," * * * * * * * * * * *")') ii
          write(output_unit_write,'("Particles per MPI process: ",I16)') spp(ii)%ppp
          write(output_unit_write,'("Energy distribution is: ",A20)') &
               TRIM(spp(ii)%energy_distribution)
          write(output_unit_write,'("Pitch-angle distribution is: ",A20)') &
               TRIM(spp(ii)%pitch_distribution)
          write(output_unit_write,'("Spatial distribution is: ",A20)') &
               TRIM(spp(ii)%spatial_distribution)
          write(output_unit_write,'("* * * * * * * * * * * * * * * * * * * * * *",/)')
       end if
       
       SELECT CASE (TRIM(spp(ii)%energy_distribution))
       CASE ('MONOENERGETIC')
          spp(ii)%go = (spp(ii)%Eo + spp(ii)%m*C_C**2)/(spp(ii)%m*C_C**2)

          spp(ii)%vars%g = spp(ii)%go ! Monoenergetic
          spp(ii)%Eo_lims = (/spp(ii)%Eo, spp(ii)%Eo /)
          
       CASE ('THERMAL')
          call thermal_distribution(params,spp(ii))

          spp(ii)%Eo_lims = (/spp(ii)%m*C_C**2*MINVAL(spp(ii)%vars%g) - &
               spp(ii)%m*C_C**2, &
               spp(ii)%m*C_C**2*MAXVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2 /)
       CASE ('AVALANCHE')
          call get_avalanche_distribution(params,spp(ii)%vars%g, &
               spp(ii)%vars%eta,spp(ii)%go,spp(ii)%etao)
          spp(ii)%Eo = spp(ii)%m*C_C**2*spp(ii)%go - spp(ii)%m*C_C**2
          spp(ii)%Eo_lims = (/spp(ii)%m*C_C**2*MINVAL(spp(ii)%vars%g) &
               - spp(ii)%m*C_C**2, &
               spp(ii)%m*C_C**2*MAXVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2 /)
       CASE ('HOLLMANN')
          call get_Hollmann_distribution(params,spp(ii))          
!          spp(ii)%Eo = spp(ii)%m*C_C**2*spp(ii)%go - spp(ii)%m*C_C**2
          spp(ii)%go = (spp(ii)%Eo + spp(ii)%m*C_C**2)/(spp(ii)%m*C_C**2)
          spp(ii)%Eo_lims = (/spp(ii)%m*C_C**2*MINVAL(spp(ii)%vars%g) &
               - spp(ii)%m*C_C**2, &
               spp(ii)%m*C_C**2*MAXVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2 /)
       CASE ('EXPERIMENTAL-GAMMA')
          call get_experimentalG_distribution(params,spp(ii)%vars%g, &
               spp(ii)%vars%eta, &
               spp(ii)%go,spp(ii)%etao)
          spp(ii)%Eo = spp(ii)%m*C_C**2*spp(ii)%go - spp(ii)%m*C_C**2
          spp(ii)%Eo_lims = (/spp(ii)%m*C_C**2*MINVAL(spp(ii)%vars%g) &
               - spp(ii)%m*C_C**2, &
               spp(ii)%m*C_C**2*MAXVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2 /)
       CASE ('GAMMA')
          call get_gamma_distribution(params,spp(ii)%vars%g,spp(ii)%go)

          spp(ii)%Eo = spp(ii)%m*C_C**2*spp(ii)%go - spp(ii)%m*C_C**2
          spp(ii)%Eo_lims = (/spp(ii)%m*C_C**2*MINVAL(spp(ii)%vars%g) &
               - spp(ii)%m*C_C**2, &
               spp(ii)%m*C_C**2*MAXVAL(spp(ii)%vars%g) - spp(ii)%m*C_C**2 /)
       CASE ('UNIFORM')
          spp(ii)%Eo = spp(ii)%Eo_lims(1)
          spp(ii)%go = (spp(ii)%Eo + spp(ii)%m*C_C**2)/(spp(ii)%m*C_C**2)

          call generate_2D_hammersley_sequence(params%mpi_params%rank, &
               params%mpi_params%nmpi,spp(ii)%vars%g,spp(ii)%vars%eta)

          spp(ii)%vars%g = (spp(ii)%Eo_lims(2) - & 
               spp(ii)%Eo_lims(1))*spp(ii)%vars%g/(spp(ii)%m*C_C**2) + &
               (spp(ii)%Eo_lims(1) + spp(ii)%m*C_C**2)/(spp(ii)%m*C_C**2)
       CASE ('AVALANCHE-4D')
          spp(ii)%go = (spp(ii)%Eo + spp(ii)%m*C_C**2)/(spp(ii)%m*C_C**2)
          spp(ii)%vars%g = spp(ii)%go
          ! Monoenergy from input file until sampled in Avalanche_4D
       CASE ('HOLLMANN-3D')
          spp(ii)%go = (spp(ii)%Eo + spp(ii)%m*C_C**2)/(spp(ii)%m*C_C**2)
          spp(ii)%vars%g = spp(ii)%go
          ! Monoenergy from input file until sampled in Hollmann_3D
       CASE ('HOLLMANN-3D-PSI')
          spp(ii)%go = (spp(ii)%Eo + spp(ii)%m*C_C**2)/(spp(ii)%m*C_C**2)
          spp(ii)%vars%g = spp(ii)%go
          ! Monoenergy from input file until sampled in Hollmann_3D
       CASE ('HOLLMANN-1DTRANSPORT')
          spp(ii)%go = (spp(ii)%Eo + spp(ii)%m*C_C**2)/(spp(ii)%m*C_C**2)
          spp(ii)%vars%g = spp(ii)%go
          ! Monoenergy from input file until sampled in Hollmann_3D
       CASE ('FIO_therm')
          spp(ii)%go = (spp(ii)%Eo + spp(ii)%m*C_C**2)/(spp(ii)%m*C_C**2)
          spp(ii)%vars%g = spp(ii)%go
       CASE DEFAULT
          ! Something to be done
       END SELECT

       call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

       SELECT CASE (TRIM(spp(ii)%pitch_distribution))
       CASE ('MONOPITCH')
          spp(ii)%vars%eta = spp(ii)%etao ! Mono-pitch-angle

          if(spp(ii)%spatial_distribution.ne.'SPONG-3D') then
             spp(ii)%etao_lims = (/spp(ii)%etao , spp(ii)%etao/)
          end if
          
       CASE ('THERMAL')
          spp(ii)%etao_lims = (/MINVAL(spp(ii)%vars%eta), &
               MAXVAL(spp(ii)%vars%eta)/)
       CASE ('AVALANCHE')
          spp(ii)%etao_lims = (/MINVAL(spp(ii)%vars%eta), &
               MAXVAL(spp(ii)%vars%eta)/)
       CASE ('HOLLMANN')
!          spp(ii)%vars%eta = spp(ii)%etao
!          spp(ii)%etao_lims = (/MINVAL(spp(ii)%vars%eta), &
!               MAXVAL(spp(ii)%vars%eta)/)
       CASE ('EXPERIMENTAL-GAMMA')
          spp(ii)%etao_lims = (/MINVAL(spp(ii)%vars%eta), &
               MAXVAL(spp(ii)%vars%eta)/)
       CASE ('UNIFORM')
          spp(ii)%etao = spp(ii)%etao_lims(1)

          spp(ii)%vars%eta = (spp(ii)%etao_lims(2) - &
               spp(ii)%etao_lims(1))*spp(ii)%vars%eta + spp(ii)%etao_lims(1)
       CASE ('SIMPLE-EQUILIBRIUM')
          call get_equilibrium_distribution(params,spp(ii)%vars%eta, &
               spp(ii)%go,spp(ii)%etao)

          spp(ii)%etao_lims = (/MINVAL(spp(ii)%vars%eta),&
               MAXVAL(spp(ii)%vars%eta)/)
       CASE ('AVALANCHE-4D')
          spp(ii)%vars%eta = spp(ii)%etao
          !Monopitch from input file until sampled in Avalanche_4D
       CASE ('HOLLMANN-3D')
          spp(ii)%vars%eta = spp(ii)%etao
          !Monopitch from input file until sampled in Hollmann_3D
       CASE ('HOLLMANN-3D-PSI')
          spp(ii)%vars%eta = spp(ii)%etao
          !Monopitch from input file until sampled in Hollmann_3D
       CASE ('HOLLMANN-1DTRANSPORT')
          spp(ii)%vars%eta = spp(ii)%etao
          !Monopitch from input file until sampled in Hollmann_3D
       CASE ('SPONG-3D')
          spp(ii)%vars%eta = spp(ii)%etao
          !Monopitch from input file until sampled in Spong_3D
       CASE ('FIO_therm')
          spp(ii)%vars%eta = spp(ii)%etao
       CASE DEFAULT
          ! Something to be done
       END SELECT



       call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
    end do
  end subroutine initial_energy_pitch_dist


  subroutine gyro_distribution(params,F,spp)

    USE, INTRINSIC :: iso_c_binding
    
    !! @Note Subroutine that initializes the gyro-angle distribution 
    !! of the particles. @endnote
    !! When evolving the particles in the 6-D phase space, in addition to 
    !! the position (3 degrees of freedom), energy (one degree of freedom), 
    !! pitch angle (one degree of freedom), we need to define the gyro-angle 
    !! of the particle (one degree of freedom), which is given by the pitch 
    !! angle and the direction of the local magnetic field. By default, this 
    !! subroutine generates a uniform gyro-angle distribution.
    !! @note Notice that all the simulation variables are normalized
    !! here. @endnote
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN) 		:: F
    !! An instance of the KORC derived type FIELDS. This structure 
    !! has the information of the magnetic field.
    TYPE(SPECIES), INTENT(INOUT) 		:: spp
    !! An instance of the derived type SPECIES containing all the
    !! parameters and 
    !! simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE	:: b1
    !! Basis vector pointing along the local magnetic field, that is,
    !!  along \(\mathbf{b} = \mathbf{B}/B\).
    REAL(rp), DIMENSION(:,:), ALLOCATABLE	:: b2
    !! Basis vector perpendicular to b1
    REAL(rp), DIMENSION(:,:), ALLOCATABLE	:: b3
    !! Basis vector perpendicular to b1 and b2.
    REAL(rp), DIMENSION(:), ALLOCATABLE 	:: Vo
    !! Initial particle speed.
    REAL(rp), DIMENSION(:), ALLOCATABLE 	:: V1
    !!  Velocity component along b1.
    REAL(rp), DIMENSION(:), ALLOCATABLE 	:: V2
    !! Velocity component along b2.
    REAL(rp), DIMENSION(:), ALLOCATABLE 	:: V3
    !! Velocity component along b3.
    REAL(rp), DIMENSION(:), ALLOCATABLE 	:: theta
    !! Uniform random number in the interval \([0,2\pi]\) 
    !! representing the gyro-angle.
    INTEGER 				:: jj
    !! Particle iterator.
    TYPE(C_PTR), DIMENSION(:), ALLOCATABLE    :: hint
    
    ALLOCATE(Vo(spp%ppp))
    ALLOCATE(V1(spp%ppp))
    ALLOCATE(V2(spp%ppp))
    ALLOCATE(V3(spp%ppp))
    ALLOCATE(b1(spp%ppp,3))
    ALLOCATE(b2(spp%ppp,3))
    ALLOCATE(b3(spp%ppp,3))
    ALLOCATE(hint(spp%ppp))

    hint=C_NULL_PTR
#ifdef FIO
    hint=spp%vars%hint
#endif
    
    ALLOCATE( theta(spp%ppp) )

    ! * * * * INITIALIZE VELOCITY * * * *

    call init_random_seed()
    call RANDOM_NUMBER(theta)
    theta = 2.0_rp*C_PI*theta
    
    if (spp%spatial_distribution.eq.'TRACER') theta=2.0*C_PI


    Vo = SQRT( 1.0_rp - 1.0_rp/(spp%vars%g(:)**2) )
    V1 = Vo*COS(C_PI*spp%vars%eta/180.0_rp)
    V2 = Vo*SIN(C_PI*spp%vars%eta/180.0_rp)*COS(theta)
    V3 = Vo*SIN(C_PI*spp%vars%eta/180.0_rp)*SIN(theta)

    !write(6,*) 'V123',V1,V2,V3
    
    !do jj=1_idef,spp%ppp
    !   write(6,*) 'MPI',params%mpi_params%rank,'X', &
    !        spp%vars%X(jj,:)*params%cpp%length
    !end do
    call unitVectors(params,spp%vars%X,F,b1,b2,b3,spp%vars%flagCon, &
         spp%vars%cart,hint)
    !! Call to subroutine [[unitVectors]] in [[korc_fields]].

    !
    !write(6,*) 'X',spp%vars%X
    !write(6,*) 'b-hat',b1
    !write(6,*) 'b-1',b2
    !write(6,*) 'b-2',b3

    
    do jj=1_idef,spp%ppp
       if ( spp%vars%flagCon(jj) .EQ. 1_idef ) then    
          spp%vars%V(jj,1) = V1(jj)*b1(jj,1) + V2(jj)*b2(jj,1) + V3(jj)*b3(jj,1)
          spp%vars%V(jj,2) = V1(jj)*b1(jj,2) + V2(jj)*b2(jj,2) + V3(jj)*b3(jj,2)
          spp%vars%V(jj,3) = V1(jj)*b1(jj,3) + V2(jj)*b2(jj,3) + V3(jj)*b3(jj,3)
       end if
    end do
    
    !write(6,'("Vx: ",E17.10)') spp%vars%V(:,1)
    !write(6,'("Vy: ",E17.10)') spp%vars%V(:,2)
    !write(6,'("Vz: ",E17.10)') spp%vars%V(:,3)
    
    DEALLOCATE(theta)
    DEALLOCATE(Vo)
    DEALLOCATE(V1)
    DEALLOCATE(V2)
    DEALLOCATE(V3)
    DEALLOCATE(b1)
    DEALLOCATE(b2)
    DEALLOCATE(b3)
    DEALLOCATE(hint)
  end subroutine gyro_distribution


  subroutine initial_gyro_distribution(params,F,spp)
    !! @note Subroutine that works as an interface for initializing various 
    !! gyro-angle distributions for the different simulated particle
    !! species. @endnote
    !! @todo At this moment this subroutine only calls the subroutine
    !! to generate 
    !! a uniform gyro-angle distribution. This will be modified later. @endtodo
    TYPE(KORC_PARAMS), INTENT(IN) 				:: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN) 					:: F
    !! An instance of the KORC derived type FIELDS. This structure has 
    !! the information of the magnetic field.
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: spp
    !! An instance of the derived type SPECIES containing all the parameters 
    !! and simulation variables of the different species in the simulation.
    INTEGER 							:: ss
    !! Species iterator.

    do ss=1_idef,params%num_species
       SELECT CASE (TRIM(spp(ss)%energy_distribution))
       CASE ('THERMAL')
          !Nothing, all was done in initialize_particles through
          !thermal_distribution
       CASE DEFAULT            
          call gyro_distribution(params,F,spp(ss))
       END SELECT
    end do
  end subroutine initial_gyro_distribution

END MODULE korc_velocity_distribution
