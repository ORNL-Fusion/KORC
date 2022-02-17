MODULE korc_spatial_distribution
  !! @note Module with subroutines for generating the initial spatial distribution 
  !! of the different partciles' species in the simulation. @endnote
  USE korc_types
  USE korc_constants
  USE korc_HDF5
  USE korc_hpc
  use korc_fields
  use korc_profiles
  use korc_rnd_numbers
  use korc_random
  use korc_hammersley_generator
  use korc_avalanche
  use korc_experimental_pdf
  USE korc_input

  IMPLICIT NONE

  REAL(rp), PRIVATE, PARAMETER :: minmax_buffer_size = 10.0_rp
  
  PUBLIC :: intitial_spatial_distribution
  PRIVATE :: uniform,&
       disk,&
       torus,&
       elliptic_torus,&
       exponential_elliptic_torus,&
       gaussian_elliptic_torus,&
       exponential_torus,&
       gaussian_torus,&
       fzero,&
       MH_gaussian_elliptic_torus,&
       indicator,&
       PSI_ROT,&
       Spong_3D,&
       Spong_2D, &
       MH_psi


CONTAINS

subroutine uniform(spp)
  !! @note Initializing to zero the particles' position when 
  !! simulating a 'UNIFORM' plasma. @endnote
  !! Even though in a simulation of a uniform plasma the particles' 
  !! position is not advanced, we initialize their position to zero.
  !! @todo Modify KORC for not allocating the particles' position 
  !! spp%vars%X and to do not use it along the simulation.
  TYPE(SPECIES), INTENT(INOUT) :: spp
    !! An instance of the derived type SPECIES containing all the 
    !! parameters and simulation variables of the different
    !!species in the simulation.

  spp%vars%X = 0.0_rp
end subroutine uniform


subroutine disk(params,spp)
  !! @note Subrotuine for generating a uniform disk/ring as the 
  !! initial spatial condition of a given species of particles  
  !! in the simulation. @endnote
  !! This uniform disk/ring distribution is generated using the 
  !! Inverse Transform Sampling method. In this case, the (toroidal) 
  !! radial distribution function of the particles is:
  !!
  !! $$f(r) = \left\{ \begin{array}{ll} 0 & r<r_{min} \\
  !! \frac{1}{2\pi^2(r_{max}^2-r_{min}^2)R_0} 
  !! & r_{min}<r<r_{max} \\ 0 & r>r_{max} \end{array} \right.,$$
  !!
  !! where \(r_{min}\) and \(r_{max}\) are the inner and outer 
  !! radius of the uniform ring distribution, and \(R_0\) is the 
  !! cylindrical radial position of the center of the disk/ring distribution.
  !! This distribution is so that \(\int_0^{2\pi}\int_{r_{min}}^{r_{max}} f(r)
  !! J(r,\theta) drd\theta = 1 \), where \(\theta\) is the poloidal angle,
  !! and \(J(r,\theta)=r(R_0 + r\cos\theta)\) is the Jacobian of the 
  !! transformation of Cartesian coordinates to toroidal coordinates.
  !! Notice that in the case of a disk \(r_{min}=0\). As a convention, 
  !! this spatial distribution will be generated on the \(xz\)-plane.
  !! Using the Inverse Transform Sampling method we sample \(f(r)\), 
  !! and obtain the radial position of the particles as \(r = \sqrt{(r_{max}^2 
  !! - r_{min}^2)U + r_{min}^2}\), where \(U\) is a uniform deviate in \([0,1]\).
  TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    !! Core KORC simulation parameters.
  TYPE(SPECIES), INTENT(INOUT) 		:: spp
    !! An instance of the derived type SPECIES containing all 
    !! the parameters and simulation variables of the different 
    !! species in the simulation.
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: r
    !! Radial position of the particles \(r\).
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: theta
    !! Uniform deviates in the range \([0,2\pi]\) representing 
    !! the uniform poloidal angle \(\theta\) distribution of the particles.

  ALLOCATE( theta(spp%ppp) )
  ALLOCATE( r(spp%ppp) )

  ! Initial condition of uniformly distributed particles on a disk in the xz-plane
  ! A unique velocity direction
  call init_u_random(10986546_8)

  call init_random_seed()
  call RANDOM_NUMBER(theta)
  
  theta = 2.0_rp*C_PI*theta

  ! Uniform distribution on a disk at a fixed azimuthal theta
  call init_random_seed()
  call RANDOM_NUMBER(r)

  r = SQRT((spp%r_outter**2 - spp%r_inner**2)*r + spp%r_inner**2)
  spp%vars%X(:,1) = ( spp%Ro + r*COS(theta) )*COS(spp%PHIo)
  spp%vars%X(:,2) = ( spp%Ro + r*COS(theta) )*SIN(spp%PHIo)
  spp%vars%X(:,3) = spp%Zo + r*SIN(theta)

  DEALLOCATE(theta)
  DEALLOCATE(r)
end subroutine disk

subroutine torus(params,spp)
    TYPE(KORC_PARAMS), INTENT(IN)       :: params
    TYPE(SPECIES), INTENT(INOUT)        :: spp
    REAL(rp), DIMENSION(:), ALLOCATABLE :: r
    REAL(rp), DIMENSION(:), ALLOCATABLE :: theta
    REAL(rp), DIMENSION(:), ALLOCATABLE :: zeta

    ALLOCATE( theta(spp%ppp) )
    ALLOCATE( zeta(spp%ppp) )
    ALLOCATE( r(spp%ppp) )

    ! Initial condition of uniformly distributed particles on a disk in the xz-plane
    ! A unique velocity direction
!    call init_u_random(10986546_8)

    call init_random_seed()
!    call RANDOM_NUMBER(theta)
!    theta = 2.0_rp*C_PI*theta
    
    call set_random_dist(0.0_rp, 2.0_rp*C_PI)
    call get_randoms(theta)

!    call init_random_seed()
!    call RANDOM_NUMBER(zeta)
!    zeta = 2.0_rp*C_PI*zeta

    call get_randoms(zeta)
!
    ! Uniform distribution on a disk at a fixed azimuthal theta
!    call init_random_seed()
!    call RANDOM_NUMBER(r)

    call set_random_dist(0.0_rp, 1.0_rp)
    call get_randoms(r)

    r = SQRT((spp%r_outter**2 - spp%r_inner**2)*r + spp%r_inner**2)

!$OMP PARALLEL WORKSHARE
    spp%vars%X(:,1) = ( spp%Ro + r*COS(theta) )*SIN(zeta)
    spp%vars%X(:,2) = ( spp%Ro + r*COS(theta) )*COS(zeta)
    spp%vars%X(:,3) = spp%Zo + r*SIN(theta)
!$OMP END PARALLEL WORKSHARE

    DEALLOCATE(theta)
    DEALLOCATE(zeta)
    DEALLOCATE(r)
end subroutine torus

!subroutine torus(params,spp)
  !! @note Subrotuine for generating a uniform torus/torus 
  !! shell as the initial spatial condition of a given species 
  !! of particles in the simulation.@endnote
  !! This distribution is generated using the Inverse Transform 
  !! Sampling method. This distribution follows the same radial 
  !! distribution of a uniform disk/ring distribution, see the 
  !! documentation of the [[disk]] subroutine.
!  TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    !! Core KORC simulation parameters.
!  TYPE(SPECIES), INTENT(INOUT) 		:: spp
    !! An instance of the derived type SPECIES 
    !! containing all the parameters and simulation variables of the 
    !! different species in the simulation.
!  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: r
    !! Radial position of the particles \(r\).
!  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: theta
  !! Uniform deviates in the range \([0,2\pi]\) 
  !! representing the uniform poloidal angle \(\theta\)
  !! distribution of the particles.
!  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: zeta
  !! Uniform deviates in the range \([0,2\pi]\) representing 
  !! the uniform toroidal angle \(\zeta\) distribution of the particles.
!  INTEGER,DIMENSION(33) :: seed=(/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)

!  ALLOCATE( theta(spp%ppp) )
!  ALLOCATE( zeta(spp%ppp) )
!  ALLOCATE( r(spp%ppp) )

  ! Initial condition of uniformly distributed particles on a disk in the xz-plane
  ! A unique velocity direction
!  call init_u_random(10986546_8)

!  if (.not.params%SameRandSeed) then
!     call init_random_seed()
!  else
!     call random_seed(put=seed)
!  end if
!  call RANDOM_NUMBER(theta)
!  theta = 2.0_rp*C_PI*theta

!  if (.not.params%SameRandSeed) then
!     call init_random_seed()
!  else
!     call random_seed(put=seed)
!  end if
!  call RANDOM_NUMBER(zeta)
!  zeta = 2.0_rp*C_PI*zeta

  !write(6,*) 'Ro',spp%Ro*params%cpp%length
  !write(6,*) 'Zo',spp%Zo*params%cpp%length
  !write(6,*) 'r_in',spp%r_inner*params%cpp%length
  !write(6,*) 'r_out',spp%r_outter*params%cpp%length
  
  ! Uniform distribution on a disk at a fixed azimuthal theta
!  if (.not.params%SameRandSeed) then
!     call init_random_seed()
!  else
!     call random_seed(put=seed)
!  end if
!  call RANDOM_NUMBER(r)

!  r = SQRT((spp%r_outter**2 - spp%r_inner**2)*r + spp%r_inner**2)
!  spp%vars%X(:,1) = ( spp%Ro + r*COS(theta) )*SIN(zeta)
!  spp%vars%X(:,2) = ( spp%Ro + r*COS(theta) )*COS(zeta)
!  spp%vars%X(:,3) = spp%Zo + r*SIN(theta)

  !write(6,*) 'r_sam',r*params%cpp%length
  !write(6,*) 'R_sam',sqrt(spp%vars%X(:,1)**2+spp%vars%X(:,2)**2)*params%cpp%length

!  DEALLOCATE(theta)
!  DEALLOCATE(zeta)
!  DEALLOCATE(r)
  
!end subroutine torus


subroutine elliptic_torus(params,spp)
  !! @note Subroutine for generating a uniform elliptic torus as the initial 
  !! spatial condition of a given particle species in the simulation. @endnote
  !! An initial spatial distribution following the uniform distribution of 
  !! [[torus]] is modified through a shear transformation and a rotation to 
  !! generate a uniform spatial distribution on tori with elliptic cross sections. 
  !! First, we obtain the uniform spatial distribution in a torus of minor radius 
  !! \(r_0\), see [[torus]]. Then, we perform a shear transformation that changes 
  !! the cross section of the torus from circular to a tilted ellipse. In 
  !! cylindrical coordinates this shear transformation is given by:
  !!
  !! $$R' = R + \alpha Z,$$
  !! $$Z' = Z,$$
  !!
  !! where \(\alpha\) is the shear factor of the transformation. 
  !! Here, \(R\) and \(Z\) are the radial and vertical position of the particles
  !! uniformly distributed in a circular torus, \(R'\) and \(Z'\) are their
  !! new positions when following a uniform distribution in a torus with
  !! elliptic circular cross section. The center of the ellipse is 
  !! \(R_0' = R_0 + \alpha Z_0\), and \(Z_0 = Z_0\), where \(R_0\) and \(Z_0\)
  !! is the center of the initial circular torus. The major and minor semi-axes 
  !! of the tilted ellipse cross section is:
  !!
  !! $$a' = \left[ - \frac{2r_0^2}{\alpha \sqrt{\alpha^2 + 4} - (2+\alpha^2)}
  !!  \right]^{1/2},$$
  !! $$b' = \left[ \frac{2r_0^2}{\alpha \sqrt{\alpha^2 + 4} + (2+\alpha^2)}
  !!  \right]^{1/2}.$$
  !!
  !! Finally, we rotate the ellipse cross section anticlockwise along 
  !! \((R_0',Z_0')\) by \(\Theta = \cot^{-1}(\alpha/2)/2\), so the major semi-axis is
  !! parallel to the \(Z\)-axis.
  TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    !! Core KORC simulation parameters.
  TYPE(SPECIES), INTENT(INOUT) 		:: spp
    !! An instance of the derived type SPECIES containing all the parameters 
    !! and simulation variables of the different species in the simulation.
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: rotation_angle
    !! This is the angle \(\Theta\) in the equations above.
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: theta
    !! Uniform deviates in the range \([0,2\pi]\) representing the uniform 
    !! poloidal angle \(\theta\) distribution of the particles.
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: r
    !! Radial position of the particles \(r\).
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: zeta
    !! Uniform deviates in the range \([0,2\pi]\) representing 
    !! the uniform toroidal angle \(\zeta\) distribution of the particles.
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: X
    !! Auxiliary vector used in the coordinate transformations.
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: Y
    !! Auxiliary vector used in the coordinate transformations.
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: X1
    !! Auxiliary vector used in the coordinate transformations.
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: Y1
    !! Auxiliary vector used in the coordinate transformations.

  ALLOCATE(X1(spp%ppp))
  ALLOCATE(Y1(spp%ppp))
  ALLOCATE(X(spp%ppp))
  ALLOCATE(Y(spp%ppp))
  ALLOCATE(rotation_angle(spp%ppp))
  ALLOCATE(theta(spp%ppp))
  ALLOCATE(zeta(spp%ppp))
  ALLOCATE(r(spp%ppp))

  ! Initial condition of uniformly distributed particles on a disk in the xz-plane
  ! A unique velocity direction
  call init_u_random(10986546_8)

  call init_random_seed()
  call RANDOM_NUMBER(theta)
  theta = 2.0_rp*C_PI*theta

  call init_random_seed()
  call RANDOM_NUMBER(zeta)
  zeta = 2.0_rp*C_PI*zeta

  ! Uniform distribution on a disk at a fixed azimuthal theta
  call init_random_seed()
  call RANDOM_NUMBER(r)

  r = SQRT((spp%r_outter**2 - spp%r_inner**2)*r + spp%r_inner**2)

  Y = r*SIN(theta)
  X = r*COS(theta) + spp%shear_factor*Y

  !> @todo Modify this approximation.
  rotation_angle = 0.5_rp*C_PI - ATAN(1.0_rp,1.0_rp + spp%shear_factor); 

  X1 = X*COS(rotation_angle) - Y*SIN(rotation_angle) + spp%Ro
  Y1 = X*SIN(rotation_angle) + Y*COS(rotation_angle) + spp%Zo

  spp%vars%X(:,1) = X1*SIN(zeta)
  spp%vars%X(:,2) = X1*COS(zeta)
  spp%vars%X(:,3) = Y1

  DEALLOCATE(X1)
  DEALLOCATE(Y1)
  DEALLOCATE(X)
  DEALLOCATE(Y)
  DEALLOCATE(rotation_angle)
  DEALLOCATE(theta)
  DEALLOCATE(zeta)
  DEALLOCATE(r)
end subroutine elliptic_torus

!> @note Function used to find the zeros of \(f(r)\) of \ref
!! korc_spatial_distribution.exponential_torus. @endnote
!!
!! @param f Value of function.
!! @param r Guess value of radial position of the particles.
!! @param a Minor radius of the toroidal distribution \(r_0\)
!! @param ko Decay rate of radial distribution, see \(f(r)\) of \ref korc_spatial_distribution.exponential_torus.
!! @param P Deviate of a random uniform distribution in the interval \([0,1]\).
FUNCTION fzero(r,a,ko,P) RESULT(f)
	REAL(rp) 				:: f
	REAL(rp), INTENT(IN) 	:: r
	REAL(rp), INTENT(IN) 	:: a
	REAL(rp), INTENT(IN) 	:: ko
	REAL(rp), INTENT(IN) 	:: P

	f = EXP(-ko*r)*(1.0_rp + r*ko) + ( 1.0_rp - EXP(-ko*a)*(1.0_rp + a*ko) )*P - 1.0_rp
END FUNCTION fzero


!> @brief Subroutine that generates a exponentially decaying radial distribution of particles in a circular cross-section torus of
!! major and minor radi \(R_0\) and \(r_0\), respectively.
!! @details We generate this exponentially decaying radial distribution \(f(r)\) following the same approach as in
!! \ref korc_spatial_distribution.disk, but this time, the radial distribution is given by:
!!
!!
!! $$f(r) = \left\{ \begin{array}{ll} \frac{k_0^2}{4\pi^2 R_0}\frac{\exp{\left( -k_0 r \right)}}{1 - \exp{\left( -k_0r_0\right)}\left[ 1 + k_0 r_0 \right]} &\
!! r<r_0 \\ 0 & r>r_0 \end{array} \right.$$
!!
!!
!! The radial position of the particles \(r\) is obtained using the Inverse Trasnform Sampling method, finding \(r\) numerically
!! through the Newton-Raphson method. First, we calculate the particles' radial distribution in a disk centered at \((R,Z) = (0,0)\).
!! Then, we transfor to a new set of coordinates where the disk is centered at \((R,Z) = (R_0,Z_0)\). Finally, we generate the
!! toroidal distribution by givin each particle a toroidal angle \(\zeta\) which follows a uniform distribution in the interval
!! \([0,2\pi]\).
!!
!! @param[in] params Core KORC simulation parameters.
!! @param[in,out] spp An instance of the derived type SPECIES containing all the parameters and simulation variables of the different
!! species in the simulation.
!! @param fl Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param fr Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param fm Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param rl Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param rr Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param rm Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param relerr Tolerance used to determine when to stop iterating the Newton-Raphson method for finding \(r\).
!! @param r Radial position of the particles \(r\).
!! @param theta Uniform deviates in the range \([0,2\pi]\) representing the uniform poloidal angle \(\theta\) distribution of the particles.
!! @param zeta Uniform deviates in the range \([0,2\pi]\) representing the uniform toroidal angle \(\zeta\) distribution of the particles.
!! @param pp Particle iterator.
subroutine exponential_torus(params,spp)
  TYPE(KORC_PARAMS), INTENT(IN) 		:: params
  TYPE(SPECIES), INTENT(INOUT) 		:: spp
  REAL(rp) 				:: fl
  REAL(rp) 				:: fr
  REAL(rp) 				:: fm
  REAL(rp) 				:: rl
  REAL(rp) 				:: rr
  REAL(rp) 				:: rm
  REAL(rp) 				:: relerr
  REAL(rp), DIMENSION(:), ALLOCATABLE :: r
  REAL(rp), DIMENSION(:), ALLOCATABLE :: theta
  REAL(rp), DIMENSION(:), ALLOCATABLE :: zeta
  INTEGER 				:: pp

  ALLOCATE( theta(spp%ppp) )
  ALLOCATE( zeta(spp%ppp) )
  ALLOCATE( r(spp%ppp) )

  ! Initial condition of uniformly distributed particles on a
  ! disk in the xz-plane
  ! A unique velocity direction
  call init_u_random(10986546_8)

  call init_random_seed()
  call RANDOM_NUMBER(theta)
  theta = 2.0_rp*C_PI*theta

  call init_random_seed()
  call RANDOM_NUMBER(zeta)
  zeta = 2.0_rp*C_PI*zeta

  ! Uniform distribution on a disk at a fixed azimuthal theta
  call init_random_seed()
  call RANDOM_NUMBER(r)

  ! Newton-Raphson applied here for finding the radial distribution
  do pp=1_idef,spp%ppp 
     rl = 0.0_rp
     rr = spp%r_outter

     fl = fzero(rl,spp%r_outter,spp%falloff_rate,r(pp))
     fr = fzero(rr,spp%r_outter,spp%falloff_rate,r(pp))
     if (fl.GT.korc_zero) then
        relerr = 100*ABS(fl - fr)/fl
     else
        relerr = 100*ABS(fl - fr)/fr
     end if

     do while(relerr.GT.1.0_rp)
        rm = 0.5_rp*(rr - rl) + rl
        fm = fzero(rm,spp%r_outter,spp%falloff_rate,r(pp))

        if (SIGN(1.0_rp,fm).EQ.SIGN(1.0_rp,fr)) then
           rr = rm
        else
           rl = rm
        end if

        fl = fzero(rl,spp%r_outter,spp%falloff_rate,r(pp))
        fr = fzero(rr,spp%r_outter,spp%falloff_rate,r(pp))

        if (fl.GT.korc_zero) then
           relerr = 100*ABS(fl - fr)/fl
        else
           relerr = 100*ABS(fl - fr)/fr
        end if
     end do
     r(pp) = rm
  end do

  spp%vars%X(:,1) = ( spp%Ro + r*COS(theta) )*SIN(zeta)
  spp%vars%X(:,2) = ( spp%Ro + r*COS(theta) )*COS(zeta)
  spp%vars%X(:,3) = spp%Zo + r*SIN(theta)

  DEALLOCATE(theta)
  DEALLOCATE(zeta)
  DEALLOCATE(r)
end subroutine exponential_torus


!> @brief Subroutine that generates an exponentially decaying radial distribution in an elliptic torus as the initial spatial
!! condition of a given particle species in the simulation.
!! @details As a first step, we generate an exponentially decaying radial distribution in a circular cross-section torus as in
!! \ref korc_spatial_distribution.exponential_torus. Then we transform this spatial distribution to a one in an torus with an
!! elliptic cross section, this following the same approach as in \ref korc_spatial_distribution.elliptic_torus.
!!
!! @param[in] params Core KORC simulation parameters.
!! @param[in,out] spp An instance of the derived type SPECIES containing all the parameters and simulation variables of the different
!! species in the simulation.
!! @param fl Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param fr Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param fm Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param rl Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param rr Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param rm Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param relerr Tolerance used to determine when to stop iterating the Newton-Raphson method for finding \(r\).
!! @param rotation_angle This is the angle \(\Theta\) in \ref korc_spatial_distribution.elliptic_torus.
!! @param r Radial position of the particles \(r\).
!! @param theta Uniform deviates in the range \([0,2\pi]\) representing the uniform poloidal angle \(\theta\) distribution of the particles.
!! @param zeta Uniform deviates in the range \([0,2\pi]\) representing the uniform toroidal angle \(\zeta\) distribution of the particles.
!! @param X Auxiliary vector used in the coordinate transformations.
!! @param Y Auxiliary vector used in the coordinate transformations.
!! @param X1 Auxiliary vector used in the coordinate transformations.
!! @param Y1 Auxiliary vector used in the coordinate transformations.
!! @param pp Particle iterator.
subroutine exponential_elliptic_torus(params,spp)
  TYPE(KORC_PARAMS), INTENT(IN) 		:: params
  TYPE(SPECIES), INTENT(INOUT) 			:: spp
  REAL(rp) 					:: fl
  REAL(rp) 					:: fr
  REAL(rp) 					:: fm
  REAL(rp) 					:: rl
  REAL(rp) 					:: rr
  REAL(rp) 					:: rm
  REAL(rp) 					:: relerr
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: rotation_angle
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: r
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: theta
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: zeta
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: X
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: Y
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: X1
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: Y1
  INTEGER 				:: pp

  ALLOCATE(X1(spp%ppp))
  ALLOCATE(Y1(spp%ppp))
  ALLOCATE(X(spp%ppp))
  ALLOCATE(Y(spp%ppp))
  ALLOCATE( rotation_angle(spp%ppp) )
  ALLOCATE( theta(spp%ppp) )
  ALLOCATE( zeta(spp%ppp) )
  ALLOCATE( r(spp%ppp) )

  ! Initial condition of uniformly distributed particles on a
  ! disk in the xz-plane
  ! A unique velocity direction
  call init_u_random(10986546_8)

  call init_random_seed()
  call RANDOM_NUMBER(theta)
  theta = 2.0_rp*C_PI*theta

  call init_random_seed()
  call RANDOM_NUMBER(zeta)
  zeta = 2.0_rp*C_PI*zeta

  ! Uniform distribution on a disk at a fixed azimuthal theta
  call init_random_seed()
  call RANDOM_NUMBER(r)

  do pp=1_idef,spp%ppp
     rl = 0.0_rp
     rr = spp%r_outter

     fl = fzero(rl,spp%r_outter,spp%falloff_rate,r(pp))
     fr = fzero(rr,spp%r_outter,spp%falloff_rate,r(pp))
     if (fl.GT.korc_zero) then
        relerr = 100*ABS(fl - fr)/fl
     else
        relerr = 100*ABS(fl - fr)/fr
     end if

     do while(relerr.GT.1.0_rp)
        rm = 0.5_rp*(rr - rl) + rl
        fm = fzero(rm,spp%r_outter,spp%falloff_rate,r(pp))

        if (SIGN(1.0_rp,fm).EQ.SIGN(1.0_rp,fr)) then
           rr = rm
        else
           rl = rm
        end if

        fl = fzero(rl,spp%r_outter,spp%falloff_rate,r(pp))
        fr = fzero(rr,spp%r_outter,spp%falloff_rate,r(pp))

        if (fl.GT.korc_zero) then
           relerr = 100*ABS(fl - fr)/fl
        else
           relerr = 100*ABS(fl - fr)/fr
        end if
     end do
     r(pp) = rm
  end do

  Y = r*SIN(theta)
  X = r*COS(theta) + spp%shear_factor*Y

  rotation_angle = 0.5_rp*C_PI - ATAN(1.0_rp,1.0_rp + spp%shear_factor);

  X1 = X*COS(rotation_angle) - Y*SIN(rotation_angle) + spp%Ro
  Y1 = X*SIN(rotation_angle) + Y*COS(rotation_angle) + spp%Zo

  spp%vars%X(:,1) = X1*SIN(zeta)
  spp%vars%X(:,2) = X1*COS(zeta)
  spp%vars%X(:,3) = Y1

  DEALLOCATE(X1)
  DEALLOCATE(Y1)
  DEALLOCATE(X)
  DEALLOCATE(Y)
  DEALLOCATE(rotation_angle)
  DEALLOCATE(theta)
  DEALLOCATE(zeta)
  DEALLOCATE(r)
end subroutine exponential_elliptic_torus


!> @brief Subroutine that generates a Gaussian radial distribution in an elliptic torus as the initial spatial
!! condition of a given particle species in the simulation.
!! @details As a first step, we generate an Gaussian radial distribution in a circular cross-section torus as in
!! \ref korc_spatial_distribution.gaussian_torus. Then we transform this spatial distribution to a one in an torus with an
!! elliptic cross section, this following the same approach as in \ref korc_spatial_distribution.elliptic_torus.
!!
!! @param[in] params Core KORC simulation parameters.
!! @param[in,out] spp An instance of the derived type SPECIES containing all the parameters and simulation variables of the different
!! species in the simulation.
!! @param rotation_angle This is the angle \(\Theta\) in \ref korc_spatial_distribution.elliptic_torus.
!! @param r Radial position of the particles \(r\).
!! @param theta Uniform deviates in the range \([0,2\pi]\) representing the uniform poloidal angle \(\theta\) distribution of the particles.
!! @param zeta Uniform deviates in the range \([0,2\pi]\) representing the uniform toroidal angle \(\zeta\) distribution of the particles.
!! @param X Auxiliary vector used in the coordinate transformations.
!! @param Y Auxiliary vector used in the coordinate transformations.
!! @param X1 Auxiliary vector used in the coordinate transformations.
!! @param Y1 Auxiliary vector used in the coordinate transformations.
!! @param sigma Standard deviation \(\sigma\) of the radial distribution function.
!! @param pp Particle iterator.
subroutine gaussian_elliptic_torus(params,spp)
  TYPE(KORC_PARAMS), INTENT(IN) 	:: params
  TYPE(SPECIES), INTENT(INOUT) 		:: spp
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: rotation_angle
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: r
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: theta
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: zeta
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: X
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: Y
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: X1
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: Y1
  REAL(rp) 				:: sigma
  INTEGER 				:: pp

  ALLOCATE(X1(spp%ppp))
  ALLOCATE(Y1(spp%ppp))
  ALLOCATE(X(spp%ppp))
  ALLOCATE(Y(spp%ppp))
  ALLOCATE( rotation_angle(spp%ppp) )
  ALLOCATE( theta(spp%ppp) )
  ALLOCATE( zeta(spp%ppp) )
  ALLOCATE( r(spp%ppp) )

  ! Initial condition of uniformly distributed particles on a
  ! disk in the xz-plane
  ! A unique velocity direction
  call init_u_random(10986546_8)

  call init_random_seed()
  call RANDOM_NUMBER(theta)
  theta = 2.0_rp*C_PI*theta

  call init_random_seed()
  call RANDOM_NUMBER(zeta)
  zeta = 2.0_rp*C_PI*zeta

  ! Uniform distribution on a disk at a fixed azimuthal theta
  call init_random_seed()
  call RANDOM_NUMBER(r)

  sigma = 1.0_rp/SQRT(2.0_rp*(spp%falloff_rate/params%cpp%length))
  sigma = sigma/params%cpp%length

  r = sigma*SQRT(-2.0_rp*LOG(1.0_rp - (1.0_rp - &
       EXP(-0.5_rp*spp%r_outter**2/sigma**2))*r))
!  spp%vars%X(:,1) = ( spp%Ro + r*COS(theta) )*SIN(zeta)
!  spp%vars%X(:,2) = ( spp%Ro + r*COS(theta) )*COS(zeta)
!  spp%vars%X(:,3) = spp%Zo + r*SIN(theta)

  Y = r*SIN(theta)
  X = r*COS(theta) + spp%shear_factor*Y

  rotation_angle = 0.5_rp*C_PI - ATAN(1.0_rp,1.0_rp + spp%shear_factor);

  X1 = X*COS(rotation_angle) - Y*SIN(rotation_angle) + spp%Ro
  Y1 = X*SIN(rotation_angle) + Y*COS(rotation_angle) + spp%Zo

  spp%vars%X(:,1) = X1*SIN(zeta)
  spp%vars%X(:,2) = X1*COS(zeta)
  spp%vars%X(:,3) = Y1

  DEALLOCATE(X1)
  DEALLOCATE(Y1)
  DEALLOCATE(X)
  DEALLOCATE(Y)
  DEALLOCATE(rotation_angle)
  DEALLOCATE(theta)
  DEALLOCATE(zeta)
  DEALLOCATE(r)
end subroutine gaussian_elliptic_torus


FUNCTION PSI_ROT(R,R0,sigR,Z,Z0,sigZ,theta)
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

  PSI_ROT=(R-R0)**2*((cos(theta))**2/(2*sigR**2)+(sin(theta))**2/(2*sigZ**2))+ &
       2*(R-R0)*(Z-Z0)*cos(theta)*sin(theta)*(1/(2*sigR**2)-1/(2*sigZ**2))+ &
       (Z-Z0)**2*((sin(theta))**2/(2*sigR**2)+(cos(theta))**2/(2*sigZ**2))
  
END FUNCTION PSI_ROT


FUNCTION indicator(psi,psi_max)
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

subroutine MH_gaussian_elliptic_torus(params,spp)
  !! @note Subroutine that generates a 2D Gaussian distribution in an 
  !! elliptic torus as the initial spatial condition of a given particle 
  !! species in the simulation. @endnote
  TYPE(KORC_PARAMS), INTENT(IN) 	:: params
  !! Core KORC simulation parameters.
  TYPE(SPECIES), INTENT(INOUT) 		:: spp
  !! An instance of the derived type SPECIES containing all the parameters
  !! and simulation variables of the different species in the simulation.
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: R_samples
  !! 
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: Z_samples
  !!
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: ZETA_samples
  !! 
  REAL(rp) 				:: psi_max_buff
  !!
  REAL(rp) 				:: theta_rad
  !!
  REAL(rp) 				:: R_buffer
  !!
  REAL(rp) 				:: Z_buffer
  !!
  REAL(rp) 				:: R_test
  !!
  REAL(rp) 				:: Z_test
  !!
  REAL(rp) 				:: psi0
  !!
  REAL(rp) 				:: psi1
  !!
  REAL(rp) 				:: rand_unif
  !!
  REAL(rp) 				:: ratio
  !!
  INTEGER				:: nsamples
  !!
  INTEGER 				:: ii
  !! Particle iterator.
  INTEGER 				:: mpierr
  LOGICAL :: accepted

  nsamples = spp%ppp*params%mpi_params%nmpi

  psi_max_buff = spp%psi_max*1.1_rp

  theta_rad=C_PI*spp%theta_gauss/180.0_rp
  
  if (params%mpi_params%rank.EQ.0_idef) then
     ALLOCATE(R_samples(nsamples))
     ! Number of samples to distribute among all MPI processes
     ALLOCATE(Z_samples(nsamples))
     ! Number of samples to distribute among all MPI processes
     ALLOCATE(ZETA_samples(nsamples))
     ! Number of samples to distribute among all MPI processes

     ! Transient !

     R_buffer = spp%Ro
     Z_buffer = spp%Zo

     ii=2_idef
     do while (ii .LE. 1000_idef)
        R_test = R_buffer + random_norm(0.0_rp,spp%sigmaR)
        Z_test = Z_buffer + random_norm(0.0_rp,spp%sigmaZ)

        psi0=PSI_ROT(R_buffer,spp%Ro,spp%sigmaR,Z_buffer,spp%Zo, &
             spp%sigmaZ,theta_rad)
        psi1=PSI_ROT(R_test,spp%Ro,spp%sigmaR,Z_test,spp%Zo,spp%sigmaZ,theta_rad)

        ratio = indicator(psi1,spp%psi_max)*R_test*EXP(-psi1)/(R_buffer*EXP(-psi0))


        if (ratio .GE. 1.0_rp) then
           R_buffer = R_test
           Z_buffer = Z_test
           ii = ii + 1_idef
        else
           call RANDOM_NUMBER(rand_unif)
           if (rand_unif .LT. ratio) then
              R_buffer = R_test
              Z_buffer = Z_test
              ii = ii + 1_idef
           end if
        end if
     end do
     ! Transient !

     ii=1_idef
     do while (ii .LE. nsamples)

        R_test = R_buffer + random_norm(0.0_rp,spp%sigmaR)
        Z_test = Z_buffer + random_norm(0.0_rp,spp%sigmaZ)

        psi0=PSI_ROT(R_buffer,spp%Ro,spp%sigmaR,Z_buffer,spp%Zo, &
             spp%sigmaZ,theta_rad)
        psi1=PSI_ROT(R_test,spp%Ro,spp%sigmaR,Z_test,spp%Zo,spp%sigmaZ,theta_rad)

        ratio = indicator(psi1,psi_max_buff)*R_test*EXP(-psi1)/(R_buffer*EXP(-psi0))

        accepted=.false.
        if (ratio .GE. 1.0_rp) then
           accepted=.true.
           R_buffer = R_test
           Z_buffer = Z_test
        else
           call RANDOM_NUMBER(rand_unif)
           if (rand_unif .LT. ratio) then
              accepted=.true.
              R_buffer = R_test
              Z_buffer = Z_test
           end if
        end if

        IF (INT(indicator(psi1,spp%psi_max)).EQ.1.and.accepted) THEN
           R_samples(ii) = R_buffer
           Z_samples(ii) = Z_buffer
!           call RANDOM_NUMBER(rand_unif)
           ZETA_samples(ii) = 2.0_rp*C_PI!!*rand_unif
           ii = ii + 1_idef 
        END IF
        
     end do

  end if

  CALL MPI_SCATTER(R_samples*sin(ZETA_samples),spp%ppp,MPI_REAL8, &
       spp%vars%X(:,1),spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  CALL MPI_SCATTER(R_samples*cos(ZETA_samples),spp%ppp,MPI_REAL8, &
       spp%vars%X(:,2),spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  CALL MPI_SCATTER(Z_samples,spp%ppp,MPI_REAL8,spp%vars%X(:,3), &
       spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

  
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

  if (params%mpi_params%rank.EQ.0_idef) then
     DEALLOCATE(R_samples)
     DEALLOCATE(Z_samples)
     DEALLOCATE(ZETA_samples)
  end if

end subroutine MH_gaussian_elliptic_torus

!> @brief Subroutine that generates a Gaussian radial distribution of particles in a circular cross-section torus of
!! major and minor radi \(R_0\) and \(r_0\), respectively.
!! @details We generate this exponentially decaying radial distribution \(f(r)\) following the same approach as in
!! \ref korc_spatial_distribution.disk, but this time, the radial distribution is given by:
!!
!!
!! $$f(r) = \left\{ \begin{array}{ll} \frac{1}{4\pi^2 \sigma^2 R_0}\frac{\exp{\left( -\frac{r^2}{2\sigma^2} \right)}}{1 - \exp{\left( -\frac{r_0^2}{2\sigma^2} \right)}} &\
!! r<r_0 \\ 0 & r>r_0 \end{array} \right. $$
!!
!!
!! The radial position of the particles \(r\) is obtained using the Inverse Trasnform Sampling method, finding \(r\) numerically
!! through the Newton-Raphson method. First, we calculate the particles' radial distribution in a disk centered at \((R,Z) = (0,0)\).
!! Then, we transfor to a new set of coordinates where the disk is centered at \((R,Z) = (R_0,Z_0)\). Finally, we generate the
!! toroidal distribution by givin each particle a toroidal angle \(\zeta\) which follows a uniform distribution in the interval
!! \([0,2\pi]\).
!!
!! @param[in] params Core KORC simulation parameters.
!! @param[in,out] spp An instance of the derived type SPECIES containing all the parameters and simulation variables of the different
!! species in the simulation.
!! @param fl Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param fr Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param fm Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param rl Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param rr Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param rm Variable used in the Newton-Raphson method for finding the radial position of each particle.
!! @param relerr Tolerance used to determine when to stop iterating the Newton-Raphson method for finding \(r\).
!! @param r Radial position of the particles \(r\).
!! @param theta Uniform deviates in the range \([0,2\pi]\) representing the uniform poloidal angle \(\theta\) distribution of the particles.
!! @param zeta Uniform deviates in the range \([0,2\pi]\) representing the uniform toroidal angle \(\zeta\) distribution of the particles.
!! @param pp Particle iterator.
subroutine gaussian_torus(params,spp)
  TYPE(KORC_PARAMS), INTENT(IN) 			:: params
  TYPE(SPECIES), INTENT(INOUT) 			:: spp
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: theta
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: zeta
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: r ! temporary vars
  REAL(rp) 				:: sigma

  ALLOCATE( theta(spp%ppp) )
  ALLOCATE( zeta(spp%ppp) )
  ALLOCATE( r(spp%ppp) )

  ! Initial condition of uniformly distributed particles on a disk in the xz-plane
  ! A unique velocity direction
  call init_u_random(10986546_8)

  call init_random_seed()
  call RANDOM_NUMBER(theta)
  theta = 2.0_rp*C_PI*theta

  call init_random_seed()
  call RANDOM_NUMBER(zeta)
  zeta = 2.0_rp*C_PI*zeta

  ! Uniform distribution on a disk at a fixed azimuthal theta
  call init_random_seed()
  call RANDOM_NUMBER(r)

  sigma = 1.0_rp/SQRT(2.0_rp*(spp%falloff_rate/params%cpp%length))
  sigma = sigma/params%cpp%length

  r = sigma*SQRT(-2.0_rp*LOG(1.0_rp - (1.0_rp - &
       EXP(-0.5_rp*spp%r_outter**2/sigma**2))*r))
  spp%vars%X(:,1) = ( spp%Ro + r*COS(theta) )*SIN(zeta)
  spp%vars%X(:,2) = ( spp%Ro + r*COS(theta) )*COS(zeta)
  spp%vars%X(:,3) = spp%Zo + r*SIN(theta)

  DEALLOCATE(theta)
  DEALLOCATE(zeta)
  DEALLOCATE(r)
end subroutine gaussian_torus

function Spong_2D(R0,b,w,dlam,R,Z,T)
  REAL(rp), INTENT(IN) :: R0
  REAL(rp), INTENT(IN) :: b
  REAL(rp), INTENT(IN) :: w
  REAL(rp), INTENT(IN) :: dlam
  REAL(rp), INTENT(IN) :: R
  REAL(rp), INTENT(IN) :: Z
  REAL(rp), INTENT(IN) :: T

  Real(rp) :: rm
  Real(rp) :: lam
  
  REAL(rp) :: Spong_2D

  rm=sqrt((R-R0)**2+Z**2)
  lam=(sin(deg2rad(T)))**2
  
  Spong_2D=(1-tanh((rm-b)/w))/(1-tanh(-b/w))*exp(-(lam/dlam)**2)
  
end function Spong_2D

subroutine Spong_3D(params,spp)
  !! @note Subroutine that generates a 2D Gaussian distribution in an 
  !! elliptic torus as the initial spatial condition of a given particle 
  !! species in the simulation. @endnote
  TYPE(KORC_PARAMS), INTENT(IN) 	:: params
  !! Core KORC simulation parameters.
  TYPE(SPECIES), INTENT(INOUT) 		:: spp
  !! An instance of the derived type SPECIES containing all the parameters
  !! and simulation variables of the different species in the simulation.

  !! An instance of the KORC derived type FIELDS.
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: R_samples
  !! Major radial location of all samples
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: PHI_samples
  !! Azimuithal angle of all samples
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: Z_samples
  !! Vertical location of all samples

  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: T_samples
  !! Pitch angle of all samples

  REAL(rp) 				:: psi_max_buff
  !! Value of buffer above desired maximum argument of 2D Gaussian spatial
  !! profile
  REAL(rp) 				:: minmax
  !! Temporary variable used for setting buffers

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
 
  REAL(rp) 				:: T_buffer
  !! Previous sample of pitch angle
  REAL(rp) 				:: R_test
  !! Present sample of R location
  REAL(rp) 				:: Z_test
  !! Present sample of Z location

  REAL(rp) 				:: T_test
  !! Present sample of pitch angle
  REAL(rp) 				:: psi0
  !! Previous value of 2D Gaussian argument based on R_buffer, Z_buffer
  REAL(rp) 				:: psi1
  !! Present value of 2D Gaussian argument based on R_test, Z_test
  REAL(rp) 				:: f0
  !! Evaluation of Avalanche distribution with previous sample
  REAL(rp) 				:: f1
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
  
  
  nsamples = spp%ppp*params%mpi_params%nmpi

  psi_max_buff = spp%psi_max*1.1_rp

  theta_rad=C_PI*spp%theta_gauss/180.0_rp


  ! buffer at minimum pitch angle boundary  
  if (spp%etao_lims(1).GE.korc_zero) then
     do ii=1_idef,INT(minmax_buffer_size,idef)
        minmax = spp%etao_lims(1) - REAL(ii,rp)* &
             (spp%etao_lims(2)-spp%etao_lims(1))/100_rp
        if (minmax.GT.0.0_rp) then
           min_pitch_angle = minmax
        end if
     end do
  else
     min_pitch_angle = spp%etao_lims(1)
  end if

  ! buffer at maximum pitch angle boundary
  do ii=1_idef,INT(minmax_buffer_size,idef)
     minmax = spp%etao_lims(2) + REAL(ii,rp)* &
          (spp%etao_lims(2)-spp%etao_lims(1))/100_rp
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
     ALLOCATE(T_samples(nsamples))
     ! Number of samples to distribute among all MPI processes
     
     ! Transient !

     R_buffer = spp%Ro
     Z_buffer = spp%Zo

     
     call RANDOM_NUMBER(rand_unif)
     T_buffer = min_pitch_angle + (max_pitch_angle  &
          - min_pitch_angle)*rand_unif

!     write(output_unit_write,'("length norm: ",E17.10)') params%cpp%length
     
     ii=1_idef
     do while (ii .LE. 1000_idef)

!        write(output_unit_write,'("burn:",I15)') ii
        
        R_test = R_buffer + random_norm(0.0_rp,spp%dR)
        Z_test = Z_buffer + random_norm(0.0_rp,spp%dZ)
        T_test = T_buffer + random_norm(0.0_rp,spp%dth)


        ! Test that pitch angle and momentum are within chosen boundary
        do while ((T_test .GT. spp%etao_lims(2)).OR. &
             (T_test .LT. spp%etao_lims(1)))
           T_test = T_buffer + random_norm(0.0_rp,spp%dth)
        end do

        ! initialize 2D gaussian argument and distribution function, or
        ! copy from previous sample
        if (ii==1) then
           psi0=PSI_ROT(R_buffer,spp%Ro,spp%sigmaR,Z_buffer,spp%Zo, &
                spp%sigmaZ,theta_rad)
           
           f0=Spong_2D(spp%Ro,spp%Spong_b,spp%Spong_w,spp%Spong_dlam, &
                R_buffer,Z_buffer,T_buffer)           
        else
           psi0=psi1
           f0=f1
        end if
        
        psi1=PSI_ROT(R_test,spp%Ro,spp%sigmaR,Z_test,spp%Zo, &
             spp%sigmaZ,theta_rad)       
        
                
        f1=Spong_2D(spp%Ro,spp%Spong_b,spp%Spong_w,spp%Spong_dlam, &
             R_test,Z_test,T_test)    

!        write(output_unit_write,'("psi0: ",E17.10)') psi0
!        write(output_unit_write,'("psi1: ",E17.10)') psi1

!        write(output_unit_write,'("f0: ",E17.10)') f0
!        write(output_unit_write,'("f1: ",E17.10)') f1

        
        ! Calculate acceptance ratio for MH algorithm. fRE function
        ! incorporates p^2 factor of spherical coordinate Jacobian
        ! for velocity phase space, factors of sin(pitch angle) for velocity
        ! phase space and cylindrical coordinate Jacobian R for spatial
        ! phase space incorporated here.
        ratio = indicator(psi1,spp%psi_max)*R_test*EXP(-psi1)*f1/ &
             (R_buffer*EXP(-psi0)*f0)
        
        if (ratio .GE. 1.0_rp) then
           R_buffer = R_test
           Z_buffer = Z_test
           T_buffer = T_test
           ii = ii + 1_idef
        else
           call RANDOM_NUMBER(rand_unif)
           if (rand_unif .LT. ratio) then
              R_buffer = R_test
              Z_buffer = Z_test
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
        
        R_test = R_buffer + random_norm(0.0_rp,spp%dR)
        Z_test = Z_buffer + random_norm(0.0_rp,spp%dZ)
        T_test = T_buffer + random_norm(0.0_rp,spp%dth)

        ! Selection boundary is set with buffer region
        do while ((T_test .GT. max_pitch_angle).OR. &
             (T_test .LT. min_pitch_angle))
           if (T_test.lt.0) then
              T_test=abs(T_test)
              exit
           end if
           T_test = T_buffer + random_norm(0.0_rp,spp%dth)
        end do


        psi0=psi1
        f0=f1
        
        psi1=PSI_ROT(R_test,spp%Ro,spp%sigmaR,Z_test,spp%Zo, &
             spp%sigmaZ,theta_rad)

        f1=Spong_2D(spp%Ro,spp%Spong_b,spp%Spong_w,spp%Spong_dlam, &
             R_test,Z_test,T_test)
        
        ratio = indicator(psi1,psi_max_buff)*R_test*EXP(-psi1)*f1/ &
             (R_buffer*EXP(-psi0)*f0)

        if (ratio .GE. 1.0_rp) then
           R_buffer = R_test
           Z_buffer = Z_test
           T_buffer = T_test
        else
           call RANDOM_NUMBER(rand_unif)
           if (rand_unif .LT. ratio) then
              R_buffer = R_test
              Z_buffer = Z_test
              T_buffer = T_test
           end if
        end if
        
        ! Only accept sample if it is within desired boundary, but
        ! add to MC above if within buffer. This helps make the boundary
        ! more defined.
        IF ((INT(indicator(psi1,spp%psi_max)).EQ.1).AND. &
             (T_buffer.LE.spp%etao_lims(2)).AND. &
             (T_buffer.GE.spp%etao_lims(1))) THEN
           R_samples(ii) = R_buffer
           Z_samples(ii) = Z_buffer
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
!  CALL MPI_SCATTER(T_samples,spp%ppp,MPI_REAL8, &
!       spp%vars%eta,spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  
  
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

!  write(output_unit_write,'("X_X: ",E17.10)') spp%vars%X(:,1)*params%cpp%length
  
  ! gamma is kept for each particle, not the momentum

  if (params%orbit_model(1:2).eq.'GC') call cart_to_cyl(spp%vars%X,spp%vars%Y)

!  write(output_unit_write,'("Y_R: ",E17.10)') spp%vars%Y(:,1)*params%cpp%length
  
!  if (minval(spp%vars%Y(:,1)).lt.1._rp/params%cpp%length) stop 'error with avalanche'
  
  if (params%mpi_params%rank.EQ.0_idef) then
     DEALLOCATE(R_samples)
     DEALLOCATE(Z_samples)
     DEALLOCATE(PHI_samples)
     DEALLOCATE(T_samples)
  end if
  
  
end subroutine Spong_3D

subroutine MH_psi(params,spp,F)
  !! @note Subroutine that generates a 2D Gaussian distribution in an 
  !! elliptic torus as the initial spatial condition of a given particle 
  !! species in the simulation. @endnote
  TYPE(KORC_PARAMS), INTENT(INOUT) 	:: params
  !! Core KORC simulation parameters.
  TYPE(SPECIES), INTENT(INOUT) 		:: spp
  !! An instance of the derived type SPECIES containing all the parameters
  !! and simulation variables of the different species in the simulation.
  TYPE(FIELDS), INTENT(IN)                                   :: F
  !! An instance of the KORC derived type FIELDS.
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: R_samples,X_samples,Y_samples
  !! Major radial location of all samples
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: PHI_samples
  !! Azimuithal angle of all samples
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: Z_samples
  !! Vertical location of all samples

  REAL(rp) 				:: min_R,max_R
  REAL(rp) 				:: min_Z,max_Z
 
  REAL(rp) 				:: R_buffer
  !! Previous sample of R location
  REAL(rp) 				:: Z_buffer
  !! Previous sample of Z location
  REAL(rp) 				:: PHI_buffer
 

  REAL(rp) 				:: R_test
  !! Present sample of R location
  REAL(rp) 				:: Z_test
  !! Present sample of Z location
  REAL(rp) 				:: PHI_test

  REAL(rp) 				:: psi_max,psi_max_buff
  REAL(rp)  :: PSIp_lim,PSIP0,PSIN,PSIN0,PSIN1,sigma,psi0,psi1

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

  LOGICAL :: accepted
  INTEGER,DIMENSION(33) :: seed=(/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)

  if (params%mpi_params%rank.EQ.0_idef) then
     write(output_unit_write,*) '*** START SAMPLING ***'
  end if
  
  nsamples = spp%ppp*params%mpi_params%nmpi

  params%GC_coords=.TRUE.
  PSIp_lim=F%PSIp_lim

  if (params%field_model.eq.'M3D_C1') then
     min_R=params%rmin/params%cpp%length
     max_R=params%rmax/params%cpp%length
     min_Z=params%Zmin/params%cpp%length
     max_Z=params%Zmax/params%cpp%length

     PSIp0=F%PSIp_0
     psi_max = spp%psi_max
     psi_max_buff = spp%psi_max
  else
     min_R=minval(F%X%R)
     max_R=maxval(F%X%R)
     min_Z=minval(F%X%Z)
     max_Z=maxval(F%X%Z)

     PSIp0=F%PSIP_min
     psi_max = spp%psi_max
     psi_max_buff = spp%psi_max*2._rp
  end if  

  sigma=spp%sigmaR*params%cpp%length
  
  !write(output_unit_write,*) min_R,max_R
  !write(output_unit_write,*) min_Z,max_Z
 
  

  ALLOCATE(R_samples(nsamples))
  ALLOCATE(X_samples(nsamples))
  ALLOCATE(Y_samples(nsamples))
  ! Number of samples to distribute among all MPI processes
  ALLOCATE(PHI_samples(nsamples))
  ! Number of samples to distribute among all MPI processes
  ALLOCATE(Z_samples(nsamples))
  ! Number of samples to distribute among all MPI processes

  if (params%mpi_params%rank.EQ.0_idef) then
        ! Transient !

     
     R_buffer = spp%Ro
     Z_buffer = spp%Zo
     PHI_buffer = 0._rp


     if (.not.params%SameRandSeed) then
        call init_random_seed()
     else
        call random_seed(put=seed)
     end if

     write(output_unit_write,'("Begin burn: ",I10)')
     accepted=.false.
     ii=1_idef
     do while (ii .LE. 1000_idef)

        if (modulo(ii,100).eq.0) then
           write(output_unit_write,'("Burn: ",I10)') ii
        end if
        !write(6,'("Burn: ",I10)') ii
        
        !R_test = R_buffer + random_norm(0.0_rp,spp%dR)
        !R_test = R_buffer + get_random_mkl_N(0.0_rp,spp%dR)
        R_test = R_buffer + get_random_N()*spp%dR
        
        !Z_test = Z_buffer + random_norm(0.0_rp,spp%dZ)
        !Z_test = Z_buffer + get_random_mkl_N(0.0_rp,spp%dZ)
        Z_test = Z_buffer + get_random_N()*spp%dZ
        




        do while ((R_test.GT.max_R).OR.(R_test .LT. min_R))
           !eta_test = eta_buffer + random_norm(0.0_rp,spp%dth)
           !eta_test = eta_buffer + get_random_mkl_N(0.0_rp,spp%dth)
           R_test = R_buffer + get_random_N()*spp%dR
        end do

        do while ((Z_test.GT.max_Z).OR.(Z_test .LT. min_Z))
           !eta_test = eta_buffer + random_norm(0.0_rp,spp%dth)
           !eta_test = eta_buffer + get_random_mkl_N(0.0_rp,spp%dth)
           Z_test = Z_buffer + get_random_N()*spp%dZ
        end do
        PHI_test = 2.0_rp*C_PI*get_random_U()

        
        ! initialize 2D gaussian argument and distribution function, or
        ! copy from previous sample
        if (ii==1) then

           spp%vars%Y(1,1)=R_buffer
           spp%vars%Y(1,2)=PHI_buffer
           spp%vars%Y(1,3)=Z_buffer

           !write(6,*) 'R',R_buffer
           !write(6,*) 'Z',Z_buffer

           spp%vars%flagCon(1)=1_is
           spp%vars%initLCFS(1)=1_is
           if (params%field_model.eq.'M3D_C1'.or.params%field_model.eq.'NIMROD') then
#ifdef FIO
              call get_fio_vector_potential(spp%vars,F,params)
#endif
           else
              call get_fields(params,spp%vars,F)
           end if

           !write(6,*) 'may have crashed'

           !write(6,*) 'R',R_buffer
           !write(6,*) 'Z',Z_buffer
           !write(6,*) 'PSIlim',PSIp_lim
           !write(6,*) 'PSI0',PSIp0
           !write(output_unit_write,*) 'PSI1',psi1
           !write(6,*) 'PSI0',psi0
           !write(output_unit_write,*) 'PSIN1',PSIN1
           !write(6,*) 'PSIN0',PSIN0
           
           psi0=spp%vars%PSI_P(1)
           PSIN0=(psi0-PSIp0)/(PSIp_lim-PSIp0)
           
        end if
        
        if (accepted) then
           PSIN0=PSIN1
        end if
        
!        psi1=PSI_ROT_exp(R_test,spp%Ro,spp%sigmaR,Z_test,spp%Zo, &
        !             spp%sigmaZ,theta_rad)
        spp%vars%Y(1,1)=R_test
        spp%vars%Y(1,2)=PHI_test
        spp%vars%Y(1,3)=Z_test

        spp%vars%flagCon(1)=1_is
        spp%vars%initLCFS(1)=1_is
        if (params%field_model.eq.'M3D_C1'.or.params%field_model.eq.'NIMROD') then
#ifdef FIO
           call get_fio_vector_potential(spp%vars,F,params)
#endif
        else
           call get_fields(params,spp%vars,F)
        end if
        
        psi1=spp%vars%PSI_P(1)

        !write(output_unit_write,*) 'PSIlim',PSIp_lim
        !write(output_unit_write,*) 'PSI0',PSIp0
        !write(output_unit_write,*) 'PSI',psi1

        PSIN1=(psi1-PSIp0)/(PSIp_lim-PSIp0)

        !write(output_unit_write,*) 'R',R_test
        !write(output_unit_write,*) 'Z',Z_test
        !write(output_unit_write,*) 'PSIlim',PSIp_lim
        !write(output_unit_write,*) 'PSI0',PSIp0
        !write(output_unit_write,*) 'PSI1',psi1
        !write(output_unit_write,*) 'PSI0',psi0
        !write(output_unit_write,*) 'PSIN',PSIN1
        !write(output_unit_write,*) 'PSIN0',PSIN0
        
        ! Calculate acceptance ratio for MH algorithm. fRE function
        ! incorporates p^2 factor of spherical coordinate Jacobian
        ! for velocity phase space, factors of sin(pitch angle) for velocity
        ! phase space and cylindrical coordinate Jacobian R for spatial
        ! phase space incorporated here.

        if (.not.F%useLCFS) spp%vars%initLCFS(1)=1_is
        ratio = real(spp%vars%flagCon(1))*real(spp%vars%initLCFS(1))* &
             indicator(PSIN1,psi_max)* &
             R_test*EXP(-PSIN1/sigma)/ &
             (R_buffer*EXP(-PSIN0/sigma))

!        ratio = f1*sin(deg2rad(eta_test))/(f0*sin(deg2rad(eta_buffer)))

        accepted=.false.
        if (ratio .GE. 1.0_rp) then
           accepted=.true.
           R_buffer = R_test
           Z_buffer = Z_test
           PHI_buffer = PHI_test
           ii = ii + 1_idef

           !write(output_unit_write,*) 'PSIN',PSIN1
        else
!           call RANDOM_NUMBER(rand_unif)
!           if (rand_unif .LT. ratio) then
           !if (get_random_mkl_U() .LT. ratio) then
           if (get_random_U() .LT. ratio) then
              accepted=.true.
              R_buffer = R_test
              Z_buffer = Z_test
              PHI_buffer = PHI_test
              ii = ii + 1_idef

              !write(output_unit_write,*) 'PSIN',PSIN1
           end if
        end if
     end do
     ! Transient !

     write(output_unit_write,'("Begin sample: ",I10)')
     ii=1_idef
     do while (ii .LE. nsamples)

!        write(output_unit_write,'("sample:",I15)') ii

#if DBG_CHECK        
#else
        if (modulo(ii,nsamples/10).eq.0) then
           write(output_unit_write,'("Sample: ",I10)') ii
        end if
#endif
        
        !R_test = R_buffer + random_norm(0.0_rp,spp%dR)
        !R_test = R_buffer + get_random_mkl_N(0.0_rp,spp%dR)
        R_test = R_buffer + get_random_N()*spp%dR
        
        !Z_test = Z_buffer + random_norm(0.0_rp,spp%dZ)
        !Z_test = Z_buffer + get_random_mkl_N(0.0_rp,spp%dZ)
        Z_test = Z_buffer + get_random_N()*spp%dZ

        PHI_test = 2.0_rp*C_PI*get_random_U()

        do while ((R_test.GT.max_R).OR.(R_test .LT. min_R))
           !eta_test = eta_buffer + random_norm(0.0_rp,spp%dth)
           !eta_test = eta_buffer + get_random_mkl_N(0.0_rp,spp%dth)
           R_test = R_buffer + get_random_N()*spp%dR
        end do

        do while ((Z_test.GT.max_Z).OR.(Z_test .LT. min_Z))
           !eta_test = eta_buffer + random_norm(0.0_rp,spp%dth)
           !eta_test = eta_buffer + get_random_mkl_N(0.0_rp,spp%dth)
           Z_test = Z_buffer + get_random_N()*spp%dZ
        end do
        
        if (accepted) then
           PSIN0=PSIN1
        end if
        
!        psi1=PSI_ROT_exp(R_test,spp%Ro,spp%sigmaR,Z_test,spp%Zo, &
!             spp%sigmaZ,theta_rad)
        spp%vars%Y(1,1)=R_test
        spp%vars%Y(1,2)=PHI_test
        spp%vars%Y(1,3)=Z_test

        spp%vars%flagCon(1)=1_is
        spp%vars%initLCFS(1)=1_is
        if (params%field_model.eq.'M3D_C1'.or.params%field_model.eq.'NIMROD') then
#ifdef FIO
           call get_fio_vector_potential(spp%vars,F,params)
#endif
        else
           call get_fields(params,spp%vars,F)
        end if

        
        psi1=spp%vars%PSI_P(1)

        !write(output_unit_write,*) 'PSIlim',PSIp_lim
        !write(output_unit_write,*) 'PSI0',PSIp0
        !write(output_unit_write,*) 'PSI',psi1
        
        PSIN1=(psi1-PSIp0)/(PSIp_lim-PSIp0)


        !write(output_unit_write,*) 'R',R_test
        !write(output_unit_write,*) 'Z',Z_test
        !write(output_unit_write,*) 'PSIlim',PSIp_lim
        !write(output_unit_write,*) 'PSI0',PSIp0
        !write(output_unit_write,*) 'PSI1',psi1
        !write(output_unit_write,*) 'PSI0',psi0
        !write(output_unit_write,*) 'PSIN',PSIN1
        !write(output_unit_write,*) 'PSIN0',PSIN0


        if (.not.F%useLCFS) spp%vars%initLCFS(1)=1_is
        ratio = real(spp%vars%flagCon(1))*real(spp%vars%initLCFS(1))* &
             indicator(PSIN1,psi_max)* &
             R_test*EXP(-PSIN1/sigma)/ &
             (R_buffer*EXP(-PSIN0/sigma))

!        ratio = f1*sin(deg2rad(eta_test))/(f0*sin(deg2rad(eta_buffer)))
        
        accepted=.false.
        if (ratio .GE. 1.0_rp) then
           accepted=.true.
           R_buffer = R_test
           Z_buffer = Z_test
           PHI_buffer = PHI_test                     
        else
           !call RANDOM_NUMBER(rand_unif)
           !if (rand_unif .LT. ratio) then
           !if (get_random_mkl_U() .LT. ratio) then
           if (get_random_U() .LT. ratio) then
              accepted=.true.
              R_buffer = R_test
              Z_buffer = Z_test
              PHI_buffer = PHI_test
           end if
        end if

!        write(output_unit_write,'("R: ",E17.10)') R_buffer
!        write(output_unit_write,'("Z: ",E17.10)') Z_buffer
        
        ! Only accept sample if it is within desired boundary, but
        ! add to MC above if within buffer. This helps make the boundary
        ! more defined.
        IF ((INT(indicator(PSIN1,psi_max)).EQ.1).AND. &
             ACCEPTED) THEN
           R_samples(ii) = R_buffer
           Z_samples(ii) = Z_buffer
           PHI_samples(ii) = PHI_buffer


           !write(output_unit_write,*) 'PSIN',PSIN1

           
!           write(output_unit_write,*) 'RS',R_buffer
           
           ! Sample phi location uniformly
           !call RANDOM_NUMBER(rand_unif)
           !PHI_samples(ii) = 2.0_rp*C_PI*rand_unif
           !PHI_samples(ii) = 2.0_rp*C_PI*get_random_mkl_U()
           ii = ii + 1_idef 
        END IF
        
     end do

!  if (minval(R_samples(:)).lt.1._rp/params%cpp%length) stop 'error with sample'
!  write(output_unit_write,'("R_sample: ",E17.10)') R_samples(:)*params%cpp%length

     X_samples=R_samples*cos(PHI_samples)
     Y_samples=R_samples*sin(PHI_samples)

!     write(output_unit_write,*) 'R_samples',R_samples
!     write(output_unit_write,*) 'PHI_samples',PHI_samples
!     write(output_unit_write,*) 'Z_samples',Z_samples
!     write(output_unit_write,*) 'G_samples',G_samples
!     write(output_unit_write,*) 'eta_samples',eta_samples
     
  end if

  params%GC_coords=.FALSE.
  
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)


  CALL MPI_SCATTER(X_samples,spp%ppp,MPI_REAL8, &
       spp%vars%X(:,1),spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  CALL MPI_SCATTER(Y_samples,spp%ppp,MPI_REAL8, &
       spp%vars%X(:,2),spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  CALL MPI_SCATTER(Z_samples,spp%ppp,MPI_REAL8, &
       spp%vars%X(:,3),spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

  
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

  if (params%orbit_model(1:2).eq.'GC') call cart_to_cyl(spp%vars%X,spp%vars%Y)

  
  
  DEALLOCATE(R_samples)
  DEALLOCATE(X_samples)
  DEALLOCATE(Y_samples)
  DEALLOCATE(Z_samples)
  DEALLOCATE(PHI_samples)
  
  
end subroutine MH_psi

subroutine FIO_therm(params,spp,F,P)
  !! @note Subroutine that generates a 2D Gaussian distribution in an 
  !! elliptic torus as the initial spatial condition of a given particle 
  !! species in the simulation. @endnote
  TYPE(KORC_PARAMS), INTENT(INOUT) 	:: params
  !! Core KORC simulation parameters.
  TYPE(SPECIES), INTENT(INOUT) 		:: spp
  !! An instance of the derived type SPECIES containing all the parameters
  !! and simulation variables of the different species in the simulation.
  TYPE(FIELDS), INTENT(IN)                                   :: F
  !! An instance of the KORC derived type FIELDS.
  TYPE(PROFILES), INTENT(IN)            :: P
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: R_samples,X_samples,Y_samples
  !! Major radial location of all samples
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: PHI_samples
  !! Azimuithal angle of all samples
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: Z_samples
  !! Vertical location of all samples
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: V_samples,G_samples
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: XI_samples,ETA_samples

  REAL(rp) 				:: min_R,max_R
  REAL(rp) 				:: min_Z,max_Z
    REAL(rp) 				:: min_V,max_V
 
  REAL(rp) 				:: R_buffer
  !! Previous sample of R location
  REAL(rp) 				:: Z_buffer
  !! Previous sample of Z location
  REAL(rp) 				:: PHI_buffer
  REAL(rp) 				:: V_buffer,XI_buffer
  REAL(rp) 				:: V_test,XI_test

  REAL(rp) 				:: R_test
  !! Present sample of R location
  REAL(rp) 				:: Z_test
  !! Present sample of Z location
  REAL(rp) 				:: PHI_test

  REAL(rp) 				:: psi_max,psi_max_buff
  REAL(rp)  :: PSIp_lim,PSIP0,PSIN,PSIN0,PSIN1,sigma,psi0,psi1
  REAL(rp)  :: vth0,vth1

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

  LOGICAL :: accepted
  INTEGER,DIMENSION(33) :: seed=(/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)

  if (params%mpi_params%rank.EQ.0_idef) then
     write(output_unit_write,*) '*** START SAMPLING ***'
  end if
  
  nsamples = spp%ppp*params%mpi_params%nmpi

  params%GC_coords=.TRUE.
  PSIp_lim=F%PSIp_lim

  if (params%field_model.eq.'M3D_C1') then
     min_R=params%rmin/params%cpp%length
     max_R=params%rmax/params%cpp%length
     min_Z=params%Zmin/params%cpp%length
     max_Z=params%Zmax/params%cpp%length

     PSIp0=F%PSIp_0
     psi_max = spp%psi_max
     psi_max_buff = spp%psi_max
  else
     min_R=minval(F%X%R)
     max_R=maxval(F%X%R)
     min_Z=minval(F%X%Z)
     max_Z=maxval(F%X%Z)

     PSIp0=F%PSIP_min
     psi_max = spp%psi_max
     psi_max_buff = spp%psi_max*2._rp
  end if  

  sigma=spp%sigmaR*params%cpp%length
  
  !write(output_unit_write,*) min_R,max_R
  !write(output_unit_write,*) min_Z,max_Z
 
  

  ALLOCATE(R_samples(nsamples))
  ALLOCATE(X_samples(nsamples))
  ALLOCATE(Y_samples(nsamples))
  ! Number of samples to distribute among all MPI processes
  ALLOCATE(PHI_samples(nsamples))
  ! Number of samples to distribute among all MPI processes
  ALLOCATE(Z_samples(nsamples))
  ! Number of samples to distribute among all MPI processes
  ALLOCATE(V_samples(nsamples))
  ALLOCATE(G_samples(nsamples))
  ALLOCATE(XI_samples(nsamples))
  ALLOCATE(ETA_samples(nsamples))
  
  if (params%mpi_params%rank.EQ.0_idef) then
        ! Transient !

     
     R_buffer = spp%Ro
     Z_buffer = spp%Zo
     PHI_buffer = 0._rp

     
     spp%vars%Y(1,1)=R_buffer
     spp%vars%Y(1,2)=PHI_buffer
     spp%vars%Y(1,3)=Z_buffer

     !write(6,*) 'R',R_buffer
     !write(6,*) 'Z',Z_buffer

     spp%vars%flagCon(1)=1_is
     if (params%field_model.eq.'M3D_C1'.or.params%field_model.eq.'NIMROD') then
#ifdef FIO
        call get_fio_vector_potential(spp%vars,F,params)
        call get_fio_profile(spp%vars,P,params)
#endif
     end if

     !write(6,*) 'may have crashed'


     psi0=spp%vars%PSI_P(1)
     PSIN0=(psi0-PSIp0)/(PSIp_lim-PSIp0)

     vth0=sqrt(2*spp%vars%te(1))

     V_buffer = vth0
     min_V=6.8e-2*V_buffer
     max_V=V_buffer*5._rp
     XI_buffer=0._rp

     !write(6,*) 'R',R_buffer
     !write(6,*) 'Z',Z_buffer
     !write(6,*) 'PSIlim',PSIp_lim
     !write(6,*) 'PSI0',PSIp0
     !write(output_unit_write,*) 'PSI1',psi1
     !write(6,*) 'PSI0',psi0
     !write(output_unit_write,*) 'PSIN1',PSIN1
     !write(6,*) 'PSIN0',PSIN0
     
     if (.not.params%SameRandSeed) then
        call init_random_seed()
     else
        call random_seed(put=seed)
     end if
     
     write(output_unit_write,'("Begin burn: ",I10)')
     accepted=.false.
     ii=1_idef
     do while (ii .LE. 1000_idef)

        if (modulo(ii,100).eq.0) then
           write(output_unit_write,'("Burn: ",I10)') ii
        end if

        
        !R_test = R_buffer + random_norm(0.0_rp,spp%dR)
        !R_test = R_buffer + get_random_mkl_N(0.0_rp,spp%dR)
        R_test = R_buffer + get_random_N()*spp%dR
        
        !Z_test = Z_buffer + random_norm(0.0_rp,spp%dZ)
        !Z_test = Z_buffer + get_random_mkl_N(0.0_rp,spp%dZ)
        Z_test = Z_buffer + get_random_N()*spp%dZ
        
        V_test=V_buffer+get_random_N()*spp%dgam


        do while ((R_test.GT.max_R).OR.(R_test .LT. min_R))
           !eta_test = eta_buffer + random_norm(0.0_rp,spp%dth)
           !eta_test = eta_buffer + get_random_mkl_N(0.0_rp,spp%dth)
           R_test = R_buffer + get_random_N()*spp%dR
        end do

        do while ((Z_test.GT.max_Z).OR.(Z_test .LT. min_Z))
           !eta_test = eta_buffer + random_norm(0.0_rp,spp%dth)
           !eta_test = eta_buffer + get_random_mkl_N(0.0_rp,spp%dth)
           Z_test = Z_buffer + get_random_N()*spp%dZ
        end do

        do while ((V_test.GT.max_V).OR.(V_test .LT. min_V))
           !eta_test = eta_buffer + random_norm(0.0_rp,spp%dth)
           !eta_test = eta_buffer + get_random_mkl_N(0.0_rp,spp%dth)
           V_test = V_buffer + get_random_N()*spp%dgam
        end do
        
        PHI_test = 2.0_rp*C_PI*get_random_U()

        XI_test = -1+2.0_rp*get_random_U()

               
        if (accepted) then
           PSIN0=PSIN1
           vth0=vth1
        end if
        
        spp%vars%Y(1,1)=R_test
        spp%vars%Y(1,2)=PHI_test
        spp%vars%Y(1,3)=Z_test

        spp%vars%flagCon(1)=1_is
        if (params%field_model.eq.'M3D_C1'.or.params%field_model.eq.'NIMROD') then
#ifdef FIO
           call get_fio_vector_potential(spp%vars,F,params)
           call get_fio_profile(spp%vars,P,params)
#endif
        end if
        
        psi1=spp%vars%PSI_P(1)

        !write(output_unit_write,*) 'PSIlim',PSIp_lim
        !write(output_unit_write,*) 'PSI0',PSIp0
        !write(output_unit_write,*) 'PSI',psi1

        PSIN1=(psi1-PSIp0)/(PSIp_lim-PSIp0)

        !write(6,*) 'R',R_test*params%cpp%length
        !write(6,*) 'Z',Z_test*params%cpp%length
        !write(6,*) 'PSIlim',PSIp_lim
        !write(6,*) 'PSI0',PSIp0
        !write(6,*) 'PSI1',psi1
        !write(6,*) 'PSI0',psi0
        !write(6,*) 'PSIN',PSIN1
        !write(6,*) 'PSIN0',PSIN0

        !write(6,*) 'te',spp%vars%te(1)
        !write(6,*) 'vth0,v_buffer',vth0,V_buffer
        
        vth1=sqrt(2*spp%vars%te(1))


        !write(6,*) 'vth1,v_test',vth1,V_test
        
        ! Calculate acceptance ratio for MH algorithm. fRE function
        ! incorporates p^2 factor of spherical coordinate Jacobian
        ! for velocity phase space, factors of sin(pitch angle) for velocity
        ! phase space and cylindrical coordinate Jacobian R for spatial
        ! phase space incorporated here.

        
        ratio = real(spp%vars%flagCon(1))*indicator(PSIN1,psi_max)* &
             R_test*EXP(-PSIN1/sigma)*V_test**2/vth1**3* &
             EXP(-(V_test/vth1)**2)/ &
             (R_buffer*EXP(-PSIN0/sigma)*V_buffer**2/vth0**3* &
             EXP(-(V_buffer/vth0)**2))

        accepted=.false.
        if (ratio .GE. 1.0_rp) then
           accepted=.true.
           R_buffer = R_test
           Z_buffer = Z_test
           PHI_buffer = PHI_test
           V_buffer = V_test
           XI_buffer = XI_test
           ii = ii + 1_idef

           !write(output_unit_write,*) 'PSIN',PSIN1
        else
!           call RANDOM_NUMBER(rand_unif)
!           if (rand_unif .LT. ratio) then
           !if (get_random_mkl_U() .LT. ratio) then
           if (get_random_U() .LT. ratio) then
              accepted=.true.
              R_buffer = R_test
              Z_buffer = Z_test
              PHI_buffer = PHI_test
              V_buffer = V_test
              XI_buffer = XI_test
              ii = ii + 1_idef

              !write(output_unit_write,*) 'PSIN',PSIN1
           end if
        end if
     end do
     ! Transient !

     write(output_unit_write,'("Begin sample: ",I10)')
     ii=1_idef
     do while (ii .LE. nsamples)

!        write(output_unit_write,'("sample:",I15)') ii
        
       if (modulo(ii,nsamples/10).eq.0) then
           write(output_unit_write,'("Sample: ",I10)') ii
        end if
        
        !R_test = R_buffer + random_norm(0.0_rp,spp%dR)
        !R_test = R_buffer + get_random_mkl_N(0.0_rp,spp%dR)
        R_test = R_buffer + get_random_N()*spp%dR
        
        !Z_test = Z_buffer + random_norm(0.0_rp,spp%dZ)
        !Z_test = Z_buffer + get_random_mkl_N(0.0_rp,spp%dZ)
        Z_test = Z_buffer + get_random_N()*spp%dZ

        PHI_test = 2.0_rp*C_PI*get_random_U()

        do while ((R_test.GT.max_R).OR.(R_test .LT. min_R))
           !eta_test = eta_buffer + random_norm(0.0_rp,spp%dth)
           !eta_test = eta_buffer + get_random_mkl_N(0.0_rp,spp%dth)
           R_test = R_buffer + get_random_N()*spp%dR
        end do

        do while ((Z_test.GT.max_Z).OR.(Z_test .LT. min_Z))
           !eta_test = eta_buffer + random_norm(0.0_rp,spp%dth)
           !eta_test = eta_buffer + get_random_mkl_N(0.0_rp,spp%dth)
           Z_test = Z_buffer + get_random_N()*spp%dZ
        end do

        V_test=V_buffer+get_random_N()*spp%dgam

        do while ((V_test.GT.max_V).OR.(V_test .LT. min_V))
           !eta_test = eta_buffer + random_norm(0.0_rp,spp%dth)
           !eta_test = eta_buffer + get_random_mkl_N(0.0_rp,spp%dth)
           V_test = V_buffer + get_random_N()*spp%dgam
        end do

        XI_test = -1+2.0_rp*get_random_U()
        
        if (accepted) then
           PSIN0=PSIN1
           vth0=vth1
        end if
        
        spp%vars%Y(1,1)=R_test
        spp%vars%Y(1,2)=PHI_test
        spp%vars%Y(1,3)=Z_test

        spp%vars%flagCon(1)=1_is
        if (params%field_model.eq.'M3D_C1'.or.params%field_model.eq.'NIMROD') then
#ifdef FIO
           call get_fio_vector_potential(spp%vars,F,params)
           call get_fio_profile(spp%vars,P,params)
#endif
        end if

        
        psi1=spp%vars%PSI_P(1)

        !write(output_unit_write,*) 'PSIlim',PSIp_lim
        !write(output_unit_write,*) 'PSI0',PSIp0
        !write(output_unit_write,*) 'PSI',psi1
        
        PSIN1=(psi1-PSIp0)/(PSIp_lim-PSIp0)


        !write(output_unit_write,*) 'R',R_test
        !write(output_unit_write,*) 'Z',Z_test
        !write(output_unit_write,*) 'PSIlim',PSIp_lim
        !write(output_unit_write,*) 'PSI0',PSIp0
        !write(output_unit_write,*) 'PSI1',psi1
        !write(output_unit_write,*) 'PSI0',psi0
        !write(output_unit_write,*) 'PSIN',PSIN1
        !write(output_unit_write,*) 'PSIN0',PSIN0

        vth1=sqrt(2*spp%vars%te(1))
        
        ratio = real(spp%vars%flagCon(1))*indicator(PSIN1,psi_max)* &
             R_test*EXP(-PSIN1/sigma)*V_test**2/vth1**3* &
             EXP(-(V_test/vth1)**2)/ &
             (R_buffer*EXP(-PSIN0/sigma)*V_buffer**2/vth0**3* &
             EXP(-(V_buffer/vth0)**2))
        
        accepted=.false.
        if (ratio .GE. 1.0_rp) then
           accepted=.true.
           R_buffer = R_test
           Z_buffer = Z_test
           PHI_buffer = PHI_test
           V_buffer = V_test
           XI_buffer = XI_test
        else
           !call RANDOM_NUMBER(rand_unif)
           !if (rand_unif .LT. ratio) then
           !if (get_random_mkl_U() .LT. ratio) then
           if (get_random_U() .LT. ratio) then
              accepted=.true.
              R_buffer = R_test
              Z_buffer = Z_test
              PHI_buffer = PHI_test
              V_buffer = V_test
              XI_buffer = XI_test
           end if
        end if

!        write(output_unit_write,'("R: ",E17.10)') R_buffer
!        write(output_unit_write,'("Z: ",E17.10)') Z_buffer
        
        ! Only accept sample if it is within desired boundary, but
        ! add to MC above if within buffer. This helps make the boundary
        ! more defined.
        IF ((INT(indicator(PSIN1,psi_max)).EQ.1).AND. &
             ACCEPTED) THEN
           R_samples(ii) = R_buffer
           Z_samples(ii) = Z_buffer
           PHI_samples(ii) = PHI_buffer
           V_samples(ii) = V_buffer
           XI_samples(ii) = XI_buffer

           !write(output_unit_write,*) 'PSIN',PSIN1

           
!           write(output_unit_write,*) 'RS',R_buffer
           
           ! Sample phi location uniformly
           !call RANDOM_NUMBER(rand_unif)
           !PHI_samples(ii) = 2.0_rp*C_PI*rand_unif
           !PHI_samples(ii) = 2.0_rp*C_PI*get_random_mkl_U()
           ii = ii + 1_idef 
        END IF
        
     end do

!  if (minval(R_samples(:)).lt.1._rp/params%cpp%length) stop 'error with sample'
!  write(output_unit_write,'("R_sample: ",E17.10)') R_samples(:)*params%cpp%length

     X_samples=R_samples*cos(PHI_samples)
     Y_samples=R_samples*sin(PHI_samples)

!     write(output_unit_write,*) 'R_samples',R_samples
!     write(output_unit_write,*) 'PHI_samples',PHI_samples
!     write(output_unit_write,*) 'Z_samples',Z_samples
!     write(output_unit_write,*) 'G_samples',G_samples
!     write(output_unit_write,*) 'eta_samples',eta_samples

     G_samples=1/sqrt(1-V_samples**2)
     ETA_samples=acos(XI_samples)*180/C_PI
     
  end if

  params%GC_coords=.FALSE.
  
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)


  CALL MPI_SCATTER(X_samples,spp%ppp,MPI_REAL8, &
       spp%vars%X(:,1),spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  CALL MPI_SCATTER(Y_samples,spp%ppp,MPI_REAL8, &
       spp%vars%X(:,2),spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  CALL MPI_SCATTER(Z_samples,spp%ppp,MPI_REAL8, &
       spp%vars%X(:,3),spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  CALL MPI_SCATTER(G_samples,spp%ppp,MPI_REAL8, &
       spp%vars%g,spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  CALL MPI_SCATTER(ETA_samples,spp%ppp,MPI_REAL8, &
       spp%vars%eta,spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

  
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

  if (params%orbit_model(1:2).eq.'GC') call cart_to_cyl(spp%vars%X,spp%vars%Y)

  
  
  DEALLOCATE(R_samples)
  DEALLOCATE(X_samples)
  DEALLOCATE(Y_samples)
  DEALLOCATE(Z_samples)
  DEALLOCATE(PHI_samples)
  DEALLOCATE(V_samples)
  DEALLOCATE(G_samples)
  DEALLOCATE(XI_samples)
  DEALLOCATE(ETA_samples)

end subroutine FIO_therm

SUBROUTINE load_data_from_hdf5_BMC(params,Nr_a,r_a,nRE)
  TYPE(KORC_PARAMS), INTENT(IN) :: params
  CHARACTER(MAX_STRING_LENGTH) :: filename
  CHARACTER(MAX_STRING_LENGTH) :: gname
  CHARACTER(MAX_STRING_LENGTH) :: subgname
  CHARACTER(MAX_STRING_LENGTH) :: dset
  INTEGER(HID_T) :: h5file_id
  INTEGER(HID_T) :: group_id
  INTEGER(HID_T) :: subgroup_id
  REAL(rp) :: rdatum
  INTEGER :: h5error
  INTEGER, intent(out) :: Nr_a
  REAL(rp), INTENT(out), ALLOCATABLE,DIMENSION(:) :: r_a,nRE

  filename = TRIM(filename_exp)
  call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
  if (h5error .EQ. -1) then
     write(output_unit_write,'("KORC ERROR: Something went wrong in: load_data_from_hdf5 (korc_experimental) --> h5fopen_f")')
     call korc_abort(20)
  end if

  dset = "/Nr_a"
  call load_from_hdf5(h5file_id,dset,rdatum)
  Nr_a = INT(rdatum)

  ALLOCATE(r_a(Nr_a))
  ALLOCATE(nRE(Nr_a))


  dset = "/r_a"
  call load_array_from_hdf5(h5file_id,dset,r_a)

  dset = "/nRE"
  call load_array_from_hdf5(h5file_id,dset,nRE)



  call h5fclose_f(h5file_id, h5error)
  if (h5error .EQ. -1) then
     write(output_unit_write,'("KORC ERROR: Something went wrong in: load_data_from_hdf5 (korc_experimental) --> h5fclose_f")')
  end if

END SUBROUTINE load_data_from_hdf5_BMC

FUNCTION fRE_BMC(Nr_a,r_a,nRE,rm)
  REAL(rp), INTENT(IN) 	:: rm
  INTEGER :: Nr_a
  REAL(rp), INTENT(IN),dimension(Nr_a) 	:: r_a,nRE
  REAL(rp) 				:: fRE_BMC
  REAL(rp) 				:: D
  REAL(rp) 				:: g0
  REAL(rp) 				:: g1
  REAL(rp) 				:: f0
  REAL(rp) 				:: f1
  REAL(rp) 				:: m
  INTEGER 				:: index

  !write(6,*) r_a(Nr_a),rm
  
  index = MINLOC(ABS(r_a - rm),1)
  ! index of gamma supplied to function in Hollmann input gamma range
  D = r_a(index) - rm

  !write(6,*) index
  !write(6,*) ''
  
  ! linear interpolation of Hollmann input gamma range to gamma supplied
  ! to function
  if (D.GT.0) then
     f0 = nRE(index-1)
     g0 = r_a(index-1)

     f1 = nRE(index)
     g1 = r_a(index)
  else
     f0 = nRE(index)
     g0 = r_a(index)

     f1 = nRE(index+1)
     g1 = r_a(index+1)
  end if

  m = (f1-f0)/(g1-g0)

  fRE_BMC = f0 + m*(rm - g0)
  ! end of linear interpolation, fRE_H is evaluation of input Hollmann energy
  ! distribution PDF at gamma supplied to function

END FUNCTION fRE_BMC

subroutine BMC_radial(params,spp,F,P)
  !! @note Subroutine that generates a 2D Gaussian distribution in an 
  !! elliptic torus as the initial spatial condition of a given particle 
  !! species in the simulation. @endnote
  TYPE(KORC_PARAMS), INTENT(INOUT) 	:: params
  !! Core KORC simulation parameters.
  TYPE(SPECIES), INTENT(INOUT) 		:: spp
  !! An instance of the derived type SPECIES containing all the parameters
  !! and simulation variables of the different species in the simulation.
  TYPE(FIELDS), INTENT(IN)                                   :: F
  !! An instance of the KORC derived type FIELDS.
  TYPE(PROFILES), INTENT(IN)            :: P
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: R_samples,X_samples,Y_samples
  !! Major radial location of all samples
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: PHI_samples
  !! Azimuithal angle of all samples
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: Z_samples
  !! Vertical location of all samples
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: V_samples,G_samples
  REAL(rp), DIMENSION(:), ALLOCATABLE 	:: XI_samples,ETA_samples

  REAL(rp) 				:: min_R,max_R
  REAL(rp) 				:: min_Z,max_Z
    REAL(rp) 				:: min_V,max_V
 
  REAL(rp) 				:: R_buffer
  !! Previous sample of R location
  REAL(rp) 				:: Z_buffer
  !! Previous sample of Z location
  REAL(rp) 				:: PHI_buffer
  REAL(rp) 				:: V_buffer,XI_buffer
  REAL(rp) 				:: V_test,XI_test

  REAL(rp) 				:: R_test
  !! Present sample of R location
  REAL(rp) 				:: Z_test
  !! Present sample of Z location
  REAL(rp) 				:: PHI_test

  REAL(rp) 				:: psi_max,psi_max_buff
  REAL(rp)  :: PSIp_lim,PSIP0,PSIN,PSIN0,PSIN1,sigma,psi0,psi1
  REAL(rp)  :: vth0,vth1

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

  LOGICAL :: accepted
  INTEGER,DIMENSION(33) :: seed=(/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)
  REAL(rp) :: rm_buffer,rm_test
  INTEGER :: Nr_a
  REAL(rp), ALLOCATABLE,DIMENSION(:) :: r_a,nRE

  if (params%mpi_params%rank.EQ.0_idef) then
     write(output_unit_write,*) '*** START SAMPLING ***'
  end if
  
  nsamples = spp%ppp*params%mpi_params%nmpi

  params%GC_coords=.TRUE.
  PSIp_lim=F%PSIp_lim

  if (params%field_model.eq.'M3D_C1') then
     min_R=params%rmin/params%cpp%length
     max_R=params%rmax/params%cpp%length
     min_Z=params%Zmin/params%cpp%length
     max_Z=params%Zmax/params%cpp%length

     PSIp0=F%PSIp_0
     psi_max = spp%psi_max
     psi_max_buff = spp%psi_max
  else
     min_R=minval(F%X%R)
     max_R=maxval(F%X%R)
     min_Z=minval(F%X%Z)
     max_Z=maxval(F%X%Z)

     PSIp0=F%PSIP_min
     psi_max = spp%psi_max
     psi_max_buff = spp%psi_max*2._rp
  end if  

  sigma=spp%sigmaR*params%cpp%length
  
  !write(output_unit_write,*) min_R,max_R
  !write(output_unit_write,*) min_Z,max_Z
 
  ALLOCATE(R_samples(nsamples))
  ALLOCATE(X_samples(nsamples))
  ALLOCATE(Y_samples(nsamples))
  ! Number of samples to distribute among all MPI processes
  ALLOCATE(PHI_samples(nsamples))
  ! Number of samples to distribute among all MPI processes
  ALLOCATE(Z_samples(nsamples))
  ! Number of samples to distribute among all MPI processes

  call load_data_from_hdf5_BMC(params,Nr_a,r_a,nRE)
  
  if (params%mpi_params%rank.EQ.0_idef) then
        ! Transient !

     
     R_buffer = spp%Ro
     Z_buffer = spp%Zo
     PHI_buffer = 0._rp

     

     write(6,*) 'R_buffer',R_buffer*params%cpp%length
     write(6,*) 'Z_buffer',Z_buffer*params%cpp%length
     !write(6,*) 'PSIlim',PSIp_lim
     !write(6,*) 'PSI0',PSIp0
     !write(output_unit_write,*) 'PSI1',psi1
     !write(6,*) 'PSI0',psi0
     !write(output_unit_write,*) 'PSIN1',PSIN1
     !write(6,*) 'PSIN0',PSIN0
     
     if (.not.params%SameRandSeed) then
        call init_random_seed()
     else
        call random_seed(put=seed)
     end if
     
     write(output_unit_write,'("Begin burn: ",I10)')
     flush(output_unit_write)
     accepted=.false.
     ii=1_idef
     do while (ii .LE. 1000_idef)

        if (modulo(ii,100).eq.0) then
           write(output_unit_write,'("Burn: ",I10)') ii        
        end if
        write(6,'("Burn: ",I10)') ii

        
        !R_test = R_buffer + random_norm(0.0_rp,spp%dR)
        !R_test = R_buffer + get_random_mkl_N(0.0_rp,spp%dR)
        R_test = R_buffer + get_random_N()*spp%dR
        
        !Z_test = Z_buffer + random_norm(0.0_rp,spp%dZ)
        !Z_test = Z_buffer + get_random_mkl_N(0.0_rp,spp%dZ)
        Z_test = Z_buffer + get_random_N()*spp%dZ
        


        do while ((R_test.GT.max_R).OR.(R_test .LT. min_R))
           !eta_test = eta_buffer + random_norm(0.0_rp,spp%dth)
           !eta_test = eta_buffer + get_random_mkl_N(0.0_rp,spp%dth)
           R_test = R_buffer + get_random_N()*spp%dR
        end do

        do while ((Z_test.GT.max_Z).OR.(Z_test .LT. min_Z))
           !eta_test = eta_buffer + random_norm(0.0_rp,spp%dth)
           !eta_test = eta_buffer + get_random_mkl_N(0.0_rp,spp%dth)
           Z_test = Z_buffer + get_random_N()*spp%dZ
        end do

        
        PHI_test = 2.0_rp*C_PI*get_random_U()

        rm_buffer=sqrt((R_buffer-F%AB%Ro)**2+(Z_buffer)**2)/F%AB%a
        rm_test=sqrt((R_test-F%AB%Ro)**2+(Z_test)**2)/F%AB%a

        write(6,*) 'Ro',F%AB%Ro*params%cpp%length
        write(6,*) 'a',F%AB%a*params%cpp%length
        
        write(6,*) 'R bounds',min_R*params%cpp%length,max_R*params%cpp%length
        write(6,*) 'Z bounds',min_Z*params%cpp%length,max_Z*params%cpp%length
        
        write(6,*) 'R_buffer',R_test*params%cpp%length
        write(6,*) 'Z_buffer',Z_test*params%cpp%length
        write(6,*) 'rm_buffer',rm_buffer*params%cpp%length
        
        write(6,*) 'R_test',R_test*params%cpp%length
        write(6,*) 'Z_test',Z_test*params%cpp%length
        write(6,*) 'rm_test',rm_test*params%cpp%length
        
        if (rm_test.gt.1._rp) cycle

        !write(6,*) 'R',R_test*params%cpp%length
        !write(6,*) 'Z',Z_test*params%cpp%length
        !write(6,*) 'PSIlim',PSIp_lim
        !write(6,*) 'PSI0',PSIp0
        !write(6,*) 'PSI1',psi1
        !write(6,*) 'PSI0',psi0
        !write(6,*) 'PSIN',PSIN1
        !write(6,*) 'PSIN0',PSIN0

        !write(6,*) 'te',spp%vars%te(1)
        !write(6,*) 'vth0,v_buffer',vth0,V_buffer



        !write(6,*) 'vth1,v_test',vth1,V_test
        
        ! Calculate acceptance ratio for MH algorithm. fRE function
        ! incorporates p^2 factor of spherical coordinate Jacobian
        ! for velocity phase space, factors of sin(pitch angle) for velocity
        ! phase space and cylindrical coordinate Jacobian R for spatial
        ! phase space incorporated here.

        
        ratio = indicator(rm_test,1._rp)*R_test*fRE_BMC(Nr_a,r_a,nRE,rm_test)/ &
             R_buffer*fRE_BMC(Nr_a,r_a,nRE,rm_buffer)

        accepted=.false.
        if (ratio .GE. 1.0_rp) then
           accepted=.true.
           R_buffer = R_test
           Z_buffer = Z_test
           PHI_buffer = PHI_test
           ii = ii + 1_idef

           !write(output_unit_write,*) 'PSIN',PSIN1
           write(6,*) 'accepted'
        else
!           call RANDOM_NUMBER(rand_unif)
!           if (rand_unif .LT. ratio) then
           !if (get_random_mkl_U() .LT. ratio) then
           if (get_random_U() .LT. ratio) then
              accepted=.true.
              R_buffer = R_test
              Z_buffer = Z_test
              PHI_buffer = PHI_test
              ii = ii + 1_idef

              !write(output_unit_write,*) 'PSIN',PSIN1
              write(6,*) 'accepted'
           end if
        end if
     end do
     ! Transient !

     write(output_unit_write,'("Begin sample: ",I10)')
     ii=1_idef
     do while (ii .LE. nsamples)

!        write(output_unit_write,'("sample:",I15)') ii
        
       if (modulo(ii,nsamples/10).eq.0) then
           write(output_unit_write,'("Sample: ",I10)') ii
        end if
        
        !R_test = R_buffer + random_norm(0.0_rp,spp%dR)
        !R_test = R_buffer + get_random_mkl_N(0.0_rp,spp%dR)
        R_test = R_buffer + get_random_N()*spp%dR
        
        !Z_test = Z_buffer + random_norm(0.0_rp,spp%dZ)
        !Z_test = Z_buffer + get_random_mkl_N(0.0_rp,spp%dZ)
        Z_test = Z_buffer + get_random_N()*spp%dZ

        PHI_test = 2.0_rp*C_PI*get_random_U()

        do while ((R_test.GT.max_R).OR.(R_test .LT. min_R))
           !eta_test = eta_buffer + random_norm(0.0_rp,spp%dth)
           !eta_test = eta_buffer + get_random_mkl_N(0.0_rp,spp%dth)
           R_test = R_buffer + get_random_N()*spp%dR
        end do

        do while ((Z_test.GT.max_Z).OR.(Z_test .LT. min_Z))
           !eta_test = eta_buffer + random_norm(0.0_rp,spp%dth)
           !eta_test = eta_buffer + get_random_mkl_N(0.0_rp,spp%dth)
           Z_test = Z_buffer + get_random_N()*spp%dZ
        end do

        
        rm_buffer=sqrt((R_buffer-F%AB%Ro)**2+(Z_buffer)**2)/F%AB%a
        rm_test=sqrt((R_test-F%AB%Ro)**2+(Z_test)**2)/F%AB%a

        !write(6,*) 'rm_buffer',rm_buffer
        !write(6,*) 'rm_test',rm_test
        
        if (rm_test.gt.1._rp) cycle

        !write(6,*) 'R',R_test*params%cpp%length
        !write(6,*) 'Z',Z_test*params%cpp%length
        
        ratio = indicator(rm_test,1._rp)*R_test*fRE_BMC(Nr_a,r_a,nRE,rm_test)/ &
             R_buffer*fRE_BMC(Nr_a,r_a,nRE,rm_buffer)
        
        accepted=.false.
        if (ratio .GE. 1.0_rp) then
           accepted=.true.
           R_buffer = R_test
           Z_buffer = Z_test
           PHI_buffer = PHI_test
        else
           !call RANDOM_NUMBER(rand_unif)
           !if (rand_unif .LT. ratio) then
           !if (get_random_mkl_U() .LT. ratio) then
           if (get_random_U() .LT. ratio) then
              accepted=.true.
              R_buffer = R_test
              Z_buffer = Z_test
              PHI_buffer = PHI_test
           end if
        end if

!        write(output_unit_write,'("R: ",E17.10)') R_buffer
!        write(output_unit_write,'("Z: ",E17.10)') Z_buffer
        
        ! Only accept sample if it is within desired boundary, but
        ! add to MC above if within buffer. This helps make the boundary
        ! more defined.
        IF (ACCEPTED) THEN
           R_samples(ii) = R_buffer
           Z_samples(ii) = Z_buffer
           PHI_samples(ii) = PHI_buffer

           !write(output_unit_write,*) 'PSIN',PSIN1

           
!           write(output_unit_write,*) 'RS',R_buffer
           
           ! Sample phi location uniformly
           !call RANDOM_NUMBER(rand_unif)
           !PHI_samples(ii) = 2.0_rp*C_PI*rand_unif
           !PHI_samples(ii) = 2.0_rp*C_PI*get_random_mkl_U()
           ii = ii + 1_idef
           !write(6,*) 'accepted'
        END IF
        
     end do

!  if (minval(R_samples(:)).lt.1._rp/params%cpp%length) stop 'error with sample'
!  write(output_unit_write,'("R_sample: ",E17.10)') R_samples(:)*params%cpp%length

     X_samples=R_samples*cos(PHI_samples)
     Y_samples=R_samples*sin(PHI_samples)

!     write(output_unit_write,*) 'R_samples',R_samples
!     write(output_unit_write,*) 'PHI_samples',PHI_samples
!     write(output_unit_write,*) 'Z_samples',Z_samples
!     write(output_unit_write,*) 'G_samples',G_samples
!     write(output_unit_write,*) 'eta_samples',eta_samples
     
  end if

  params%GC_coords=.FALSE.
  
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)


  CALL MPI_SCATTER(X_samples,spp%ppp,MPI_REAL8, &
       spp%vars%X(:,1),spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  CALL MPI_SCATTER(Y_samples,spp%ppp,MPI_REAL8, &
       spp%vars%X(:,2),spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  CALL MPI_SCATTER(Z_samples,spp%ppp,MPI_REAL8, &
       spp%vars%X(:,3),spp%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

  
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

  if (params%orbit_model(1:2).eq.'GC') call cart_to_cyl(spp%vars%X,spp%vars%Y)
  
  
  DEALLOCATE(R_samples)
  DEALLOCATE(X_samples)
  DEALLOCATE(Y_samples)
  DEALLOCATE(Z_samples)
  DEALLOCATE(PHI_samples)
  DEALLOCATE(r_a)
  DEALLOCATE(nRE)
  
  
end subroutine BMC_radial

subroutine intitial_spatial_distribution(params,spp,P,F)
  !! @note Subroutine that contains calls to the different subroutines 
  !! for initializing the simulated particles with various
  !! spatial distribution functions. @endnote
  TYPE(KORC_PARAMS), INTENT(INOUT) 			  :: params
  !! Core KORC simulation parameters.
  TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
  !! An instance of the derived type SPECIES containing all the parameters and 
  !! simulation variables of the different species in the simulation.
  TYPE(PROFILES), INTENT(IN)                              :: P
  !! An instance of the KORC derived type PROFILES.
  TYPE(FIELDS), INTENT(IN)                                   :: F
  !! An instance of the KORC derived type FIELDS.
  INTEGER 						  :: ss
  !! Species iterator.
  INTEGER 				:: mpierr

  do ss=1_idef,params%num_species
     SELECT CASE (TRIM(spp(ss)%spatial_distribution))
     CASE ('UNIFORM')
        call uniform(spp(ss))
     CASE ('DISK')
        call disk(params,spp(ss))
     CASE ('TORUS')
        call torus(params,spp(ss))
     CASE ('EXPONENTIAL-TORUS')
        call exponential_torus(params,spp(ss))
     CASE ('GAUSSIAN-TORUS')
        call gaussian_torus(params,spp(ss))
     CASE ('ELLIPTIC-TORUS')
        call elliptic_torus(params,spp(ss))
     CASE ('EXPONENTIAL-ELLIPTIC-TORUS')
        call exponential_elliptic_torus(params,spp(ss))
     CASE ('GAUSSIAN-ELLIPTIC-TORUS')
        call gaussian_elliptic_torus(params,spp(ss))
     CASE ('2D-GAUSSIAN-ELLIPTIC-TORUS-MH')
        call MH_gaussian_elliptic_torus(params,spp(ss))
     CASE ('AVALANCHE-4D')
        call get_Avalanche_4D(params,spp(ss),P,F)
        !! In addition to spatial distribution function, [[Avalanche_4D]]
        !! samples the avalanche distribution function used to initialize
        !! the components of velocity for all particles.
     CASE ('TRACER')
        spp(ss)%vars%X(:,1)=spp(ss)%Xtrace(1)
        spp(ss)%vars%X(:,2)=spp(ss)%Xtrace(2)
        spp(ss)%vars%X(:,3)=spp(ss)%Xtrace(3)
     CASE ('SPONG-3D')
        call Spong_3D(params,spp(ss))
     CASE ('HOLLMANN-3D')
        call get_Hollmann_distribution_3D(params,spp(ss),F)
     CASE ('HOLLMANN-3D-PSI')
        call get_Hollmann_distribution_3D_psi(params,spp(ss),F)
     CASE ('HOLLMANN-1DTRANSPORT')
        call get_Hollmann_distribution_1Dtransport(params,spp(ss),F)
     CASE('MH_psi')

#if DBG_CHECK        
#else
        if (spp(ss)%ppp*params%mpi_params%nmpi.lt.10) then
           if(params%mpi_params%rank.eq.0) then
              write(6,*) &
                   'num_samples need to be atleast 10 but is only: ', &
                   spp(ss)%ppp*params%mpi_params%nmpi
           end if
           call korc_abort(19)
        end if
#endif
        
        call MH_psi(params,spp(ss),F)
     CASE('FIO_therm')

        if (spp(ss)%ppp*params%mpi_params%nmpi.lt.10) then
           if(params%mpi_params%rank.eq.0) then
              write(6,*) &
                   'num_samples need to be atleast 10 but is only: ', &
                   spp(ss)%ppp*params%mpi_params%nmpi
           end if
           call korc_abort(19)
        end if
        
        call FIO_therm(params,spp(ss),F,P)
     CASE('BMC_radial')

        if (spp(ss)%ppp*params%mpi_params%nmpi.lt.10) then
           if(params%mpi_params%rank.eq.0) then
              write(6,*) &
                   'num_samples need to be atleast 10 but is only: ', &
                   spp(ss)%ppp*params%mpi_params%nmpi
           end if
           call korc_abort(19)
        end if
        
        call BMC_radial(params,spp(ss),F,P)
     CASE DEFAULT
        call torus(params,spp(ss))
     END SELECT
  end do
end subroutine intitial_spatial_distribution


END MODULE korc_spatial_distribution
