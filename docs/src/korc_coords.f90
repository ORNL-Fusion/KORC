module korc_coords
  !! @note Module containing subroutines to calculate the position of
  !! the simulated particles in toroidal and cylindrical coordinates. @endnote
  use korc_types
  use korc_constants

  IMPLICIT NONE

  PUBLIC :: cart_to_cyl,&
       cart_to_cyl_p,&
       cart_to_tor_check_if_confined,&
       cart_to_tor_p,&
       cyl_to_cart,&
       cyl_check_if_confined,&
       cyl_check_if_confined_p

CONTAINS


  subroutine cart_to_cyl(X,Xcyl)
    !! @note  Subroutine that converts the position of simulated particles
    !! from Cartesian \((x,y,z)\) to cylindrical \((R,\phi,Z)\) coordinates.
    !! @endnote
    !! Here, the coordinate transformation is:
    !!
    !! $$R = \sqrt{x^2 + y^2},$$
    !! $$\phi = \arctan{\left( \frac{y}{x} \right)},$$
    !! $$Z = z.$$
    implicit none
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: X
    !! Particles' position in Cartesian coordinates. X(1,:) = \(x\), X(2,:)
    !! = \(y\), X(3,:) = \(z\)
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: Xcyl
    !! Particles' position in cylindrical coordinates. Xcyl(1,:) = \(R\),
    !! Xcyl(2,:) = \(\phi\), Xcyl(3,:) = \(Z\)
    INTEGER                                                :: pp
    !! Iterator.
    INTEGER                                                :: ss
    !! Iterator.

!    write(output_unit_write,'("X_X: ",E17.10)') X(1:10,1)
!    write(output_unit_write,'("X_Y: ",E17.10)') X(1:10,2)
!    write(output_unit_write,'("X_Z: ",E17.10)') X(1:10,3)

    if (size(X,1).eq.1) then
       ss = size(X,1)
    else
       if (X(2,1).eq.0) then
          ss=1_idef
       else
          ss = size(X,1)
       end if
    endif
    

!    write(output_unit_write,*) 'varX',X(:,1)
!    write(output_unit_write,*) 'varY',X(:,2)
!    write(output_unit_write,*) 'varR',Xcyl(:,1)
!    write(output_unit_write,*) 'varPHI',Xcyl(:,2)

!    !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp) SHARED(X,Xcyl)
    do pp=1_idef,ss
!       write(output_unit_write,*) 'pp',pp
       Xcyl(pp,1) = SQRT(X(pp,1)**2 + X(pp,2)**2)
       Xcyl(pp,2) = ATAN2(X(pp,2), X(pp,1))
       Xcyl(pp,2) = MODULO(Xcyl(pp,2), 2.0_rp*C_PI)
       Xcyl(pp,3) = X(pp,3)
    end do
!    !$OMP END PARALLEL DO

!    write(output_unit_write,*) 'varX',X(:,1)
!    write(output_unit_write,*) 'varY',X(:,2)
!    write(output_unit_write,*) 'varR',Xcyl(:,1)
!    write(output_unit_write,*) 'varPHI',Xcyl(:,2)


    
  end subroutine cart_to_cyl

  subroutine cart_to_cyl_p(pchunk,X_X,X_Y,X_Z,Y_R,Y_PHI,Y_Z)
    implicit none
    INTEGER, INTENT(IN)  :: pchunk
    REAL(rp), DIMENSION(pchunk), INTENT(IN)      :: X_X
    REAL(rp), DIMENSION(pchunk), INTENT(IN)      :: X_Y
    REAL(rp), DIMENSION(pchunk), INTENT(IN)      :: X_Z

    REAL(rp), DIMENSION(pchunk), INTENT(OUT)   :: Y_R
    REAL(rp), DIMENSION(pchunk), INTENT(OUT)   :: Y_PHI
    REAL(rp), DIMENSION(pchunk), INTENT(OUT)   :: Y_Z

    INTEGER                                                :: pp
    !! Iterator.

    !$OMP SIMD
!    !$OMP& aligned(Y_R,Y_PHI,Y_Z,X_X,X_Y,X_Z)
    do pp=1_idef,pchunk
       Y_R(pp) = SQRT(X_X(pp)*X_X(pp) + X_Y(pp)*X_Y(pp))
       Y_PHI(pp) = ATAN2(X_Y(pp), X_X(pp))
       Y_PHI(pp) = MODULO(Y_PHI(pp), 2.0_rp*C_PI)
       Y_Z(pp) = X_Z(pp)
    end do
    !$OMP END SIMD

  end subroutine cart_to_cyl_p

  subroutine cyl_to_cart(Xcyl,X)
    !! @note  Subroutine that converts the position of simulated particles
    !! from cylindrical \((R,\phi,Z)\) to Cartesian \((x,y,z)\ coordinates.
    !! @endnote
    !! Here, the coordinate transformation is:
    !!
    !! $$x=R\cos(\phi),$$
    !! $$y=R\sin(\phi),,$$
    !! $$Z = z.$$
    implicit none
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT):: X
    !! Particles' position in Cartesian coordinates. X(1,:) = \(x\), X(2,:)
    !! = \(y\), X(3,:) = \(z\)
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)   :: Xcyl
    !! Particles' position in cylindrical coordinates. Xcyl(1,:) = \(R\),
    !! Xcyl(2,:) = \(\phi\), Xcyl(3,:) = \(Z\)
    INTEGER                                                :: pp
    !! Iterator.
    INTEGER                                                :: ss
    !! Iterator.
    
    if (size(Xcyl,1).eq.1) then
       ss = size(Xcyl,1)
    else
       if (Xcyl(2,1).eq.0) then
          ss=1_idef
       else
          ss = size(Xcyl,1)
       end if
    endif
    

    !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp) SHARED(X,Xcyl)
    do pp=1_idef,ss
       X(pp,1) = Xcyl(pp,1)*cos(Xcyl(pp,2))
       X(pp,2) = Xcyl(pp,1)*sin(Xcyl(pp,2))
       X(pp,3) = Xcyl(pp,3)
    end do
    !$OMP END PARALLEL DO
  end subroutine cyl_to_cart

  subroutine cyl_check_if_confined(F,Xcyl,flag)
    implicit none
    TYPE(FIELDS), INTENT(IN)                               :: F
    !! An instance of the KORC derived type FIELDS.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Xcyl
    !! Particles' position in cylindrical coordinates. Xcyl(1,:) = \(R\),
    !! Xcyl(2,:) = \(\phi\), Xcyl(3,:) = \(Z\)
    INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)   :: flag
    REAL(rp)                                               :: a
    !! Distance to plasma edge as measured from the magnetic axis.
    REAL(rp)                                               :: Ro
    !! Radial position of the magnetic axis.
    INTEGER                                                :: pp
    !! Iterator.
    INTEGER                                                :: ss
    !! Iterator.
 

    if (size(Xcyl,1).eq.1) then
       ss = size(Xcyl,1)
    else
       if (Xcyl(2,1).eq.0) then
          ss=1_idef
       else
          ss = size(Xcyl,1)
       end if
    endif

    
    a = F%AB%a
    Ro = F%AB%Ro
    
    !$OMP PARALLEL DO FIRSTPRIVATE(ss,a,Ro) PRIVATE(pp) SHARED(Xcyl,flag)
    do pp=1_idef,ss
       if (sqrt((Xcyl(pp,1)-Ro)**2+Xcyl(pp,3)**2) .gt. a) then
          flag(pp)=0_is
       endif       
    end do
    !$OMP END PARALLEL DO
  end subroutine cyl_check_if_confined

  subroutine cyl_check_if_confined_p(pchunk,a,R0,Xcyl_R,Xcyl_Z,flag)
    implicit none
    INTEGER, INTENT(IN)  :: pchunk
    REAL(rp),DIMENSION(pchunk),  INTENT(IN)      :: Xcyl_R
    REAL(rp),DIMENSION(pchunk),  INTENT(IN)      :: Xcyl_Z
    INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flag
    REAL(rp),  INTENT(IN)                            :: a,R0
    !! Distance to plasma edge as measured from the magnetic axis.
    INTEGER                                  :: cc
    
    !$OMP SIMD
!    !$OMP& aligned(Xcyl_R,Xcyl_Z,flag)
    do cc=1_idef,pchunk
       if (sqrt((Xcyl_R(cc)-R0)**2+Xcyl_Z(cc)**2) .gt. a) flag(cc)=0_is
    end do
    !$OMP END SIMD

  end subroutine cyl_check_if_confined_p

  subroutine cart_to_tor_check_if_confined(X,F,Xtor,flag)
    !! @note Subroutine that converts the position of simulated particles
    !! from Cartesian \((x,y,z)\) to toroidal \((r,\theta, \zeta)\) coordinates.
    !! In addition to performing the coordinate transformation, this
    !! subroutine checks whether a given particle is within the plasma or not.
    !! A particle is not longer considered to be within the plasma if its
    !! minor radius \(r > r_{edge}\), where \(r_{edge}\) is the radial
    !! distance to the plasma edge as measured from the magnetic axis. For
    !! more details see the analytical model of the magnetic field in
    !! [[korc_types]] and [[korc_fields]].
    !!
    !! The coordinate transformation is given by:
    !!
    !! $$r = \sqrt{ \left[\sqrt{x^2 + y^2}-R_0\right]^2 + z^2 },$$
    !! $$\theta = \arctan{\left( \frac{z}{\sqrt{x^2 + y^2}-Ro} \right)}.$$
    !! $$\zeta = \arctan{\left( \frac{x}{y} \right)},$$
    !!
    !! where \(R_0\) is the radial position of the magnetic axis.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: X
    !! Particles' position in Cartesian coordinates. X(1,:) = \(x\), X(2,:)
    !! = \(y\), X(3,:) = \(z\)
    TYPE(FIELDS), INTENT(IN)                               :: F
    !! An instance of the KORC derived type FIELDS.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: Xtor
    !! Particles' position in cylindrical coordinates. Xtor(1,:) =
    !! \(r\), Xtor(2,:) = \(\theta\), Xtor(3,:) = \(\zeta\)
    INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
    REAL(rp)                                               :: a
    !! Distance to plasma edge as measured from the magnetic axis.
    REAL(rp)                                               :: Ro
    !! Radial position of the magnetic axis.
    INTEGER                                                :: pp
    !! Iterator.
    INTEGER                                                :: ss
    !! Iterator.

    if (size(X,1).eq.1) then
       ss = size(X,1)
    else
       if (X(2,1).eq.0) then
          ss=1_idef
       else
          ss = size(X,1)
       end if
    endif
    
    a = F%AB%a
    Ro = F%AB%Ro

!    write(output_unit_write,'("X c2tor: ",E17.10)') X(1,:)
    
    !$OMP PARALLEL DO FIRSTPRIVATE(ss,a,Ro) PRIVATE(pp) SHARED(X,Xtor,flag)
    do pp=1_idef,ss
       if ( flag(pp) .EQ. 1_is ) then
          Xtor(pp,1) = SQRT( (SQRT(X(pp,1)**2 + X(pp,2)**2) - Ro)**2 + &
               X(pp,3)**2 )
          Xtor(pp,2) = ATAN2(X(pp,3), SQRT(X(pp,1)**2 + X(pp,2)**2) - Ro)
          Xtor(pp,2) = MODULO(Xtor(pp,2),2.0_rp*C_PI)
          Xtor(pp,3) = ATAN2(X(pp,1),X(pp,2))
          Xtor(pp,3) = MODULO(Xtor(pp,3),2.0_rp*C_PI)

!          write(output_unit_write,'("r: ",E17.10)') Xtor(1,1)
!          write(output_unit_write,'("a: ",E17.10)') a
!          write(output_unit_write,'("Ro: ",E17.10)') Ro
          
          if (Xtor(pp,1) .GT. F%AB%a) then
             flag(pp) = 0_is

!             stop 'error in dist init'
          end if
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine cart_to_tor_check_if_confined

  subroutine cart_to_tor_p(pchunk,R0,X_X,X_Y,X_Z,T_R,T_T,T_Z)
    INTEGER, INTENT(IN)  :: pchunk
    REAL(rp),  INTENT(IN)      :: R0
    REAL(rp),  INTENT(IN),DIMENSION(pchunk)      :: X_X,X_Y,X_Z
    REAL(rp),  INTENT(OUT),DIMENSION(pchunk)      :: T_R,T_T,T_Z
    REAL(rp),DIMENSION(pchunk)   :: RR
    INTEGER                                      :: cc
    !! Particle chunk iterator.

    !$OMP SIMD
!    !$OMP& aligned(RR,X_X,X_Y,T_R,T_T,T_Z,X_Z)
    do cc=1_idef,pchunk
       RR(cc)=SQRT(X_X(cc)*X_X(cc) + X_Y(cc)*X_Y(cc)) - R0


       T_R(cc) = SQRT( RR(cc)*RR(cc) + X_Z(cc)*X_Z(cc) )
       T_T(cc) = ATAN2(X_Z(cc), RR(cc))
       T_T(cc) = MODULO(T_T(cc),2.0_rp*C_PI)
       T_Z(cc) = ATAN2(X_X(cc),X_Y(cc))
       T_Z(cc) = MODULO(T_Z(cc),2.0_rp*C_PI)
    end do
    !$OMP END SIMD

  end subroutine cart_to_tor_p

  subroutine cart_to_tor_check_if_confined_p(pchunk,ar,R0,X_X,X_Y,X_Z, &
       T_R,T_T,T_Z,flag_cache)
    INTEGER, INTENT(IN)  :: pchunk
    REAL(rp),  INTENT(IN)      :: R0,ar
    REAL(rp),  INTENT(IN),DIMENSION(pchunk)      :: X_X,X_Y,X_Z
    REAL(rp),  INTENT(OUT),DIMENSION(pchunk)      :: T_R,T_T,T_Z
    INTEGER(is),  INTENT(INOUT),DIMENSION(pchunk)      :: flag_cache
    REAL(rp),DIMENSION(pchunk)   :: RR
    INTEGER                                      :: cc
    !! Particle chunk iterator.

    !$OMP SIMD
!    !$OMP& aligned(RR,X_X,X_Y,T_R,T_T,T_Z,X_Z)
    do cc=1_idef,pchunk
       RR(cc)=SQRT(X_X(cc)*X_X(cc) + X_Y(cc)*X_Y(cc)) - R0


       T_R(cc) = SQRT( RR(cc)*RR(cc) + X_Z(cc)*X_Z(cc) )
       T_T(cc) = ATAN2(X_Z(cc), RR(cc))
       T_T(cc) = MODULO(T_T(cc),2.0_rp*C_PI)
       T_Z(cc) = ATAN2(X_X(cc),X_Y(cc))
       T_Z(cc) = MODULO(T_Z(cc),2.0_rp*C_PI)
    end do
    !$OMP END SIMD

    !$OMP SIMD
    do cc=1_idef,pchunk
       if (T_R(cc) .GT. ar) then
          flag_cache(cc) = 0_is
       end if
    end do
    !$OMP END SIMD
    
  end subroutine cart_to_tor_check_if_confined_p
  
end module korc_coords
