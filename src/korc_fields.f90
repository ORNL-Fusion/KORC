module korc_fields
  !! @note Module containing subroutines to initialize externally
  !! generated fields, and to calculate the electric and magnetic
  !! fields when using an analytical model. @endnote
  use korc_types
  use korc_hpc
  use korc_coords
  use korc_interp
  use korc_HDF5
  use korc_input

  IMPLICIT NONE

  PUBLIC :: get_fields,&
       add_analytical_E_p,&
       analytical_fields_GC_p,&
       analytical_fields_Bmag_p,&
       analytical_fields_p,&
       initialize_fields,&
       load_field_data_from_hdf5,&
       load_dim_data_from_hdf5,&
       ALLOCATE_2D_FIELDS_ARRAYS,&
       ALLOCATE_3D_FIELDS_ARRAYS,&
       DEALLOCATE_FIELDS_ARRAYS,&
       uniform_fields_p
  PRIVATE :: get_analytical_fields,&
       analytical_fields,&
       analytical_fields_GC_init,&
       analytical_fields_GC,&
       uniform_magnetic_field,&
       uniform_electric_field,&
       uniform_fields,&
       cross,&
       analytical_electric_field_cyl,&
       ALLOCATE_V_FIELD_2D,&
       ALLOCATE_V_FIELD_3D

CONTAINS

  subroutine analytical_fields(F,Y,E,B,flag)
    !! @note Subroutine that calculates and returns the analytic electric and
    !! magnetic field for each particle in the simulation. @endnote
    !! The analytical magnetic field is given by:
    !!
    !! $$\mathbf{B}(r,\vartheta) = \frac{1}{1 + \eta \cos{\vartheta}}
    !! \left[ B_0 \hat{e}_\zeta + B_\vartheta(r) \hat{e}_\vartheta \right],$$
    !!
    !! where \(\eta = r/R_0\) is the aspect ratio, the constant \(B_0\)
    !! denotes the magnitude of the toroidal magnetic field,
    !! and \(B_\vartheta(r) = \eta B_0/q(r)\) is the poloidal magnetic
    !! field with
    !! safety factor \(q(r) = q_0\left( 1 + \frac{r^2}{\lambda^2} \right)\).
    !! The constant \(q_0\) is the safety factor at the magnetic axis and
    !! the constant \(\lambda\) is obtained from the values of \(q_0\)
    !! and \(q(r)\) at the plasma edge \(r=r_{edge}\). On the other hand,
    !! the analytical electric fields is given by:
    !!
    !! $$\mathbf{E}(r,\vartheta) = \frac{1}{1 + \eta \cos{\vartheta}}
    !! E_0 \hat{e}_\zeta,$$
    !!
    !! where \(E_0\) is the electric field as measured at the mangetic axis.
    TYPE(FIELDS), INTENT(IN)                               :: F
    !! An instance of the KORC derived type FIELDS.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
    !! Toroidal coordinates of each particle in the simulation;
    !! Y(1,:) = \(r\), Y(2,:) = \(\theta\), Y(3,:) = \(\zeta\).
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: B
    !! Magnetic field components in Cartesian coordinates;
    !! B(1,:) = \(B_x\), B(2,:) = \(B_y\), B(3,:) = \(B_z\)
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: E
    !! Electric field components in Cartesian coordinates;
    !! E(1,:) = \(E_x\), E(2,:) = \(E_y\), E(3,:) = \(E_z\)
    INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(IN)     :: flag
    !! Flag for each particle to decide whether it is being followed (flag=T)
    !! or not (flag=F).
    REAL(rp)                                               :: Ezeta
    !! Toroidal electric field \(E_\zeta\).
    REAL(rp)                                               :: Bzeta
    !! Toroidal magnetic field \(B_\zeta\).
    REAL(rp)                                               :: Bp
    !! Poloidal magnetic field \(B_\theta(r)\).
    REAL(rp)                                               :: eta
    !! Aspect ratio \(\eta\).
    REAL(rp)                                               :: q
    !! Safety profile \(q(r)\).
    INTEGER(ip)                                            :: pp ! Iterator(s)
    !! Particle iterator.
    INTEGER(ip)                                            :: ss
    !! Particle species iterator.

    ss = SIZE(Y,1)

    !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,Ezeta,Bp,Bzeta,eta,q) &
    !$OMP& SHARED(F,Y,E,B,flag)
    do pp=1_idef,ss
       if ( flag(pp) .EQ. 1_is ) then
          eta = Y(pp,1)/F%Ro
          q = F%AB%qo*(1.0_rp + (Y(pp,1)/F%AB%lambda)**2)
          Bp = -eta*F%AB%Bo/(q*(1.0_rp + eta*COS(Y(pp,2))))
          Bzeta = F%AB%Bo/( 1.0_rp + eta*COS(Y(pp,2)) )


          B(pp,1) =  Bzeta*COS(Y(pp,3)) - Bp*SIN(Y(pp,2))*SIN(Y(pp,3))
          B(pp,2) = -Bzeta*SIN(Y(pp,3)) - Bp*SIN(Y(pp,2))*COS(Y(pp,3))
          B(pp,3) = Bp*COS(Y(pp,2))

          if (abs(F%Eo) > 0) then
             Ezeta = -F%Eo/( 1.0_rp + eta*COS(Y(pp,2)) )

             E(pp,1) = Ezeta*COS(Y(pp,3))
             E(pp,2) = -Ezeta*SIN(Y(pp,3))
             E(pp,3) = 0.0_rp
          end if
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine analytical_fields

  subroutine analytical_fields_p(params,pchunk,F,X_X,X_Y,X_Z, &
       B_X,B_Y,B_Z,E_X,E_Y,E_Z,flag_cache)
    TYPE(KORC_PARAMS), INTENT(IN)                              :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                                   :: F
    INTEGER, INTENT(IN)  :: pchunk
    REAL(rp)      :: R0,B0,lam,q0,E0,ar
    REAL(rp),  INTENT(IN),DIMENSION(pchunk)      :: X_X,X_Y,X_Z
    REAL(rp),  INTENT(OUT),DIMENSION(pchunk)     :: B_X,B_Y,B_Z
    REAL(rp),  INTENT(OUT),DIMENSION(pchunk)     :: E_X,E_Y,E_Z
    INTEGER(is),  INTENT(INOUT),DIMENSION(pchunk)     :: flag_cache
    REAL(rp),DIMENSION(pchunk)     :: T_R,T_T,T_Z
    REAL(rp),DIMENSION(pchunk)                               :: Ezeta
    !! Toroidal electric field \(E_\zeta\).
    REAL(rp),DIMENSION(pchunk)                               :: Er
    !! changed YG Radial electric field \(E_\r\).
    REAL(rp),DIMENSION(pchunk)                               :: Bzeta
    !! Toroidal magnetic field \(B_\zeta\).
    REAL(rp),DIMENSION(pchunk)                              :: Bp
    !! Poloidal magnetic field \(B_\theta(r)\).
    REAL(rp),DIMENSION(pchunk)                               :: eta
    !! Aspect ratio \(\eta\).
    REAL(rp),DIMENSION(pchunk)                                :: q
    !! Safety profile \(q(r)\).
    REAL(rp),DIMENSION(pchunk)                             :: cT,sT,cZ,sZ
    INTEGER                                      :: cc
    !! Particle chunk iterator.
    REAL(rp) :: Er0,rrmn,sigmaamn

    B0=F%Bo
    E0=F%Eo
    lam=F%AB%lambda
    R0=F%AB%Ro
    q0=F%AB%qo
    ar=F%AB%a

    Er0=F%AB%Ero
    rrmn=F%AB%rmn
    sigmaamn=F%AB%sigmamn

    call cart_to_tor_check_if_confined_p(pchunk,ar,R0,X_X,X_Y,X_Z, &
         T_R,T_T,T_Z,flag_cache)

    !$OMP SIMD
    !    !$OMP& aligned(cT,sT,cZ,sZ,eta,q,Bp,Bzeta,B_X,B_Y,B_Z, &
    !    !$OMP& Ezeta,E_X,E_Y,E_Z,T_T,T_Z,T_R)
    do cc=1_idef,pchunk
       cT(cc)=cos(T_T(cc))
       sT(cc)=sin(T_T(cc))
       cZ(cc)=cos(T_Z(cc))
       sZ(cc)=sin(T_Z(cc))

       eta(cc) = T_R(cc)/R0
       q(cc) = q0*(1.0_rp + (T_R(cc)*T_R(cc)/(lam*lam)))
       Bp(cc) = -eta(cc)*B0/(q(cc)*(1.0_rp + eta(cc)*cT(cc)))
       !changed kappa  YG
       Bzeta(cc) = B0/( 1.0_rp + eta(cc)*cT(cc))


       B_X(cc) =  Bzeta(cc)*cZ(cc) - Bp(cc)*sT(cc)*sZ(cc)
       B_Y(cc) = -Bzeta(cc)*sZ(cc) - Bp(cc)*sT(cc)*cZ(cc)
       B_Z(cc) = Bp(cc)*cT(cc)

       !write(6,*) 'Ero ',Ero,'Er0 ',Er0
       !write(6,*) 'rmn ',rmn,'rrmn ',rrmn
       !write(6,*) 'sigmamn ',sigmamn,'sigmaamn ',sigmaamn
       !write(6,*) 'T_R ',T_R(cc)*params%cpp%length

       Ezeta(cc) = -E0/( 1.0_rp + eta(cc)*cT(cc))
       Er(cc) =Er0*(1/cosh((T_R(cc)-rrmn)/sigmaamn))

       E_X(cc) = Ezeta(cc)*cZ(cc)+Er(cc)*cT(cc)*sZ(cc)
       E_Y(cc) = -Ezeta(cc)*sZ(cc)+Er(cc)*cT(cc)*cZ(cc)
       E_Z(cc) = Er(cc)*sT(cc)

       !write(6,*) 'Er ',Er(cc)
    end do
    !$OMP END SIMD

  end subroutine analytical_fields_p

  subroutine analytical_fields_GC_init(params,F,Y,E,B,gradB,curlb,flag,PSIp)
    TYPE(KORC_PARAMS), INTENT(IN)      :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                               :: F
    !! An instance of the KORC derived type FIELDS.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
    !! Cylindrical coordinates of each particle in the simulation;
    !! Y(1,:) = \(r\), Y(2,:) = \(\theta\), Y(3,:) = \(\zeta\).
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: B
    !! Magnetic field components in cylindrical coordinates;
    !! B(1,:) = \(B_R\), B(2,:) = \(B_\phi\), B(3,:) = \(B_Z\)
    REAL(rp), DIMENSION(3)   :: Btmp
    !! Placeholder for magnetic field components in cylindrical coordinates;
    !! B(1,:) = \(B_R\), B(2,:) = \(B_\phi\), B(3,:) = \(B_Z\)
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: gradB
    !! Gradient of magnitude of magnetic field in cylindrical coordinates;
    !! gradB(1,:) = \(\nabla_R B\), B(2,:) = \(\nabla_\phi B_\),
    !! B(3,:) = \(\nabla_Z B\)
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: curlB
    !! Curl of magnetic field unit vector in cylindrical coordinates
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: E
    !! Electric field components in cylindricalcoordinates;
    !! E(1,:) = \(E_R\), E(2,:) = \(E_\phi\), E(3,:) = \(E_Z\)
    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)   :: PSIp
    INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(IN)     :: flag
    !! Flag for each particle to decide whether it is being followed (flag=T)
    !! or not (flag=F).
    REAL(rp)                                               :: Ezeta
    !! Toroidal electric field \(E_\zeta\).
    REAL(rp)                                               :: Bzeta
    !! Toroidal magnetic field \(B_\zeta\).
    REAL(rp)                                               :: Bp
    !! Poloidal magnetic field \(B_\theta(r)\).
    REAL(rp)                                               :: eta
    !! Aspect ratio \(\eta\).
    REAL(rp)                                               :: q
    !! Safety profile \(q(r)\).
    INTEGER(ip)                                            :: pp ! Iterator(s)
    !! Particle iterator.
    INTEGER(ip)                                            :: ss
    !! Particle species iterator.
    REAL(rp)    ::  dRBR
    REAL(rp)    ::  dRBPHI
    REAL(rp)    ::  dRBZ
    REAL(rp)    ::  dZBR
    REAL(rp)    ::  dZBPHI
    REAL(rp)    ::  dZBZ
    REAL(rp)    ::  Bmag
    REAL(rp)    ::  dRbhatPHI
    REAL(rp)    ::  dRbhatZ
    REAL(rp)    ::  dZbhatR
    REAL(rp)    ::  dZbhatPHI
    REAL(rp)    ::  qprof
    REAL(rp)    ::  rm,theta

    !write(output_unit_write,'("Y: ",E17.10)') Y

    ss = SIZE(Y,1)

    !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,rm,Btmp,qprof,dRBR,dRBPHI, &
    !$OMP dRBZ,dZBR,dZBPHI,dZBZ,Bmag,dRbhatPHI,dRbhatZ,dZbhatR,dZbhatPHI, &
    !$OMP theta,PSIp) &
    !$OMP& SHARED(F,Y,E,B,gradB,curlb,flag)
    do pp=1_idef,ss
       !       if ( flag(pp) .EQ. 1_is ) then

       rm=sqrt((Y(pp,1)-F%AB%Ro)**2+Y(pp,3)**2)
       theta=atan2(Y(pp,3),(Y(pp,1)-F%AB%Ro))
       qprof = 1.0_rp + (rm/F%AB%lambda)**2

       !write(output_unit_write,*) 'rm: ',rm
       !write(output_unit_write,*) 'R0: ',F%AB%Ro
       !write(output_unit_write,*) 'Y_R: ',Y(pp,1)
       !write(output_unit_write,*) 'theta: ',theta

       PSIp(pp)=Y(pp,1)*F%AB%lambda**2*F%Bo/ &
            (2*F%AB%qo*(F%AB%Ro+rm*cos(theta)))* &
            log(1+(rm/F%AB%lambda)**2)

       Btmp(1)=F%AB%Bo*Y(pp,3)/(F%AB%qo*qprof*Y(pp,1))
       Btmp(2)=-F%AB%Bo*F%AB%Ro/Y(pp,1)
       Btmp(3)=-F%AB%Bo*(Y(pp,1)-F%AB%Ro)/(F%AB%qo*qprof*Y(pp,1))


       B(pp,1) = Btmp(1)*COS(Y(pp,2)) - Btmp(2)*SIN(Y(pp,2))
       B(pp,2) = Btmp(1)*SIN(Y(pp,2)) + Btmp(2)*COS(Y(pp,2))
       B(pp,3) = Btmp(3)

       dRBR=-F%AB%Bo*Y(pp,3)/(F%AB%qo*qprof*Y(pp,1))*(1./Y(pp,1)+ &
            2*(Y(pp,1)-F%AB%Ro)/(F%AB%lambda**2*qprof))
       dRBPHI=F%AB%Bo*F%AB%Ro/Y(pp,1)**2
       dRBZ=F%AB%Bo/(F%AB%qo*qprof*Y(pp,1))*(-F%AB%Ro/Y(pp,1)+2*(Y(pp,1)- &
            F%AB%Ro)**2/(F%AB%lambda**2*qprof))
       dZBR=F%AB%Bo/(F%AB%qo*qprof*Y(pp,1))*(1-2*Y(pp,3)*Y(pp,3)/ &
            (F%AB%lambda**2*qprof))
       dZBPHI=0._rp
       dZBZ=F%AB%Bo*(Y(pp,1)-F%AB%Ro)/(F%AB%qo*Y(pp,1))*2*Y(pp,3)/ &
            ((F%AB%lambda*qprof)**2)

       Bmag=sqrt(B(pp,1)*B(pp,1)+B(pp,2)*B(pp,2)+B(pp,3)*B(pp,3))

       gradB(pp,1)=(B(pp,1)*dRBR+B(pp,2)*dRBPHI+B(pp,3)*dRBZ)/Bmag
       gradB(pp,2)=0._rp
       gradB(pp,3)=(B(pp,1)*dZBR+B(pp,2)*dZBPHI+B(pp,3)*dZBZ)/Bmag

       dRbhatPHI=(Bmag*dRBPHI-B(pp,2)*gradB(pp,1))/Bmag**2
       dRbhatZ=(Bmag*dRBZ-B(pp,3)*gradB(pp,1))/Bmag**2
       dZbhatR=(Bmag*dZBR-B(pp,1)*gradB(pp,3))/Bmag**2
       dZbhatPHI=(Bmag*dZBPHI-B(pp,2)*gradB(pp,3))/Bmag**2

       curlb(pp,1)=-dZbhatPHI
       curlb(pp,2)=dZbhatR-dRbhatZ
       curlb(pp,3)=B(pp,2)/(Bmag*Y(pp,1))+dRbhatPHI



       !          if (abs(F%Eo) > 0) then
       E(pp,1) = 0.0_rp
       E(pp,2) = F%Eo*F%AB%Ro/Y(pp,1)
       E(pp,3) = 0.0_rp
       !         end if
       !      end if
    end do
    !$OMP END PARALLEL DO

    !write(output_unit_write,'("B: ",E17.10)') B

  end subroutine analytical_fields_GC_init

  subroutine analytical_fields_GC(params,F,Y,E,B,gradB,curlb,flag,PSIp)
    TYPE(KORC_PARAMS), INTENT(IN)      :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                               :: F
    !! An instance of the KORC derived type FIELDS.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
    !! Cylindrical coordinates of each particle in the simulation;
    !! Y(1,:) = \(r\), Y(2,:) = \(\theta\), Y(3,:) = \(\zeta\).
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: B
    !! Magnetic field components in cylindrical coordinates;
    !! B(1,:) = \(B_R\), B(2,:) = \(B_\phi\), B(3,:) = \(B_Z\)
    REAL(rp), DIMENSION(3)   :: Btmp
    !! Placeholder for magnetic field components in cylindrical coordinates;
    !! B(1,:) = \(B_R\), B(2,:) = \(B_\phi\), B(3,:) = \(B_Z\)
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: gradB
    !! Gradient of magnitude of magnetic field in cylindrical coordinates;
    !! gradB(1,:) = \(\nabla_R B\), B(2,:) = \(\nabla_\phi B_\),
    !! B(3,:) = \(\nabla_Z B\)
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: curlB
    !! Curl of magnetic field unit vector in cylindrical coordinates
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: E
    !! Electric field components in cylindricalcoordinates;
    !! E(1,:) = \(E_R\), E(2,:) = \(E_\phi\), E(3,:) = \(E_Z\)
    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)   :: PSIp
    INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(IN)     :: flag
    !! Flag for each particle to decide whether it is being followed (flag=T)
    !! or not (flag=F).
    REAL(rp)                                               :: Ezeta
    !! Toroidal electric field \(E_\zeta\).
    REAL(rp)                                               :: Bzeta
    !! Toroidal magnetic field \(B_\zeta\).
    REAL(rp)                                               :: Bp
    !! Poloidal magnetic field \(B_\theta(r)\).
    REAL(rp)                                               :: eta
    !! Aspect ratio \(\eta\).

    REAL(rp)                                               :: q
    !! Safety profile \(q(r)\).
    INTEGER(ip)                                            :: pp ! Iterator(s)
    !! Particle iterator.
    INTEGER(ip)                                            :: ss
    !! Particle species iterator.
    REAL(rp)    ::  dRBR
    REAL(rp)    ::  dRBPHI
    REAL(rp)    ::  dRBZ
    REAL(rp)    ::  dZBR
    REAL(rp)    ::  dZBPHI
    REAL(rp)    ::  dZBZ
    REAL(rp)    ::  Bmag
    REAL(rp)    ::  dRbhatPHI
    REAL(rp)    ::  dRbhatZ
    REAL(rp)    ::  dZbhatR
    REAL(rp)    ::  dZbhatPHI
    REAL(rp)    ::  qprof
    REAL(rp)    ::  rm,theta


    ss = SIZE(Y,1)

    !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,rm,Btmp,qprof,dRBR,dRBPHI, &
    !$OMP dRBZ,dZBR,dZBPHI,dZBZ,Bmag,dRbhatPHI,dRbhatZ,dZbhatR,dZbhatPHI, &
    !$OMP theta) &
    !$OMP& SHARED(F,Y,E,B,gradB,curlb,flag,PSIp)
    do pp=1_idef,ss
       !       if ( flag(pp) .EQ. 1_is ) then

       rm=sqrt((Y(pp,1)-F%AB%Ro)**2+Y(pp,3)**2)
       theta=atan2(Y(pp,3),(Y(pp,1)-F%AB%Ro))
       qprof = 1.0_rp + (rm/F%AB%lambda)**2

       !write(6,*) 'rm: ',rm*params%cpp%length
       !write(6,*) 'R0: ',F%AB%Ro*params%cpp%length
       !write(6,*) 'Y_R: ',Y(pp,1)*params%cpp%length
       !write(6,*) 'theta: ',theta

       PSIp(pp)=Y(pp,1)*F%AB%lambda**2*F%Bo/ &
            (2*F%AB%qo*(F%AB%Ro+rm*cos(theta)))* &
            log(1+(rm/F%AB%lambda)**2)

       !write(output_unit_write,*) 'PSIp: ',PSIp(pp)

       Btmp(1)=F%AB%Bo*Y(pp,3)/(F%AB%qo*qprof*Y(pp,1))
       Btmp(2)=-F%AB%Bo*F%AB%Ro/Y(pp,1)
       Btmp(3)=-F%AB%Bo*(Y(pp,1)-F%AB%Ro)/(F%AB%qo*qprof*Y(pp,1))

       B(pp,1) = Btmp(1)
       B(pp,2) = Btmp(2)
       B(pp,3) = Btmp(3)

       dRBR=-F%AB%Bo*Y(pp,3)/(F%AB%qo*qprof*Y(pp,1))*(1./Y(pp,1)+ &
            2*(Y(pp,1)-F%AB%Ro)/(F%AB%lambda**2*qprof))
       dRBPHI=F%AB%Bo*F%AB%Ro/Y(pp,1)**2
       dRBZ=F%AB%Bo/(F%AB%qo*qprof*Y(pp,1))*(-F%AB%Ro/Y(pp,1)+2*(Y(pp,1)- &
            F%AB%Ro)**2/(F%AB%lambda**2*qprof))
       dZBR=F%AB%Bo/(F%AB%qo*qprof*Y(pp,1))*(1-2*Y(pp,3)*Y(pp,3)/ &
            (F%AB%lambda**2*qprof))
       dZBPHI=0._rp
       dZBZ=F%AB%Bo*(Y(pp,1)-F%AB%Ro)/(F%AB%qo*Y(pp,1))*2*Y(pp,3)/ &
            ((F%AB%lambda*qprof)**2)

       Bmag=sqrt(B(pp,1)*B(pp,1)+B(pp,2)*B(pp,2)+B(pp,3)*B(pp,3))

       !write(6,*) pp,B(pp,:),Bmag,dRBR,dRBPHI,dRBZ
       !write(6,*) gradB(pp,:)

       gradB(pp,1)=(B(pp,1)*dRBR+B(pp,2)*dRBPHI+B(pp,3)*dRBZ)/Bmag
       gradB(pp,2)=0._rp
       gradB(pp,3)=(B(pp,1)*dZBR+B(pp,2)*dZBPHI+B(pp,3)*dZBZ)/Bmag

       dRbhatPHI=(Bmag*dRBPHI-B(pp,2)*gradB(pp,1))/Bmag**2
       dRbhatZ=(Bmag*dRBZ-B(pp,3)*gradB(pp,1))/Bmag**2
       dZbhatR=(Bmag*dZBR-B(pp,1)*gradB(pp,3))/Bmag**2
       dZbhatPHI=(Bmag*dZBPHI-B(pp,2)*gradB(pp,3))/Bmag**2

       curlb(pp,1)=-dZbhatPHI
       curlb(pp,2)=dZbhatR-dRbhatZ
       curlb(pp,3)=B(pp,2)/(Bmag*Y(pp,1))+dRbhatPHI

       !          if (abs(F%Eo) > 0) then
       E(pp,1) = 0.0_rp
       E(pp,2) = F%Eo*F%AB%Ro/Y(pp,1)
       E(pp,3) = 0.0_rp
       !         end if
       !      end if
    end do
    !$OMP END PARALLEL DO

    !write(output_unit_write,*) 'PSIp: ',PSIp(:)
    !write(output_unit_write,*) 'B_PHI: ',B(:,2)

  end subroutine analytical_fields_GC

  subroutine analytical_fields_Bmag_p(pchunk,F,Y_R,Y_PHI,Y_Z,Bmag,E_PHI)
    INTEGER, INTENT(IN)  :: pchunk
    TYPE(FIELDS), INTENT(IN)                                   :: F
    REAL(rp)  :: R0,B0,lam,q0,EF0
    REAL(rp),DIMENSION(pchunk),INTENT(IN)  :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(pchunk) :: B_R,B_PHI,B_Z,rm,qprof
    REAL(rp),DIMENSION(pchunk),INTENT(OUT) :: Bmag,E_PHI
    integer(ip) :: cc


    B0=F%Bo
    EF0=F%Eo
    lam=F%AB%lambda
    R0=F%AB%Ro
    q0=F%AB%qo

    !$OMP SIMD
    !    !$OMP& aligned(Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,rm,qprof,Bmag)
    do cc=1_idef,pchunk
       rm(cc)=sqrt((Y_R(cc)-R0)*(Y_R(cc)-R0)+Y_Z(cc)*Y_Z(cc))
       qprof(cc) = 1.0_rp + (rm(cc)*rm(cc)/(lam*lam))

       B_R(cc)=B0*Y_Z(cc)/(q0*qprof(cc)*Y_R(cc))
       B_PHI(cc)=-B0*R0/Y_R(cc)
       B_Z(cc)=-B0*(Y_R(cc)-R0)/(q0*qprof(cc)*Y_R(cc))

       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       E_PHI(cc)=EF0*R0/Y_R(cc)
    end do
    !$OMP END SIMD

  end subroutine analytical_fields_Bmag_p

  subroutine add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    TYPE(FIELDS), INTENT(IN)                                   :: F
    INTEGER(ip),INTENT(IN)  :: tt
    REAL(rp)  :: E_dyn,E_pulse,E_width,time,arg,arg1,R0,Z0,a,E_edge
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: E_PHI
    REAL(rp),DIMENSION(params%pchunk),INTENT(IN) :: Y_R,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: rm,r_a
    integer(ip) :: cc,pchunk

    pchunk=params%pchunk

    SELECT CASE (TRIM(F%E_profile))
    CASE('D3D_LOOP')

       time=params%init_time+(params%it-1+tt)*params%dt

       E_dyn=F%E_dyn
       E_pulse=F%E_pulse
       E_width=F%E_width
       R0=F%Ro

       !write(output_unit_write,*) E_dyn,E_pulse,E_width,R0

       !$OMP SIMD
       !    !$OMP& aligned(E_PHI)
       do cc=1_idef,pchunk

          arg=(time-E_pulse)**2/(2._rp*E_width**2)
          arg1=10._rp*(time-E_pulse)/(sqrt(2._rp)*E_width)

          E_PHI(cc)=E_PHI(cc)+R0*E_dyn/Y_R(cc)*exp(-arg)*(1._rp+erf(-arg1))/2._rp
       end do
       !$OMP END SIMD
    CASE('MST_FSA')

       R0=F%AB%Ro
       Z0=F%Zo
       a=F%AB%a
       E_dyn=F%E_dyn

       !$OMP SIMD
       do cc=1_idef,pchunk

          !write(6,*) 'E_dyn',E_dyn,'E_PHI_in',E_PHI(cc)

          rm(cc)=sqrt((Y_R(cc)-R0)**2+(Y_Z(cc)-Z0)**2)
          r_a(cc)=rm(cc)/a
          E_PHI(cc) = E_PHI(cc)+E_dyn-sign((2._rp*r_a(cc)**3._rp- &
               3._rp*r_a(cc)**2._rp+1._rp)*0.05/params%cpp%Eo,E_dyn)

          !write(6,*) 'r/a',r_a,'E_PHI_out',E_PHI(cc)

       end do
       !$OMP END SIMD
    CASE('MST_FSA1')

       R0=F%AB%Ro
       Z0=F%Zo
       a=F%AB%a
       E_dyn=F%E_dyn
       E_edge=F%E_edge

       !$OMP SIMD
       do cc=1_idef,pchunk

          !write(6,*) 'E_dyn',E_dyn*params%cpp%Eo,'E_edge',E_edge*params%cpp%Eo,'E_PHI_in',E_PHI(cc)*params%cpp%Eo

          rm(cc)=sqrt((Y_R(cc)-R0)**2+(Y_Z(cc)-Z0)**2)
          r_a(cc)=rm(cc)/a
          E_PHI(cc) = E_PHI(cc)+(E_dyn-E_edge)*(1._rp-r_a(cc)**4._rp)**4._rp+E_edge

          !write(6,*) 'r/a',r_a,'E_PHI_out',E_PHI(cc)*params%cpp%Eo

       end do
       !$OMP END SIMD
    CASE('NONE')
       !$OMP SIMD
       do cc=1_idef,pchunk
          E_PHI(cc) = E_PHI(cc)
       end do
       !$OMP END SIMD
    CASE DEFAULT
       !$OMP SIMD
       do cc=1_idef,pchunk
          E_PHI(cc)=E_PHI(cc)
       end do
       !$OMP END SIMD
    END SELECT

    !write(output_unit_write,*) arg,arg1

  end subroutine add_analytical_E_p

  subroutine analytical_fields_GC_p(pchunk,F,Y_R,Y_PHI, &
       Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
       gradB_PHI,gradB_Z,PSIp)
    INTEGER, INTENT(IN) :: pchunk
    TYPE(FIELDS), INTENT(IN)                                   :: F
    REAL(rp),DIMENSION(pchunk),INTENT(IN)  :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(pchunk),INTENT(OUT) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(pchunk),INTENT(OUT) :: gradB_R,gradB_PHI,gradB_Z
    REAL(rp),DIMENSION(pchunk),INTENT(OUT) :: curlB_R,curlB_PHI,curlB_Z
    REAL(rp),DIMENSION(pchunk),INTENT(OUT) :: E_R,E_PHI,E_Z
    REAL(rp),DIMENSION(pchunk),INTENT(OUT) :: PSIp
    REAL(rp),DIMENSION(pchunk)  :: dRBR,dRBPHI,dRBZ,dZBR,dZBPHI,dZBZ,Bmag,dRbhatPHI
    REAL(rp),DIMENSION(pchunk)  :: dRbhatZ,dZbhatR,dZbhatPHI,qprof,rm,theta
    REAL(rp)  :: B0,E0,lam,R0,q0
    integer(ip) :: cc

    B0=F%Bo
    E0=F%Eo
    lam=F%AB%lambda
    R0=F%AB%Ro
    q0=F%AB%qo

    !$OMP SIMD
    !!$OMP& aligned(Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,gradB_R,gradB_PHI,gradB_Z, &
    !!$OMP& curlB_R,curlB_PHI,curlB_Z,E_R,E_PHI,E_Z,PSIp)
    do cc=1_idef,pchunk
       rm(cc)=sqrt((Y_R(cc)-R0)*(Y_R(cc)-R0)+Y_Z(cc)*Y_Z(cc))
       theta(cc)=atan2(Y_Z(cc),(Y_R(cc)-R0))
       qprof(cc) = 1.0_rp + (rm(cc)*rm(cc)/(lam*lam))

       PSIp(cc)=Y_R(cc)*lam**2*B0/ &
            (2*q0*(R0+rm(cc)*cos(theta(cc))))* &
            log(1+(rm(cc)/lam)**2)

       B_R(cc)=B0*Y_Z(cc)/(q0*qprof(cc)*Y_R(cc))
       B_PHI(cc)=-B0*R0/Y_R(cc)
       B_Z(cc)=-B0*(Y_R(cc)-R0)/(q0*qprof(cc)*Y_R(cc))

       dRBR(cc)=-B0*Y_Z(cc)/(q0*qprof(cc)*Y_R(cc))*(1./Y_R(cc)+ &
            2*(Y_R(cc)-R0)/(lam*lam*qprof(cc)))
       dRBPHI(cc)=B0*R0/(Y_R(cc)*Y_R(cc))
       dRBZ(cc)=B0/(q0*qprof(cc)*Y_R(cc))*(-R0/Y_R(cc)+2*(Y_R(cc)- &
            R0)*(Y_R(cc)-R0)/(lam*lam*qprof(cc)))
       dZBR(cc)=B0/(q0*qprof(cc)*Y_R(cc))*(1-2*Y_Z(cc)*Y_Z(cc)/ &
            (lam*lam*qprof(cc)))
       dZBPHI(cc)=0._rp
       dZBZ(cc)=B0*(Y_R(cc)-R0)/(q0*Y_R(cc))*2*Y_Z(cc)/ &
            (lam*lam*qprof(cc)*qprof(cc))

       Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       gradB_R(cc)=(B_R(cc)*dRBR(cc)+B_PHI(cc)*dRBPHI(cc)+B_Z(cc)*dRBZ(cc))/ &
            Bmag(cc)
       gradB_PHI(cc)=0._rp
       gradB_Z(cc)=(B_R(cc)*dZBR(cc)+B_PHI(cc)*dZBPHI(cc)+B_Z(cc)*dZBZ(cc))/ &
            Bmag(cc)

       dRbhatPHI(cc)=(Bmag(cc)*dRBPHI(cc)-B_PHI(cc)*gradB_R(cc))/ &
            (Bmag(cc)*Bmag(cc))
       dRbhatZ(cc)=(Bmag(cc)*dRBZ(cc)-B_Z(cc)*gradB_R(cc))/(Bmag(cc)*Bmag(cc))
       dZbhatR(cc)=(Bmag(cc)*dZBR(cc)-B_R(cc)*gradB_Z(cc))/(Bmag(cc)*Bmag(cc))
       dZbhatPHI(cc)=(Bmag(cc)*dZBPHI(cc)-B_PHI(cc)*gradB_Z(cc))/ &
            (Bmag(cc)*Bmag(cc))

       curlb_R(cc)=-dZbhatPHI(cc)
       curlb_PHI(cc)=dZbhatR(cc)-dRbhatZ(cc)
       curlb_Z(cc)=B_PHI(cc)/(Bmag(cc)*Y_R(cc))+dRbhatPHI(cc)


       E_R(cc) = 0.0_rp
       E_PHI(cc) = E0*R0/Y_R(cc)
       E_Z(cc) = 0.0_rp



    end do
    !$OMP END SIMD

  end subroutine analytical_fields_GC_p

  subroutine uniform_magnetic_field(F,B)
    !! @note Subroutine that returns the value of a uniform magnetic
    !! field. @endnote
    !! This subroutine is used only when the simulation is ran for a
    !! 'UNIFORM' plasma. As a convention, in a uniform plasma we
    !! set \(\mathbf{B} = B_0 \hat{x}\).
    TYPE(FIELDS), INTENT(IN)                               :: F
    !! An instance of the KORC derived type FIELDS.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: B
    !! Magnetic field components in Cartesian coordinates;
    !! B(1,:) = \(B_x\), B(2,:) = \(B_y\), B(3,:) = \(B_z\)
    B(:,1) = F%Bo
    B(:,2:3) = 0.0_rp
  end subroutine uniform_magnetic_field

  subroutine uniform_fields_p(pchunk,F,B_X,B_Y,B_Z,E_X,E_Y,E_Z)
    INTEGER, INTENT(IN) :: pchunk
    !! @note Subroutine that returns the value of a uniform magnetic
    !! field. @endnote
    !! This subroutine is used only when the simulation is ran for a
    !! 'UNIFORM' plasma. As a convention, in a uniform plasma we
    !! set \(\mathbf{B} = B_0 \hat{x}\).
    TYPE(FIELDS), INTENT(IN)                               :: F
    !! An instance of the KORC derived type FIELDS.
    REAL(rp),DIMENSION(pchunk), INTENT(OUT)   :: B_X,B_Y,B_Z
    REAL(rp),DIMENSION(pchunk), INTENT(OUT)   :: E_X,E_Y,E_Z
    !! Magnetic field components in Cartesian coordinates;
    !! B(1,:) = \(B_x\), B(2,:) = \(B_y\), B(3,:) = \(B_z\)
    integer(ip) :: cc

    !$OMP SIMD
    do cc=1_idef,pchunk
       B_X(cc) = F%Bo
       B_Y(cc) = 0._rp
       B_Z(cc) = 0._rp

       E_X(cc) = F%Eo
       E_Y(cc) = 0._rp
       E_Z(cc) = 0._rp
    end do
    !$OMP END SIMD

  end subroutine uniform_fields_p

  subroutine uniform_electric_field(F,E)
    !! @note Subroutine that returns the value of a uniform electric
    !! field. @endnote
    !! This subroutie is used only when the simulation is ran for a
    !! 'UNIFORM' plasma. As a convention, in a uniform plasma we set
    !! \(\mathbf{E} = E_0 \hat{x}\).
    TYPE(FIELDS), INTENT(IN)                               :: F
    !! An instance of the KORC derived type FIELDS.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: E
    !! Electric field components in Cartesian coordinates;
    !! E(1,:) = \(E_x\), E(2,:) = \(E_y\), E(3,:) = \(E_z\)

    E(:,1) = F%Eo
    E(:,2:3) = 0.0_rp
  end subroutine uniform_electric_field

  subroutine analytical_electric_field_cyl(F,Y,E,flag)
    !! @note Subrotuine that calculates and returns the electric field using the
    !! same analytical model of the 'analytical_fields' subroutine. @endnote
    TYPE(FIELDS), INTENT(IN)                               :: F
    !! An instance of the KORC derived type FIELDS.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
    !! Cylindrical coordinates of each particle in the simulation;
    !! Y(1,:) = \(R\), Y(2,:) = \(\phi\), Y(3,:) = \(Z\).
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: E
    !! Electric field components in Cartesian coordinates;
    !!  E(1,:) = \(E_x\), E(2,:) = \(E_y\), E(3,:) = \(E_z\)
    INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(IN)     :: flag
    !! Flag for each particle to decide whether it is being followed (flag=T)
    !! or not (flag=F).
    REAL(rp)                                               :: Ephi
    !! Azimuthal electric field.
    INTEGER(ip)                                            :: pp
    !! Particle iterator.
    INTEGER(ip)                                            :: ss
    !! Particle species iterator.

    if (abs(F%Eo) > 0) then
       ss = SIZE(Y,1)
       !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,Ephi) SHARED(F,Y,E,flag)
       do pp=1_idef,ss
          if ( flag(pp) .EQ. 1_is ) then
             Ephi = F%Eo*F%Ro/Y(pp,1)

             E(pp,1) = -Ephi*SIN(Y(pp,2))
             E(pp,2) = Ephi*COS(Y(pp,2))
             E(pp,3) = 0.0_rp
          end if
       end do
       !$OMP END PARALLEL DO
    end if
  end subroutine analytical_electric_field_cyl

  subroutine get_analytical_fields(params,vars,F)
    !! @note Interface for calculating the analytical electric and magnetic
    !! fields for each particle in the simulation. @endnote
    TYPE(KORC_PARAMS), INTENT(IN)      :: params
    !! Core KORC simulation parameters.
    TYPE(PARTICLES), INTENT(INOUT) :: vars
    !! An instance of the KORC derived type PARTICLES.
    TYPE(FIELDS), INTENT(IN)       :: F
    !! An instance of the KORC derived type FIELDS.

    if (params%orbit_model(1:2).eq.'FO') then

       call cart_to_tor_check_if_confined(vars%X,F,vars%Y,vars%flagCon)

       call analytical_fields(F,vars%Y, vars%E, vars%B, vars%flagCon)

       !       call cart_to_cyl(vars%X,vars%Y)

    elseif (params%orbit_model(1:2).eq.'GC') then

       if (.not.params%GC_coords) then

          call cart_to_cyl(vars%X,vars%Y)

          call cyl_check_if_confined(F,vars%Y,vars%flagCon)

          call analytical_fields_GC_init(params,F,vars%Y, vars%E, vars%B, &
               vars%gradB,vars%curlb, vars%flagCon, vars%PSI_P)

       else

          call cyl_check_if_confined(F,vars%Y,vars%flagCon)

          call analytical_fields_GC(params,F,vars%Y, vars%E, vars%B, &
               vars%gradB,vars%curlb, vars%flagCon,vars%PSI_P)

       end if

    endif

  end subroutine get_analytical_fields

  subroutine uniform_fields(vars,F)
    !! @note Interface for calculating the uniform electric and magnetic
    !! fields for each particle in the simulation. @endnote
    TYPE(PARTICLES), INTENT(INOUT) :: vars
    !! An instance of the KORC derived type PARTICLES.
    TYPE(FIELDS), INTENT(IN)       :: F
    !! An instance of the KORC derived type FIELDS.

    call uniform_magnetic_field(F, vars%B)

    call uniform_electric_field(F, vars%E)
  end subroutine uniform_fields

  pure function cross(a,b)
    !! @note Function that calculates the cross product of the two
    !! vectors \(\mathbf{a}\) and \(\mathbf{b}\). @endnote
    REAL(rp), DIMENSION(3)             :: cross
    !! Cross product \(\mathbf{a}\times \mathbf{b}\)
    REAL(rp), DIMENSION(3), INTENT(IN) :: a
    !!  Vector \(\mathbf{a}\).
    REAL(rp), DIMENSION(3), INTENT(IN) :: b
    !!  Vector \(\mathbf{b}\).

    cross(1) = a(2)*b(3) - a(3)*b(2)
    cross(2) = a(3)*b(1) - a(1)*b(3)
    cross(3) = a(1)*b(2) - a(2)*b(1)
  end function cross


  subroutine unitVectors(params,Xo,F,b1,b2,b3,flag,cart,hint)
    !! @note Subrotuine that calculates an orthonormal basis using information
    !! of the (local) magnetic field at position \(\mathbf{X}_0\). @endnote
    TYPE(KORC_PARAMS), INTENT(IN)                                      :: params
    !! Core KORC simulation parameters.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)                  :: Xo
    !! Array with the position of the simulated particles.
    TYPE(FIELDS), INTENT(IN)                                           :: F
    !! F An instance of the KORC derived type FIELDS.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)               :: b1
    !! Basis vector pointing along the local magnetic field,
    !! that is, along \(\mathbf{b} = \mathbf{B}/B\).
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)               :: b2
    !!  Basis vector perpendicular to b1
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)               :: b3
    !! Basis vector perpendicular to b1 and b2.
    INTEGER(is), DIMENSION(:), ALLOCATABLE, OPTIONAL, INTENT(INOUT)    :: flag
    !! Flag for each particle to decide whether it is being
    !! followed (flag=T) or not (flag=F).
    TYPE(C_PTR), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: hint
    !! Flag for each particle to decide whether it is being
    !! followed (flag=T) or not (flag=F).
    TYPE(PARTICLES)                                                    :: vars
    !! A temporary instance of the KORC derived type PARTICLES.
    INTEGER                                                            :: ii
    !! Iterator.
    INTEGER                                                            :: ppp
    !! Number of particles.
    LOGICAL :: cart
    REAL(rp), DIMENSION(3) ::b1tmp,b2tmp,b3tmp,tmpvec

    !write(output_unit_write,*) 'in unitVector'

    ppp = SIZE(Xo,1) ! Number of particles

    ALLOCATE( vars%X(ppp,3) )
    ALLOCATE( vars%Y(ppp,3) )
    ALLOCATE( vars%B(ppp,3) )
    ALLOCATE( vars%gradB(ppp,3) )
    ALLOCATE( vars%curlb(ppp,3) )
    ALLOCATE( vars%PSI_P(ppp) )
    ALLOCATE( vars%E(ppp,3) )
    ALLOCATE( vars%flagCon(ppp) )
    ALLOCATE( vars%initLCFS(ppp) )

#ifdef FIO
    ALLOCATE( vars%hint(ppp) )
#endif

    vars%X = Xo
#ifdef FIO
    vars%hint = hint
#endif
    vars%flagCon = flag
    vars%initLCFS = 0_is
    vars%B=0._rp
    vars%PSI_P=0._rp
    vars%cart=.false.

    !write(output_unit_write,*) 'before init_random_seed'

    call init_random_seed(params)

   ! write(output_unit_write,*) 'before get_fields'

    !write(6,*) 'before first get fields'
    call get_fields(params,vars,F)
    !write(6,*) 'before second get fields'

    !write(6,'("Bx: ",E17.10)') vars%B(:,1)*params%cpp%Bo
    !write(6,'("By: ",E17.10)') vars%B(:,2)*params%cpp%Bo
    !write(6,'("Bz: ",E17.10)') vars%B(:,3)*params%cpp%Bo

        !write(output_unit_write,*) 'before b1,b2,b3 calculation'

    tmpvec=(/1.0_rp,1.0_rp,1.0_rp/)

    do ii=1_idef,ppp
       !write(6,*) 'ii',ii
       if ( vars%flagCon(ii) .EQ. 1_idef ) then
          b1tmp = vars%B(ii,:)/sqrt(vars%B(ii,1)*vars%B(ii,1)+ &
               vars%B(ii,2)*vars%B(ii,2)+vars%B(ii,3)*vars%B(ii,3))

          b2tmp = cross(b1tmp,tmpvec)
          b2tmp = b2tmp/sqrt(b2tmp(1)*b2tmp(1)+b2tmp(2)*b2tmp(2)+ &
               b2tmp(3)*b2tmp(3))

          b3tmp = cross(b1tmp,b2tmp)
          b3tmp = b3tmp/sqrt(b3tmp(1)*b3tmp(1)+b3tmp(2)*b3tmp(2)+ &
               b3tmp(3)*b3tmp(3))
       end if
       b1(ii,:)=b1tmp
       b2(ii,:)=b2tmp
       b3(ii,:)=b3tmp
    end do

    !write(output_unit_write,*) 'before copying hint and flag'
#ifdef FIO
    hint = vars%hint
#endif

    if (PRESENT(flag)) then
       flag = vars%flagCon
    end if

    DEALLOCATE( vars%X )
    DEALLOCATE( vars%Y )
    DEALLOCATE( vars%B )
    DEALLOCATE( vars%PSI_P )
    DEALLOCATE( vars%gradB )
    DEALLOCATE( vars%curlb )
    DEALLOCATE( vars%E )
    DEALLOCATE( vars%flagCon )
#ifdef FIO
    DEALLOCATE( vars%hint)
#endif

    !write(output_unit_write,*) 'out unitVectors'

  end subroutine unitVectors


  subroutine get_fields(params,vars,F)
    !! @note Inferface with calls to subroutines for calculating the electric
    !! and magnetic field for each particle in the simulation. @endnote
    TYPE(KORC_PARAMS), INTENT(IN)      :: params
    !! Core KORC simulation parameters.
    TYPE(PARTICLES), INTENT(INOUT)     :: vars
    !!  An instance of the KORC derived type PARTICLES.
    TYPE(FIELDS), INTENT(IN)           :: F
    !! An instance of the KORC derived type FIELDS.

    !write(6,*) params%field_model

    if (params%field_model(1:10).eq.'ANALYTICAL') then
    !SELECT CASE (TRIM(params%field_model))
    !CASE('ANALYTICAL')
       if (params%field_eval.eq.'eqn') then
          call get_analytical_fields(params,vars, F)
       else
          call interp_fields(params,vars, F)
       end if
    else if (params%field_model(1:8).eq.'EXTERNAL') then

       !       write(output_unit_write,'("2 size of PSI_P: ",I16)') size(vars%PSI_P)

       call interp_fields(params,vars, F)

       !write(output_unit_write,'("get_fields")')
       !write(output_unit_write,'("B_X: ",E17.10)') vars%B(:,1)
       !write(output_unit_write,'("B_Z: ",E17.10)') vars%B(:,2)
       !write(output_unit_write,'("B_Y: ",E17.10)') vars%B(:,3)

       !if (F%Efield.AND..NOT.F%Efield_in_file) then
       !   call analytical_electric_field_cyl(F,vars%Y,vars%E,vars%flagCon)
       !end if
    else if (TRIM(params%field_model).eq.'M3D_C1'.or. &
         TRIM(params%field_model).eq.'NIMROD') then
       !write(6,*) 'get_fields'

       call interp_fields(params,vars, F)

    else if (params%field_model.eq.'UNIFORM') then

       call uniform_fields(vars, F)
    end if
  end subroutine get_fields

  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !
  ! * * *  SUBROUTINES FOR INITIALIZING FIELDS   * * * !
  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !

  subroutine initialize_fields(params,F)
    !! @note Subroutine that initializes the analytical or externally
    !! calculated electric and magnetic fields. @endnote
    !! In this subroutine we load the parameters of the electric and
    !! magnetic fields from the namelists 'analytical_fields_params' and
    !! 'externalPlasmaModel' in the input file.
    TYPE(KORC_PARAMS), INTENT(INOUT)  :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(OUT)      :: F
    !! An instance of the KORC derived type FIELDS.
    !REAL(rp)                       :: Bo
    !! Magnetic field at magnetic axis for an 'ANALYTICAL' magnetic field,
    !! or the magnitude of the magnetic field for a 'UNFIROM' plasma.
    !REAL(rp)                       :: minor_radius
    !! Plasma edge \(r_{edge}\) as measured from the magnetic axis.
    !REAL(rp)                       :: major_radius
    !! Radial position of the magnetic axis \(R_0\)
    !REAL(rp)                       :: qa
    !! Safety factor at the plasma edge.
    !REAL(rp)                       :: qo
    !! Safety factor at the magnetic axis \(q_0\).
    !CHARACTER(MAX_STRING_LENGTH)   :: current_direction
    !! String with information about the direction of the plasma current,
    !! 'PARALLEL'  or 'ANTI-PARALLEL' to the toroidal magnetic field.
    !CHARACTER(MAX_STRING_LENGTH)   :: E_model
    !REAL(rp)                       :: Eo,E_dyn,E_pulse,E_width
    !! Electric field at the magnetic axis.
    !LOGICAL                        :: Efield
    !! Logical variable that specifies if the electric field is
    !! going to be used on in a given simulation.
    !LOGICAL                        :: dBfield
    !LOGICAL                        :: Bfield
    !! Logical variable that specifies if the magnetic field is
    !! going to be used on in a given simulation.
    !LOGICAL                        :: Bflux
    !LOGICAL                        :: Bflux3D
    !LOGICAL                        :: Dim2x1t
    !LOGICAL                        :: E_2x1t,ReInterp_2x1t
    !! Logical variable that specifies if the poloidal magnetic
    !! flux is going to be used on in a given simulation.
    !LOGICAL                        :: axisymmetric_fields
    !! Logical variable that specifies if the plasma is axisymmetric.
    INTEGER                        :: ii
    !! Iterators for creating mesh for GC model with analytic fields
    INTEGER                        :: kk
    !! Iterators for creating mesh for GC model with analytic fields
    !INTEGER                        :: nR
    !! Number of mesh points in R for grid in GC model of analytical field
    !INTEGER                        :: nZ,nPHI
    !! Number of mesh points in Z for grid in GC model of analytical field
    real(rp)                       :: rm
    !! Minor radius at each position in the grid for
    !! GC model of analytical field
    real(rp)                       :: qr
    !! Safety factor at each position in the grid for
    !! GC model of analytical field
    real(rp)                       :: theta
    !! Poloidal angle at each position in the grid for
    !! GC model of analytical field
    logical :: test
    !integer :: res_double
    real(rp) :: RMAX,RMIN,ZMAX,ZMIN
    !integer :: dim_1D,ind0_2x1t
    !real(rp) :: dt_E_SC,Ip_exp,PSIp_lim,PSIp_0
    !real(rp) :: t0_2x1t


    !NAMELIST /analytical_fields_params/ Bo,minor_radius,major_radius,&
    !     qa,qo,Eo,current_direction,nR,nZ,nPHI,dim_1D,dt_E_SC,Ip_exp, &
    !     E_dyn,E_pulse,E_width
    !NAMELIST /externalPlasmaModel/ Efield, Bfield, Bflux,Bflux3D,dBfield, &
    !     axisymmetric_fields, Eo,E_dyn,E_pulse,E_width,res_double, &
    !     dim_1D,dt_E_SC,Ip_exp,PSIp_lim,Dim2x1t,t0_2x1t,E_2x1t,ReInterp_2x1t, &
    !     ind0_2x1t,PSIp_0

#ifdef FIO
    F%FIO_B = -1
    F%FIO_E = -1
#endif

    if (params%mpi_params%rank .EQ. 0) then
       write(output_unit_write,'(/,"* * * * * * * * INITIALIZING FIELDS * * * * * * * *")')
    end if

!    SELECT CASE (TRIM(params%field_model))
    if (params%field_model(1:10).eq.'ANALYTICAL') then
!    CASE('ANALYTICAL')
       ! Load the parameters of the analytical magnetic field
       !open(unit=default_unit_open,file=TRIM(params%path_to_inputs), &
       !     status='OLD',form='formatted')
       !read(default_unit_open,nml=analytical_fields_params)
       !close(default_unit_open)

       F%AB%Bo = Bo
       F%AB%a = minor_radius
       F%AB%Ro = major_radius
       F%Ro = major_radius
       F%Zo = 0.0_rp
       F%AB%qa = qa
       F%AB%qo = qo
       F%AB%lambda = F%AB%a/SQRT(qa/qo - 1.0_rp)
       F%AB%Bpo = F%AB%lambda*F%AB%Bo/(F%AB%qo*F%AB%Ro)
       F%AB%current_direction = TRIM(current_direction)
       SELECT CASE (TRIM(F%AB%current_direction))
       CASE('PARALLEL')
          F%AB%Bp_sign = 1.0_rp
       CASE('ANTI-PARALLEL')
          F%AB%Bp_sign = -1.0_rp
       CASE DEFAULT
       END SELECT
       F%Eo = Eo
       F%Bo = F%AB%Bo

       F%E_dyn = E_dyn
       F%E_edge = E_edge
       F%E_pulse = E_pulse
       F%E_width = E_width

       F%PSIp_lim=PSIp_lim

       F%AB%Ero=Ero
       F%AB%rmn=rmn
       F%AB%sigmamn=sigmamn

       !write(output_unit_write,*) E_dyn,E_pulse,E_width

       F%res_double=res_double

       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'("ANALYTIC")')
          write(output_unit_write,'("Magnetic field: ",E17.10)') F%Bo
          write(output_unit_write,'("Electric field: ",E17.10)') F%Eo

       end if


       if (params%field_eval.eq.'interp') then
          F%dims(1) = nR
          F%dims(2) = nPHI
          F%dims(3) = nZ

          if (params%field_model(12:14).eq.'PSI') then

             F%axisymmetric_fields = .TRUE.
             F%Bfield=.TRUE.
             F%Efield=.TRUE.
             F%Bflux=.TRUE.

             call ALLOCATE_2D_FIELDS_ARRAYS(params,F,F%Bfield,F%Bflux, &
                  .false.,F%Efield,.FALSE.,.FALSE.)

             do ii=1_idef,F%dims(1)
                F%X%R(ii)=(F%Ro-F%AB%a)+(ii-1)*2*F%AB%a/(F%dims(1)-1)
             end do
             do ii=1_idef,F%dims(3)
                F%X%Z(ii)=(F%Zo-F%AB%a)+(ii-1)*2*F%AB%a/(F%dims(3)-1)
             end do

             !write(6,*) F%X%R
             !write(6,*) F%X%Z

             do ii=1_idef,F%dims(1)
                do kk=1_idef,F%dims(3)
                   rm=sqrt((F%X%R(ii)-F%Ro)**2+(F%X%Z(kk)-F%Zo)**2)
                   qr=F%AB%qo*(1+(rm/F%AB%lambda)**2)
                   theta=atan2(F%X%Z(kk)-F%Zo,F%X%R(ii)-F%Ro)
                   F%B_2D%R(ii,kk)=(rm/F%X%R(ii))* &
                        (F%AB%Bo/qr)*sin(theta)
                   F%B_2D%PHI(ii,kk)=-(F%Ro/F%X%R(ii))*F%AB%Bo
                   F%B_2D%Z(ii,kk)=-(rm/F%X%R(ii))* &
                        (F%AB%Bo/qr)*cos(theta)
                   F%E_2D%R(ii,kk)=0.0_rp
                   F%E_2D%PHI(ii,kk)=-(F%Ro/F%X%R(ii))*F%Eo
                   F%E_2D%Z(ii,kk)=0.0_rp

                   F%PSIp(ii,kk)=F%X%R(ii)*F%AB%lambda**2*F%Bo/ &
                        (2*F%AB%qo*(F%Ro+rm*cos(theta)))* &
                        log(1+(rm/F%AB%lambda)**2)

                   !! Sign convention in analytical fields corresponds to
                   !! DIII-D fields with \(B_\phi<0\) and \(B_\theta<0\).
                   F%FLAG2D=1.
                end do
             end do

             F%FLAG2D(1:2,:)=0.
             F%FLAG2D(F%dims(1)-1:F%dims(1),:)=0.
             F%FLAG2D(:,1:2)=0.
             F%FLAG2D(:,F%dims(3)-1:F%dims(3))=0.

             if (F%Bflux) F%PSIP_min=minval(F%PSIp)

             F%Bfield=.FALSE.

          else if (params%field_model(12:13).eq.'3D') then

             F%axisymmetric_fields = .FALSE.
             F%Bfield=.TRUE.
             F%Efield=.TRUE.

             call ALLOCATE_3D_FIELDS_ARRAYS(params,F,F%Bfield,F%Efield,.false.)

             do ii=1_idef,F%dims(1)
                F%X%R(ii)=(F%Ro-F%AB%a)+(ii-1)*2*F%AB%a/(F%dims(1)-1)
             end do
             do ii=1_idef,F%dims(2)
                F%X%PHI(ii)=0._rp+(ii-1)*2*C_PI/(F%dims(1)-1)
             end do
             do ii=1_idef,F%dims(3)
                F%X%Z(ii)=(F%Zo-F%AB%a)+(ii-1)*2*F%AB%a/(F%dims(3)-1)
             end do

             !write(output_unit_write,*) size(F%B_3D%R)

             do ii=1_idef,F%dims(1)
                do kk=1_idef,F%dims(3)

                   !write(output_unit_write,*) ii,kk

                   rm=sqrt((F%X%R(ii)-F%Ro)**2+(F%X%Z(kk)-F%Zo)**2)
                   qr=F%AB%qo*(1+(rm/F%AB%lambda)**2)
                   theta=atan2(F%X%Z(kk)-F%Zo,F%X%R(ii)-F%Ro)
                   F%B_3D%R(ii,:,kk)=(rm/F%X%R(ii))* &
                        (F%AB%Bo/qr)*sin(theta)
                   F%B_3D%PHI(ii,:,kk)=-(F%Ro/F%X%R(ii))*F%AB%Bo

                   !write(output_unit_write,*) F%B_3D%PHI(ii,1,kk)

                   F%B_3D%Z(ii,:,kk)=-(rm/F%X%R(ii))* &
                        (F%AB%Bo/qr)*cos(theta)
                   F%E_3D%R(ii,:,kk)=0.0_rp
                   F%E_3D%PHI(ii,:,kk)=-(F%Ro/F%X%R(ii))*F%Eo
                   F%E_3D%Z(ii,:,kk)=0.0_rp


                   !! Sign convention in analytical fields corresponds to
                   !! DIII-D fields with \(B_\phi<0\) and \(B_\theta<0\).
                   F%FLAG3D=1.
                end do
             end do

             F%FLAG3D(1:2,:,:)=0.
             F%FLAG3D(F%dims(1)-1:F%dims(1),:,:)=0.
             F%FLAG3D(:,:,1:2)=0.
             F%FLAG3D(:,:,F%dims(3)-1:F%dims(3))=0.




          end if


          if (params%orbit_model(3:5).eq.'pre') then
             if (params%mpi_params%rank .EQ. 0) then
                write(output_unit_write,'("Initializing GC fields from analytic EM fields")')
             end if

             if (params%field_model(12:13).eq.'2D') then
                call initialize_GC_fields(F)
             else if (params%field_model(12:13).eq.'3D') then
                call initialize_GC_fields_3D(F)
             end if

          end if

          !F%Bfield= .FALSE.
          !F%axisymmetric_fields = .TRUE.
          !F%Bflux=.TRUE.
          !F%Efield=.FALSE.

       end if

!    CASE('EXTERNAL')
    else if (params%field_model(1:8).eq.'EXTERNAL') then
       ! Load the magnetic field from an external HDF5 file
       !open(unit=default_unit_open,file=TRIM(params%path_to_inputs), &
       !     status='OLD',form='formatted')
       !read(default_unit_open,nml=externalPlasmaModel)
       !close(default_unit_open)

       F%Bfield = Bfield
       F%B1field = B1field
       F%dBfield = dBfield
       F%Bflux = Bflux
       F%Bflux3D = Bflux3D
       F%Efield = Efield
       F%E1field = E1field
       F%axisymmetric_fields = axisymmetric_fields
       F%Dim2x1t = Dim2x1t
       F%ReInterp_2x1t = ReInterp_2x1t
       F%t0_2x1t = t0_2x1t
       F%ind0_2x1t = ind0_2x1t
       F%psip_conv = psip_conv
       F%MARS_AMP_Scale = MARS_AMP_Scale
       F%MARS_phase = MARS_phase
       F%Analytic_IWL=Analytic_IWL
       F%ntiles=ntiles
       F%circumradius=circumradius
       F%AORSA_AMP_Scale=AORSA_AMP_Scale
       F%AORSA_freq=AORSA_freq
       F%AORSA_nmode=AORSA_nmode
       F%useLCFS = useLCFS
       F%useDiMES = useDiMES
       F%DiMESloc = DiMESloc
       F%DiMESdims = DiMESdims

       if (params%proceed.and.F%ReInterp_2x1t) then
          call load_prev_iter(params)
          F%ind_2x1t=params%prev_iter_2x1t
       else
          F%ind_2x1t=F%ind0_2x1t
       end if


       F%E_2x1t = E_2x1t

       F%E_profile = E_profile
       F%E_dyn = E_dyn
       F%E_edge = E_edge
       F%E_pulse = E_pulse
       F%E_width = E_width
       F%AB%a = minor_radius
       F%AB%Ro = major_radius

       F%PSIp_lim=PSIp_lim

       F%res_double=res_double

       !write(output_unit_write,'("E_dyn: ",E17.10)') E_dyn
!       write(output_unit_write,'("E_pulse: ",E17.10)') E_pulse
!       write(output_unit_write,'("E_width: ",E17.10)') E_width

       call load_dim_data_from_hdf5(params,F)
       !sets F%dims for 2D or 3D data

       !write(6,*) F%dims

       call which_fields_in_file(params,F%Bfield_in_file,F%Efield_in_file, &
            F%Bflux_in_file,F%dBfield_in_file,F%B1field_in_file, &
            F%E1field_in_file)


       if (F%Bflux.AND..NOT.F%Bflux_in_file) then
          write(output_unit_write,'("ERROR: Magnetic flux to be used but no data in file!")')
          call KORC_ABORT(18)
       end if

       if (F%Bfield.AND..NOT.F%Bfield_in_file) then
          write(output_unit_write,'("ERROR: Magnetic field to be used but no data in file!")')
          call KORC_ABORT(18)
       end if

       if (F%B1field.AND..NOT.F%B1field_in_file) then
          write(output_unit_write,'("ERROR: Magnetic perturbation field to be used but no data in file!")')
          call KORC_ABORT(18)
       end if

       if (F%E1field.AND..NOT.F%E1field_in_file) then
          write(output_unit_write,'("ERROR: Electric perturbation field to be used but no data in file!")')
          call KORC_ABORT(18)
       end if

       if (F%dBfield.AND..NOT.F%dBfield_in_file) then
          write(output_unit_write,'("ERROR: differential Magnetic field to be used &
               but no data in file!")')
          call KORC_ABORT(18)
       end if

       if (F%Efield.AND..NOT.F%Efield_in_file) then
          if (params%mpi_params%rank.EQ.0_idef) then
             write(output_unit_write,'(/,"* * * * * * * * * *  FIELDS  * * * * * * * * * *")')
             write(output_unit_write,'("MESSAGE: Analytical electric field will be used.")')
             write(output_unit_write,'("* * * * * * * * * * * * ** * * * * * * * * * * *",/)')
             flush(output_unit_write)
          end if
       end if

       if (F%axisymmetric_fields) then

          if (F%Dim2x1t) then

             call ALLOCATE_2D_FIELDS_ARRAYS(params,F,F%Bfield, &
                  F%Bflux,F%dBfield,F%Efield.AND.F%Efield_in_file,F%B1field, &
                  F%E1field)

             call ALLOCATE_3D_FIELDS_ARRAYS(params,F,F%Bfield, &
                  F%Efield,F%dBfield)

          else if (params%orbit_model(1:2).eq.'FO') then

             call ALLOCATE_2D_FIELDS_ARRAYS(params,F,F%Bfield, &
                  .TRUE.,F%dBfield,F%Efield,F%B1field,F%E1field)

          else

             call ALLOCATE_2D_FIELDS_ARRAYS(params,F,F%Bfield, &
                  F%Bflux,F%dBfield,F%Efield,F%B1field,F%E1field)

          end if

       else if ((params%field_model(10:13).eq.'MARS').OR. &
            (params%field_model(10:14).eq.'AORSA')) then

          call ALLOCATE_2D_FIELDS_ARRAYS(params,F,F%Bfield, &
               F%Bflux,F%dBfield,F%Efield,F%B1field,F%E1field)

       else
          call ALLOCATE_3D_FIELDS_ARRAYS(params,F,F%Bfield,F%Efield,F%dBfield)

       end if


       !allocates 2D or 3D data arrays (fields and spatial)


       call load_field_data_from_hdf5(params,F)

              !write(output_unit_write,*) F%PSIp

       !write(output_unit_write,*) F%E_3D%PHI(:,F%ind0_2x1t,:)

       !write(6,*) F%B1Re_2D%R(:,200)
       !flush(6)

!       end if

       if (F%Bflux.and..not.(params%field_model(10:13).eq.'MARS')) then
          F%PSIP_min=minval(F%PSIp)

       else if(F%Bflux3D) then
          F%PSIP_min=minval(F%PSIp3D(:,1,:))
       end if

       if ((.not.F%Efield_in_file).and.(.not.F%Dim2x1t).and.F%Efield) then

          if (F%Eo.eq.0) F%Eo = Eo

          if (F%axisymmetric_fields) then
             F%E_2D%R=0._rp
             do ii=1_idef,F%dims(1)
                F%E_2D%PHI(ii,:)=F%Eo*F%Ro/F%X%R(ii)
             end do
             F%E_2D%Z=0._rp

          else
             F%E_3D%R=0._rp
             do ii=1_idef,F%dims(1)
                F%E_3D%PHI(ii,:,:)=F%Eo*F%Ro/F%X%R(ii)
             end do
             F%E_3D%Z=0._rp
          end if
       end if

       if(F%dBfield.and..not.F%dBfield_in_file) then
          if (F%axisymmetric_fields) then
             F%dBdR_2D%R=0._rp
             F%dBdR_2D%PHI=0._rp
             F%dBdR_2D%Z=0._rp

             F%dBdPHI_2D%R=0._rp
             F%dBdPHI_2D%PHI=0._rp
             F%dBdPHI_2D%Z=0._rp

             F%dBdZ_2D%R=0._rp
             F%dBdZ_2D%PHI=0._rp
             F%dBdZ_2D%Z=0._rp
          else
             F%dBdR_3D%R=0._rp
             F%dBdR_3D%PHI=0._rp
             F%dBdR_3D%Z=0._rp

             F%dBdPHI_3D%R=0._rp
             F%dBdPHI_3D%PHI=0._rp
             F%dBdPHI_3D%Z=0._rp

             F%dBdZ_3D%R=0._rp
             F%dBdZ_3D%PHI=0._rp
             F%dBdZ_3D%Z=0._rp
          end if
       end if

       if (params%mpi_params%rank .EQ. 0) then

          write(output_unit_write,'("EXTERNAL")')
          write(output_unit_write,'("Magnetic field: ",E17.10)') F%Bo
          write(output_unit_write,'("Electric field: ",E17.10)') F%Eo

       end if

       if (params%mpi_params%rank.EQ.0) then

          if (F%axisymmetric_fields) then
             if (F%Bflux) then
                write(output_unit_write,'("PSIp(r=0)",E17.10)') F%PSIp(F%dims(1)/2,F%dims(3)/2)
                write(output_unit_write,'("BPHI(r=0)",E17.10)') F%Bo
                write(output_unit_write,'("EPHI(r=0)",E17.10)') F%Eo
             else if (F%Bflux3D) then
                write(output_unit_write,'("PSIp(r=0)",E17.10)') F%PSIp3D(F%dims(1)/2,1,F%dims(3)/2)
             else
                write(output_unit_write,'("BR(r=0)",E17.10)') F%B_2D%R(F%dims(1)/2,F%dims(3)/2)
                write(output_unit_write,'("BPHI(r=0)",E17.10)') &
                     F%B_2D%PHI(F%dims(1)/2,F%dims(3)/2)
                write(output_unit_write,'("BZ(r=0)",E17.10)') F%B_2D%Z(F%dims(1)/2,F%dims(3)/2)
                write(output_unit_write,'("EPHI(r=0)",E17.10)') &
                     F%E_2D%PHI(F%dims(1)/2,F%dims(3)/2)
             end if
          end if


       end if

    else
       F%Bo=Bo
       F%Eo=Eo

       if (params%mpi_params%rank.EQ.0) then
          write(output_unit_write,'("UNIFORM")')
          write(output_unit_write,'("BPHI(r=0)",E17.10)') F%Bo
          write(output_unit_write,'("EPHI(r=0)",E17.10)') F%Eo
       end if
    end if

    if (params%mpi_params%rank.eq.0) then
       write(output_unit_write,'("* * * * * * * * * * * * * * * * * * * * * * * * *",/)')
    end if

  end subroutine initialize_field

  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  ! Subroutines for getting the fields data from HDF5 files
  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  !> @brief Subroutine that loads the size of the arrays having the electric and magnetic field data.
  !! @details All the information of externally calculated fields must be given in a rectangular, equally spaced mesh in the \((R,\phi,Z)\) space of cylindrical coordinates.
  !! If the fields are axisymmetric, then the fields must be in a rectangular mesh on the \(RZ\)-plane.
  !!
  !! @param[in] params Core KORC simulation parameters.
  !! @param[in,out] F An instance of the KORC derived type FIELDS.
  !! @param filename String containing the name of the HDF5 file.
  !! @param gname String containing the group name of a parameter in the HDF5 file.
  !! @param subgname String containing the subgroup name of a parameter in the HDF5 file.
  !! @param dset Name of data set to read from file.
  !! @param h5file_id HDF5 file identifier.
  !! @param group_id HDF5 group identifier.
  !! @param subgroup_id HDF5 subgroup identifier.
  !! @dims Array containing the size of the mesh with the data of the electric and magnetic fields. dims(1) = dimension along the \(R\) coordinate,
  !! dims(2) = dimension along the \(\phi\) coordinate, and dims(3) = dimension along the \(Z\) coordinate.
  !! @param h5error HDF5 error status.
  !! @param rdamum Temporary variable keeping the read data.
  subroutine load_dim_data_from_hdf5(params,F)
    TYPE(KORC_PARAMS), INTENT(IN)                  :: params
    TYPE(FIELDS), INTENT(INOUT)                    :: F
    CHARACTER(MAX_STRING_LENGTH)                   :: filename
    CHARACTER(MAX_STRING_LENGTH)                   :: gname
    CHARACTER(MAX_STRING_LENGTH)                   :: subgname
    CHARACTER(MAX_STRING_LENGTH)                   :: dset
    INTEGER(HID_T)                                 :: h5file_id
    INTEGER(HID_T)                                 :: group_id
    INTEGER(HID_T)                                 :: subgroup_id
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE    :: dims
    INTEGER                                        :: h5error
    REAL(rp)                                       :: rdatum

    filename = TRIM(params%magnetic_field_filename)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: load_dim_data_from_hdf5 --> h5fopen_f")')
    end if

    if (F%Bflux.OR.F%axisymmetric_fields) then
       dset = "/NR"
       call load_from_hdf5(h5file_id,dset,rdatum)
       F%dims(1) = INT(rdatum)

       F%dims(2) = 0

       dset = "/NZ"
       call load_from_hdf5(h5file_id,dset,rdatum)
       F%dims(3) = INT(rdatum)

       if(params%SC_E) then

          dset = "/OSNR"
          call load_from_hdf5(h5file_id,dset,rdatum)
          F%dims(1) = INT(rdatum)

          dset = "/OSNZ"
          call load_from_hdf5(h5file_id,dset,rdatum)
          F%dims(3) = INT(rdatum)

       end if

       if(F%Dim2x1t) then

          dset = "/NPHI"
          call load_from_hdf5(h5file_id,dset,rdatum)
          F%dims(2) = INT(rdatum)

       end if

    else
       dset = "/NR"
       call load_from_hdf5(h5file_id,dset,rdatum)
       F%dims(1) = INT(rdatum)

       dset = "/NPHI"
       call load_from_hdf5(h5file_id,dset,rdatum)
       F%dims(2) = INT(rdatum)

       dset = "/NZ"
       call load_from_hdf5(h5file_id,dset,rdatum)
       F%dims(3) = INT(rdatum)
    end if


    call h5fclose_f(h5file_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: load_dim_data_from_hdf5 --> h5fclose_f")')
    end if
  end subroutine load_dim_data_from_hdf5


  !> @brief Subroutine that queries the HDF5 file what data are present in the HDF5 input file (sanity check).
  !!
  !! @param[in] params Core KORC simulation parameters.
  !! @param Bfield Logical variable that specifies if the magnetic field is present in the HDF5 file.
  !! @param Efield Logical variable that specifies if the electric field is present in the HDF5 file.
  !! @param Bflux Logical variable that specifies if the poloidal magnetic flux is present in the HDF5 file.
  !! @param filename String containing the name of the HDF5 file.
  !! @param gname String containing the group name of a parameter in the HDF5 file.
  !! @param subgname String containing the subgroup name of a parameter in the HDF5 file.
  !! @param dset Name of data set to read from file.
  !! @param h5file_id HDF5 file identifier.
  !! @param group_id HDF5 group identifier.
  !! @param subgroup_id HDF5 subgroup identifier.
  !! @param h5error HDF5 error status.
  subroutine which_fields_in_file(params,Bfield,Efield,Bflux,B1field, &
       E1field)
    TYPE(KORC_PARAMS), INTENT(IN)      :: params
    LOGICAL, INTENT(OUT)               :: Bfield
    LOGICAL, INTENT(OUT)               :: B1field
    LOGICAL, INTENT(OUT)               :: E1field
    LOGICAL, INTENT(OUT)               :: Efield
    LOGICAL, INTENT(OUT)               :: Bflux
    CHARACTER(MAX_STRING_LENGTH)       :: filename
    CHARACTER(MAX_STRING_LENGTH)       :: gname
    CHARACTER(MAX_STRING_LENGTH)       :: subgname
    CHARACTER(MAX_STRING_LENGTH)       :: dset
    INTEGER(HID_T)                     :: h5file_id
    INTEGER(HID_T)                     :: group_id
    INTEGER(HID_T)                     :: subgroup_id
    INTEGER                            :: h5error

    filename = TRIM(params%magnetic_field_filename)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: load_field_data_from_hdf5 --> h5fopen_f")')
    end if

    gname = "BR"
    call h5lexists_f(h5file_id,TRIM(gname),Bfield,h5error)

    gname = "ER"
    call h5lexists_f(h5file_id,TRIM(gname),Efield,h5error)

    gname = "PSIp"
    call h5lexists_f(h5file_id,TRIM(gname),Bflux,h5error)

    gname = "ReBZ"
    call h5lexists_f(h5file_id,TRIM(gname),B1field,h5error)

    gname = "ReEZ"
    call h5lexists_f(h5file_id,TRIM(gname),E1field,h5error)


    call h5fclose_f(h5file_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: load_field_data_from_hdf5 --> h5fclose_f")')
    end if
  end subroutine which_fields_in_file


  !> @brief Subroutine that loads the fields data from the HDF5 input file.
  !!
  !! @param[in] params Core KORC simulation parameters.
  !! @param[in,out] F An instance of the KORC derived type FIELDS. In this variable we keep the loaded data.
  !! @param filename String containing the name of the HDF5 file.
  !! @param gname String containing the group name of a parameter in the HDF5 file.
  !! @param subgname String containing the subgroup name of a parameter in the HDF5 file.
  !! @param dset Name of data set to read from file.
  !! @param h5file_id HDF5 file identifier.
  !! @param group_id HDF5 group identifier.
  !! @param subgroup_id HDF5 subgroup identifier.
  !! @param h5error HDF5 error status.
   subroutine load_field_data_from_hdf5(params,F)
      TYPE(KORC_PARAMS), INTENT(IN)          :: params
      TYPE(FIELDS), INTENT(INOUT)            :: F
      CHARACTER(MAX_STRING_LENGTH)           :: filename
      CHARACTER(MAX_STRING_LENGTH)           :: gname
      CHARACTER(MAX_STRING_LENGTH)           :: subgname
      CHARACTER(MAX_STRING_LENGTH)           :: dset
      INTEGER(HID_T)                         :: h5file_id
      INTEGER(HID_T)                         :: group_id
      INTEGER(HID_T)                         :: subgroup_id
      INTEGER                                :: h5error
      LOGICAL :: Efield

      filename = TRIM(params%magnetic_field_filename)
      call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
      if (h5error .EQ. -1) then
         write(output_unit_write,'("KORC ERROR: Something went wrong in: load_field_data_from_hdf5 --> h5fopen_f")')
      end if


      if (((.NOT.F%Bflux).AND.(.NOT.F%axisymmetric_fields)).OR. &
            F%Dim2x1t) then
         dset = "/PHI"
         call load_array_from_hdf5(h5file_id,dset,F%X%PHI)
      end if



      dset = "/R"
      call load_array_from_hdf5(h5file_id,dset,F%X%R)

      dset = "/Z"
      call load_array_from_hdf5(h5file_id,dset,F%X%Z)


      dset = '/Bo'
      call load_from_hdf5(h5file_id,dset,F%Bo)

      if (F%Efield) then
         dset = '/Eo'
         gname = 'Eo'

         call h5lexists_f(h5file_id,TRIM(gname),Efield,h5error)

         if (Efield) then
            call load_from_hdf5(h5file_id,dset,F%Eo)
         else
            F%Eo = 0.0_rp
         end if
      else
         F%Eo = 0.0_rp
      end if

      if (params%field_model(10:13).eq.'MARS') then

         dset = '/PSIP0'
         call load_from_hdf5(h5file_id,dset,F%PSIP_min)

         dset = '/PSIPlim'
         call load_from_hdf5(h5file_id,dset,F%PSIp_lim)

         dset = '/AMP'
         call load_from_hdf5(h5file_id,dset,F%AMP)

         F%AMP=F%AMP*F%MARS_AMP_Scale

      end if

      if (params%field_model(10:14).eq.'AORSA') then

         dset = '/PSIP0'
         call load_from_hdf5(h5file_id,dset,F%PSIP_min)

         dset = '/PSIPlim'
         call load_from_hdf5(h5file_id,dset,F%PSIp_lim)

         F%AMP=F%AORSA_AMP_Scale

      end if

      dset = '/Ro'
      call load_from_hdf5(h5file_id,dset,F%Ro)

      dset = '/Zo'
      call load_from_hdf5(h5file_id,dset,F%Zo)

      if ((F%Bflux.OR.F%axisymmetric_fields).AND.(.NOT.F%Dim2x1t)) then

         dset = "/FLAG"
         call load_array_from_hdf5(h5file_id,dset,F%FLAG2D)

         if (F%useLCFS) then
            dset = "/LCFS"
            call load_array_from_hdf5(h5file_id,dset,F%LCFS2D)
         else
            F%LCFS2D = 0._rp
         end if

      else
         dset = "/FLAG"
         call load_array_from_hdf5(h5file_id,dset,F%FLAG3D)

         if (F%useLCFS) then
            dset = "/LCFS"
            call load_array_from_hdf5(h5file_id,dset,F%LCFS3D)
         else
            F%LCFS3D = 0._rp
         end if
      end if

      if (F%Bflux) then

         dset = "/PSIp"
         gname = 'PSIp'

         call h5lexists_f(h5file_id,TRIM(gname),Efield,h5error)

         if (Efield) then
            call load_array_from_hdf5(h5file_id,dset,F%PSIp)
         else
            F%PSIp = 0.0_rp
         end if

            !          F%PSIp=2*C_PI*(F%PSIp-minval(F%PSIp))

      end if

      if (F%Bflux3D) then
         dset = "/PSIp"
         gname = 'PSIp'

         call h5lexists_f(h5file_id,TRIM(gname),Efield,h5error)

         if (Efield) then
            call load_array_from_hdf5(h5file_id,dset,F%PSIp3D)
         else
            F%PSIp3D = 0.0_rp
         end if

   !       F%PSIp3D=2*C_PI*(F%PSIp3D-minval(F%PSIp3D))

      end if

      if (F%B1field) then

         if (params%field_model(10:13).eq.'MARS') then

            dset = "/ReBR"
            call load_array_from_hdf5(h5file_id,dset,F%B1Re_2D%R)

            dset = "/ReBPHI"
            call load_array_from_hdf5(h5file_id,dset,F%B1Re_2D%PHI)

            dset = "/ReBZ"
            call load_array_from_hdf5(h5file_id,dset,F%B1Re_2D%Z)

            dset = "/ImBR"
            call load_array_from_hdf5(h5file_id,dset,F%B1Im_2D%R)

            dset = "/ImBPHI"
            call load_array_from_hdf5(h5file_id,dset,F%B1Im_2D%PHI)

            dset = "/ImBZ"
            call load_array_from_hdf5(h5file_id,dset,F%B1Im_2D%Z)

         else if (params%field_model(10:14).eq.'AORSA') then

            dset = "/ReBX"
            call load_array_from_hdf5(h5file_id,dset,F%B1Re_2DX%X)

            dset = "/ReBY"
            call load_array_from_hdf5(h5file_id,dset,F%B1Re_2DX%Y)

            dset = "/ReBZ"
            call load_array_from_hdf5(h5file_id,dset,F%B1Re_2DX%Z)

            dset = "/ImBX"
            call load_array_from_hdf5(h5file_id,dset,F%B1Im_2DX%X)

            dset = "/ImBY"
            call load_array_from_hdf5(h5file_id,dset,F%B1Im_2DX%Y)

            dset = "/ImBZ"
            call load_array_from_hdf5(h5file_id,dset,F%B1Im_2DX%Z)

         endif

      end if

      if (F%E1field) then

         dset = "/ReEX"
         call load_array_from_hdf5(h5file_id,dset,F%E1Re_2DX%X)

         dset = "/ReEY"
         call load_array_from_hdf5(h5file_id,dset,F%E1Re_2DX%Y)

         dset = "/ReEZ"
         call load_array_from_hdf5(h5file_id,dset,F%E1Re_2DX%Z)

         dset = "/ImEX"
         call load_array_from_hdf5(h5file_id,dset,F%E1Im_2DX%X)

         dset = "/ImEY"
         call load_array_from_hdf5(h5file_id,dset,F%E1Im_2DX%Y)

         dset = "/ImEZ"
         call load_array_from_hdf5(h5file_id,dset,F%E1Im_2DX%Z)

      end if

      if (F%Bfield) then
         if (F%axisymmetric_fields) then
            dset = "/BR"
            call load_array_from_hdf5(h5file_id,dset,F%B_2D%R)

            dset = "/BPHI"
            call load_array_from_hdf5(h5file_id,dset,F%B_2D%PHI)

            dset = "/BZ"
            call load_array_from_hdf5(h5file_id,dset,F%B_2D%Z)
         else
            dset = "/BR"
            call load_array_from_hdf5(h5file_id,dset,F%B_3D%R)

            !write(6,*) 'BR(25,1,:)',F%B_3D%R(25,1,:)

            dset = "/BPHI"
            call load_array_from_hdf5(h5file_id,dset,F%B_3D%PHI)

            dset = "/BZ"
            call load_array_from_hdf5(h5file_id,dset,F%B_3D%Z)
         end if
      end if

      if (F%Efield.AND.F%Efield_in_file) then
         if (F%axisymmetric_fields.and.(.not.F%Dim2x1t)) then
            dset = "/ER"
            call load_array_from_hdf5(h5file_id,dset,F%E_2D%R)

            dset = "/EPHI"
            call load_array_from_hdf5(h5file_id,dset,F%E_2D%PHI)

            dset = "/EZ"
            call load_array_from_hdf5(h5file_id,dset,F%E_2D%Z)
         else

            dset = "/ER"
            call load_array_from_hdf5(h5file_id,dset,F%E_3D%R)

            dset = "/EPHI"
            call load_array_from_hdf5(h5file_id,dset,F%E_3D%PHI)

            dset = "/EZ"
            call load_array_from_hdf5(h5file_id,dset,F%E_3D%Z)
         end if
      end if

      call h5fclose_f(h5file_id, h5error)
      if (h5error .EQ. -1) then
         write(output_unit_write,'("KORC ERROR: Something went wrong in: load_field_data_from_hdf5 --> h5fclose_f")')
      end if
  end subroutine load_field_data_from_hdf5

  subroutine ALLOCATE_2D_FIELDS_ARRAYS(params,F,bfield,bflux,dbfield, &
       efield,b1field,e1field)
    !! @note Subroutine that allocates the variables keeping the axisymmetric
    !! fields data. @endnote
    TYPE (KORC_PARAMS), INTENT(IN) 	:: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(INOUT)    :: F
    !! An instance of the KORC derived type FIELDS. In this variable we keep
    !! the loaded data.
    LOGICAL, INTENT(IN)            :: bfield
    LOGICAL, INTENT(IN)            :: b1field
    LOGICAL, INTENT(IN)            :: e1field
    LOGICAL, INTENT(IN)            :: dbfield
    !! Logical variable that specifies if the variables that keep the magnetic
    !! field data is allocated (bfield=T) or not (bfield=F).
    LOGICAL, INTENT(IN)            :: bflux
    !! Logical variable that specifies if the variables that keep the poloidal
    !! magnetic flux data is allocated (bflux=T) or not (bflux=F).
    LOGICAL, INTENT(IN)            :: efield
    !! Logical variable that specifies if the variables that keep the electric
    !! field data is allocated (efield=T) or not (efield=F).

    if (bfield.and.(.not.ALLOCATED(F%B_2D%R))) then
       call ALLOCATE_V_FIELD_2D(F%B_2D,F%dims)

       if(params%orbit_model(3:5).EQ.'pre') then
          call ALLOCATE_V_FIELD_2D(F%curlb_2D,F%dims)
          call ALLOCATE_V_FIELD_2D(F%gradB_2D,F%dims)
       end if

    end if

    if (bflux.and.(.not.ALLOCATED(F%PSIp))) then
       ALLOCATE(F%PSIp(F%dims(1),F%dims(3)))
       F%PSIp=0._rp
    end if

    if (dbfield.and.(.not.ALLOCATED(F%dBdR_2D%R))) then
       call ALLOCATE_V_FIELD_2D(F%dBdR_2D,F%dims)
       call ALLOCATE_V_FIELD_2D(F%dBdPHI_2D,F%dims)
       call ALLOCATE_V_FIELD_2D(F%dBdZ_2D,F%dims)
    end if

    if (params%field_model(10:13).eq.'MARS') then

       if (B1field.and.(.not.ALLOCATED(F%B1Re_2D%R))) then
          call ALLOCATE_V_FIELD_2D(F%B1Re_2D,F%dims)
          call ALLOCATE_V_FIELD_2D(F%B1Im_2D,F%dims)
       end if

    else if (params%field_model(10:14).eq.'AORSA') then

       if (B1field.and.(.not.ALLOCATED(F%B1Re_2DX%X))) then
          call ALLOCATE_V_FIELD_2DX(F%B1Re_2DX,F%dims)
          call ALLOCATE_V_FIELD_2DX(F%B1Im_2DX,F%dims)
       end if

       if (E1field.and.(.not.ALLOCATED(F%E1Re_2DX%X))) then
          call ALLOCATE_V_FIELD_2DX(F%E1Re_2DX,F%dims)
          call ALLOCATE_V_FIELD_2DX(F%E1Im_2DX,F%dims)
       end if

    endif



    if (efield.and.(.not.ALLOCATED(F%E_2D%R))) then
       call ALLOCATE_V_FIELD_2D(F%E_2D,F%dims)

    end if

    if (.NOT.ALLOCATED(F%FLAG2D)) ALLOCATE(F%FLAG2D(F%dims(1),F%dims(3)))
    if (.NOT.ALLOCATED(F%LCFS2D)) ALLOCATE(F%LCFS2D(F%dims(1),F%dims(3)))

    if (.NOT.ALLOCATED(F%X%R)) ALLOCATE(F%X%R(F%dims(1)))
    if (.NOT.ALLOCATED(F%X%Z)) ALLOCATE(F%X%Z(F%dims(3)))
  end subroutine ALLOCATE_2D_FIELDS_ARRAYS


  !> @brief Subroutine that allocates the variables keeping the 3-D fields data.
  !!
  !! @param[in,out] F An instance of the KORC derived type FIELDS. In this variable we keep the loaded data.
  !! @param[in] bfield Logical variable that specifies if the variables that keep the magnetic field data is allocated (bfield=T) or not (bfield=F).
  !! @param[in] efield Logical variable that specifies if the variables that keep the electric field data is allocated (efield=T) or not (efield=F).
  subroutine ALLOCATE_3D_FIELDS_ARRAYS(params,F,bfield,efield,dbfield)
    TYPE (KORC_PARAMS), INTENT(IN) 	:: params
    TYPE(FIELDS), INTENT(INOUT)    :: F
    LOGICAL, INTENT(IN)            :: bfield
    LOGICAL, INTENT(IN)            :: dbfield
    LOGICAL, INTENT(IN)            :: efield

    if (bfield) then
       call ALLOCATE_V_FIELD_3D(F%B_3D,F%dims)

       if(params%orbit_model(3:5).EQ.'pre') then
          call ALLOCATE_V_FIELD_3D(F%curlb_3D,F%dims)
          call ALLOCATE_V_FIELD_3D(F%gradB_3D,F%dims)
       end if

    end if

    if (F%Bflux3D.and.(.not.ALLOCATED(F%PSIp3D))) then
       ALLOCATE(F%PSIp3D(F%dims(1),F%dims(2),F%dims(3)))
    end if

    if (dbfield.and.(.not.ALLOCATED(F%dBdR_3D%R))) then
       call ALLOCATE_V_FIELD_3D(F%dBdR_3D,F%dims)
       call ALLOCATE_V_FIELD_3D(F%dBdPHI_3D,F%dims)
       call ALLOCATE_V_FIELD_3D(F%dBdZ_3D,F%dims)
    end if

    if (efield) then
       call ALLOCATE_V_FIELD_3D(F%E_3D,F%dims)
    end if

    if (.NOT.ALLOCATED(F%FLAG3D)) ALLOCATE(F%FLAG3D(F%dims(1),F%dims(2),F%dims(3)))
    if (.NOT.ALLOCATED(F%LCFS3D)) ALLOCATE(F%LCFS3D(F%dims(1),F%dims(2),F%dims(3)))

    if (.NOT.ALLOCATED(F%X%R)) ALLOCATE(F%X%R(F%dims(1)))
    if (.NOT.ALLOCATED(F%X%PHI)) ALLOCATE(F%X%PHI(F%dims(2)))
    if (.NOT.ALLOCATED(F%X%Z)) ALLOCATE(F%X%Z(F%dims(3)))
  end subroutine ALLOCATE_3D_FIELDS_ARRAYS


  !> @brief Subroutine that allocates the cylindrical components of an axisymmetric field.
  !!
  !! @param[in,out] F Vector field to be allocated.
  !! @param[in] dims Dimension of the mesh containing the field data.
  subroutine ALLOCATE_V_FIELD_2D(F,dims)
    TYPE(V_FIELD_2D), INTENT(INOUT)    :: F
    INTEGER, DIMENSION(3), INTENT(IN)  :: dims

    ALLOCATE(F%R(dims(1),dims(3)))
    ALLOCATE(F%PHI(dims(1),dims(3)))
    ALLOCATE(F%Z(dims(1),dims(3)))
  end subroutine ALLOCATE_V_FIELD_2D

    !> @brief Subroutine that allocates the cartesian components of an axisymmetric field.
  !!
  !! @param[in,out] F Vector field to be allocated.
  !! @param[in] dims Dimension of the mesh containing the field data.
  subroutine ALLOCATE_V_FIELD_2DX(F,dims)
    TYPE(V_FIELD_2DX), INTENT(INOUT)    :: F
    INTEGER, DIMENSION(3), INTENT(IN)  :: dims

    ALLOCATE(F%X(dims(1),dims(3)))
    ALLOCATE(F%Y(dims(1),dims(3)))
    ALLOCATE(F%Z(dims(1),dims(3)))
  end subroutine ALLOCATE_V_FIELD_2DX


  !> @brief Subroutine that allocates the cylindrical components of a 3-D field.
  !!
  !! @param[in,out] F Vector field to be allocated.
  !! @param[in] dims Dimension of the mesh containing the field data.
  subroutine ALLOCATE_V_FIELD_3D(F,dims)
    TYPE(V_FIELD_3D), INTENT(INOUT)    :: F
    INTEGER, DIMENSION(3), INTENT(IN)  :: dims

    ALLOCATE(F%R(dims(1),dims(2),dims(3)))
    ALLOCATE(F%PHI(dims(1),dims(2),dims(3)))
    ALLOCATE(F%Z(dims(1),dims(2),dims(3)))
  end subroutine ALLOCATE_V_FIELD_3D

  !> @brief Subroutine that deallocates all the variables of the electric and magnetic fields.
  !!
  !! @param[in,out] F An instance of the KORC derived type FIELDS.
  subroutine DEALLOCATE_FIELDS_ARRAYS(F)
    TYPE(FIELDS), INTENT(INOUT) :: F

    if (ALLOCATED(F%PSIp)) DEALLOCATE(F%PSIp)

    if (ALLOCATED(F%B_2D%R)) DEALLOCATE(F%B_2D%R)
    if (ALLOCATED(F%B_2D%PHI)) DEALLOCATE(F%B_2D%PHI)
    if (ALLOCATED(F%B_2D%Z)) DEALLOCATE(F%B_2D%Z)

    if (ALLOCATED(F%gradB_2D%R)) DEALLOCATE(F%gradB_2D%R)
    if (ALLOCATED(F%gradB_2D%PHI)) DEALLOCATE(F%gradB_2D%PHI)
    if (ALLOCATED(F%gradB_2D%Z)) DEALLOCATE(F%gradB_2D%Z)

    if (ALLOCATED(F%curlb_2D%R)) DEALLOCATE(F%curlb_2D%R)
    if (ALLOCATED(F%curlb_2D%PHI)) DEALLOCATE(F%curlb_2D%PHI)
    if (ALLOCATED(F%curlb_2D%Z)) DEALLOCATE(F%curlb_2D%Z)

    if (ALLOCATED(F%B_3D%R)) DEALLOCATE(F%B_3D%R)
    if (ALLOCATED(F%B_3D%PHI)) DEALLOCATE(F%B_3D%PHI)
    if (ALLOCATED(F%B_3D%Z)) DEALLOCATE(F%B_3D%Z)

    if (ALLOCATED(F%E_2D%R)) DEALLOCATE(F%E_2D%R)
    if (ALLOCATED(F%E_2D%PHI)) DEALLOCATE(F%E_2D%PHI)
    if (ALLOCATED(F%E_2D%Z)) DEALLOCATE(F%E_2D%Z)

    if (ALLOCATED(F%E_3D%R)) DEALLOCATE(F%E_3D%R)
    if (ALLOCATED(F%E_3D%PHI)) DEALLOCATE(F%E_3D%PHI)
    if (ALLOCATED(F%E_3D%Z)) DEALLOCATE(F%E_3D%Z)

    if (ALLOCATED(F%X%R)) DEALLOCATE(F%X%R)
    if (ALLOCATED(F%X%PHI)) DEALLOCATE(F%X%PHI)
    if (ALLOCATED(F%X%Z)) DEALLOCATE(F%X%Z)

    if (ALLOCATED(F%FLAG2D)) DEALLOCATE(F%FLAG2D)
    if (ALLOCATED(F%FLAG3D)) DEALLOCATE(F%FLAG3D)

    if (ALLOCATED(F%LCFS2D)) DEALLOCATE(F%LCFS2D)
    if (ALLOCATED(F%LCFS3D)) DEALLOCATE(F%LCFS3D)
  end subroutine DEALLOCATE_FIELDS_ARRAYS
end module korc_fields
