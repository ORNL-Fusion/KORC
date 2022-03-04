
module korc_ppusher
  !! @note Module with subroutines for advancing the particles' position and
  !! velocity in the simulations. @endnote
  use korc_types
  use korc_constants
  use korc_fields
  use korc_profiles
  use korc_interp
  use korc_collisions
  use korc_hpc

  IMPLICIT NONE

  REAL(rp), PRIVATE :: E0
  !! Dimensionless vacuum permittivity \(\epsilon_0 \times (m_{ch}^2
  !! v_{ch}^3/q_{ch}^3 B_{ch})\), see [[korc_units]].

  PRIVATE :: cross,&
       radiation_force_p,&
       GCEoM_p,&
       GCEoM1_p,&
       aux_fields
  PUBLIC :: initialize_particle_pusher,&
       GC_init,&
       FO_init,&
#ifdef FIO
       adv_FOfio_top,&
#endif
#ifdef PSPLINE
       adv_FOinterp_top,&
       adv_FOinterp_mars_top,&
       adv_FOinterp_aorsa_top,&
       advance_FOinterp_vars,&
       advance_GCinterp_psi_vars,&
       advance_GCinterp_B_vars,&
       adv_GCinterp_psi_top,&
       adv_GCinterp_psiwE_top,&
       adv_GCinterp_psi2x1t_top,&
       adv_GCinterp_psi_top_FS,&
       adv_GCinterp_B_top,&
       adv_GCinterp_B2D_top,&
       adv_GCinterp_2DBdB_top,&
       adv_GCinterp_3DBdB_top,&
       adv_GCinterp_3DBdB1_top,&
#endif
       adv_FOeqn_top,&
       advance_FOeqn_vars,&
       adv_GCeqn_top,&
       advance_GCeqn_vars

contains



  subroutine initialize_particle_pusher(params)
    !! @note This subroutine initializes all the variables needed for advancing
    !! the particles' position and velocity. @endnote
    !! This subroutine is specially useful when we need to define or initialize
    !! values of parameters used to calculate derived quantities.
    !! The intent of this subroutine is to work as a constructor of the module.
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    !! Core KORC simulation parameters.

    E0 = C_E0*(params%cpp%mass**2*params%cpp%velocity**3)/ &
         (params%cpp%charge**3*params%cpp%Bo)

  end subroutine initialize_particle_pusher


  pure function cross(a,b)
    !! @note Function that calculates and returns the cross product
    !! \(\mathbf{a}\times \mathbf{b}\). These vectors are in Cartesian
    !! coordinates. @endnote
    !! @note Notice that all the variables in this subroutine have been
    !! normalized using the characteristic scales in [[korc_units]]. @endnote
    REAL(rp), DIMENSION(3), INTENT(IN) :: a
    !! Vector \(\mathbf{a}\).
    REAL(rp), DIMENSION(3), INTENT(IN) :: b
    !! Vector \(\mathbf{b}\).
    REAL(rp), DIMENSION(3)             :: cross
    !!Value of \(\mathbf{a}\times \mathbf{b}\)

    cross(1) = a(2)*b(3) - a(3)*b(2)
    cross(2) = a(3)*b(1) - a(1)*b(3)
    cross(3) = a(1)*b(2) - a(2)*b(1)
  end function cross



  subroutine radiation_force_p(pchunk,q_cache,m_cache,U_X,U_Y,U_Z,E_X,E_Y,E_Z, &
       B_X,B_Y,B_Z,Frad_X,Frad_Y,Frad_Z)
    INTEGER, INTENT(IN)  :: pchunk
    REAL(rp), INTENT(IN)                       :: m_cache,q_cache

    REAL(rp), DIMENSION(pchunk), INTENT(IN)     :: U_X,U_Y,U_Z
    !! \(\mathbf{u} = \gamma \mathbf{v}\), where \(\mathbf{v}\) is the
    !! particle's velocity.
    REAL(rp), DIMENSION(pchunk), INTENT(IN)     :: E_X,E_Y,E_Z
    !! Electric field \(\mathbf{E}\) seen by each particle. This is given
    !! in Cartesian coordinates.
    REAL(rp), DIMENSION(pchunk), INTENT(IN)     :: B_X,B_Y,B_Z
    !! Magnetic field \(\mathbf{B}\) seen by each particle. This is given
    !! in Cartesian coordinates.
    REAL(rp), DIMENSION(pchunk), INTENT(OUT)    :: Frad_X,Frad_Y,Frad_Z
    !! The calculated synchrotron radiation reaction force \(\mathbf{F}_R\).
    REAL(rp), DIMENSION(3)                 :: F1
    !! The component \(\mathbf{F}_1\) of \(\mathbf{F}_R\).
    REAL(rp), DIMENSION(pchunk)                 :: F2_X,F2_Y,F2_Z
    !! The component \(\mathbf{F}_2\) of \(\mathbf{F}_R\).
    REAL(rp), DIMENSION(pchunk)                 :: F3_X,F3_Y,F3_Z
    !! The component \(\mathbf{F}_3\) of \(\mathbf{F}_R\).
    REAL(rp), DIMENSION(pchunk)                 :: V_X,V_Y,V_Z
    !! The particle's velocity \(\mathbf{v}\).
    REAL(rp), DIMENSION(pchunk)                 :: vec_X,vec_Y,vec_Z
    REAL(rp), DIMENSION(pchunk)                 :: cross_EB_X,cross_EB_Y,cross_EB_Z
    REAL(rp), DIMENSION(pchunk)                 :: cross_BV_X,cross_BV_Y,cross_BV_Z
    REAL(rp), DIMENSION(pchunk)                 :: cross_BBV_X,cross_BBV_Y,cross_BBV_Z
    REAL(rp), DIMENSION(pchunk)                 :: dot_EV,dot_vecvec
    !! An auxiliary 3-D vector.
    REAL(rp),DIMENSION(pchunk)                               :: g
    !! The relativistic \(\gamma\) factor of the particle.
    REAL(rp)                               :: tmp
    INTEGER :: cc

    !$OMP SIMD
    !    !$OMP& aligned(g,U_X,U_Y,U_Z,V_X,V_Y,V_Z, &
    !    !$OMP& cross_EB_X,cross_EB_Y,cross_EB_Z,E_X,E_Y,E_Z,B_X,B_Y,B_Z, &
    !    !$OMP& dot_EV,cross_BV_X,cross_BV_Y,cross_BV_Z, &
    !    !$OMP& cross_BBV_X,cross_BBV_Y,cross_BBV_Z,F2_X,F2_Y,F2_Z, &
    !    !$OMP& vec_X,vec_Y,vec_Z,dot_vecvec,F3_X,F3_Y,F3_Z, &
    !    !$OMP& Frad_X,Frad_Y,Frad_Z)
    do cc=1_idef,pchunk
       g(cc) = SQRT(1.0_rp + U_X(cc)*U_X(cc)+ U_Y(cc)*U_Y(cc)+ U_Z(cc)*U_Z(cc))

       V_X(cc) = U_X(cc)/g(cc)
       V_Y(cc) = U_Y(cc)/g(cc)
       V_Z(cc) = U_Z(cc)/g(cc)

       tmp = q_cache**4/(6.0_rp*C_PI*E0*m_cache**2)

       cross_EB_X(cc)=E_Y(cc)*B_Z(cc)-E_Z(cc)*B_Y(cc)
       cross_EB_Y(cc)=E_Z(cc)*B_X(cc)-E_X(cc)*B_Z(cc)
       cross_EB_Z(cc)=E_X(cc)*B_Y(cc)-E_Y(cc)*B_X(cc)

       dot_EV(cc)=E_X(cc)*V_X(cc)+E_Y(cc)*V_Y(cc)+E_Z(cc)*V_Z(cc)

       cross_BV_X(cc)=B_Y(cc)*V_Z(cc)-B_Z(cc)*V_Y(cc)
       cross_BV_Y(cc)=B_Z(cc)*V_X(cc)-B_X(cc)*V_Z(cc)
       cross_BV_Z(cc)=B_X(cc)*V_Y(cc)-B_Y(cc)*V_X(cc)

       cross_BBV_X(cc)=B_Y(cc)*cross_BV_Z(cc)-B_Z(cc)*cross_BV_Y(cc)
       cross_BBV_Y(cc)=B_Z(cc)*cross_BV_X(cc)-B_X(cc)*cross_BV_Z(cc)
       cross_BBV_Z(cc)=B_X(cc)*cross_BV_Y(cc)-B_Y(cc)*cross_BV_X(cc)

       F2_X(cc) = tmp*( dot_EV(cc)*E_X(cc) + cross_EB_X(cc) + cross_BBV_X(cc) )
       F2_Y(cc) = tmp*( dot_EV(cc)*E_Y(cc) + cross_EB_Y(cc) + cross_BBV_Y(cc) )
       F2_Z(cc) = tmp*( dot_EV(cc)*E_Z(cc) + cross_EB_Z(cc) + cross_BBV_Z(cc) )

       vec_X(cc) = E_X(cc) - cross_BV_X(cc)
       vec_Y(cc) = E_Y(cc) - cross_BV_Y(cc)
       vec_Z(cc) = E_Z(cc) - cross_BV_Z(cc)

       dot_vecvec(cc)=vec_X(cc)*vec_X(cc)+vec_Y(cc)*vec_Y(cc)+vec_Z(cc)*vec_Z(cc)

       F3_X(cc) = (tmp*g(cc)**2)*( dot_EV(cc)**2 - dot_vecvec(cc) )*V_X(cc)
       F3_Y(cc) = (tmp*g(cc)**2)*( dot_EV(cc)**2 - dot_vecvec(cc) )*V_Y(cc)
       F3_Z(cc) = (tmp*g(cc)**2)*( dot_EV(cc)**2 - dot_vecvec(cc) )*V_Z(cc)

       Frad_X(cc) = F2_X(cc) + F3_X(cc)
       Frad_Y(cc) = F2_Y(cc) + F3_Y(cc)
       Frad_Z(cc) = F2_Z(cc) + F3_Z(cc)

    end do
    !$OMP END SIMD

  end subroutine radiation_force_p




  subroutine FO_init(params,F,spp,output,step)
    TYPE(KORC_PARAMS), INTENT(IN)                              :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                                   :: F
    !! An instance of the KORC derived type FIELDS.
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.

    REAL(rp)                                      :: dt
    !! Time step used in the leapfrog step (\(\Delta t\)).

    REAL(rp), DIMENSION(params%pchunk)               :: Bmag
    REAL(rp), DIMENSION(params%pchunk)               :: b_unit_X,b_unit_Y,b_unit_Z
    REAL(rp), DIMENSION(params%pchunk)               :: v,vpar,vperp
    REAL(rp), DIMENSION(params%pchunk)               :: tmp
    REAL(rp), DIMENSION(params%pchunk)               :: g
    REAL(rp), DIMENSION(params%pchunk)               :: cross_X,cross_Y,cross_Z
    REAL(rp), DIMENSION(params%pchunk)               :: vec_X,vec_Y,vec_Z
    INTEGER                                      :: ii
    !! Species iterator.
    INTEGER                                      :: pp
    !! Particles iterator.
    INTEGER                                      :: cc,pchunk
    !! Chunk iterator.

    LOGICAL,intent(in) :: output
    LOGICAL,intent(in) :: step   

    REAL(rp),DIMENSION(params%pchunk) :: X_X,X_Y,X_Z
    REAL(rp),DIMENSION(params%pchunk) :: V_X,V_Y,V_Z
    REAL(rp),DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: B_X,B_Y,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_X,E_Y,E_Z
    REAL(rp),DIMENSION(params%pchunk) :: PSIp
    REAL(rp) :: m_cache,q_cache,B0,EF0,lam,R0,q0,ar
    TYPE(C_PTR),DIMENSION(params%pchunk) :: hint
    INTEGER(is) ,DIMENSION(params%pchunk) :: flagCon,flagCol

    pchunk=params%pchunk

    B0=F%Bo
    EF0=F%Eo
    lam=F%AB%lambda
    R0=F%AB%Ro
    q0=F%AB%qo
    ar=F%AB%a

    do ii = 1_idef,params%num_species
       
       m_cache=spp(ii)%m
       q_cache=spp(ii)%q
       
       if(output) then

          !$OMP PARALLEL DO default(none) &
          !$OMP firstprivate(E0,m_cache,q_cache,B0,EF0,lam,R0,q0,ar,pchunk) &
          !$OMP& shared(params,ii,spp,F) &
          !$OMP& PRIVATE(pp,cc,X_X,X_Y,X_Z,B_X,B_Y,B_Z,V_X,V_Y,V_Z, &
          !$OMP& E_X,E_Y,E_Z,Y_R,Y_PHI,Y_Z,flagCon,flagCol,PSIp,hint,Bmag, &
          !$OMP& b_unit_X,b_unit_Y,b_unit_Z,v,vpar,vperp,tmp, &
          !$OMP& cross_X,cross_Y,cross_Z,vec_X,vec_Y,vec_Z,g)
          do pp=1_idef,spp(ii)%ppp,pchunk

             !$OMP SIMD
             do cc=1_idef,pchunk
                X_X(cc)=spp(ii)%vars%X(pp-1+cc,1)
                X_Y(cc)=spp(ii)%vars%X(pp-1+cc,2)
                X_Z(cc)=spp(ii)%vars%X(pp-1+cc,3)

                V_X(cc)=spp(ii)%vars%V(pp-1+cc,1)
                V_Y(cc)=spp(ii)%vars%V(pp-1+cc,2)
                V_Z(cc)=spp(ii)%vars%V(pp-1+cc,3)

                Y_R(cc)=0._rp
                Y_PHI(cc)=0._rp
                Y_Z(cc)=0._rp
                
                B_X(cc)=0._rp
                B_Y(cc)=0._rp
                B_Z(cc)=0._rp

                E_X(cc)=0._rp
                E_Y(cc)=0._rp
                E_Z(cc)=0._rp

                PSIp(cc)=100._rp
                
                flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
                flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)
             end do
             !$OMP END SIMD

#ifdef FIO
             if (TRIM(params%field_model).eq.'M3D_C1'.or. &
                  (TRIM(params%field_model).eq.'NIMROD')) then
                !$OMP SIMD
                do cc=1_idef,pchunk
                   hint(cc)=spp(ii)%vars%hint(pp-1+cc)
                end do
                !$OMP END SIMD
             end if
#endif

             call cart_to_cyl_p(pchunk,X_X,X_Y,X_Z,Y_R,Y_PHI,Y_Z)             

             if (params%field_model(1:3).eq.'ANA') then
                call analytical_fields_p(params,pchunk,F, &
                     X_X,X_Y,X_Z,B_X,B_Y,B_Z,E_X,E_Y,E_Z,flagCon)
             else if (params%field_model(1:3).eq.'UNI') then
                call uniform_fields_p(pchunk,F, &
                    B_X,B_Y,B_Z,E_X,E_Y,E_Z)
#ifdef PSPLINE
             else if (F%axisymmetric_fields.and. &
                  (params%orbit_model(3:3).eq.'B')) then
                call interp_FOfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_X,B_Y,B_Z, &
                     E_X,E_Y,E_Z,PSIp,flagCon)
             else if ((.not.F%axisymmetric_fields).and. &
                  (params%orbit_model(3:3).eq.'B')) then
                call interp_FO3Dfields_p(pchunk,F,Y_R,Y_PHI,Y_Z, &
                     B_X,B_Y,B_Z,E_X,E_Y,E_Z,flagCon)
             else if (params%orbit_model(3:5).eq.'psi') then
                call interp_FOfields1_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_X,B_Y,B_Z, &
                     E_X,E_Y,E_Z,PSIp,flagCon)
             else if (params%field_model(10:13).eq.'MARS') then
                call interp_FOfields_mars_p(params,pchunk,F,Y_R,Y_PHI,Y_Z, &
                     B_X,B_Y,B_Z,PSIp,flagCon)
             else if (params%field_model(10:14).eq.'AORSA') then
                call interp_FOfields_aorsa_p(0._rp,params,pchunk,F, &
                     Y_R,Y_PHI,Y_Z,B_X,B_Y,B_Z,E_X,E_Y,E_Z,PSIp,flagCon)
#endif
#ifdef FIO
             else if (TRIM(params%field_model).eq.'M3D_C1'.or. &
                  TRIM(params%field_model).eq.'NIMROD') then
                call get_fio_FOmagnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
                     B_X,B_Y,B_Z,flagCon,hint)
                if (F%FIO_E .ge. 0) then
                   call get_fio_FOelectric_fields_p(params,F, &
                        Y_R,Y_PHI,Y_Z,E_X,E_Y,E_Z,flagCon,hint)
                end if
                if (F%FIO_A .ge. 0) then
                   call get_fio_vector_potential_p(params,F,Y_R,Y_PHI,Y_Z, &
                        PSIp,flagCon,hint)
                end if
#endif                
             end if

             !write(6,'("Y_R: ",E17.10)') Y_R*params%cpp%length
             !write(6,'("Y_PHI: ",E17.10)') Y_PHI
             !write(6,'("Y_Z: ",E17.10)') Y_Z*params%cpp%length
             !write(6,*) 'r_sam',sqrt((Y_R-spp(ii)%Ro)**2+Y_Z**2)*params%cpp%length
             
             !write(6,'("B_X: ",E17.10)') B_X*params%cpp%Bo
             !write(6,'("B_Y: ",E17.10)') B_Y*params%cpp%Bo
             !write(6,'("B_Z: ",E17.10)') B_Z*params%cpp%Bo
             
             !$OMP SIMD
             do cc=1_idef,pchunk
                spp(ii)%vars%B(pp-1+cc,1) = B_X(cc)
                spp(ii)%vars%B(pp-1+cc,2) = B_Y(cc)
                spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                spp(ii)%vars%E(pp-1+cc,1) = E_X(cc)
                spp(ii)%vars%E(pp-1+cc,2) = E_Y(cc)
                spp(ii)%vars%E(pp-1+cc,3) = E_Z(cc)

                spp(ii)%vars%PSI_P(pp-1+cc) = PSIp(cc)
             end do
             !$OMP END SIMD

#ifdef FIO
             if (TRIM(params%field_model).eq.'M3D_C1'.or. &
                  TRIM(params%field_model).eq.'NIMROD') then
                !$OMP SIMD
                do cc=1_idef,pchunk
                   spp(ii)%vars%hint(pp-1+cc) = hint(cc)
                end do
                !$OMP END SIMD
             end if
#endif

             !$OMP SIMD
             !          !$OMP& aligned(Bmag,B_X,B_Y,B_Z, &
             !          !$OMP& b_unit_X,b_unit_Y,b_unit_Z,v,V_X,V_Y,V_Z,vpar, &
             !          !$OMP& vperp,tmp,cross_X,cross_Y,cross_Z, &
             !          !$OMP& vec_X,vec_Y,vec_Z,E_X,E_Y,E_Z)
             do cc=1_idef,pchunk
                !Derived output data
                Bmag(cc) = SQRT(B_X(cc)*B_X(cc)+B_Y(cc)*B_Y(cc)+B_Z(cc)*B_Z(cc))

                ! Parallel unit vector
                b_unit_X(cc) = B_X(cc)/Bmag(cc)
                b_unit_Y(cc) = B_Y(cc)/Bmag(cc)
                b_unit_Z(cc) = B_Z(cc)/Bmag(cc)

                !write(6,*) 'X',X_X,X_Y,X_Z
                !write(6,*) 'b_unit',b_unit_X,b_unit_Y,b_unit_Z
                
                v(cc) = SQRT(V_X(cc)*V_X(cc)+V_Y(cc)*V_Y(cc)+V_Z(cc)*V_Z(cc))
                if (v(cc).GT.korc_zero) then
                   ! Parallel and perpendicular components of velocity
                   vpar(cc) = (V_X(cc)*b_unit_X(cc)+V_Y(cc)*b_unit_Y(cc)+ &
                        V_Z(cc)*b_unit_Z(cc))

                   vperp(cc) =  v(cc)**2 - vpar(cc)**2
                   if ( vperp(cc) .GE. korc_zero ) then
                      vperp(cc) = SQRT( vperp(cc) )
                   else
                      vperp(cc) = 0.0_rp
                   end if

                   !write(6,*) 'v,vpar,vperp',v(cc),vpar(cc),vperp(cc)
                   
                   ! Pitch angle
                   spp(ii)%vars%eta(pp-1+cc) = 180.0_rp* &
                        MODULO(ATAN2(vperp(cc),vpar(cc)),2.0_rp*C_PI)/C_PI

                   ! Magnetic moment
                   spp(ii)%vars%mu(pp-1+cc) = 0.5_rp*m_cache* &
                        g(cc)**2*vperp(cc)**2/Bmag(cc)
                   ! See Northrop's book (The adiabatic motion of charged
                   ! particles)

                   ! Radiated power
                   tmp(cc) = q_cache**4/(6.0_rp*C_PI*E0*m_cache**2)

                   cross_X(cc) = V_Y(cc)*B_Z(cc)-V_Z(cc)*B_Y(cc)
                   cross_Y(cc) = V_Z(cc)*B_X(cc)-V_X(cc)*B_Z(cc)
                   cross_Z(cc) = V_X(cc)*B_Y(cc)-V_Y(cc)*B_X(cc)

                   vec_X(cc) = E_X(cc) + cross_X(cc)
                   vec_Y(cc) = E_Y(cc) + cross_Y(cc)
                   vec_Z(cc) = E_Z(cc) + cross_Z(cc)

                   spp(ii)%vars%Prad(pp-1+cc) = tmp(cc)* &
                        ( E_X(cc)*E_X(cc)+E_Y(cc)*E_Y(cc)+E_Z(cc)*E_Z(cc) + &
                        cross_X(cc)*E_X(cc)+cross_Y(cc)*E_Y(cc)+ &
                        cross_Z(cc)*E_Z(cc) + g(cc)**2* &
                        ((E_X(cc)*V_X(cc)+E_Y(cc)*V_Y(cc)+E_Z(cc)*V_Z(cc))**2 &
                        - vec_X(cc)*vec_X(cc)-vec_Y(cc)*vec_Y(cc)- &
                        vec_Z(cc)*vec_Z(cc)) )

                   ! Input power due to electric field
                   spp(ii)%vars%Pin(pp-1+cc) = q_cache*(E_X(cc)*V_X(cc)+ &
                        E_Y(cc)*V_Y(cc)+E_Z(cc)*V_Z(cc))
                else
                   spp(ii)%vars%eta(pp-1+cc) = 0.0_rp
                   spp(ii)%vars%mu(pp-1+cc) = 0.0_rp
                   spp(ii)%vars%Prad(pp-1+cc) = 0.0_rp
                   spp(ii)%vars%Pin(pp-1+cc) = 0.0_rp
                end if

             end do
             !$OMP END SIMD
             
          end do
          !$OMP END PARALLEL DO

       end if !(if output)

       if(step.and.(.not.params%FokPlan)) then
          dt=0.5_rp*params%dt

          !$OMP PARALLEL DO FIRSTPRIVATE(dt) PRIVATE(pp,cc) &
          !$OMP& SHARED(ii,spp,params)
          do pp=1_idef,spp(ii)%ppp,pchunk

             !$OMP SIMD
             do cc=1_idef,pchunk
                spp(ii)%vars%X(pp-1+cc,1) = spp(ii)%vars%X(pp-1+cc,1) + &
                     dt*spp(ii)%vars%V(pp-1+cc,1)
                spp(ii)%vars%X(pp-1+cc,2) = spp(ii)%vars%X(pp-1+cc,2) + &
                     dt*spp(ii)%vars%V(pp-1+cc,2)
                spp(ii)%vars%X(pp-1+cc,3) = spp(ii)%vars%X(pp-1+cc,3) + &
                     dt*spp(ii)%vars%V(pp-1+cc,3)
             end do
             !$OMP END SIMD

          end do
          !$OMP END PARALLEL DO

       end if !(if step)

    end do ! over species

  end subroutine FO_init

  subroutine adv_FOeqn_top(params,F,P,spp)

    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                                   :: F
    !! An instance of the KORC derived type FIELDS.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    !! An instance of the KORC derived type PROFILES.
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(params%pchunk)               :: Bmag
    REAL(rp), DIMENSION(params%pchunk)               :: b_unit_X,b_unit_Y,b_unit_Z
    REAL(rp), DIMENSION(params%pchunk)               :: v,vpar,vperp
    REAL(rp), DIMENSION(params%pchunk)               :: tmp
    REAL(rp), DIMENSION(params%pchunk)               :: g
    REAL(rp), DIMENSION(params%pchunk)               :: cross_X,cross_Y,cross_Z
    REAL(rp), DIMENSION(params%pchunk)               :: vec_X,vec_Y,vec_Z
    REAL(rp),DIMENSION(params%pchunk) :: X_X,X_Y,X_Z
    REAL(rp),DIMENSION(params%pchunk) :: V_X,V_Y,V_Z
    REAL(rp),DIMENSION(params%pchunk) :: B_X,B_Y,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_X,E_Y,E_Z,PSIp
    INTEGER(is),DIMENSION(params%pchunk) :: flagCon,flagCol

    REAL(rp) :: B0,EF0,R0,q0,lam,ar
    REAL(rp) :: a,m_cache,q_cache
    REAL(rp) :: ne0,Te0,Zeff0



    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    !! time iterator.


    do ii = 1_idef,params%num_species      

       pchunk=params%pchunk
       m_cache=spp(ii)%m
       q_cache=spp(ii)%q
       a = params%dt

       !write(6,*) 'q,m,dt,a',q_cache,m_cache,params%dt,a
       
       B0=F%Bo
       EF0=F%Eo
       lam=F%AB%lambda
       R0=F%AB%Ro
       q0=F%AB%qo
       ar=F%AB%a



       !$OMP PARALLEL DO default(none) &
       !$OMP& FIRSTPRIVATE(E0,a,m_cache,q_cache,B0,EF0,lam,R0,q0,ar,pchunk)&
       !$OMP& shared(params,ii,spp,P,F) &
       !$OMP& PRIVATE(pp,tt,Bmag,cc,X_X,X_Y,X_Z,V_X,V_Y,V_Z,B_X,B_Y,B_Z, &
       !$OMP& E_X,E_Y,E_Z,b_unit_X,b_unit_Y,b_unit_Z,v,vpar,vperp,tmp, &
       !$OMP& cross_X,cross_Y,cross_Z,vec_X,vec_Y,vec_Z,g,flagCon,flagCol,PSIp)
       do pp=1_idef,spp(ii)%ppp,pchunk

          !$OMP SIMD
          do cc=1_idef,pchunk
             X_X(cc)=spp(ii)%vars%X(pp-1+cc,1)
             X_Y(cc)=spp(ii)%vars%X(pp-1+cc,2)
             X_Z(cc)=spp(ii)%vars%X(pp-1+cc,3)

             V_X(cc)=spp(ii)%vars%V(pp-1+cc,1)
             V_Y(cc)=spp(ii)%vars%V(pp-1+cc,2)
             V_Z(cc)=spp(ii)%vars%V(pp-1+cc,3)

             PSIp(cc)=spp(ii)%vars%PSI_P(pp-1+cc)

             g(cc)=spp(ii)%vars%g(pp-1+cc)
             flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
             flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)
          end do
          !$OMP END SIMD

          if (.not.params%FokPlan) then
             do tt=1_ip,params%t_skip

                !write(6,*) 'tt',tt

                if (params%field_model(1:3).eq.'ANA') then
                   call analytical_fields_p(params,pchunk,F, &
                        X_X,X_Y,X_Z,B_X,B_Y,B_Z,E_X,E_Y,E_Z,flagCon)
                else if (params%field_model(1:3).eq.'UNI') then
                   call uniform_fields_p(pchunk,F, &
                        B_X,B_Y,B_Z,E_X,E_Y,E_Z)
                end if

                call advance_FOeqn_vars(tt,a,q_cache,m_cache,params, &
                     X_X,X_Y,X_Z,V_X,V_Y,V_Z,B_X,B_Y,B_Z,E_X,E_Y,E_Z, &
                     P,F,g,flagCon,flagCol,PSIp)
             end do !timestep iterator

             !$OMP SIMD
             do cc=1_idef,pchunk
                spp(ii)%vars%X(pp-1+cc,1)=X_X(cc)
                spp(ii)%vars%X(pp-1+cc,2)=X_Y(cc)
                spp(ii)%vars%X(pp-1+cc,3)=X_Z(cc)

                spp(ii)%vars%V(pp-1+cc,1)=V_X(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_Y(cc)
                spp(ii)%vars%V(pp-1+cc,3)=V_Z(cc)

                spp(ii)%vars%g(pp-1+cc) = g(cc)

                spp(ii)%vars%flagCon(pp-1+cc) = flagCon(cc)
                spp(ii)%vars%flagCol(pp-1+cc) = flagCol(cc)

                spp(ii)%vars%B(pp-1+cc,1) = B_X(cc)
                spp(ii)%vars%B(pp-1+cc,2) = B_Y(cc)
                spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                spp(ii)%vars%E(pp-1+cc,1) = E_X(cc)
                spp(ii)%vars%E(pp-1+cc,2) = E_Y(cc)
                spp(ii)%vars%E(pp-1+cc,3) = E_Z(cc)
             end do
             !$OMP END SIMD

          else

             !$OMP SIMD
             do cc=1_idef,pchunk
                B_X(cc)=spp(ii)%vars%B(pp-1+cc,1)
                B_Y(cc)=spp(ii)%vars%B(pp-1+cc,2)
                B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)

                E_X(cc)=spp(ii)%vars%E(pp-1+cc,1)
                E_Y(cc)=spp(ii)%vars%E(pp-1+cc,2)
                E_Z(cc)=spp(ii)%vars%E(pp-1+cc,3)
             end do
             !$OMP END SIMD

             call advance_FP3Deqn_vars(params,X_X,X_Y,X_Z,V_X,V_Y,V_Z, &
                  g,m_cache,B0,lam,R0,q0,EF0,B_X,B_Y,B_Z,E_X,E_Y,E_Z, &
                  P,F,flagCon,flagCol,PSIp)

             !$OMP SIMD
             do cc=1_idef,pchunk

                spp(ii)%vars%V(pp-1+cc,1)=V_X(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_Y(cc)
                spp(ii)%vars%V(pp-1+cc,3)=V_Z(cc)

                spp(ii)%vars%g(pp-1+cc) = g(cc)

                flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)
             end do
             !$OMP END SIMD

          end if

          !$OMP SIMD
          !          !$OMP& aligned(Bmag,B_X,B_Y,B_Z, &
          !          !$OMP& b_unit_X,b_unit_Y,b_unit_Z,v,V_X,V_Y,V_Z,vpar, &
          !          !$OMP& vperp,tmp,cross_X,cross_Y,cross_Z, &
          !          !$OMP& vec_X,vec_Y,vec_Z,E_X,E_Y,E_Z)
          do cc=1_idef,pchunk
             !Derived output data
             Bmag(cc) = SQRT(B_X(cc)*B_X(cc)+B_Y(cc)*B_Y(cc)+B_Z(cc)*B_Z(cc))

             ! Parallel unit vector
             b_unit_X(cc) = B_X(cc)/Bmag(cc)
             b_unit_Y(cc) = B_Y(cc)/Bmag(cc)
             b_unit_Z(cc) = B_Z(cc)/Bmag(cc)

             v(cc) = SQRT(V_X(cc)*V_X(cc)+V_Y(cc)*V_Y(cc)+V_Z(cc)*V_Z(cc))
             if (v(cc).GT.korc_zero) then
                ! Parallel and perpendicular components of velocity
                vpar(cc) = (V_X(cc)*b_unit_X(cc)+V_Y(cc)*b_unit_Y(cc)+ &
                     V_Z(cc)*b_unit_Z(cc))

                vperp(cc) =  v(cc)**2 - vpar(cc)**2
                if ( vperp(cc) .GE. korc_zero ) then
                   vperp(cc) = SQRT( vperp(cc) )
                else
                   vperp(cc) = 0.0_rp
                end if

                ! Pitch angle
                spp(ii)%vars%eta(pp-1+cc) = 180.0_rp* &
                     MODULO(ATAN2(vperp(cc),vpar(cc)),2.0_rp*C_PI)/C_PI

                ! Magnetic moment
                spp(ii)%vars%mu(pp-1+cc) = 0.5_rp*m_cache* &
                     g(cc)**2*vperp(cc)**2/Bmag(cc)
                ! See Northrop's book (The adiabatic motion of charged
                ! particles)

                ! Radiated power
                tmp(cc) = q_cache**4/(6.0_rp*C_PI*E0*m_cache**2)

                cross_X(cc) = V_Y(cc)*B_Z(cc)-V_Z(cc)*B_Y(cc)
                cross_Y(cc) = V_Z(cc)*B_X(cc)-V_X(cc)*B_Z(cc)
                cross_Z(cc) = V_X(cc)*B_Y(cc)-V_Y(cc)*B_X(cc)

                vec_X(cc) = E_X(cc) + cross_X(cc)
                vec_Y(cc) = E_Y(cc) + cross_Y(cc)
                vec_Z(cc) = E_Z(cc) + cross_Z(cc)

                spp(ii)%vars%Prad(pp-1+cc) = tmp(cc)* &
                     ( E_X(cc)*E_X(cc)+E_Y(cc)*E_Y(cc)+E_Z(cc)*E_Z(cc) + &
                     cross_X(cc)*E_X(cc)+cross_Y(cc)*E_Y(cc)+ &
                     cross_Z(cc)*E_Z(cc) + g(cc)**2* &
                     ((E_X(cc)*V_X(cc)+E_Y(cc)*V_Y(cc)+E_Z(cc)*V_Z(cc))**2 &
                     - vec_X(cc)*vec_X(cc)-vec_Y(cc)*vec_Y(cc)- &
                     vec_Z(cc)*vec_Z(cc)) )

                ! Input power due to electric field
                spp(ii)%vars%Pin(pp-1+cc) = q_cache*(E_X(cc)*V_X(cc)+ &
                     E_Y(cc)*V_Y(cc)+E_Z(cc)*V_Z(cc))
             else
                spp(ii)%vars%eta(pp-1+cc) = 0.0_rp
                spp(ii)%vars%mu(pp-1+cc) = 0.0_rp
                spp(ii)%vars%Prad(pp-1+cc) = 0.0_rp
                spp(ii)%vars%Pin(pp-1+cc) = 0.0_rp
             end if

          end do
          !$OMP END SIMD


       end do !particle chunk iterator
       !$OMP END PARALLEL DO

    end do !species iterator

  end subroutine adv_FOeqn_top

  subroutine advance_FOeqn_vars(tt,a,q_cache,m_cache,params,X_X,X_Y,X_Z, &
       V_X,V_Y,V_Z,B_X,B_Y,B_Z,E_X,E_Y,E_Z,P,F,g,flagCon,flagCol,PSIp)
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)      :: F
    TYPE(KORC_PARAMS), INTENT(IN)                              :: params
    !! Core KORC simulation parameters.

    INTEGER(ip), INTENT(IN)                                       :: tt
    !! Time step used in the leapfrog step (\(\Delta t\)).
    REAL(rp)                                      :: dt
    !! Time step used in the leapfrog step (\(\Delta t\)).
    REAL(rp), INTENT(IN)                       :: m_cache,q_cache
    !! Time step used in the leapfrog step (\(\Delta t\)).

    REAL(rp),DIMENSION(params%pchunk)                                  :: Bmag



    REAL(rp),INTENT(in)                                       :: a
    !! This variable is used to simplify notation in the code, and
    !! is given by \(a=q\Delta t/m\),
    REAL(rp),DIMENSION(params%pchunk)                                    :: sigma
    !! This variable is \(\sigma = \gamma'^2 - \tau^2\) in the above equations.
    REAL(rp),DIMENSION(params%pchunk)                               :: us
    !! This variable is \(u^{*} = p^{*}/m\) where \( p^{*} =
    !! \mathbf{p}'\cdot \mathbf{\tau}/mc\).
    !! Variable 'u^*' in Vay, J.-L. PoP (200params%pchunk).
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT)                 :: g
    REAL(rp),DIMENSION(params%pchunk) :: gp,g0
    !! Relativistic factor \(\gamma\).
    REAL(rp),DIMENSION(params%pchunk)                                 :: s
    !! This variable is \(s = 1/(1+t^2)\) in the equations above.
    !! Variable 's' in Vay, J.-L. PoP (2008).
    REAL(rp),DIMENSION(params%pchunk)                            :: U_hs_X,U_hs_Y,U_hs_Z
    !! Is \(\mathbf{u}=\mathbf{p}/m\) at half-time step (\(i+1/2\)) in
    !! the absence of radiation losses or collisions. \(\mathbf{u}^{i+1/2} =
    !! \mathbf{u}^i + \frac{q\Delta t}{2m}\left( \mathbf{E}^{i+1/2} +
    !! \mathbf{v}^i\times \mathbf{B}^{i+1/2} \right)\).
    REAL(rp),DIMENSION(params%pchunk)                           :: tau_X,tau_Y,tau_Z
    !! This variable is \(\mathbf{\tau} = (q\Delta t/2)\mathbf{B}^{i+1/2}\).
    REAL(rp),DIMENSION(params%pchunk)                            :: up_X,up_Y,up_Z
    !! This variable is \(\mathbf{u}'= \mathbf{p}'/m\), where \(\mathbf{p}'
    !! = \mathbf{p}^i + q\Delta t \left( \mathbf{E}^{i+1/2} +
    !! \frac{\mathbf{v}^i}{2} \times \mathbf{B}^{i+1/2} \right)\).
    REAL(rp),DIMENSION(params%pchunk)                                     :: t_X,t_Y,t_Z
    !! This variable is \(\mathbf{t} = {\mathbf \tau}/\gamma^{i+1}\).
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT)                     :: X_X,X_Y,X_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT)                      :: V_X,V_Y,V_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(IN)                      :: B_X,B_Y,B_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(IN)          :: E_X,E_Y,E_Z,PSIp
    REAL(rp),DIMENSION(params%pchunk)                     :: U_L_X,U_L_Y,U_L_Z
    REAL(rp),DIMENSION(params%pchunk)                     :: U_X,U_Y,U_Z
    REAL(rp),DIMENSION(params%pchunk)                     :: U_RC_X,U_RC_Y,U_RC_Z
    REAL(rp),DIMENSION(params%pchunk)                     :: U_os_X,U_os_Y,U_os_Z
    !! This variable is \(\mathbf{u}^{i+1}= \mathbf{p}^{i+1}/m\).
    REAL(rp),DIMENSION(params%pchunk)                          :: cross_X,cross_Y,cross_Z

    REAL(rp), DIMENSION(params%pchunk)                       :: Frad_X,Frad_Y,Frad_Z
    !! Synchrotron radiation reaction force of each particle.

    REAL(rp),DIMENSION(params%pchunk) :: ne,Te,Zeff,Y_R,Y_PHI,Y_Z

    INTEGER                                      :: cc,pchunk
    !! Chunk iterator.

    INTEGER(is),DIMENSION(params%pchunk),intent(inout)             :: flagCon,flagCol

    dt=params%dt

    pchunk=params%pchunk

    !write(6,*) 'X',X_X,X_Y,X_Z
    !write(6,*) 'V',V_X,V_Y,V_Z
    !write(6,*) 'B',B_X,B_Y,B_Z
    
    !$OMP SIMD
    !    !$OMP& aligned(g0,g,U_X,U_Y,U_Z,V_X,V_Y,V_Z,Bmag,B_X,B_Y,B_Z, &
    !    !$OMP& U_L_X,U_L_Y,U_L_Z,U_RC_X,U_RC_Y,U_RC_Z, &
    !    !$OMP& cross_X,cross_Y,cross_Z,U_hs_X,U_hs_Y,U_hs_Z,E_X,E_Y,E_Z, &
    !    !$OMP& tau_X,tau_Y,tau_Z,up_X,up_Y,up_Z,gp,sigma,us,t_X,t_Y,t_Z,s, &
    !    !$OMP& U_os_X,U_os_Y,U_os_Z,Frad_X,Frad_Y,Frad_Z)
    do cc=1_idef,pchunk

       g0(cc)=g(cc)

       U_X(cc) = g(cc)*V_X(cc)
       U_Y(cc) = g(cc)*V_Y(cc)
       U_Z(cc) = g(cc)*V_Z(cc)


       ! Magnitude of magnetic field
       Bmag(cc) = SQRT(B_X(cc)*B_X(cc)+B_Y(cc)*B_Y(cc)+B_Z(cc)*B_Z(cc))

       U_L_X(cc)=U_X(cc)
       U_L_Y(cc)=U_Y(cc)
       U_L_Z(cc)=U_Z(cc)

       U_RC_X(cc)=U_X(cc)
       U_RC_Y(cc)=U_Y(cc)
       U_RC_Z(cc)=U_Z(cc)

       ! LEAP-FROG SCHEME FOR LORENTZ FORCE !

       cross_X(cc)=V_Y(cc)*B_Z(cc)-V_Z(cc)*B_Y(cc)
       cross_Y(cc)=V_Z(cc)*B_X(cc)-V_X(cc)*B_Z(cc)
       cross_Z(cc)=V_X(cc)*B_Y(cc)-V_Y(cc)*B_X(cc)

       !write(6,*) 'vcrossB',cross_X,cross_Y,cross_Z

       U_hs_X(cc) = U_L_X(cc) + 0.5_rp*a*(E_X(cc) +cross_X(cc))
       U_hs_Y(cc) = U_L_Y(cc) + 0.5_rp*a*(E_Y(cc) +cross_Y(cc))
       U_hs_Z(cc) = U_L_Z(cc) + 0.5_rp*a*(E_Z(cc) +cross_Z(cc))


       !write(6,*) 'half step',0.5_rp*a*cross_X(cc),0.5_rp*a*cross_Y(cc),0.5_rp*a*cross_Z(cc)

       tau_X(cc) = 0.5_rp*a*B_X(cc)
       tau_Y(cc) = 0.5_rp*a*B_Y(cc)
       tau_Z(cc) = 0.5_rp*a*B_Z(cc)



       up_X(cc) = U_hs_X(cc) + 0.5_rp*a*E_X(cc)
       up_Y(cc) = U_hs_Y(cc) + 0.5_rp*a*E_Y(cc)
       up_Z(cc) = U_hs_Z(cc) + 0.5_rp*a*E_Z(cc)

       gp(cc) = SQRT( 1.0_rp + up_X(cc)*up_X(cc)+up_Y(cc)*up_Y(cc)+ &
            up_Z(cc)*up_Z(cc) )

       sigma(cc) = gp(cc)*gp(cc) - (tau_X(cc)*tau_X(cc)+ &
            tau_Y(cc)*tau_Y(cc)+tau_Z(cc)*tau_Z(cc))

       us(cc) = up_X(cc)*tau_X(cc)+up_Y(cc)*tau_Y(cc)+ &
            up_Z(cc)*tau_Z(cc)

       ! variable 'u^*' in Vay, J.-L. PoP (2008)
       g(cc) = SQRT( 0.5_rp*(sigma(cc) + SQRT(sigma(cc)*sigma(cc) + &
            4.0_rp*(tau_X(cc)*tau_X(cc)+tau_Y(cc)*tau_Y(cc)+ &
            tau_Z(cc)*tau_Z(cc) + us(cc)*us(cc)))) )

       t_X(cc) = tau_X(cc)/g(cc)
       t_Y(cc) = tau_Y(cc)/g(cc)
       t_Z(cc) = tau_Z(cc)/g(cc)


       s(cc) = 1.0_rp/(1.0_rp + t_X(cc)*t_X(cc)+t_Y(cc)*t_Y(cc)+ &
            t_Z(cc)*t_Z(cc))
       ! variable 's' in Vay, J.-L. PoP (2008)

       cross_X(cc)=up_Y(cc)*t_Z(cc)-up_Z(cc)*t_Y(cc)
       cross_Y(cc)=up_Z(cc)*t_X(cc)-up_X(cc)*t_Z(cc)
       cross_Z(cc)=up_X(cc)*t_Y(cc)-up_Y(cc)*t_X(cc)

       U_L_X(cc) = s(cc)*(up_X(cc) + (up_X(cc)*t_X(cc)+ &
            up_Y(cc)*t_Y(cc)+up_Z(cc)*t_Z(cc))*t_X(cc) + cross_X(cc))
       U_L_Y(cc) = s(cc)*(up_Y(cc) + (up_X(cc)*t_X(cc)+ &
            up_Y(cc)*t_Y(cc)+up_Z(cc)*t_Z(cc))*t_Y(cc) + cross_Y(cc))
       U_L_Z(cc) = s(cc)*(up_Z(cc) + (up_X(cc)*t_X(cc)+ &
            up_Y(cc)*t_Y(cc)+up_Z(cc)*t_Z(cc))*t_Z(cc) + cross_Z(cc))
       ! LEAP-FROG SCHEME FOR LORENTZ FORCE !      

       U_os_X(cc) = 0.5_rp*(U_L_X(cc) + U_X(cc))
       U_os_Y(cc) = 0.5_rp*(U_L_Y(cc) + U_Y(cc))
       U_os_Z(cc) = 0.5_rp*(U_L_Z(cc) + U_Z(cc))
       ! Splitting operator for including radiation

       if (params%radiation) then
          !! Calls [[radiation_force_p]] in [[korc_ppusher]].
          call radiation_force_p(pchunk,q_cache,m_cache,U_os_X,U_os_Y,U_os_Z, &
               E_X,E_Y,E_Z,B_Z,B_Y,B_Z,Frad_X,Frad_Y,Frad_Z)
          U_RC_X(cc) = U_RC_X(cc) + a*Frad_X(cc)/q_cache
          U_RC_Y(cc) = U_RC_Y(cc) + a*Frad_Y(cc)/q_cache
          U_RC_Z(cc) = U_RC_Z(cc) + a*Frad_Z(cc)/q_cache
       end if
       ! Splitting operator for including radiation

       U_X(cc) = U_L_X(cc) + U_RC_X(cc) - U_X(cc)
       U_Y(cc) = U_L_Y(cc) + U_RC_Y(cc) - U_Y(cc)
       U_Z(cc) = U_L_Z(cc) + U_RC_Z(cc) - U_Z(cc)

    end do
    !$OMP END SIMD


    if (params%collisions) then

       call include_CoulombCollisions_FO_p(tt,params,X_X,X_Y,X_Z, &
            U_X,U_Y,U_Z,B_X,B_Y,B_Z,m_cache,P,F,flagCon,flagCol,PSIp)

    end if

    if (params%radiation.or.params%collisions) then

       !$OMP SIMD
       !       !$OMP& aligned(g,U_X,U_Y,U_Z)
       do cc=1_idef,pchunk
          g(cc)=sqrt(1._rp+U_X(cc)*U_X(cc)+U_Y(cc)*U_Y(cc)+U_Z(cc)*U_Z(cc))
       end do
       !$OMP END SIMD

    end if

    !$OMP SIMD
    !    !$OMP& aligned(g,g0,V_X,V_Y,V_Z,U_X,U_Y,U_Z,X_X,X_Y,X_Z,flagCon,flagCol)
    do cc=1_idef,pchunk

       if ((flagCon(cc).eq.0_is).or.(flagCol(cc).eq.0_is)) then
          g(cc)=g0(cc)
       else
          V_X(cc) = U_X(cc)/g(cc)
          V_Y(cc) = U_Y(cc)/g(cc)
          V_Z(cc) = U_Z(cc)/g(cc)
       end if

       X_X(cc) = X_X(cc) + dt*V_X(cc)*REAL(flagCon(cc))*REAL(flagCol(cc))
       X_Y(cc) = X_Y(cc) + dt*V_Y(cc)*REAL(flagCon(cc))*REAL(flagCol(cc))
       X_Z(cc) = X_Z(cc) + dt*V_Z(cc)*REAL(flagCon(cc))*REAL(flagCol(cc))
    end do
    !$OMP END SIMD

  end subroutine advance_FOeqn_vars

  subroutine advance_FP3Deqn_vars(params,X_X,X_Y,X_Z,V_X,V_Y,V_Z,g, &
       m_cache,B0,lam,R0,q0,EF0,B_X,B_Y,B_Z,E_X,E_Y,E_Z, &
       P,F,flagCon,flagCol,PSIp)
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)      :: F
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    !! time iterator.
    REAL(rp),DIMENSION(params%pchunk), INTENT(IN)  :: X_X,X_Y,X_Z
    REAL(rp),DIMENSION(params%pchunk), INTENT(IN)  :: E_X,E_Y,E_Z,PSIp
    REAL(rp),DIMENSION(params%pchunk), INTENT(IN)  :: B_X,B_Y,B_Z
    INTEGER(is),DIMENSION(params%pchunk), INTENT(INOUT)  :: flagCon,flagCol
    REAL(rp),DIMENSION(params%pchunk) :: U_X,U_Y,U_Z
    REAL(rp),DIMENSION(params%pchunk), INTENT(INOUT)  :: V_X,V_Y,V_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: g
    REAL(rp),intent(in) :: B0,EF0,R0,q0,lam,m_cache

    pchunk=params%pchunk

    !    call analytical_fields_p(B0,EF0,R0,q0,lam,X_X,X_Y,X_Z, &
    !         B_X,B_Y,B_Z,E_X,E_Y,E_Z)

    !$OMP SIMD
    !    !$OMP& aligned(U_X,U_Y,U_Z,V_X,V_Y,V_Z,g)
    do cc=1_idef,pchunk
       U_X(cc)=V_X(cc)*g(cc)
       U_Y(cc)=V_Y(cc)*g(cc)
       U_Z(cc)=V_Z(cc)*g(cc)
    end do
    !$OMP END SIMD

    do tt=1_ip,params%t_skip

       call include_CoulombCollisions_FO_p(tt,params,X_X,X_Y,X_Z, &
            U_X,U_Y,U_Z,B_X,B_Y,B_Z,m_cache,P,F,flagCon,flagCol,PSIp)

    end do

    !$OMP SIMD
    !    !$OMP& aligned(U_X,U_Y,U_Z,V_X,V_Y,V_Z,g)
    do cc=1_idef,pchunk

       g(cc)=sqrt(1._rp+U_X(cc)*U_X(cc)+U_Y(cc)*U_Y(cc)+U_Z(cc)*U_Z(cc))

       V_X(cc)=U_X(cc)/g(cc)
       V_Y(cc)=U_Y(cc)/g(cc)
       V_Z(cc)=U_Z(cc)/g(cc)
    end do
    !$OMP END SIMD

  end subroutine advance_FP3Deqn_vars

#ifdef FIO
  subroutine adv_FOfio_top(params,F,P,spp)  
    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                                   :: F
    !! An instance of the KORC derived type FIELDS.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    !! An instance of the KORC derived type PROFILES.
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(params%pchunk)               :: Bmag
    REAL(rp), DIMENSION(params%pchunk)               :: b_unit_X,b_unit_Y,b_unit_Z
    REAL(rp), DIMENSION(params%pchunk)               :: v,vpar,vperp
    REAL(rp), DIMENSION(params%pchunk)               :: tmp
    REAL(rp), DIMENSION(params%pchunk)               :: g
    REAL(rp), DIMENSION(params%pchunk)               :: cross_X,cross_Y,cross_Z
    REAL(rp), DIMENSION(params%pchunk)               :: vec_X,vec_Y,vec_Z
    REAL(rp),DIMENSION(params%pchunk) :: X_X,X_Y,X_Z
    REAL(rp),DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: V_X,V_Y,V_Z
    REAL(rp),DIMENSION(params%pchunk) :: B_X,B_Y,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_X,E_Y,E_Z
    REAL(rp),DIMENSION(params%pchunk) :: PSIp
    TYPE(C_PTR),DIMENSION(params%pchunk) :: hint
    INTEGER(is),DIMENSION(params%pchunk) :: flagCon,flagCol
    REAL(rp) :: a,m_cache,q_cache    
    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp,pchunk
    !! Particles iterator.
    INTEGER                                                    :: cc
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    !! time iterator.


    do ii = 1_idef,params%num_species      

       pchunk=params%pchunk
       m_cache=spp(ii)%m
       q_cache=spp(ii)%q
       a = params%dt


       !$OMP PARALLEL DO default(none) &
       !$OMP& FIRSTPRIVATE(a,m_cache,q_cache,E0,pchunk) &
       !$OMP& shared(params,ii,spp,P,F) &
       !$OMP& PRIVATE(pp,tt,Bmag,cc,X_X,X_Y,X_Z,V_X,V_Y,V_Z,B_X,B_Y,B_Z, &
       !$OMP& E_X,E_Y,E_Z,b_unit_X,b_unit_Y,b_unit_Z,v,vpar,vperp,tmp, &
       !$OMP& cross_X,cross_Y,cross_Z,vec_X,vec_Y,vec_Z,g, &
       !$OMP& Y_R,Y_PHI,Y_Z,flagCon,flagCol,PSIp,hint)
       do pp=1_idef,spp(ii)%ppp,pchunk

          !$OMP SIMD
          do cc=1_idef,pchunk
             X_X(cc)=spp(ii)%vars%X(pp-1+cc,1)
             X_Y(cc)=spp(ii)%vars%X(pp-1+cc,2)
             X_Z(cc)=spp(ii)%vars%X(pp-1+cc,3)

             Y_R(cc)=0._rp
             Y_PHI(cc)=0._rp
             Y_Z(cc)=0._rp

             V_X(cc)=spp(ii)%vars%V(pp-1+cc,1)
             V_Y(cc)=spp(ii)%vars%V(pp-1+cc,2)
             V_Z(cc)=spp(ii)%vars%V(pp-1+cc,3)

             B_X(cc)=spp(ii)%vars%B(pp-1+cc,1)
             B_Y(cc)=spp(ii)%vars%B(pp-1+cc,2)
             B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)

             E_X(cc)=spp(ii)%vars%E(pp-1+cc,1)
             E_Y(cc)=spp(ii)%vars%E(pp-1+cc,2)
             E_Z(cc)=spp(ii)%vars%E(pp-1+cc,3)

             PSIp(cc)=spp(ii)%vars%PSI_P(pp-1+cc)

             hint(cc)=spp(ii)%vars%hint(pp-1+cc)

             g(cc)=spp(ii)%vars%g(pp-1+cc)

             flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
             flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)
          end do
          !$OMP END SIMD

          !write(output_unit_write,*) 'Yin: ',Y_R,Y_PHI,Y_Z
          !write(output_unit_write,*) 'Bin: ',B_X,B_Y,B_Z


          if (.not.params%FokPlan) then
             do tt=1_ip,params%t_skip

                call cart_to_cyl_p(pchunk,X_X,X_Y,X_Z,Y_R,Y_PHI,Y_Z)

                call get_fio_FOmagnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
                     B_X,B_Y,B_Z,flagCon,hint)
                if (F%FIO_E .ge. 0) then
                   call get_fio_FOelectric_fields_p(params,F, &
                        Y_R,Y_PHI,Y_Z,E_X,E_Y,E_Z,flagCon,hint)
                end if
                if (F%FIO_A .ge. 0) then
                   call get_fio_vector_potential_p(params,F,Y_R,Y_PHI,Y_Z, &
                        PSIp,flagCon,hint)
                end if

                call advance_FOfio_vars(tt,a,q_cache,m_cache,params, &
                     X_X,X_Y,X_Z,V_X,V_Y,V_Z,B_X,B_Y,B_Z,E_X,E_Y,E_Z, &
                     g,flagCon,flagCol,P,F,PSIp,hint)

             end do !timestep iterator

             call cart_to_cyl_p(pchunk,X_X,X_Y,X_Z,Y_R,Y_PHI,Y_Z)

             if (F%FIO_A .ge. 0) then
                call get_fio_vector_potential_p(params,F,Y_R,Y_PHI,Y_Z, &
                     PSIp,flagCon,hint)
             endif
                
             !write(output_unit_write,*) 'Yout: ',Y_R,Y_PHI,Y_Z
             !write(output_unit_write,*) 'Bout: ',B_X,B_Y,B_Z



             !$OMP SIMD
             do cc=1_idef,pchunk
                spp(ii)%vars%X(pp-1+cc,1)=X_X(cc)
                spp(ii)%vars%X(pp-1+cc,2)=X_Y(cc)
                spp(ii)%vars%X(pp-1+cc,3)=X_Z(cc)

                spp(ii)%vars%Y(pp-1+cc,1)=Y_R(cc)
                spp(ii)%vars%Y(pp-1+cc,2)=Y_PHI(cc)
                spp(ii)%vars%Y(pp-1+cc,3)=Y_Z(cc)

                spp(ii)%vars%V(pp-1+cc,1)=V_X(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_Y(cc)
                spp(ii)%vars%V(pp-1+cc,3)=V_Z(cc)

                spp(ii)%vars%g(pp-1+cc) = g(cc)
                spp(ii)%vars%flagCon(pp-1+cc) = flagCon(cc)
                spp(ii)%vars%flagCol(pp-1+cc) = flagCol(cc)

                spp(ii)%vars%B(pp-1+cc,1) = B_X(cc)
                spp(ii)%vars%B(pp-1+cc,2) = B_Y(cc)
                spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                spp(ii)%vars%E(pp-1+cc,1) = E_X(cc)
                spp(ii)%vars%E(pp-1+cc,2) = E_Y(cc)
                spp(ii)%vars%E(pp-1+cc,3) = E_Z(cc)

                spp(ii)%vars%PSI_P(pp-1+cc) = PSIp(cc)

                spp(ii)%vars%hint(pp-1+cc) = hint(cc)
             end do
             !$OMP END SIMD

          else
             !$OMP SIMD
             do cc=1_idef,pchunk
                B_X(cc)=spp(ii)%vars%B(pp-1+cc,1)
                B_Y(cc)=spp(ii)%vars%B(pp-1+cc,2)
                B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)

                E_X(cc)=spp(ii)%vars%E(pp-1+cc,1)
                E_Y(cc)=spp(ii)%vars%E(pp-1+cc,2)
                E_Z(cc)=spp(ii)%vars%E(pp-1+cc,3)
             end do
             !$OMP END SIMD

             call advance_FP3Dinterp_vars(params,X_X,X_Y,X_Z,V_X,V_Y,V_Z, &
                  g,m_cache,B_X,B_Y,B_Z,E_X,E_Y,E_Z,flagCon,flagCol,P,F,PSIp)

             !$OMP SIMD
             do cc=1_idef,pchunk

                spp(ii)%vars%V(pp-1+cc,1)=V_X(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_Y(cc)
                spp(ii)%vars%V(pp-1+cc,3)=V_Z(cc)

                spp(ii)%vars%g(pp-1+cc) = g(cc)
                spp(ii)%vars%flagCol(pp-1+cc) = flagCol(cc)

             end do
             !$OMP END SIMD
          end if

          !$OMP SIMD
          !          !$OMP& aligned(Bmag,B_X,B_Y,B_Z, &
          !          !$OMP& b_unit_X,b_unit_Y,b_unit_Z,v,V_X,V_Y,V_Z,vpar, &
          !          !$OMP& vperp,tmp,cross_X,cross_Y,cross_Z, &
          !          !$OMP& vec_X,vec_Y,vec_Z,E_X,E_Y,E_Z)
          do cc=1_idef,pchunk
             !Derived output data
             Bmag(cc) = SQRT(B_X(cc)*B_X(cc)+B_Y(cc)*B_Y(cc)+B_Z(cc)*B_Z(cc))

             ! Parallel unit vector
             b_unit_X(cc) = B_X(cc)/Bmag(cc)
             b_unit_Y(cc) = B_Y(cc)/Bmag(cc)
             b_unit_Z(cc) = B_Z(cc)/Bmag(cc)

             v(cc) = SQRT(V_X(cc)*V_X(cc)+V_Y(cc)*V_Y(cc)+V_Z(cc)*V_Z(cc))
             if (v(cc).GT.korc_zero) then
                ! Parallel and perpendicular components of velocity
                vpar(cc) = (V_X(cc)*b_unit_X(cc)+V_Y(cc)*b_unit_Y(cc)+ &
                     V_Z(cc)*b_unit_Z(cc))

                vperp(cc) =  v(cc)**2 - vpar(cc)**2
                if ( vperp(cc) .GE. korc_zero ) then
                   vperp(cc) = SQRT( vperp(cc) )
                else
                   vperp(cc) = 0.0_rp
                end if

                ! Pitch angle
                spp(ii)%vars%eta(pp-1+cc) = 180.0_rp* &
                     MODULO(ATAN2(vperp(cc),vpar(cc)),2.0_rp*C_PI)/C_PI

                ! Magnetic moment
                spp(ii)%vars%mu(pp-1+cc) = 0.5_rp*m_cache* &
                     g(cc)**2*vperp(cc)**2/Bmag(cc)
                ! See Northrop's book (The adiabatic motion of charged
                ! particles)

                ! Radiated power
                tmp(cc) = q_cache**4/(6.0_rp*C_PI*E0*m_cache**2)

                cross_X(cc) = V_Y(cc)*B_Z(cc)-V_Z(cc)*B_Y(cc)
                cross_Y(cc) = V_Z(cc)*B_X(cc)-V_X(cc)*B_Z(cc)
                cross_Z(cc) = V_X(cc)*B_Y(cc)-V_Y(cc)*B_X(cc)

                vec_X(cc) = E_X(cc) + cross_X(cc)
                vec_Y(cc) = E_Y(cc) + cross_Y(cc)
                vec_Z(cc) = E_Z(cc) + cross_Z(cc)

                spp(ii)%vars%Prad(pp-1+cc) = tmp(cc)* &
                     ( E_X(cc)*E_X(cc)+E_Y(cc)*E_Y(cc)+E_Z(cc)*E_Z(cc) + &
                     cross_X(cc)*E_X(cc)+cross_Y(cc)*E_Y(cc)+ &
                     cross_Z(cc)*E_Z(cc) + g(cc)**2* &
                     ((E_X(cc)*V_X(cc)+E_Y(cc)*V_Y(cc)+E_Z(cc)*V_Z(cc))**2 &
                     - vec_X(cc)*vec_X(cc)+vec_Y(cc)*vec_Y(cc)+ &
                     vec_Z(cc)*vec_Z(cc)) )

                ! Input power due to electric field
                spp(ii)%vars%Pin(pp-1+cc) = q_cache*(E_X(cc)*V_X(cc)+ &
                     E_Y(cc)*V_Y(cc)+E_Z(cc)*V_Z(cc))
             else
                spp(ii)%vars%eta(pp-1+cc) = 0.0_rp
                spp(ii)%vars%mu(pp-1+cc) = 0.0_rp
                spp(ii)%vars%Prad(pp-1+cc) = 0.0_rp
                spp(ii)%vars%Pin(pp-1+cc) = 0.0_rp
             end if

          end do
          !$OMP END SIMD


       end do !particle chunk iterator
       !$OMP END PARALLEL DO

    end do !species iterator

  end subroutine adv_FOfio_top
#endif 

#ifdef PSPLINE
  subroutine adv_FOinterp_top(params,F,P,spp)  
    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                                   :: F
    !! An instance of the KORC derived type FIELDS.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    !! An instance of the KORC derived type PROFILES.
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(params%pchunk)               :: Bmag
    REAL(rp), DIMENSION(params%pchunk)               :: b_unit_X,b_unit_Y,b_unit_Z
    REAL(rp), DIMENSION(params%pchunk)               :: v,vpar,vperp
    REAL(rp), DIMENSION(params%pchunk)               :: tmp
    REAL(rp), DIMENSION(params%pchunk)               :: g
    REAL(rp), DIMENSION(params%pchunk)               :: cross_X,cross_Y,cross_Z
    REAL(rp), DIMENSION(params%pchunk)               :: vec_X,vec_Y,vec_Z
    REAL(rp),DIMENSION(params%pchunk) :: X_X,X_Y,X_Z
    REAL(rp),DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: V_X,V_Y,V_Z
    REAL(rp),DIMENSION(params%pchunk) :: B_X,B_Y,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_X,E_Y,E_Z
    REAL(rp),DIMENSION(params%pchunk) :: PSIp
    INTEGER(is),DIMENSION(params%pchunk) :: flagCon,flagCol
    REAL(rp) :: a,m_cache,q_cache    
    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp,pchunk
    !! Particles iterator.
    INTEGER                                                    :: cc
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    !! time iterator.     

    do ii = 1_idef,params%num_species      

       pchunk=params%pchunk
       m_cache=spp(ii)%m
       q_cache=spp(ii)%q
       a = params%dt


       !$OMP PARALLEL DO default(none) &
       !$OMP& FIRSTPRIVATE(a,m_cache,q_cache,pchunk,E0) &
       !$OMP& shared(params,ii,spp,P,F) &
       !$OMP& PRIVATE(pp,tt,Bmag,cc,X_X,X_Y,X_Z,V_X,V_Y,V_Z,B_X,B_Y,B_Z, &
       !$OMP& E_X,E_Y,E_Z,b_unit_X,b_unit_Y,b_unit_Z,v,vpar,vperp,tmp, &
       !$OMP& cross_X,cross_Y,cross_Z,vec_X,vec_Y,vec_Z,g, &
       !$OMP& Y_R,Y_PHI,Y_Z,flagCon,flagCol,PSIp)
       do pp=1_idef,spp(ii)%ppp,pchunk

          !$OMP SIMD
          do cc=1_idef,pchunk
             X_X(cc)=spp(ii)%vars%X(pp-1+cc,1)
             X_Y(cc)=spp(ii)%vars%X(pp-1+cc,2)
             X_Z(cc)=spp(ii)%vars%X(pp-1+cc,3)

             Y_R(cc)=0._rp
             Y_PHI(cc)=0._rp
             Y_Z(cc)=0._rp

             V_X(cc)=spp(ii)%vars%V(pp-1+cc,1)
             V_Y(cc)=spp(ii)%vars%V(pp-1+cc,2)
             V_Z(cc)=spp(ii)%vars%V(pp-1+cc,3)

             B_X(cc)=spp(ii)%vars%B(pp-1+cc,1)
             B_Y(cc)=spp(ii)%vars%B(pp-1+cc,2)
             B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)

             E_X(cc)=spp(ii)%vars%E(pp-1+cc,1)
             E_Y(cc)=spp(ii)%vars%E(pp-1+cc,2)
             E_Z(cc)=spp(ii)%vars%E(pp-1+cc,3)
             
             PSIp(cc)=spp(ii)%vars%PSI_P(pp-1+cc)

             g(cc)=spp(ii)%vars%g(pp-1+cc)

             flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
             flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)
          end do
          !$OMP END SIMD

          if (.not.params%FokPlan) then
             do tt=1_ip,params%t_skip

                call cart_to_cyl_p(pchunk,X_X,X_Y,X_Z,Y_R,Y_PHI,Y_Z)

                if (F%axisymmetric_fields.and. &
                     (params%orbit_model(3:3).eq.'B')) then
                   call interp_FOfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_X,B_Y,B_Z, &
                        E_X,E_Y,E_Z,PSIp,flagCon)
                else if ((.not.F%axisymmetric_fields).and. &
                     (params%orbit_model(3:3).eq.'B')) then
                   call interp_FO3Dfields_p(pchunk,F,Y_R,Y_PHI,Y_Z, &
                        B_X,B_Y,B_Z,E_X,E_Y,E_Z,flagCon)
                else if (params%orbit_model(3:5).eq.'psi') then
                   call interp_FOfields1_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_X,B_Y,B_Z, &
                        E_X,E_Y,E_Z,PSIp,flagCon)
                end if


                !               write(output_unit_write,'("B_X: ",E17.10)') B_X(1)
                !               write(output_unit_write,'("B_Y: ",E17.10)') B_Y(1)
                !               write(output_unit_write,'("B_Z: ",E17.10)') B_Z(1)

                call advance_FOinterp_vars(tt,a,q_cache,m_cache,params, &
                     X_X,X_Y,X_Z,V_X,V_Y,V_Z,B_X,B_Y,B_Z,E_X,E_Y,E_Z, &
                     g,flagCon,flagCol,P,F,PSIp)
             end do !timestep iterator

             !$OMP SIMD
             do cc=1_idef,pchunk
                spp(ii)%vars%X(pp-1+cc,1)=X_X(cc)
                spp(ii)%vars%X(pp-1+cc,2)=X_Y(cc)
                spp(ii)%vars%X(pp-1+cc,3)=X_Z(cc)

                spp(ii)%vars%Y(pp-1+cc,1)=Y_R(cc)
                spp(ii)%vars%Y(pp-1+cc,2)=Y_PHI(cc)
                spp(ii)%vars%Y(pp-1+cc,3)=Y_Z(cc)
                
                spp(ii)%vars%V(pp-1+cc,1)=V_X(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_Y(cc)
                spp(ii)%vars%V(pp-1+cc,3)=V_Z(cc)

                spp(ii)%vars%g(pp-1+cc) = g(cc)
                spp(ii)%vars%flagCon(pp-1+cc) = flagCon(cc)
                spp(ii)%vars%flagCol(pp-1+cc) = flagCol(cc)

                spp(ii)%vars%B(pp-1+cc,1) = B_X(cc)
                spp(ii)%vars%B(pp-1+cc,2) = B_Y(cc)
                spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                spp(ii)%vars%E(pp-1+cc,1) = E_X(cc)
                spp(ii)%vars%E(pp-1+cc,2) = E_Y(cc)
                spp(ii)%vars%E(pp-1+cc,3) = E_Z(cc)

                spp(ii)%vars%PSI_P(pp-1+cc) = PSIp(cc)
             end do
             !$OMP END SIMD

          else
             !$OMP SIMD
             do cc=1_idef,pchunk
                B_X(cc)=spp(ii)%vars%B(pp-1+cc,1)
                B_Y(cc)=spp(ii)%vars%B(pp-1+cc,2)
                B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)

                E_X(cc)=spp(ii)%vars%E(pp-1+cc,1)
                E_Y(cc)=spp(ii)%vars%E(pp-1+cc,2)
                E_Z(cc)=spp(ii)%vars%E(pp-1+cc,3)
             end do
             !$OMP END SIMD

             call advance_FP3Dinterp_vars(params,X_X,X_Y,X_Z,V_X,V_Y,V_Z, &
                  g,m_cache,B_X,B_Y,B_Z,E_X,E_Y,E_Z,flagCon,flagCol,P,F,PSIp)

             !$OMP SIMD
             do cc=1_idef,pchunk

                spp(ii)%vars%V(pp-1+cc,1)=V_X(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_Y(cc)
                spp(ii)%vars%V(pp-1+cc,3)=V_Z(cc)

                spp(ii)%vars%g(pp-1+cc) = g(cc)
                spp(ii)%vars%flagCol(pp-1+cc) = flagCol(cc)

             end do
             !$OMP END SIMD
          end if

          !$OMP SIMD
          !          !$OMP& aligned(Bmag,B_X,B_Y,B_Z, &
          !          !$OMP& b_unit_X,b_unit_Y,b_unit_Z,v,V_X,V_Y,V_Z,vpar, &
          !          !$OMP& vperp,tmp,cross_X,cross_Y,cross_Z, &
          !          !$OMP& vec_X,vec_Y,vec_Z,E_X,E_Y,E_Z)
          do cc=1_idef,pchunk
             !Derived output data
             Bmag(cc) = SQRT(B_X(cc)*B_X(cc)+B_Y(cc)*B_Y(cc)+B_Z(cc)*B_Z(cc))

             ! Parallel unit vector
             b_unit_X(cc) = B_X(cc)/Bmag(cc)
             b_unit_Y(cc) = B_Y(cc)/Bmag(cc)
             b_unit_Z(cc) = B_Z(cc)/Bmag(cc)

             v(cc) = SQRT(V_X(cc)*V_X(cc)+V_Y(cc)*V_Y(cc)+V_Z(cc)*V_Z(cc))
             if (v(cc).GT.korc_zero) then
                ! Parallel and perpendicular components of velocity
                vpar(cc) = (V_X(cc)*b_unit_X(cc)+V_Y(cc)*b_unit_Y(cc)+ &
                     V_Z(cc)*b_unit_Z(cc))

                vperp(cc) =  v(cc)**2 - vpar(cc)**2
                if ( vperp(cc) .GE. korc_zero ) then
                   vperp(cc) = SQRT( vperp(cc) )
                else
                   vperp(cc) = 0.0_rp
                end if

                ! Pitch angle
                spp(ii)%vars%eta(pp-1+cc) = 180.0_rp* &
                     MODULO(ATAN2(vperp(cc),vpar(cc)),2.0_rp*C_PI)/C_PI

                ! Magnetic moment
                spp(ii)%vars%mu(pp-1+cc) = 0.5_rp*m_cache* &
                     g(cc)**2*vperp(cc)**2/Bmag(cc)
                ! See Northrop's book (The adiabatic motion of charged
                ! particles)

                ! Radiated power
                tmp(cc) = q_cache**4/(6.0_rp*C_PI*E0*m_cache**2)

                cross_X(cc) = V_Y(cc)*B_Z(cc)-V_Z(cc)*B_Y(cc)
                cross_Y(cc) = V_Z(cc)*B_X(cc)-V_X(cc)*B_Z(cc)
                cross_Z(cc) = V_X(cc)*B_Y(cc)-V_Y(cc)*B_X(cc)

                vec_X(cc) = E_X(cc) + cross_X(cc)
                vec_Y(cc) = E_Y(cc) + cross_Y(cc)
                vec_Z(cc) = E_Z(cc) + cross_Z(cc)

                spp(ii)%vars%Prad(pp-1+cc) = tmp(cc)* &
                     ( E_X(cc)*E_X(cc)+E_Y(cc)*E_Y(cc)+E_Z(cc)*E_Z(cc) + &
                     cross_X(cc)*E_X(cc)+cross_Y(cc)*E_Y(cc)+ &
                     cross_Z(cc)*E_Z(cc) + g(cc)**2* &
                     ((E_X(cc)*V_X(cc)+E_Y(cc)*V_Y(cc)+E_Z(cc)*V_Z(cc))**2 &
                     - vec_X(cc)*vec_X(cc)+vec_Y(cc)*vec_Y(cc)+ &
                     vec_Z(cc)*vec_Z(cc)) )

                ! Input power due to electric field
                spp(ii)%vars%Pin(pp-1+cc) = q_cache*(E_X(cc)*V_X(cc)+ &
                     E_Y(cc)*V_Y(cc)+E_Z(cc)*V_Z(cc))
             else
                spp(ii)%vars%eta(pp-1+cc) = 0.0_rp
                spp(ii)%vars%mu(pp-1+cc) = 0.0_rp
                spp(ii)%vars%Prad(pp-1+cc) = 0.0_rp
                spp(ii)%vars%Pin(pp-1+cc) = 0.0_rp
             end if

          end do
          !$OMP END SIMD


       end do !particle chunk iterator
       !$OMP END PARALLEL DO

    end do !species iterator

  end subroutine adv_FOinterp_top

  subroutine adv_FOinterp_mars_top(params,F,P,spp)  
    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                                   :: F
    !! An instance of the KORC derived type FIELDS.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    !! An instance of the KORC derived type PROFILES.
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(params%pchunk) :: Bmag
    REAL(rp), DIMENSION(params%pchunk) :: b_unit_X,b_unit_Y,b_unit_Z
    REAL(rp), DIMENSION(params%pchunk) :: v,vpar,vperp
    REAL(rp), DIMENSION(params%pchunk) :: tmp
    REAL(rp), DIMENSION(params%pchunk) :: g
    REAL(rp), DIMENSION(params%pchunk) :: cross_X,cross_Y,cross_Z
    REAL(rp), DIMENSION(params%pchunk) :: vec_X,vec_Y,vec_Z
    REAL(rp), DIMENSION(params%pchunk) :: X_X,X_Y,X_Z
    REAL(rp), DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp), DIMENSION(params%pchunk) :: V_X,V_Y,V_Z
    REAL(rp), DIMENSION(params%pchunk) :: B_X,B_Y,B_Z
    REAL(rp), DIMENSION(params%pchunk) :: E_X,E_Y,E_Z
    REAL(rp), DIMENSION(params%pchunk) :: PSIp
    INTEGER(is),DIMENSION(params%pchunk) :: flagCon,flagCol
    REAL(rp) :: a,m_cache,q_cache    
    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp,pchunk
    !! Particles iterator.
    INTEGER                                                    :: cc
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    !! time iterator.     

    do ii = 1_idef,params%num_species      

       pchunk=params%pchunk
       m_cache=spp(ii)%m
       q_cache=spp(ii)%q
       a =params%dt


       !$OMP PARALLEL DO default(none) &
       !$OMP& FIRSTPRIVATE(a,m_cache,q_cache,pchunk,E0) &
       !$OMP& shared(params,ii,spp,P,F) &
       !$OMP& PRIVATE(pp,tt,Bmag,cc,X_X,X_Y,X_Z,V_X,V_Y,V_Z,B_X,B_Y,B_Z, &
       !$OMP& E_X,E_Y,E_Z,b_unit_X,b_unit_Y,b_unit_Z,v,vpar,vperp,tmp, &
       !$OMP& cross_X,cross_Y,cross_Z,vec_X,vec_Y,vec_Z,g, &
       !$OMP& Y_R,Y_PHI,Y_Z,flagCon,flagCol,PSIp)
       do pp=1_idef,spp(ii)%ppp,pchunk

          !$OMP SIMD
          do cc=1_idef,pchunk
             X_X(cc)=spp(ii)%vars%X(pp-1+cc,1)
             X_Y(cc)=spp(ii)%vars%X(pp-1+cc,2)
             X_Z(cc)=spp(ii)%vars%X(pp-1+cc,3)

             Y_R(cc)=0._rp
             Y_PHI(cc)=0._rp
             Y_Z(cc)=0._rp
             
             V_X(cc)=spp(ii)%vars%V(pp-1+cc,1)
             V_Y(cc)=spp(ii)%vars%V(pp-1+cc,2)
             V_Z(cc)=spp(ii)%vars%V(pp-1+cc,3)

             B_X(cc)=spp(ii)%vars%B(pp-1+cc,1)
             B_Y(cc)=spp(ii)%vars%B(pp-1+cc,2)
             B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)
             
             E_X(cc)=spp(ii)%vars%E(pp-1+cc,1)
             E_Y(cc)=spp(ii)%vars%E(pp-1+cc,2)
             E_Z(cc)=spp(ii)%vars%E(pp-1+cc,3)

             PSIp(cc)=spp(ii)%vars%PSI_P(pp-1+cc)

             g(cc)=spp(ii)%vars%g(pp-1+cc)

             flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
             flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)
          end do
          !$OMP END SIMD


          do tt=1_ip,params%t_skip

             call cart_to_cyl_p(pchunk,X_X,X_Y,X_Z,Y_R,Y_PHI,Y_Z)

             call interp_FOfields_mars_p(params,pchunk,F,Y_R,Y_PHI,Y_Z, &
                  B_X,B_Y,B_Z,PSIp,flagCon)



             !               write(output_unit_write,'("B_X: ",E17.10)') B_X(1)
             !               write(output_unit_write,'("B_Y: ",E17.10)') B_Y(1)
             !               write(output_unit_write,'("B_Z: ",E17.10)') B_Z(1)

             call advance_FOinterp_vars(tt,a,q_cache,m_cache,params, &
                  X_X,X_Y,X_Z,V_X,V_Y,V_Z,B_X,B_Y,B_Z,E_X,E_Y,E_Z, &
                  g,flagCon,flagCol,P,F,PSIp)
          end do !timestep iterator

          !$OMP SIMD
          do cc=1_idef,pchunk
             spp(ii)%vars%X(pp-1+cc,1)=X_X(cc)
             spp(ii)%vars%X(pp-1+cc,2)=X_Y(cc)
             spp(ii)%vars%X(pp-1+cc,3)=X_Z(cc)

             spp(ii)%vars%V(pp-1+cc,1)=V_X(cc)
             spp(ii)%vars%V(pp-1+cc,2)=V_Y(cc)
             spp(ii)%vars%V(pp-1+cc,3)=V_Z(cc)

             spp(ii)%vars%g(pp-1+cc) = g(cc)
             spp(ii)%vars%flagCon(pp-1+cc) = flagCon(cc)
             spp(ii)%vars%flagCol(pp-1+cc) = flagCol(cc)

             spp(ii)%vars%B(pp-1+cc,1) = B_X(cc)
             spp(ii)%vars%B(pp-1+cc,2) = B_Y(cc)
             spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

             spp(ii)%vars%E(pp-1+cc,1) = E_X(cc)
             spp(ii)%vars%E(pp-1+cc,2) = E_Y(cc)
             spp(ii)%vars%E(pp-1+cc,3) = E_Z(cc)

             spp(ii)%vars%PSI_P(pp-1+cc) = PSIp(cc)
          end do
          !$OMP END SIMD


          !$OMP SIMD
          !          !$OMP& aligned(Bmag,B_X,B_Y,B_Z, &
          !          !$OMP& b_unit_X,b_unit_Y,b_unit_Z,v,V_X,V_Y,V_Z,vpar, &
          !          !$OMP& vperp,tmp,cross_X,cross_Y,cross_Z, &
          !          !$OMP& vec_X,vec_Y,vec_Z,E_X,E_Y,E_Z)
          do cc=1_idef,pchunk
             !Derived output data
             Bmag(cc) = SQRT(B_X(cc)*B_X(cc)+B_Y(cc)*B_Y(cc)+B_Z(cc)*B_Z(cc))

             ! Parallel unit vector
             b_unit_X(cc) = B_X(cc)/Bmag(cc)
             b_unit_Y(cc) = B_Y(cc)/Bmag(cc)
             b_unit_Z(cc) = B_Z(cc)/Bmag(cc)

             v(cc) = SQRT(V_X(cc)*V_X(cc)+V_Y(cc)*V_Y(cc)+V_Z(cc)*V_Z(cc))
             if (v(cc).GT.korc_zero) then
                ! Parallel and perpendicular components of velocity
                vpar(cc) = (V_X(cc)*b_unit_X(cc)+V_Y(cc)*b_unit_Y(cc)+ &
                     V_Z(cc)*b_unit_Z(cc))

                vperp(cc) =  v(cc)**2 - vpar(cc)**2
                if ( vperp(cc) .GE. korc_zero ) then
                   vperp(cc) = SQRT( vperp(cc) )
                else
                   vperp(cc) = 0.0_rp
                end if

                ! Pitch angle
                spp(ii)%vars%eta(pp-1+cc) = 180.0_rp* &
                     MODULO(ATAN2(vperp(cc),vpar(cc)),2.0_rp*C_PI)/C_PI

                ! Magnetic moment
                spp(ii)%vars%mu(pp-1+cc) = 0.5_rp*m_cache* &
                     g(cc)**2*vperp(cc)**2/Bmag(cc)
                ! See Northrop's book (The adiabatic motion of charged
                ! particles)


                
                ! Radiated power
                tmp(cc) = q_cache**4/(6.0_rp*C_PI*E0*m_cache**2)

                cross_X(cc) = V_Y(cc)*B_Z(cc)-V_Z(cc)*B_Y(cc)
                cross_Y(cc) = V_Z(cc)*B_X(cc)-V_X(cc)*B_Z(cc)
                cross_Z(cc) = V_X(cc)*B_Y(cc)-V_Y(cc)*B_X(cc)

                vec_X(cc) = E_X(cc) + cross_X(cc)
                vec_Y(cc) = E_Y(cc) + cross_Y(cc)
                vec_Z(cc) = E_Z(cc) + cross_Z(cc)

                spp(ii)%vars%Prad(pp-1+cc) = tmp(cc)* &
                     ( E_X(cc)*E_X(cc)+E_Y(cc)*E_Y(cc)+E_Z(cc)*E_Z(cc) + &
                     cross_X(cc)*E_X(cc)+cross_Y(cc)*E_Y(cc)+ &
                     cross_Z(cc)*E_Z(cc) + g(cc)**2* &
                     ((E_X(cc)*V_X(cc)+E_Y(cc)*V_Y(cc)+E_Z(cc)*V_Z(cc))**2 &
                     - vec_X(cc)*vec_X(cc)+vec_Y(cc)*vec_Y(cc)+ &
                     vec_Z(cc)*vec_Z(cc)) )

                ! Input power due to electric field
                spp(ii)%vars%Pin(pp-1+cc) = q_cache*(E_X(cc)*V_X(cc)+ &
                     E_Y(cc)*V_Y(cc)+E_Z(cc)*V_Z(cc))
             else
                spp(ii)%vars%eta(pp-1+cc) = 0.0_rp
                spp(ii)%vars%mu(pp-1+cc) = 0.0_rp
                spp(ii)%vars%Prad(pp-1+cc) = 0.0_rp
                spp(ii)%vars%Pin(pp-1+cc) = 0.0_rp
             end if

          end do
          !$OMP END SIMD


       end do !particle chunk iterator
       !$OMP END PARALLEL DO

    end do !species iterator

  end subroutine adv_FOinterp_mars_top
  
  subroutine adv_FOinterp_aorsa_top(params,F,P,spp)  
    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)                                   :: F
    !! An instance of the KORC derived type FIELDS.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    !! An instance of the KORC derived type PROFILES.
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(params%pchunk) :: Bmag
    REAL(rp), DIMENSION(params%pchunk) :: b_unit_X,b_unit_Y,b_unit_Z
    REAL(rp), DIMENSION(params%pchunk) :: v,vpar,vperp
    REAL(rp), DIMENSION(params%pchunk) :: tmp
    REAL(rp), DIMENSION(params%pchunk) :: g
    REAL(rp), DIMENSION(params%pchunk) :: cross_X,cross_Y,cross_Z
    REAL(rp), DIMENSION(params%pchunk) :: vec_X,vec_Y,vec_Z
    REAL(rp), DIMENSION(params%pchunk) :: X_X,X_Y,X_Z
    REAL(rp), DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp), DIMENSION(params%pchunk) :: V_X,V_Y,V_Z
    REAL(rp), DIMENSION(params%pchunk) :: B_X,B_Y,B_Z
    REAL(rp), DIMENSION(params%pchunk) :: E_X,E_Y,E_Z
    REAL(rp), DIMENSION(params%pchunk) :: PSIp
    INTEGER(is),DIMENSION(params%pchunk) :: flagCon,flagCol
    REAL(rp) :: a,m_cache,q_cache,time    
    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp,pchunk
    !! Particles iterator.
    INTEGER                                                    :: cc
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    !! time iterator.     

    do ii = 1_idef,params%num_species      

       pchunk=params%pchunk
       m_cache=spp(ii)%m
       q_cache=spp(ii)%q
       a = params%dt


       !$OMP PARALLEL DO default(none) &
       !$OMP& FIRSTPRIVATE(a,m_cache,q_cache,pchunk,E0) &
       !$OMP& shared(params,ii,spp,P,F) &
       !$OMP& PRIVATE(pp,tt,Bmag,cc,X_X,X_Y,X_Z,V_X,V_Y,V_Z,B_X,B_Y,B_Z, &
       !$OMP& E_X,E_Y,E_Z,b_unit_X,b_unit_Y,b_unit_Z,v,vpar,vperp,tmp, &
       !$OMP& cross_X,cross_Y,cross_Z,vec_X,vec_Y,vec_Z,g, &
       !$OMP& Y_R,Y_PHI,Y_Z,flagCon,flagCol,PSIp,time)
       do pp=1_idef,spp(ii)%ppp,pchunk

          !$OMP SIMD
          do cc=1_idef,pchunk
             X_X(cc)=spp(ii)%vars%X(pp-1+cc,1)
             X_Y(cc)=spp(ii)%vars%X(pp-1+cc,2)
             X_Z(cc)=spp(ii)%vars%X(pp-1+cc,3)

             Y_R(cc)=0._rp
             Y_PHI(cc)=0._rp
             Y_Z(cc)=0._rp
             
             V_X(cc)=spp(ii)%vars%V(pp-1+cc,1)
             V_Y(cc)=spp(ii)%vars%V(pp-1+cc,2)
             V_Z(cc)=spp(ii)%vars%V(pp-1+cc,3)

             B_X(cc)=spp(ii)%vars%B(pp-1+cc,1)
             B_Y(cc)=spp(ii)%vars%B(pp-1+cc,2)
             B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)
             
             E_X(cc)=spp(ii)%vars%E(pp-1+cc,1)
             E_Y(cc)=spp(ii)%vars%E(pp-1+cc,2)
             E_Z(cc)=spp(ii)%vars%E(pp-1+cc,3)

             PSIp(cc)=spp(ii)%vars%PSI_P(pp-1+cc)

             g(cc)=spp(ii)%vars%g(pp-1+cc)

             flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
             flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)
          end do
          !$OMP END SIMD


          do tt=1_ip,params%t_skip

             time=(params%time+params%dt*tt)*params%cpp%time
             
             call cart_to_cyl_p(pchunk,X_X,X_Y,X_Z,Y_R,Y_PHI,Y_Z)

             call interp_FOfields_aorsa_p(time,params,pchunk,F,Y_R,Y_PHI,Y_Z, &
                  B_X,B_Y,B_Z,E_X,E_Y,E_Z,PSIp,flagCon)



             !               write(output_unit_write,'("B_X: ",E17.10)') B_X(1)
             !               write(output_unit_write,'("B_Y: ",E17.10)') B_Y(1)
             !               write(output_unit_write,'("B_Z: ",E17.10)') B_Z(1)

             call advance_FOinterp_vars(tt,a,q_cache,m_cache,params, &
                  X_X,X_Y,X_Z,V_X,V_Y,V_Z,B_X,B_Y,B_Z,E_X,E_Y,E_Z, &
                  g,flagCon,flagCol,P,F,PSIp)
          end do !timestep iterator

          !$OMP SIMD
          do cc=1_idef,pchunk
             spp(ii)%vars%X(pp-1+cc,1)=X_X(cc)
             spp(ii)%vars%X(pp-1+cc,2)=X_Y(cc)
             spp(ii)%vars%X(pp-1+cc,3)=X_Z(cc)

             spp(ii)%vars%V(pp-1+cc,1)=V_X(cc)
             spp(ii)%vars%V(pp-1+cc,2)=V_Y(cc)
             spp(ii)%vars%V(pp-1+cc,3)=V_Z(cc)

             spp(ii)%vars%g(pp-1+cc) = g(cc)
             spp(ii)%vars%flagCon(pp-1+cc) = flagCon(cc)
             spp(ii)%vars%flagCol(pp-1+cc) = flagCol(cc)

             spp(ii)%vars%B(pp-1+cc,1) = B_X(cc)
             spp(ii)%vars%B(pp-1+cc,2) = B_Y(cc)
             spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

             spp(ii)%vars%E(pp-1+cc,1) = E_X(cc)
             spp(ii)%vars%E(pp-1+cc,2) = E_Y(cc)
             spp(ii)%vars%E(pp-1+cc,3) = E_Z(cc)

             spp(ii)%vars%PSI_P(pp-1+cc) = PSIp(cc)
          end do
          !$OMP END SIMD


          !$OMP SIMD
          !          !$OMP& aligned(Bmag,B_X,B_Y,B_Z, &
          !          !$OMP& b_unit_X,b_unit_Y,b_unit_Z,v,V_X,V_Y,V_Z,vpar, &
          !          !$OMP& vperp,tmp,cross_X,cross_Y,cross_Z, &
          !          !$OMP& vec_X,vec_Y,vec_Z,E_X,E_Y,E_Z)
          do cc=1_idef,pchunk
             !Derived output data
             Bmag(cc) = SQRT(B_X(cc)*B_X(cc)+B_Y(cc)*B_Y(cc)+B_Z(cc)*B_Z(cc))

             ! Parallel unit vector
             b_unit_X(cc) = B_X(cc)/Bmag(cc)
             b_unit_Y(cc) = B_Y(cc)/Bmag(cc)
             b_unit_Z(cc) = B_Z(cc)/Bmag(cc)

             v(cc) = SQRT(V_X(cc)*V_X(cc)+V_Y(cc)*V_Y(cc)+V_Z(cc)*V_Z(cc))
             if (v(cc).GT.korc_zero) then
                ! Parallel and perpendicular components of velocity
                vpar(cc) = (V_X(cc)*b_unit_X(cc)+V_Y(cc)*b_unit_Y(cc)+ &
                     V_Z(cc)*b_unit_Z(cc))

                vperp(cc) =  v(cc)**2 - vpar(cc)**2
                if ( vperp(cc) .GE. korc_zero ) then
                   vperp(cc) = SQRT( vperp(cc) )
                else
                   vperp(cc) = 0.0_rp
                end if

                ! Pitch angle
                spp(ii)%vars%eta(pp-1+cc) = 180.0_rp* &
                     MODULO(ATAN2(vperp(cc),vpar(cc)),2.0_rp*C_PI)/C_PI

                ! Magnetic moment
                spp(ii)%vars%mu(pp-1+cc) = 0.5_rp*m_cache* &
                     g(cc)**2*vperp(cc)**2/Bmag(cc)
                ! See Northrop's book (The adiabatic motion of charged
                ! particles)


                
                ! Radiated power
                tmp(cc) = q_cache**4/(6.0_rp*C_PI*E0*m_cache**2)

                cross_X(cc) = V_Y(cc)*B_Z(cc)-V_Z(cc)*B_Y(cc)
                cross_Y(cc) = V_Z(cc)*B_X(cc)-V_X(cc)*B_Z(cc)
                cross_Z(cc) = V_X(cc)*B_Y(cc)-V_Y(cc)*B_X(cc)

                vec_X(cc) = E_X(cc) + cross_X(cc)
                vec_Y(cc) = E_Y(cc) + cross_Y(cc)
                vec_Z(cc) = E_Z(cc) + cross_Z(cc)

                spp(ii)%vars%Prad(pp-1+cc) = tmp(cc)* &
                     ( E_X(cc)*E_X(cc)+E_Y(cc)*E_Y(cc)+E_Z(cc)*E_Z(cc) + &
                     cross_X(cc)*E_X(cc)+cross_Y(cc)*E_Y(cc)+ &
                     cross_Z(cc)*E_Z(cc) + g(cc)**2* &
                     ((E_X(cc)*V_X(cc)+E_Y(cc)*V_Y(cc)+E_Z(cc)*V_Z(cc))**2 &
                     - vec_X(cc)*vec_X(cc)+vec_Y(cc)*vec_Y(cc)+ &
                     vec_Z(cc)*vec_Z(cc)) )

                ! Input power due to electric field
                spp(ii)%vars%Pin(pp-1+cc) = q_cache*(E_X(cc)*V_X(cc)+ &
                     E_Y(cc)*V_Y(cc)+E_Z(cc)*V_Z(cc))
             else
                spp(ii)%vars%eta(pp-1+cc) = 0.0_rp
                spp(ii)%vars%mu(pp-1+cc) = 0.0_rp
                spp(ii)%vars%Prad(pp-1+cc) = 0.0_rp
                spp(ii)%vars%Pin(pp-1+cc) = 0.0_rp
             end if

          end do
          !$OMP END SIMD


       end do !particle chunk iterator
       !$OMP END PARALLEL DO

    end do !species iterator

  end subroutine adv_FOinterp_aorsa_top
  
  subroutine advance_FOinterp_vars(tt,a,q_cache,m_cache,params,X_X,X_Y,X_Z, &
       V_X,V_Y,V_Z,B_X,B_Y,B_Z,E_X,E_Y,E_Z,g,flagCon,flagCol,P,F,PSIp)
    TYPE(KORC_PARAMS), INTENT(IN)                              :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)      :: F
    INTEGER(ip), INTENT(IN)                                       :: tt
    !! Time step used in the leapfrog step (\(\Delta t\)).
    REAL(rp)                                      :: dt
    !! Time step used in the leapfrog step (\(\Delta t\)).
    REAL(rp), INTENT(IN)                                       :: m_cache,q_cache
    !! Time step used in the leapfrog step (\(\Delta t\)).

    REAL(rp),DIMENSION(params%pchunk)                                  :: Bmag

    REAL(rp),INTENT(in)                                       :: a
    !! This variable is used to simplify notation in the code, and
    !! is given by \(a=q\Delta t/m\),
    REAL(rp),DIMENSION(params%pchunk)                                    :: sigma
    !! This variable is \(\sigma = \gamma'^2 - \tau^2\) in the above equations.
    REAL(rp),DIMENSION(params%pchunk)                               :: us
    !! This variable is \(u^{*} = p^{*}/m\) where \( p^{*} =
    !! \mathbf{p}'\cdot \mathbf{\tau}/mc\).
    !! Variable 'u^*' in Vay, J.-L. PoP (2008).
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT)                 :: g
    REAL(rp),DIMENSION(params%pchunk) :: gp,g0
    !! Relativistic factor \(\gamma\).
    REAL(rp),DIMENSION(params%pchunk)                                 :: s
    !! This variable is \(s = 1/(1+t^2)\) in the equations above.
    !! Variable 's' in Vay, J.-L. PoP (200params%pchunk).
    REAL(rp),DIMENSION(params%pchunk)                            :: U_hs_X,U_hs_Y,U_hs_Z
    !! Is \(\mathbf{u}=\mathbf{p}/m\) at half-time step (\(i+1/2\)) in
    !! the absence of radiation losses or collisions. \(\mathbf{u}^{i+1/2} =
    !! \mathbf{u}^i + \frac{q\Delta t}{2m}\left( \mathbf{E}^{i+1/2} +
    !! \mathbf{v}^i\times \mathbf{B}^{i+1/2} \right)\).
    REAL(rp),DIMENSION(params%pchunk)                           :: tau_X,tau_Y,tau_Z
    !! This variable is \(\mathbf{\tau} = (q\Delta t/2)\mathbf{B}^{i+1/2}\).
    REAL(rp),DIMENSION(params%pchunk)                            :: up_X,up_Y,up_Z
    !! This variable is \(\mathbf{u}'= \mathbf{p}'/m\), where \(\mathbf{p}'
    !! = \mathbf{p}^i + q\Delta t \left( \mathbf{E}^{i+1/2} +
    !! \frac{\mathbf{v}^i}{2} \times \mathbf{B}^{i+1/2} \right)\).
    REAL(rp),DIMENSION(params%pchunk)                                     :: t_X,t_Y,t_Z
    !! This variable is \(\mathbf{t} = {\mathbf \tau}/\gamma^{i+1}\).
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT)                     :: X_X,X_Y,X_Z
    REAL(rp),DIMENSION(params%pchunk)                    :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT)                      :: V_X,V_Y,V_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(IN)                      :: B_X,B_Y,B_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(IN)                     :: E_X,E_Y,E_Z,PSIp
    REAL(rp),DIMENSION(params%pchunk)                     :: U_L_X,U_L_Y,U_L_Z
    REAL(rp),DIMENSION(params%pchunk)                     :: U_X,U_Y,U_Z
    REAL(rp),DIMENSION(params%pchunk)                     :: U_RC_X,U_RC_Y,U_RC_Z
    REAL(rp),DIMENSION(params%pchunk)                     :: U_os_X,U_os_Y,U_os_Z
    !! This variable is \(\mathbf{u}^{i+1}= \mathbf{p}^{i+1}/m\).
    REAL(rp),DIMENSION(params%pchunk)                          :: cross_X,cross_Y,cross_Z

    REAL(rp), DIMENSION(params%pchunk)                       :: Frad_X,Frad_Y,Frad_Z
    !! Synchrotron radiation reaction force of each particle.

    REAL(rp),DIMENSION(params%pchunk) :: ne,Te,Zeff

    INTEGER                                      :: cc,pchunk
    !! Chunk iterator.

    INTEGER(is) ,DIMENSION(params%pchunk),intent(inout)                   :: flagCon,flagCol

    dt=params%dt
    pchunk=params%pchunk


    !$OMP SIMD
    !    !$OMP& aligned(g0,g,U_X,U_Y,U_Z,V_X,V_Y,V_Z,Bmag,B_X,B_Y,B_Z, &
    !    !$OMP& U_L_X,U_L_Y,U_L_Z,U_RC_X,U_RC_Y,U_RC_Z, &
    !    !$OMP& cross_X,cross_Y,cross_Z,U_hs_X,U_hs_Y,U_hs_Z,E_X,E_Y,E_Z, &
    !    !$OMP& tau_X,tau_Y,tau_Z,up_X,up_Y,up_Z,gp,sigma,us,t_X,t_Y,t_Z,s, &
    !    !$OMP& U_os_X,U_os_Y,U_os_Z,Frad_X,Frad_Y,Frad_Z)
    do cc=1_idef,pchunk

       g0(cc)=g(cc)

       U_X(cc) = g(cc)*V_X(cc)
       U_Y(cc) = g(cc)*V_Y(cc)
       U_Z(cc) = g(cc)*V_Z(cc)


       ! Magnitude of magnetic field
       Bmag(cc) = SQRT(B_X(cc)*B_X(cc)+B_Y(cc)*B_Y(cc)+B_Z(cc)*B_Z(cc))

       U_L_X(cc)=U_X(cc)
       U_L_Y(cc)=U_Y(cc)
       U_L_Z(cc)=U_Z(cc)

       U_RC_X(cc)=U_X(cc)
       U_RC_Y(cc)=U_Y(cc)
       U_RC_Z(cc)=U_Z(cc)

       ! LEAP-FROG SCHEME FOR LORENTZ FORCE !

       cross_X(cc)=V_Y(cc)*B_Z(cc)-V_Z(cc)*B_Y(cc)
       cross_Y(cc)=V_Z(cc)*B_X(cc)-V_X(cc)*B_Z(cc)
       cross_Z(cc)=V_X(cc)*B_Y(cc)-V_Y(cc)*B_X(cc)



       U_hs_X(cc) = U_L_X(cc) + 0.5_rp*a*(E_X(cc) +cross_X(cc))
       U_hs_Y(cc) = U_L_Y(cc) + 0.5_rp*a*(E_Y(cc) +cross_Y(cc))
       U_hs_Z(cc) = U_L_Z(cc) + 0.5_rp*a*(E_Z(cc) +cross_Z(cc))



       tau_X(cc) = 0.5_rp*a*B_X(cc)
       tau_Y(cc) = 0.5_rp*a*B_Y(cc)
       tau_Z(cc) = 0.5_rp*a*B_Z(cc)



       up_X(cc) = U_hs_X(cc) + 0.5_rp*a*E_X(cc)
       up_Y(cc) = U_hs_Y(cc) + 0.5_rp*a*E_Y(cc)
       up_Z(cc) = U_hs_Z(cc) + 0.5_rp*a*E_Z(cc)

       gp(cc) = SQRT( 1.0_rp + up_X(cc)*up_X(cc)+up_Y(cc)*up_Y(cc)+ &
            up_Z(cc)*up_Z(cc) )

       sigma(cc) = gp(cc)*gp(cc) - (tau_X(cc)*tau_X(cc)+ &
            tau_Y(cc)*tau_Y(cc)+tau_Z(cc)*tau_Z(cc))

       us(cc) = up_X(cc)*tau_X(cc)+up_Y(cc)*tau_Y(cc)+ &
            up_Z(cc)*tau_Z(cc)

       ! variable 'u^*' in Vay, J.-L. PoP (2008)
       g(cc) = SQRT( 0.5_rp*(sigma(cc) + SQRT(sigma(cc)*sigma(cc) + &
            4.0_rp*(tau_X(cc)*tau_X(cc)+tau_Y(cc)*tau_Y(cc)+ &
            tau_Z(cc)*tau_Z(cc) + us(cc)*us(cc)))) )

       t_X(cc) = tau_X(cc)/g(cc)
       t_Y(cc) = tau_Y(cc)/g(cc)
       t_Z(cc) = tau_Z(cc)/g(cc)


       s(cc) = 1.0_rp/(1.0_rp + t_X(cc)*t_X(cc)+t_Y(cc)*t_Y(cc)+ &
            t_Z(cc)*t_Z(cc))
       ! variable 's' in Vay, J.-L. PoP (2008)

       cross_X(cc)=up_Y(cc)*t_Z(cc)-up_Z(cc)*t_Y(cc)
       cross_Y(cc)=up_Z(cc)*t_X(cc)-up_X(cc)*t_Z(cc)
       cross_Z(cc)=up_X(cc)*t_Y(cc)-up_Y(cc)*t_X(cc)

       U_L_X(cc) = s(cc)*(up_X(cc) + (up_X(cc)*t_X(cc)+ &
            up_Y(cc)*t_Y(cc)+up_Z(cc)*t_Z(cc))*t_X(cc) + cross_X(cc))
       U_L_Y(cc) = s(cc)*(up_Y(cc) + (up_X(cc)*t_X(cc)+ &
            up_Y(cc)*t_Y(cc)+up_Z(cc)*t_Z(cc))*t_Y(cc) + cross_Y(cc))
       U_L_Z(cc) = s(cc)*(up_Z(cc) + (up_X(cc)*t_X(cc)+ &
            up_Y(cc)*t_Y(cc)+up_Z(cc)*t_Z(cc))*t_Z(cc) + cross_Z(cc))
       ! LEAP-FROG SCHEME FOR LORENTZ FORCE !

       U_os_X(cc) = 0.5_rp*(U_L_X(cc) + U_X(cc))
       U_os_Y(cc) = 0.5_rp*(U_L_Y(cc) + U_Y(cc))
       U_os_Z(cc) = 0.5_rp*(U_L_Z(cc) + U_Z(cc))
       ! Splitting operator for including radiation

       if (params%radiation) then
          !! Calls [[radiation_force_p]] in [[korc_ppusher]].
          call radiation_force_p(pchunk,q_cache,m_cache,U_os_X,U_os_Y,U_os_Z, &
               E_X,E_Y,E_Z,B_Z,B_Y,B_Z,Frad_X,Frad_Y,Frad_Z)
          U_RC_X(cc) = U_RC_X(cc) + a*Frad_X(cc)/q_cache
          U_RC_Y(cc) = U_RC_Y(cc) + a*Frad_Y(cc)/q_cache
          U_RC_Z(cc) = U_RC_Z(cc) + a*Frad_Z(cc)/q_cache
       end if
       ! Splitting operator for including radiation

       U_X(cc) = U_L_X(cc) + U_RC_X(cc) - U_X(cc)
       U_Y(cc) = U_L_Y(cc) + U_RC_Y(cc) - U_Y(cc)
       U_Z(cc) = U_L_Z(cc) + U_RC_Z(cc) - U_Z(cc)

    end do
    !$OMP END SIMD

    if (params%collisions) then

       call include_CoulombCollisions_FO_p(tt,params,X_X,X_Y,X_Z, &
            U_X,U_Y,U_Z,B_X,B_Y,B_Z,m_cache,P,F,flagCon,flagCol,PSIp)

    end if

    if (params%radiation.or.params%collisions) then

       !$OMP SIMD
       !       !$OMP& aligned(g,U_X,U_Y,U_Z)
       do cc=1_idef,pchunk
          g(cc)=sqrt(1._rp+U_X(cc)*U_X(cc)+U_Y(cc)*U_Y(cc)+U_Z(cc)*U_Z(cc))
       end do
       !$OMP END SIMD

    end if

    !$OMP SIMD
    !    !$OMP& aligned(g,g0,V_X,V_Y,V_Z,U_X,U_Y,U_Z,X_X,X_Y,X_Z,flagCon,flagCol)
    do cc=1_idef,pchunk

       if ((flagCon(cc).eq.0_is).or.(flagCol(cc).eq.0_is)) then
          g(cc)=g0(cc)
       else
          V_X(cc) = U_X(cc)/g(cc)
          V_Y(cc) = U_Y(cc)/g(cc)
          V_Z(cc) = U_Z(cc)/g(cc)
       end if

       X_X(cc) = X_X(cc) + dt*V_X(cc)*REAL(flagCon(cc))*REAL(flagCol(cc))
       X_Y(cc) = X_Y(cc) + dt*V_Y(cc)*REAL(flagCon(cc))*REAL(flagCol(cc))
       X_Z(cc) = X_Z(cc) + dt*V_Z(cc)*REAL(flagCon(cc))*REAL(flagCol(cc))

    end do
    !$OMP END SIMD

  end subroutine advance_FOinterp_vars
#endif
  
#ifdef FIO
  subroutine advance_FOfio_vars(tt,a,q_cache,m_cache,params,X_X,X_Y,X_Z, &
       V_X,V_Y,V_Z,B_X,B_Y,B_Z,E_X,E_Y,E_Z,g,flagCon,flagCol,P,F,PSIp,hint)
    TYPE(KORC_PARAMS), INTENT(IN)                              :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)      :: F
    INTEGER(ip), INTENT(IN)                                       :: tt
    !! Time step used in the leapfrog step (\(\Delta t\)).
    REAL(rp)                                      :: dt
    !! Time step used in the leapfrog step (\(\Delta t\)).
    REAL(rp), INTENT(IN)                                       :: m_cache,q_cache
    !! Time step used in the leapfrog step (\(\Delta t\)).

    REAL(rp),DIMENSION(params%pchunk)                                  :: Bmag

    REAL(rp),INTENT(in)                                       :: a
    !! This variable is used to simplify notation in the code, and
    !! is given by \(a=q\Delta t/m\),
    REAL(rp),DIMENSION(params%pchunk)                                    :: sigma
    !! This variable is \(\sigma = \gamma'^2 - \tau^2\) in the above equations.
    REAL(rp),DIMENSION(params%pchunk)                               :: us
    !! This variable is \(u^{*} = p^{*}/m\) where \( p^{*} =
    !! \mathbf{p}'\cdot \mathbf{\tau}/mc\).
    !! Variable 'u^*' in Vay, J.-L. PoP (2008).
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT)                 :: g
    REAL(rp),DIMENSION(params%pchunk) :: gp,g0
    !! Relativistic factor \(\gamma\).
    REAL(rp),DIMENSION(params%pchunk)                                 :: s
    !! This variable is \(s = 1/(1+t^2)\) in the equations above.
    !! Variable 's' in Vay, J.-L. PoP (2008).
    REAL(rp),DIMENSION(params%pchunk)                            :: U_hs_X,U_hs_Y,U_hs_Z
    !! Is \(\mathbf{u}=\mathbf{p}/m\) at half-time step (\(i+1/2\)) in
    !! the absence of radiation losses or collisions. \(\mathbf{u}^{i+1/2} =
    !! \mathbf{u}^i + \frac{q\Delta t}{2m}\left( \mathbf{E}^{i+1/2} +
    !! \mathbf{v}^i\times \mathbf{B}^{i+1/2} \right)\).
    REAL(rp),DIMENSION(params%pchunk)                           :: tau_X,tau_Y,tau_Z
    !! This variable is \(\mathbf{\tau} = (q\Delta t/2)\mathbf{B}^{i+1/2}\).
    REAL(rp),DIMENSION(params%pchunk)                            :: up_X,up_Y,up_Z
    !! This variable is \(\mathbf{u}'= \mathbf{p}'/m\), where \(\mathbf{p}'
    !! = \mathbf{p}^i + q\Delta t \left( \mathbf{E}^{i+1/2} +
    !! \frac{\mathbf{v}^i}{2} \times \mathbf{B}^{i+1/2} \right)\).
    REAL(rp),DIMENSION(params%pchunk)                                     :: t_X,t_Y,t_Z
    !! This variable is \(\mathbf{t} = {\mathbf \tau}/\gamma^{i+1}\).
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT)                     :: X_X,X_Y,X_Z
    REAL(rp),DIMENSION(params%pchunk)                    :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT)                      :: V_X,V_Y,V_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(IN)                      :: B_X,B_Y,B_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(IN)                     :: E_X,E_Y,E_Z,PSIp
    REAL(rp),DIMENSION(params%pchunk)                     :: U_L_X,U_L_Y,U_L_Z
    REAL(rp),DIMENSION(params%pchunk)                     :: U_X,U_Y,U_Z
    REAL(rp),DIMENSION(params%pchunk)                     :: U_RC_X,U_RC_Y,U_RC_Z
    REAL(rp),DIMENSION(params%pchunk)                     :: U_os_X,U_os_Y,U_os_Z
    !! This variable is \(\mathbf{u}^{i+1}= \mathbf{p}^{i+1}/m\).
    REAL(rp),DIMENSION(params%pchunk)                          :: cross_X,cross_Y,cross_Z

    REAL(rp), DIMENSION(params%pchunk)                       :: Frad_X,Frad_Y,Frad_Z
    !! Synchrotron radiation reaction force of each particle.

    REAL(rp),DIMENSION(params%pchunk) :: ne,Te,Zeff

    INTEGER                                      :: cc,pchunk
    !! Chunk iterator.

    INTEGER(is) ,DIMENSION(params%pchunk),intent(inout)                   :: flagCon,flagCol
    TYPE(C_PTR),DIMENSION(params%pchunk) :: hint

    dt=params%dt
    pchunk=params%pchunk


    !$OMP SIMD
    !    !$OMP& aligned(g0,g,U_X,U_Y,U_Z,V_X,V_Y,V_Z,Bmag,B_X,B_Y,B_Z, &
    !    !$OMP& U_L_X,U_L_Y,U_L_Z,U_RC_X,U_RC_Y,U_RC_Z, &
    !    !$OMP& cross_X,cross_Y,cross_Z,U_hs_X,U_hs_Y,U_hs_Z,E_X,E_Y,E_Z, &
    !    !$OMP& tau_X,tau_Y,tau_Z,up_X,up_Y,up_Z,gp,sigma,us,t_X,t_Y,t_Z,s, &
    !    !$OMP& U_os_X,U_os_Y,U_os_Z,Frad_X,Frad_Y,Frad_Z)
    do cc=1_idef,pchunk

       g0(cc)=g(cc)

       U_X(cc) = g(cc)*V_X(cc)
       U_Y(cc) = g(cc)*V_Y(cc)
       U_Z(cc) = g(cc)*V_Z(cc)


       ! Magnitude of magnetic field
       Bmag(cc) = SQRT(B_X(cc)*B_X(cc)+B_Y(cc)*B_Y(cc)+B_Z(cc)*B_Z(cc))

       U_L_X(cc)=U_X(cc)
       U_L_Y(cc)=U_Y(cc)
       U_L_Z(cc)=U_Z(cc)

       U_RC_X(cc)=U_X(cc)
       U_RC_Y(cc)=U_Y(cc)
       U_RC_Z(cc)=U_Z(cc)

       ! LEAP-FROG SCHEME FOR LORENTZ FORCE !

       cross_X(cc)=V_Y(cc)*B_Z(cc)-V_Z(cc)*B_Y(cc)
       cross_Y(cc)=V_Z(cc)*B_X(cc)-V_X(cc)*B_Z(cc)
       cross_Z(cc)=V_X(cc)*B_Y(cc)-V_Y(cc)*B_X(cc)



       U_hs_X(cc) = U_L_X(cc) + 0.5_rp*a*(E_X(cc) +cross_X(cc))
       U_hs_Y(cc) = U_L_Y(cc) + 0.5_rp*a*(E_Y(cc) +cross_Y(cc))
       U_hs_Z(cc) = U_L_Z(cc) + 0.5_rp*a*(E_Z(cc) +cross_Z(cc))



       tau_X(cc) = 0.5_rp*a*B_X(cc)
       tau_Y(cc) = 0.5_rp*a*B_Y(cc)
       tau_Z(cc) = 0.5_rp*a*B_Z(cc)



       up_X(cc) = U_hs_X(cc) + 0.5_rp*a*E_X(cc)
       up_Y(cc) = U_hs_Y(cc) + 0.5_rp*a*E_Y(cc)
       up_Z(cc) = U_hs_Z(cc) + 0.5_rp*a*E_Z(cc)

       gp(cc) = SQRT( 1.0_rp + up_X(cc)*up_X(cc)+up_Y(cc)*up_Y(cc)+ &
            up_Z(cc)*up_Z(cc) )

       sigma(cc) = gp(cc)*gp(cc) - (tau_X(cc)*tau_X(cc)+ &
            tau_Y(cc)*tau_Y(cc)+tau_Z(cc)*tau_Z(cc))

       us(cc) = up_X(cc)*tau_X(cc)+up_Y(cc)*tau_Y(cc)+ &
            up_Z(cc)*tau_Z(cc)

       ! variable 'u^*' in Vay, J.-L. PoP (2008)
       g(cc) = SQRT( 0.5_rp*(sigma(cc) + SQRT(sigma(cc)*sigma(cc) + &
            4.0_rp*(tau_X(cc)*tau_X(cc)+tau_Y(cc)*tau_Y(cc)+ &
            tau_Z(cc)*tau_Z(cc) + us(cc)*us(cc)))) )

       t_X(cc) = tau_X(cc)/g(cc)
       t_Y(cc) = tau_Y(cc)/g(cc)
       t_Z(cc) = tau_Z(cc)/g(cc)


       s(cc) = 1.0_rp/(1.0_rp + t_X(cc)*t_X(cc)+t_Y(cc)*t_Y(cc)+ &
            t_Z(cc)*t_Z(cc))
       ! variable 's' in Vay, J.-L. PoP (2008)

       cross_X(cc)=up_Y(cc)*t_Z(cc)-up_Z(cc)*t_Y(cc)
       cross_Y(cc)=up_Z(cc)*t_X(cc)-up_X(cc)*t_Z(cc)
       cross_Z(cc)=up_X(cc)*t_Y(cc)-up_Y(cc)*t_X(cc)

       U_L_X(cc) = s(cc)*(up_X(cc) + (up_X(cc)*t_X(cc)+ &
            up_Y(cc)*t_Y(cc)+up_Z(cc)*t_Z(cc))*t_X(cc) + cross_X(cc))
       U_L_Y(cc) = s(cc)*(up_Y(cc) + (up_X(cc)*t_X(cc)+ &
            up_Y(cc)*t_Y(cc)+up_Z(cc)*t_Z(cc))*t_Y(cc) + cross_Y(cc))
       U_L_Z(cc) = s(cc)*(up_Z(cc) + (up_X(cc)*t_X(cc)+ &
            up_Y(cc)*t_Y(cc)+up_Z(cc)*t_Z(cc))*t_Z(cc) + cross_Z(cc))
       ! LEAP-FROG SCHEME FOR LORENTZ FORCE !

       U_os_X(cc) = 0.5_rp*(U_L_X(cc) + U_X(cc))
       U_os_Y(cc) = 0.5_rp*(U_L_Y(cc) + U_Y(cc))
       U_os_Z(cc) = 0.5_rp*(U_L_Z(cc) + U_Z(cc))
       ! Splitting operator for including radiation

       if (params%radiation) then
          !! Calls [[radiation_force_p]] in [[korc_ppusher]].
          call radiation_force_p(pchunk,q_cache,m_cache,U_os_X,U_os_Y,U_os_Z, &
               E_X,E_Y,E_Z,B_Z,B_Y,B_Z,Frad_X,Frad_Y,Frad_Z)
          U_RC_X(cc) = U_RC_X(cc) + a*Frad_X(cc)/q_cache
          U_RC_Y(cc) = U_RC_Y(cc) + a*Frad_Y(cc)/q_cache
          U_RC_Z(cc) = U_RC_Z(cc) + a*Frad_Z(cc)/q_cache
       end if
       ! Splitting operator for including radiation

       U_X(cc) = U_L_X(cc) + U_RC_X(cc) - U_X(cc)
       U_Y(cc) = U_L_Y(cc) + U_RC_Y(cc) - U_Y(cc)
       U_Z(cc) = U_L_Z(cc) + U_RC_Z(cc) - U_Z(cc)

    end do
    !$OMP END SIMD

    if (params%collisions) then

       call include_CoulombCollisions_FOfio_p(tt,params,X_X,X_Y,X_Z, &
            U_X,U_Y,U_Z,B_X,B_Y,B_Z,m_cache,P,F,flagCon,flagCol,PSIp,hint)

    end if

    if (params%radiation.or.params%collisions) then

       !$OMP SIMD
       !       !$OMP& aligned(g,U_X,U_Y,U_Z)
       do cc=1_idef,pchunk
          g(cc)=sqrt(1._rp+U_X(cc)*U_X(cc)+U_Y(cc)*U_Y(cc)+U_Z(cc)*U_Z(cc))
       end do
       !$OMP END SIMD

    end if

    !$OMP SIMD
    !    !$OMP& aligned(g,g0,V_X,V_Y,V_Z,U_X,U_Y,U_Z,X_X,X_Y,X_Z,flagCon,flagCol)
    do cc=1_idef,pchunk

       if ((flagCon(cc).eq.0_is).or.(flagCol(cc).eq.0_is)) then
          g(cc)=g0(cc)
       else
          V_X(cc) = U_X(cc)/g(cc)
          V_Y(cc) = U_Y(cc)/g(cc)
          V_Z(cc) = U_Z(cc)/g(cc)
       end if

       X_X(cc) = X_X(cc) + dt*V_X(cc)*REAL(flagCon(cc))*REAL(flagCol(cc))
       X_Y(cc) = X_Y(cc) + dt*V_Y(cc)*REAL(flagCon(cc))*REAL(flagCol(cc))
       X_Z(cc) = X_Z(cc) + dt*V_Z(cc)*REAL(flagCon(cc))*REAL(flagCol(cc))

    end do
    !$OMP END SIMD

  end subroutine advance_FOfio_vars
#endif

  subroutine advance_FP3Dinterp_vars(params,X_X,X_Y,X_Z,V_X,V_Y,V_Z,g, &
       m_cache,B_X,B_Y,B_Z,E_X,E_Y,E_Z,flagCon,flagCol,P,F,PSIp)    
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)      :: F
    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    !! time iterator.
    REAL(rp),DIMENSION(params%pchunk), INTENT(IN)  :: X_X,X_Y,X_Z
    REAL(rp),DIMENSION(params%pchunk)  :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk), INTENT(IN)  :: E_X,E_Y,E_Z,PSIp
    REAL(rp),DIMENSION(params%pchunk), INTENT(IN)  :: B_X,B_Y,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: U_X,U_Y,U_Z
    REAL(rp),DIMENSION(params%pchunk), INTENT(INOUT)  :: V_X,V_Y,V_Z
    REAL(rp),DIMENSION(params%pchunk) :: ne,Te,Zeff
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: g
    INTEGER(is),DIMENSION(params%pchunk),INTENT(INOUT) :: flagCon,flagCol
    REAL(rp),intent(in) :: m_cache

    pchunk=params%pchunk  

    !$OMP SIMD
    !    !$OMP& aligned(U_X,U_Y,U_Z,V_X,V_Y,V_Z,g)
    do cc=1_idef,pchunk
       U_X(cc)=V_X(cc)*g(cc)
       U_Y(cc)=V_Y(cc)*g(cc)
       U_Z(cc)=V_Z(cc)*g(cc)
    end do
    !$OMP END SIMD

    do tt=1_ip,params%t_skip

       call include_CoulombCollisions_FO_p(tt,params,X_X,X_Y,X_Z, &
            U_X,U_Y,U_Z,B_X,B_Y,B_Z,m_cache,P,F,flagCon,flagCol,PSIp)

    end do

    !$OMP SIMD
    !    !$OMP& aligned(U_X,U_Y,U_Z,V_X,V_Y,V_Z,g)
    do cc=1_idef,pchunk

       g(cc)=sqrt(1._rp+U_X(cc)*U_X(cc)+U_Y(cc)*U_Y(cc)+U_Z(cc)*U_Z(cc))

       V_X(cc)=U_X(cc)/g(cc)
       V_Y(cc)=U_Y(cc)/g(cc)
       V_Z(cc)=U_Z(cc)/g(cc)
    end do
    !$OMP END SIMD

  end subroutine advance_FP3Dinterp_vars


  subroutine GC_init(params,F,spp)
    !! @note Subroutine to advance GC variables \(({\bf X},p_\parallel)\)
    !! @endnote
    !! Comment this section further with evolution equations, numerical
    !! methods, and descriptions of both.
    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(INOUT)                                   :: F
    !! An instance of the KORC derived type FIELDS.

    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.

    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    REAL(rp)               :: Bmag1,pmag
    REAL(rp)               :: Bmagc
    REAL(rp)               :: rm
    REAL(rp),DIMENSION(:,:),ALLOCATABLE               :: RAphi
    REAL(rp), DIMENSION(3) :: bhat
    REAL(rp), DIMENSION(3) :: bhatc
    REAL(rp), DIMENSION(params%pchunk) :: E_PHI
    REAL(rp),DIMENSION(:),ALLOCATABLE               :: RVphi

    REAL(rp),DIMENSION(params%pchunk) :: rm8,Y_R,Y_Z,V_PLL,vpll,gam
    real(rp),dimension(F%dim_1D) :: Vpart,Vpartave,VpartOMP
    real(rp) :: dr
    integer :: rind

    !    write(output_unit_write,'("eta",E17.10)') spp(ii)%vars%eta(pp)
    !    write(output_unit_write,'("gam",E17.10)') spp(ii)%vars%g(pp)

    do ii = 1_idef,params%num_species

       pchunk=params%pchunk

       if (spp(ii)%spatial_distribution.eq.'TRACER'.and. &
            params%FO_GC_compare) then
          call get_fields(params,spp(ii)%vars,F)
          !! Calls [[get_fields]] in [[korc_fields]].
          ! Interpolates fields at local particles' position and keeps in
          ! spp%vars. Fields in (R,\(\phi\),Z) coordinates. 

          ALLOCATE(RAphi(spp(ii)%ppp,2))
          ALLOCATE(RVphi(spp(ii)%ppp))
          RAphi=0.0_rp

          call cart_to_cyl(spp(ii)%vars%X,spp(ii)%vars%Y)

          !$OMP PARALLEL DO SHARED(params,ii,spp,F,RAphi,RVphi) &
          !$OMP&  PRIVATE(pp,Bmag1,bhat,rm)
          ! Call OpenMP to calculate p_par and mu for each particle and
          ! put into spp%vars%V
          do pp=1_idef,spp(ii)%ppp
             if ( spp(ii)%vars%flagCon(pp) .EQ. 1_is ) then

                RVphi(pp)=(-sin(spp(ii)%vars%Y(pp,2))*spp(ii)%vars%V(pp,1)+ &
                     cos(spp(ii)%vars%Y(pp,2))*spp(ii)%vars%V(pp,2))* &
                     spp(ii)%vars%Y(pp,1)

                Bmag1 = SQRT(spp(ii)%vars%B(pp,1)*spp(ii)%vars%B(pp,1)+ &
                     spp(ii)%vars%B(pp,2)*spp(ii)%vars%B(pp,2)+ &
                     spp(ii)%vars%B(pp,3)*spp(ii)%vars%B(pp,3))

                !             write(output_unit_write,'("pp: ",I16)') pp
                !             write(output_unit_write,'("Bmag: ",E17.10)') Bmag


                bhat = spp(ii)%vars%B(pp,:)/Bmag1

                if (params%field_model(1:10).eq.'ANALYTICAL') then
                   rm=sqrt((spp(ii)%vars%Y(pp,1)-F%AB%Ro)**2+ &
                        (spp(ii)%vars%Y(pp,3))**2)

                   RAphi(pp,1)=-F%AB%lambda**2*F%AB%Bo/(2*F%AB%qo)* &
                        log(1+(rm/F%AB%lambda)**2)

                else if (params%field_model(1:8).eq.'EXTERNAL') then

                   RAphi(pp,1)=spp(ii)%vars%PSI_P(pp)/(2*C_PI)

                end if

                !             write(output_unit_write,'("bhat: ",E17.10)') bhat 
                !             write(output_unit_write,'("V: ",E17.10)') spp(ii)%vars%V(pp,:)


                spp(ii)%vars%X(pp,:)=spp(ii)%vars%X(pp,:)- &
                     spp(ii)%m*spp(ii)%vars%g(pp)* &
                     cross(bhat,spp(ii)%vars%V(pp,:))/(spp(ii)%q*Bmag1)

                ! transforming from particle location to associated
                ! GC location

             end if ! if particle in domain, i.e. spp%vars%flagCon==1
          end do ! loop over particles on an mpi process
          !$OMP END PARALLEL DO

          call cart_to_cyl(spp(ii)%vars%X,spp(ii)%vars%Y)
          call get_fields(params,spp(ii)%vars,F)

          !$OMP PARALLEL DO SHARED(params,ii,spp,F,RAphi,RVphi) &
          !$OMP&  PRIVATE(pp,rm)
          ! Call OpenMP to calculate p_par and mu for each particle and
          ! put into spp%vars%V
          do pp=1_idef,spp(ii)%ppp
             if ( spp(ii)%vars%flagCon(pp) .EQ. 1_is ) then

                if (params%field_model(1:10).eq.'ANALYTICAL') then
                   rm=sqrt((spp(ii)%vars%Y(pp,1)-F%AB%Ro)**2+ &
                        (spp(ii)%vars%Y(pp,3))**2)
                   RAphi(pp,2)=-F%AB%lambda**2*F%AB%Bo/(2*F%AB%qo)* &
                        log(1+(rm/F%AB%lambda)**2)

                else if (params%field_model(1:8).eq.'EXTERNAL') then

                   RAphi(pp,2)=spp(ii)%vars%PSI_P(pp)/(2*C_PI)

                end if

                write(output_unit_write,'("RAphi1: ",E17.10)') RAphi(pp,1)
                write(output_unit_write,'("RAphi2: ",E17.10)') RAphi(pp,2)

                spp(ii)%vars%V(pp,1)=(spp(ii)%m*spp(ii)%vars%g(pp)* &
                     RVphi(pp)+spp(ii)%q*(RAphi(pp,1)-RAphi(pp,2)))/ &
                     spp(ii)%vars%Y(pp,1)
                !GC ppar              

             end if ! if particle in domain, i.e. spp%vars%flagCon==1
          end do ! loop over particles on an mpi process
          !$OMP END PARALLEL DO

          !$OMP PARALLEL DO SHARED(ii,spp) PRIVATE(pp,Bmagc,bhatc)
          ! Call OpenMP to calculate p_par and mu for each particle and
          ! put into spp%vars%V
          do pp=1_idef,spp(ii)%ppp
             if ( spp(ii)%vars%flagCon(pp) .EQ. 1_is ) then

                Bmagc = SQRT( DOT_PRODUCT(spp(ii)%vars%B(pp,:), &
                     spp(ii)%vars%B(pp,:)))

                bhatc = spp(ii)%vars%B(pp,:)/Bmagc

                spp(ii)%vars%V(pp,1)=spp(ii)%vars%V(pp,1)/ &
                     bhatc(2)
                !GC ppar

                spp(ii)%vars%V(pp,2)=spp(ii)%m/(2*Bmagc)* &
                     (spp(ii)%vars%g(pp)**2- &
                     (1+(spp(ii)%vars%V(pp,1)/spp(ii)%m)**2))           
                !GC mu


             end if ! if particle in domain, i.e. spp%vars%flagCon==1
          end do ! loop over particles on an mpi process
          !$OMP END PARALLEL DO

          params%GC_coords=.TRUE.
          DEALLOCATE(RAphi)
          DEALLOCATE(RVphi)

          !Preparing Output Data
          call get_fields(params,spp(ii)%vars,F)

          !$OMP PARALLEL DO shared(F,params,spp) &
          !$OMP& PRIVATE(cc,pp,E_PHI,Y_R) firstprivate(pchunk)
          do pp=1_idef,spp(ii)%ppp,pchunk

             !$OMP SIMD
             do cc=1_idef,pchunk
                E_PHI(cc)=spp(ii)%vars%E(pp-1+cc,2)
                Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)                
             end do
             !$OMP END SIMD

             call add_analytical_E_p(params,0_ip,F,E_PHI,Y_R,Y_Z)


             !$OMP SIMD
             do cc=1_idef,pchunk
                spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
             end do
             !$OMP END SIMD

          end do
          !$OMP END PARALLEL DO


          !$OMP PARALLEL DO SHARED(ii,spp) PRIVATE(pp,Bmag1)
          ! Call OpenMP to calculate p_par and mu for each particle and
          ! put into spp%vars%V
          do pp=1_idef,spp(ii)%ppp
             if ( spp(ii)%vars%flagCon(pp) .EQ. 1_is ) then

                Bmag1 = SQRT( DOT_PRODUCT(spp(ii)%vars%B(pp,:), &
                     spp(ii)%vars%B(pp,:)))

                spp(ii)%vars%g(pp)=sqrt(1+(spp(ii)%vars%V(pp,1))**2+ &
                     2*spp(ii)%vars%V(pp,2)*Bmag1)

                !                write(output_unit_write,'("Bmag:",E17.10)') Bmag1
                !                write(output_unit_write,'("PPLL:",E17.10)') spp(ii)%vars%V(pp,1)
                !                write(output_unit_write,'("MU:",E17.10)') spp(ii)%vars%V(pp,2)

                spp(ii)%vars%eta(pp) = atan2(sqrt(2*spp(ii)%m*Bmag1* &
                     spp(ii)%vars%V(pp,2)),spp(ii)%vars%V(pp,1))*180.0_rp/C_PI

                !                             write(output_unit_write,'("BR",E17.10)') spp(ii)%vars%B(pp,1)
                !                             write(output_unit_write,'("BPHI",E17.10)') spp(ii)%vars%B(pp,2)
                !                             write(output_unit_write,'("BZ",E17.10)') spp(ii)%vars%B(pp,3)

                !             write(output_unit_write,'("ppll",E17.10)') spp(ii)%vars%V(pp,1)
                !             write(output_unit_write,'("pperp",E17.10)') sqrt(2*spp(ii)%m*Bmag1* &
                !                  spp(ii)%vars%V(pp,2))

                !                             write(output_unit_write,'("eta GCinit",E17.10)') spp(ii)%vars%eta(pp)
                !             write(output_unit_write,'("gam",E17.10)') spp(ii)%vars%g(pp)


             end if ! if particle in domain, i.e. spp%vars%flagCon==1
          end do ! loop over particles on an mpi process
          !$OMP END PARALLEL DO                
       else

          if ((spp(ii)%spatial_distribution.eq.'TRACER').or. &
               (spp(ii)%spatial_distribution.eq.'TORUS').or. &
               (spp(ii)%spatial_distribution.eq.'DISK').or. &
               (spp(ii)%spatial_distribution.eq. &
               '2D-GAUSSIAN-ELLIPTIC-TORUS-MH')) &
               call cart_to_cyl(spp(ii)%vars%X,spp(ii)%vars%Y)

          params%GC_coords=.TRUE.

          do pp=1_idef,spp(ii)%ppp
             spp(ii)%vars%E(pp,1)=0._rp
             spp(ii)%vars%E(pp,2)=0._rp
             spp(ii)%vars%E(pp,3)=0._rp
          end do


          !write(6,*) 'before second get fields'
          call get_fields(params,spp(ii)%vars,F)        
          !write(6,*) 'after second get fields'
          
          !write(output_unit_write,*) spp(1)%vars%PSI_P

          !$OMP PARALLEL DO SHARED(ii,spp) PRIVATE(pp,Bmag1)

          do pp=1_idef,spp(ii)%ppp
             !             if ( spp(ii)%vars%flagCon(pp) .EQ. 1_is ) then

             !                write(output_unit_write,'("BR: ",E17.10)') spp(ii)%vars%B(pp,1)
             !                write(output_unit_write,'("BPHI: ",E17.10)') spp(ii)%vars%B(pp,2)
             !                write(output_unit_write,'("BZ: ",E17.10)') spp(ii)%vars%B(pp,3)

             Bmag1 = SQRT( DOT_PRODUCT(spp(ii)%vars%B(pp,:), &
                  spp(ii)%vars%B(pp,:)))

             pmag=sqrt(spp(ii)%vars%g(pp)**2-1)

             spp(ii)%vars%V(pp,1)=pmag*cos(deg2rad(spp(ii)%vars%eta(pp)))

             !write(6,*)  'GC_init'
             !write(6,*) spp(ii)%m,Bmag1
             
             spp(ii)%vars%V(pp,2)=(pmag* &
                  sin(deg2rad(spp(ii)%vars%eta(pp))))**2/ &
                  (2*spp(ii)%m*Bmag1)

             !    write(output_unit_write,'("BR",E17.10)') spp(ii)%vars%B(pp,1)
             !    write(output_unit_write,'("BPHI",E17.10)') spp(ii)%vars%B(pp,2)
             !    write(output_unit_write,'("BZ",E17.10)') spp(ii)%vars%B(pp,3)

             !write(output_unit_write,'("ppll",E17.10)') spp(ii)%vars%V(pp,1)
             !write(output_unit_write,'("mu",E17.10)') spp(ii)%vars%V(pp,2)

             !     write(output_unit_write,'("eta",E17.10)') spp(ii)%vars%eta(pp)
             !     write(output_unit_write,'("gam",E17.10)') spp(ii)%vars%g(pp)


             !             end if ! if particle in domain, i.e. spp%vars%flagCon==1
          end do ! loop over particles on an mpi process
          !$OMP END PARALLEL DO  

          !$OMP PARALLEL DO shared(F,params,spp) &
          !$OMP& PRIVATE(pp,cc,E_PHI,Y_R) firstprivate(pchunk)
          do pp=1_idef,spp(ii)%ppp,pchunk

             !$OMP SIMD
             do cc=1_idef,pchunk
                E_PHI(cc)=spp(ii)%vars%E(pp-1+cc,2)
                Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)                
             end do
             !$OMP END SIMD

             if (params%field_model(1:8).eq.'EXTERNAL') then
                call add_analytical_E_p(params,0_ip,F,E_PHI,Y_R,Y_Z)
             end if

             !$OMP SIMD
             do cc=1_idef,pchunk

                spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
             end do
             !$OMP END SIMD

          end do
          !$OMP END PARALLEL DO

       end if

    end do ! loop over particle species

  end subroutine GC_init

  FUNCTION deg2rad(x)
    REAL(rp), INTENT(IN) :: x
    REAL(rp) :: deg2rad

    deg2rad = C_PI*x/180.0_rp
  END FUNCTION deg2rad

  FUNCTION rad2deg(x)
    REAL(rp), INTENT(IN) :: x
    REAL(rp) :: rad2deg

    rad2deg = x*180.0_rp/C_PI
  END FUNCTION rad2deg

  subroutine adv_GCeqn_top(params,F,P,spp)

    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(INOUT)                                   :: F
    !! An instance of the KORC derived type FIELDS.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    !! An instance of the KORC derived type PROFILES.
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(params%pchunk)               :: Bmag
    REAL(rp),DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z,E_PHI
    REAL(rp),DIMENSION(params%pchunk) :: PSIp,ne,Te
    REAL(rp),DIMENSION(params%pchunk) :: V_PLL,V_MU
    REAL(rp) :: B0,EF0,R0,q0,lam,ar,m_cache,q_cache,ne0,Te0,Zeff0
    INTEGER(is),DIMENSION(params%pchunk)  :: flagCon,flagCol,flagRE

    LOGICAL                                                    :: ss_collisions
    !! Logical variable that indicates if collisions are included in
    !! the simulation.

    INTEGER           :: ii
    !! Species iterator.
    INTEGER           :: pp
    !! Particles iterator.
    INTEGER           :: cc,pchunk,achunk
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    INTEGER(ip)                                                    :: ttt,tttt
    !! time iterator.
    real(rp),dimension(F%dim_1D) :: Vden,Vdenave,VdenOMP
    INTEGER :: newREs


    do ii = 1_idef,params%num_species      

       pchunk=params%pchunk
       q_cache=spp(ii)%q
       m_cache=spp(ii)%m


       do ttt=1_ip,params%t_it_SC

          VdenOMP=0._rp

          if(.not.params%LargeCollisions) then

             !$OMP PARALLEL DO default(none) &
             !$OMP& FIRSTPRIVATE(E0,q_cache,m_cache,pchunk) &
             !$OMP& shared(F,P,params,ii,spp) &
             !$OMP& PRIVATE(pp,tt,ttt,Bmag,cc,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
             !$OMP& flagCon,flagCol,B_R,B_PHI,B_Z,E_PHI,PSIp,ne, &
             !$OMP& Vden,Vdenave) &
             !$OMP& REDUCTION(+:VdenOMP)
             do pp=1_idef,spp(ii)%ppp,pchunk


                !$OMP SIMD
                do cc=1_idef,pchunk
                   Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)
                   Y_PHI(cc)=spp(ii)%vars%Y(pp-1+cc,2)
                   Y_Z(cc)=spp(ii)%vars%Y(pp-1+cc,3)

                   V_PLL(cc)=spp(ii)%vars%V(pp-1+cc,1)
                   V_MU(cc)=spp(ii)%vars%V(pp-1+cc,2)

                   PSIp(cc)=spp(ii)%vars%PSI_p(pp-1+cc)
                   ne(cc)=spp(ii)%vars%ne(pp-1+cc)

                   flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
                   flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)
                end do
                !$OMP END SIMD

                if (.not.params%FokPlan) then
                   Vdenave=0._rp
                   do tt=1_ip,params%t_skip

                      !                   write(output_unit_write,*) params%mpi_params%rank,'Y_R',Y_R

                      call advance_GCeqn_vars(spp(ii)%vars,pp, &
                           tt+params%t_skip*(ttt-1),params, &
                           Y_R,Y_PHI, Y_Z,V_PLL,V_MU,flagCon,flagCol,q_cache,m_cache, &
                           B_R,B_PHI,B_Z,F,P,PSIp,E_PHI)

                      !                   write(output_unit_write,*) params%mpi_params%rank,'Y_R',Y_R

                      if (params%collisions) then

                         call include_CoulombCollisions_GC_p(tt+params%t_skip*(ttt-1),params, &
                              Y_R,Y_PHI,Y_Z,V_PLL,V_MU,m_cache,flagCon,flagCol, &
                              F,P,E_PHI,ne,PSIp)

                      end if
                      
                      if (params%SC_E) then                  
                         call calculate_SC_p(params,F,B_R,B_PHI,B_Z,Y_R,Y_Z, &
                              V_PLL,V_MU,m_cache,flagCon,flagCol,Vden)
                         Vdenave=(Vdenave*REAL(tt-1_ip)+Vden)/REAL(tt)
                      end if

                   end do !timestep iterator

                   VdenOMP=VdenOMP+Vdenave


                   !$OMP SIMD
                   do cc=1_idef,pchunk
                      spp(ii)%vars%Y(pp-1+cc,1)=Y_R(cc)
                      spp(ii)%vars%Y(pp-1+cc,2)=Y_PHI(cc)
                      spp(ii)%vars%Y(pp-1+cc,3)=Y_Z(cc)

                      spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                      spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                      spp(ii)%vars%flagCon(pp-1+cc)=flagCon(cc)
                      spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                      spp(ii)%vars%B(pp-1+cc,1) = B_R(cc)
                      spp(ii)%vars%B(pp-1+cc,2) = B_PHI(cc)
                      spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                      spp(ii)%vars%PSI_P(pp-1+cc) = PSIp(cc)
                      spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
                   end do
                   !$OMP END SIMD

                else

                   do tt=1_ip,params%t_skip


                      !if (mod(tt,50).eq.0) then
                      !write(output_unit_write,*) 'iteration',tt
                      !flush(output_unit_write)
                      !endif

                      
                      call include_CoulombCollisions_GC_p(tt,params, &
                           Y_R,Y_PHI,Y_Z,V_PLL,V_MU,m_cache,flagCon,flagCol, &
                           F,P,E_PHI,ne,PSIp)

                   end do


                   !$OMP SIMD
                   do cc=1_idef,pchunk
                      spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                      spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                      spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                      spp(ii)%vars%ne(pp-1+cc) = ne(cc)
                   end do
                   !$OMP END SIMD


                end if

                call analytical_fields_Bmag_p(pchunk,F,Y_R,Y_PHI,Y_Z, &
                     Bmag,E_PHI)

                !$OMP SIMD
                do cc=1_idef,pchunk
                   spp(ii)%vars%g(pp-1+cc)=sqrt(1+V_PLL(cc)**2+ &
                        2*V_MU(cc)*Bmag(cc)*m_cache)

                   spp(ii)%vars%eta(pp-1+cc) = rad2deg(atan2(sqrt(2*m_cache* &
                        Bmag(cc)*spp(ii)%vars%V(pp-1+cc,2)), &
                        spp(ii)%vars%V(pp-1+cc,1)))
                end do
                !$OMP END SIMD

             end do !particle chunk iterator
             !$OMP END PARALLEL DO

          else if (params%LargeCollisions.and.sample_test) then

             do tt=1_ip,params%coll_per_dump

                !$OMP PARALLEL DO default(none) &
                !$OMP& FIRSTPRIVATE(m_cache,pchunk) &
                !$OMP& shared(F,P,params,ii,spp,tt) &
                !$OMP& PRIVATE(pp,Bmag,cc,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
                !$OMP& flagCon,flagCol,B_R,B_PHI,B_Z,E_PHI,PSIp,ne, &
                !$OMP& achunk,Te)
                do pp=1_idef,spp(ii)%pinit,pchunk

                   if ((spp(ii)%pRE-pp).lt.pchunk) then
                      achunk=spp(ii)%pRE-pp+1
                   else
                      achunk=pchunk
                   end if

                   !$OMP SIMD
                   do cc=1_idef,achunk
                      Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)
                      Y_PHI(cc)=spp(ii)%vars%Y(pp-1+cc,2)
                      Y_Z(cc)=spp(ii)%vars%Y(pp-1+cc,3)

                      V_PLL(cc)=spp(ii)%vars%V(pp-1+cc,1)
                      V_MU(cc)=spp(ii)%vars%V(pp-1+cc,2)

                      PSIp(cc)=spp(ii)%vars%PSI_p(pp-1+cc)
                      ne(cc)=spp(ii)%vars%ne(pp-1+cc)

                      flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
                      flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)
                   end do
                   !$OMP END SIMD

                   if (params%FokPlan) then

                      !if (mod(tt,50).eq.0) then
                      !write(6,*) 'particle',pp,'iteration',tt, &
                      !     ' of',params%t_skip
                      !flush(6)
                      !endif
                      
                      call include_CoulombCollisionsLA_GC_p(spp(ii),achunk, &
                           tt,params,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,m_cache, &
                           flagCon,flagCol,F,P,E_PHI,ne,Te,PSIp)

                   end if

                   !$OMP SIMD
                   do cc=1_idef,achunk
                      spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                      spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                      spp(ii)%vars%ne(pp-1+cc)=ne(cc)

                      spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)
                   end do
                   !$OMP END SIMD
                   
                end do !particle chunk iterator
                !$OMP END PARALLEL DO

             end do

             !$OMP PARALLEL DO default(none) &
             !$OMP& FIRSTPRIVATE(m_cache,pchunk) &
             !$OMP& shared(F,ii,spp) &
             !$OMP& PRIVATE(pp,Bmag,cc,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,achunk,E_PHI)
             do pp=1_idef,spp(ii)%pRE,pchunk

                if ((spp(ii)%pRE-pp).lt.pchunk) then
                   achunk=spp(ii)%pRE-pp+1
                else
                   achunk=pchunk
                end if

                !$OMP SIMD
                do cc=1_idef,achunk
                   Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)
                   Y_PHI(cc)=spp(ii)%vars%Y(pp-1+cc,2)
                   Y_Z(cc)=spp(ii)%vars%Y(pp-1+cc,3)
                end do
                !$OMP END SIMD

                !write(6,*) 'Y_R',Y_R
                
                call analytical_fields_Bmag_p(achunk,F,Y_R,Y_PHI,Y_Z, &
                     Bmag,E_PHI)                

                !$OMP SIMD
                do cc=1_idef,achunk
                   V_PLL(cc)=spp(ii)%vars%V(pp-1+cc,1)
                   V_MU(cc)=spp(ii)%vars%V(pp-1+cc,2)
                   
                   spp(ii)%vars%g(pp-1+cc)=sqrt(1+V_PLL(cc)**2+ &
                        2*V_MU(cc)*Bmag(cc)*m_cache)

                   spp(ii)%vars%eta(pp-1+cc) = rad2deg(atan2(sqrt(2*m_cache* &
                        Bmag(cc)*spp(ii)%vars%V(pp-1+cc,2)), &
                        spp(ii)%vars%V(pp-1+cc,1)))
                end do
                !$OMP END SIMD

             end do !particle chunk iterator
             !$OMP END PARALLEL DO

          else
             
             do tt=1_ip,params%coll_per_dump

                !$OMP PARALLEL DO default(none) &
                !$OMP& FIRSTPRIVATE(m_cache,pchunk,q_cache) &
                !$OMP& shared(F,P,params,ii,spp,tt) &
                !$OMP& PRIVATE(pp,Bmag,cc,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
                !$OMP& flagCon,flagCol,B_R,B_PHI,B_Z,E_PHI,PSIp,ne, &
                !$OMP& achunk,tttt,Te)
                do pp=1_idef,spp(ii)%pRE,pchunk

                   if ((spp(ii)%pRE-pp).lt.pchunk) then
                      achunk=spp(ii)%pRE-pp+1
                   else
                      achunk=pchunk
                   end if

                   !$OMP SIMD
                   do cc=1_idef,achunk
                      Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)
                      Y_PHI(cc)=spp(ii)%vars%Y(pp-1+cc,2)
                      Y_Z(cc)=spp(ii)%vars%Y(pp-1+cc,3)

                      V_PLL(cc)=spp(ii)%vars%V(pp-1+cc,1)
                      V_MU(cc)=spp(ii)%vars%V(pp-1+cc,2)

                      PSIp(cc)=spp(ii)%vars%PSI_p(pp-1+cc)
                      ne(cc)=spp(ii)%vars%ne(pp-1+cc)

                      flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
                      flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)
                   end do
                   !$OMP END SIMD
                   
                      !if (mod(tt,50).eq.0) then
                      !write(6,*) 'particle',pp,'iteration',tt, &
                      !     ' of',params%t_skip
                      !flush(6)
                      !endif

                   if (.not.params%FokPlan) then
                      do tttt=1_ip,params%orbits_per_coll
                         call advance_GCeqn_vars(spp(ii)%vars,pp, &
                              tttt,params, &
                              Y_R,Y_PHI,Y_Z,V_PLL,V_MU,flagCon,flagCol,q_cache,m_cache, &
                              B_R,B_PHI,B_Z,F,P,PSIp,E_PHI)
                      end do
                   endif
                   
                   call include_CoulombCollisionsLA_GC_p(spp(ii),achunk, &
                        tt,params,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,m_cache, &
                        flagCon,flagCol,F,P,E_PHI,ne,Te,PSIp)


                   !$OMP SIMD
                   do cc=1_idef,achunk
                      spp(ii)%vars%ne(pp-1+cc)=ne(cc)

                      spp(ii)%vars%Y(pp-1+cc,1)=Y_R(cc)
                      spp(ii)%vars%Y(pp-1+cc,2)=Y_PHI(cc)
                      spp(ii)%vars%Y(pp-1+cc,3)=Y_Z(cc)

                      spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                      spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                      spp(ii)%vars%flagCon(pp-1+cc)=flagCon(cc)
                      spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                      spp(ii)%vars%B(pp-1+cc,1) = B_R(cc)
                      spp(ii)%vars%B(pp-1+cc,2) = B_PHI(cc)
                      spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                      spp(ii)%vars%PSI_P(pp-1+cc) = PSIp(cc)
                      spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
                   end do
                   !$OMP END SIMD
                   
                end do !particle chunk iterator
                !$OMP END PARALLEL DO

             end do

             !$OMP PARALLEL DO default(none) &
             !$OMP& FIRSTPRIVATE(m_cache,pchunk) &
             !$OMP& shared(F,ii,spp) &
             !$OMP& PRIVATE(pp,Bmag,cc,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,achunk,E_PHI)
             do pp=1_idef,spp(ii)%pRE,pchunk

                if ((spp(ii)%pRE-pp).lt.pchunk) then
                   achunk=spp(ii)%pRE-pp+1
                else
                   achunk=pchunk
                end if

                !$OMP SIMD
                do cc=1_idef,achunk
                   Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)
                   Y_PHI(cc)=spp(ii)%vars%Y(pp-1+cc,2)
                   Y_Z(cc)=spp(ii)%vars%Y(pp-1+cc,3)
                end do
                !$OMP END SIMD

                !write(6,*) 'Y_R',Y_R
                
                call analytical_fields_Bmag_p(achunk,F,Y_R,Y_PHI,Y_Z, &
                     Bmag,E_PHI)                

                !$OMP SIMD
                do cc=1_idef,achunk
                   V_PLL(cc)=spp(ii)%vars%V(pp-1+cc,1)
                   V_MU(cc)=spp(ii)%vars%V(pp-1+cc,2)
                   
                   spp(ii)%vars%g(pp-1+cc)=sqrt(1+V_PLL(cc)**2+ &
                        2*V_MU(cc)*Bmag(cc)*m_cache)

                   spp(ii)%vars%eta(pp-1+cc) = rad2deg(atan2(sqrt(2*m_cache* &
                        Bmag(cc)*spp(ii)%vars%V(pp-1+cc,2)), &
                        spp(ii)%vars%V(pp-1+cc,1)))
                end do
                !$OMP END SIMD

             end do !particle chunk iterator
             !$OMP END PARALLEL DO
             
          endif

          

       end do

       if (params%SC_E) then
          call calculate_SC_E1D(params,F,VdenOMP)             
       end if



    end do !species iterator

  end subroutine adv_GCeqn_top

  subroutine advance_GCeqn_vars(vars,pp,tt,params,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
       flagCon,flagCol,q_cache,m_cache,B_R,B_PHI,B_Z,F,P,PSIp,E_PHI)
    !! @note Subroutine to advance GC variables \(({\bf X},p_\parallel)\)
    !! @endnote
    !! Comment this section further with evolution equations, numerical
    !! methods, and descriptions of both.
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(PARTICLES), INTENT(INOUT)     :: vars
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                 :: F
    !! An instance of the KORC derived type PROFILES.
    REAL(rp)                                      :: dt
    !! Time step used in the leapfrog step (\(\Delta t\)).

    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    INTEGER(ip),INTENT(IN)                     :: tt
    !! time iterator.
    INTEGER,INTENT(IN)                                      :: pp

    REAL(rp) :: a1 = 1./5._rp
    REAL(rp) :: a21 = 3./40._rp,a22=9./40._rp
    REAL(rp) :: a31 = 3./10._rp,a32=-9./10._rp,a33=6./5._rp
    REAL(rp) :: a41 = -11./54._rp,a42=5./2._rp,a43=-70./27._rp,a44=35./27._rp
    REAL(rp) :: a51 = 1631./55296._rp,a52=175./512._rp,a53=575./13824._rp,a54=44275./110592._rp,a55=253./4096._rp
    REAL(rp) :: b1=37./378._rp,b2=0._rp,b3=250./621._rp,b4=125./594._rp,b5=0._rp,b6=512./1771._rp

    REAL(rp),DIMENSION(params%pchunk) :: k1_R,k1_PHI,k1_Z,k1_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k2_R,k2_PHI,k2_Z,k2_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k3_R,k3_PHI,k3_Z,k3_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k4_R,k4_PHI,k4_Z,k4_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k5_R,k5_PHI,k5_Z,k5_PLL
    REAL(rp),DIMENSION(params%pchunk) :: k6_R,k6_PHI,k6_Z,k6_PLL
    REAL(rp),DIMENSION(params%pchunk) :: Y0_R,Y0_PHI,Y0_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: PSIp
    REAL(rp),DIMENSION(params%pchunk) :: curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: V_PLL,V_MU
    REAL(rp),DIMENSION(params%pchunk) :: RHS_R,RHS_PHI,RHS_Z,RHS_PLL,V0,E_Z,E_R
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: E_PHI
    REAL(rp),DIMENSION(params%pchunk) :: Bmag,ne,Te,Zeff
    INTEGER(is),dimension(params%pchunk), intent(inout) :: flagCon,flagCol

    REAL(rp) :: ar,R0
    REAL(rp),intent(IN) :: q_cache,m_cache

    ar=F%AB%a
    R0=F%AB%Ro

    pchunk=params%pchunk
    dt=params%dt

    !    write(output_unit_write,'("Y_R 0: ",E17.10)') Y_R(1)
    !    write(output_unit_write,'("Y_PHI 0: ",E17.10)') Y_PHI(1)
    !    write(output_unit_write,'("Y_Z 0: ",E17.10)') Y_Z(1)

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL)
    do cc=1_idef,pchunk
       Y0_R(cc)=Y_R(cc)
       Y0_PHI(cc)=Y_PHI(cc)
       Y0_Z(cc)=Y_Z(cc)
       V0(cc)=V_PLL(cc)
    end do
    !$OMP END SIMD

    call analytical_fields_GC_p(pchunk,F,Y_R,Y_PHI, &
         Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
         gradB_R,gradB_PHI,gradB_Z,PSIp)

    if (params%SC_E_add) then
#ifdef PSPLINE
       call add_interp_SCE_p(params,F,Y_R,Y_PHI,Y_Z,E_PHI)
#endif
    end if

    !    write(output_unit_write,'("ER:",E17.10)') E_R
    !    write(output_unit_write,'("EPHI:",E17.10)') E_PHI
    !    write(output_unit_write,'("EZ:",E17.10)') E_Z


    call GCEoM_p(params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,q_cache,m_cache) 

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,k1_R,k1_PHI,k1_Z,k1_PLL)
    do cc=1_idef,pchunk
       k1_R(cc)=dt*RHS_R(cc)              
       k1_PHI(cc)=dt*RHS_PHI(cc)    
       k1_Z(cc)=dt*RHS_Z(cc)    
       k1_PLL(cc)=dt*RHS_PLL(cc)    

       !       vars%RHS(pp-1+cc,1)=RHS_R(cc)
       !       vars%RHS(pp-1+cc,2)=RHS_PHI(cc)
       !       vars%RHS(pp-1+cc,3)=RHS_Z(cc)
       !       vars%RHS(pp-1+cc,4)=RHS_PLL(cc)

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
       V_PLL(cc)=V0(cc)   +a1*k1_PLL(cc)
    end do
    !$OMP END SIMD


    !    write(output_unit_write,'("Y_R 1: ",E17.10)') Y_R(1)
    !    write(output_unit_write,'("Y_PHI 1: ",E17.10)') Y_PHI(1)
    !    write(output_unit_write,'("Y_Z 1: ",E17.10)') Y_Z(1)

    !    write(output_unit_write,'("k1R: ",E17.10)') k1_R(1)
    !    write(output_unit_write,'("k1PHI: ",E17.10)') k1_PHI(1)
    !    write(output_unit_write,'("k1Z: ",E17.10)') k1_Z(1)
    !    write(output_unit_write,'("k1PLL: ",E17.10)') k1_PLL(1) 

    call analytical_fields_GC_p(pchunk,F,Y_R,Y_PHI, &
         Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
         gradB_R,gradB_PHI,gradB_Z,PSIp)

    if (params%SC_E_add) then
#ifdef PSPLINE
       call add_interp_SCE_p(params,F,Y_R,Y_PHI,Y_Z,E_PHI)
#endif
    end if

    call GCEoM_p(params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,q_cache,m_cache) 

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,k2_R,k2_PHI,k2_Z,k2_PLL)
    do cc=1_idef,pchunk
       k2_R(cc)=dt*RHS_R(cc)    
       k2_PHI(cc)=dt*RHS_PHI (cc)   
       k2_Z(cc)=dt*RHS_Z(cc)   
       k2_PLL(cc)=dt*RHS_PLL(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
       V_PLL(cc)=V0(cc)   +a21*k1_PLL(cc)+a22*k2_PLL(cc)
    end do
    !$OMP END SIMD


    !    write(output_unit_write,'("Y_R 2: ",E17.10)') Y_R(1)
    !    write(output_unit_write,'("Y_PHI 2: ",E17.10)') Y_PHI(1)
    !    write(output_unit_write,'("Y_Z 2: ",E17.10)') Y_Z(1)

    call analytical_fields_GC_p(pchunk,F,Y_R,Y_PHI, &
         Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
         gradB_R,gradB_PHI,gradB_Z,PSIp)

    if (params%SC_E_add) then
#ifdef PSPLINE
       call add_interp_SCE_p(params,F,Y_R,Y_PHI,Y_Z,E_PHI)
#endif
    end if

    call GCEoM_p(params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,q_cache,m_cache)

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,k3_R,k3_PHI,k3_Z,k3_PLL)
    do cc=1_idef,pchunk
       k3_R(cc)=dt*RHS_R(cc)   
       k3_PHI(cc)=dt*RHS_PHI(cc)    
       k3_Z(cc)=dt*RHS_Z(cc)    
       k3_PLL(cc)=dt*RHS_PLL(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)
       V_PLL(cc)=V0(cc)   +a31*k1_PLL(cc)+a32*k2_PLL(cc)+a33*k3_PLL(cc)
    end do
    !$OMP END SIMD

    !    write(output_unit_write,'("Y_R 3: ",E17.10)') Y_R(1)
    !    write(output_unit_write,'("Y_PHI 3: ",E17.10)') Y_PHI(1)
    !    write(output_unit_write,'("Y_Z 3: ",E17.10)') Y_Z(1)

    call analytical_fields_GC_p(pchunk,F,Y_R,Y_PHI, &
         Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
         gradB_R,gradB_PHI,gradB_Z,PSIp)

    if (params%SC_E_add) then
#ifdef PSPLINE
       call add_interp_SCE_p(params,F,Y_R,Y_PHI,Y_Z,E_PHI)
#endif
    end if

    call GCEoM_p(params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,q_cache,m_cache)     

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,k4_R,k4_PHI,k4_Z,k4_PLL)
    do cc=1_idef,pchunk
       k4_R(cc)=dt*RHS_R(cc)   
       k4_PHI(cc)=dt*RHS_PHI(cc)    
       k4_Z(cc)=dt*RHS_Z(cc)    
       k4_PLL(cc)=dt*RHS_PLL(cc)

       Y_R(cc)=Y0_R(cc)+a41*k1_R(cc)+a42*k2_R(cc)+a43*k3_R(cc)+ &
            a44*k4_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a41*k1_PHI(cc)+a42*k2_PHI(cc)+ &
            a43*k3_PHI(cc)+a44*k4_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a41*k1_Z(cc)+a42*k2_Z(cc)+a43*k3_Z(cc)+ &
            a44*k4_Z(cc)
       V_PLL(cc)=V0(cc)   +a41*k1_PLL(cc)+a42*k2_PLL(cc)+ &
            a43*k3_PLL(cc)+a44*k4_PLL(cc)
    end do
    !$OMP END SIMD

    !    write(output_unit_write,'("Y_R 4: ",E17.10)') Y_R(1)
    !    write(output_unit_write,'("Y_PHI 4: ",E17.10)') Y_PHI(1)
    !    write(output_unit_write,'("Y_Z 4: ",E17.10)') Y_Z(1)

    call analytical_fields_GC_p(pchunk,F,Y_R,Y_PHI, &
         Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
         gradB_R,gradB_PHI,gradB_Z,PSIp)

    if (params%SC_E_add) then
#ifdef PSPLINE
       call add_interp_SCE_p(params,F,Y_R,Y_PHI,Y_Z,E_PHI)
#endif
    end if

    call GCEoM_p(params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,q_cache,m_cache)   

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,k5_R,k5_PHI,k5_Z,k5_PLL)
    do cc=1_idef,pchunk
       k5_R(cc)=dt*RHS_R(cc)    
       k5_PHI(cc)=dt*RHS_PHI(cc)    
       k5_Z(cc)=dt*RHS_Z(cc)    
       k5_PLL(cc)=dt*RHS_PLL(cc)

       Y_R(cc)=Y0_R(cc)+a51*k1_R(cc)+a52*k2_R(cc)+a53*k3_R(cc)+ &
            a54*k4_R(cc)+a55*k5_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a51*k1_PHI(cc)+a52*k2_PHI(cc)+ &
            a53*k3_PHI(cc)+a54*k4_PHI(cc)+a55*k5_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a51*k1_Z(cc)+a52*k2_Z(cc)+a53*k3_Z(cc)+ &
            a54*k4_Z(cc)+a55*k5_Z(cc)
       V_PLL(cc)=V0(cc)   +a51*k1_PLL(cc)+a52*k2_PLL(cc)+ &
            a53*k3_PLL(cc)+a54*k4_PLL(cc)+a55*k5_PLL(cc)
    end do
    !$OMP END SIMD

    !    write(output_unit_write,'("Y_R 5: ",E17.10)') Y_R(1)
    !    write(output_unit_write,'("Y_PHI 5: ",E17.10)') Y_PHI(1)
    !    write(output_unit_write,'("Y_Z 5: ",E17.10)') Y_Z(1)

    call analytical_fields_GC_p(pchunk,F,Y_R,Y_PHI, &
         Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
         gradB_R,gradB_PHI,gradB_Z,PSIp)

    if (params%SC_E_add) then
#ifdef PSPLINE
       call add_interp_SCE_p(params,F,Y_R,Y_PHI,Y_Z,E_PHI)
#endif
    end if

    call GCEoM_p(params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,q_cache,m_cache)         

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0,Y_R,Y_PHI,Y_Z,V_PLL, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,k6_R,k6_PHI,k6_Z,k6_PLL)
    do cc=1_idef,pchunk
       k6_R(cc)=dt*RHS_R(cc)    
       k6_PHI(cc)=dt*RHS_PHI(cc)    
       k6_Z(cc)=dt*RHS_Z(cc)    
       k6_PLL(cc)=dt*RHS_PLL(cc)

       Y_R(cc)=Y0_R(cc)+b1*k1_R(cc)+b2*k2_R(cc)+ &
            b3*k3_R(cc)+b4*k4_R(cc)+b5*k5_R(cc)+b6*k6_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+b1*k1_PHI(cc)+b2*k2_PHI(cc)+ &
            b3*k3_PHI(cc)+b4*k4_PHI(cc)+b5*k5_PHI(cc)+b6*k6_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+b1*k1_Z(cc)+b2*k2_Z(cc)+ &
            b3*k3_Z(cc)+b4*k4_Z(cc)+b5*k5_Z(cc)+b6*k6_Z(cc)
       V_PLL(cc)=V0(cc)+b1*k1_PLL(cc)+b2*k2_PLL(cc)+ &
            b3*k3_PLL(cc)+b4*k4_PLL(cc)+b5*k5_PLL(cc)+b6*k6_PLL(cc) 
    end do
    !$OMP END SIMD

    !    write(output_unit_write,'("Y_R 6: ",E17.10)') Y_R(1)
    !    write(output_unit_write,'("Y_PHI 6: ",E17.10)') Y_PHI(1)
    !    write(output_unit_write,'("Y_Z 6: ",E17.10)') Y_Z(1)

    call cyl_check_if_confined_p(pchunk,ar,R0,Y_R,Y_Z,flagCon)

    !$OMP SIMD
    !    !$OMP& aligned(Y_R,Y_PHI,Y_Z,V_PLL,Y0_R,Y0_PHI,Y0_Z,V0)
    do cc=1_idef,pchunk

       if ((flagCon(cc).eq.0_is).or.(flagCol(cc).eq.0_is)) then
          Y_R(cc)=Y0_R(cc)
          Y_PHI(cc)=Y0_PHI(cc)
          Y_Z(cc)=Y0_Z(cc)
          V_PLL(cc)=V0(cc)
       end if

    end do
    !$OMP END SIMD

    call analytical_fields_GC_p(pchunk,F,Y_R,Y_PHI, &
         Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
         gradB_R,gradB_PHI,gradB_Z,PSIp)

    if (params%SC_E_add) then
#ifdef PSPLINE
       call add_interp_SCE_p(params,F,Y_R,Y_PHI,Y_Z,E_PHI)
#endif
    end if

  end subroutine advance_GCeqn_vars

  subroutine advance_FPeqn_vars(params,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,flagCon,flagCol, &
       m_cache,F,P,PSIp)

    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    INTEGER(ip)                                                    :: tt
    !! time iterator.
    REAL(rp),DIMENSION(params%pchunk), INTENT(INOUT)  :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk), INTENT(INOUT)  :: V_PLL,V_MU,PSIp
    REAL(rp),DIMENSION(params%pchunk)  :: E_PHI
    INTEGER(is),DIMENSION(params%pchunk), INTENT(INOUT)  :: flagCon,flagCol
    REAL(rp),intent(in) :: m_cache
    REAL(rp),DIMENSION(params%pchunk) :: ne

    do tt=1_ip,params%t_skip

       call include_CoulombCollisions_GC_p(tt,params,Y_R,Y_PHI,Y_Z, &
            V_PLL,V_MU,m_cache,flagCon,flagCol,F,P,E_PHI,ne,PSIp)

       !       write(output_unit_write,'("Collision Loop in FP")')

    end do




  end subroutine advance_FPeqn_vars

#ifdef PSPLINE
  subroutine adv_GCinterp_psi_top_FS(params,spp,P,F)

    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(INOUT)                                   :: F
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(params%pchunk)               :: Bmag
    REAL(rp),DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_R,E_PHI,E_Z
    REAL(rp),DIMENSION(params%pchunk) :: ne,Te,Zeff    
    REAL(rp),DIMENSION(params%pchunk) :: V_PLL,V_MU
    REAL(rp),DIMENSION(params%pchunk) :: PSIp
    REAL(rp),DIMENSION(params%pchunk) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp),DIMENSION(params%pchunk) :: gradB_R,gradB_PHI,gradB_Z
    INTEGER(is),DIMENSION(params%pchunk) :: flagCon,flagCol
    REAL(rp) :: m_cache,q_cache,B0,EF0,R0,q0,lam,ar


    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    INTEGER(ip)                                                    :: ttt
    !! time iterator.

    real(rp),dimension(F%dim_1D) :: Vden,Vdenave,VdenOMP
    !! time iterator.


    do ii = 1_idef,params%num_species      

       pchunk=params%pchunk
       q_cache=spp(ii)%q
       m_cache=spp(ii)%m

       do ttt=1_ip,params%t_it_SC

          VdenOMP=0._rp

          !$OMP PARALLEL DO default(none) &
          !$OMP& FIRSTPRIVATE(q_cache,m_cache,pchunk) &
          !$OMP& SHARED(params,ii,spp,P,F) &
          !$OMP& PRIVATE(pp,tt,Bmag,cc,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,B_R,B_PHI,B_Z, &
          !$OMP& flagCon,flagCol,E_PHI,PSIp,curlb_R,curlb_PHI,curlb_Z, &
          !$OMP& gradB_R,gradB_PHI,gradB_Z,ne, &
          !$OMP& Vden,Vdenave) &
          !$OMP& REDUCTION(+:VdenOMP)
          do pp=1_idef,spp(ii)%ppp,pchunk

             !          write(output_unit_write,'("pp: ",I16)') pp

             !$OMP SIMD
             do cc=1_idef,pchunk
                Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)
                Y_PHI(cc)=spp(ii)%vars%Y(pp-1+cc,2)
                Y_Z(cc)=spp(ii)%vars%Y(pp-1+cc,3)

                V_PLL(cc)=spp(ii)%vars%V(pp-1+cc,1)
                V_MU(cc)=spp(ii)%vars%V(pp-1+cc,2)

                PSIp(cc)=spp(ii)%vars%PSI_P(pp-1+cc)

                flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
                flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)           
             end do
             !$OMP END SIMD

             if (.not.params%FokPlan) then
                Vdenave=0._rp
                do tt=1_ip,params%t_skip

                   call advance_GCinterp_psi_vars_FS(spp(ii)%vars,pp,tt, &
                        params, &
                        Y_R,Y_PHI,Y_Z,V_PLL,V_MU,q_cache,m_cache,flagCon,flagCol, &
                        F,P,B_R,B_PHI,B_Z,E_PHI,PSIp,curlb_R,curlb_PHI, &
                        curlb_Z,gradB_R,gradB_PHI,gradB_Z)

                   call calculate_SC_p_FS(params,F,B_R,B_PHI,B_Z,PSIp, &
                        V_PLL,V_MU,m_cache,flagCon,flagCol,Vden)

                   !                   write(output_unit_write,*) 'pre-Vdenave',Vdenave(F%dim_1D)
                   Vdenave=(Vdenave*REAL(tt-1_ip)+Vden)/REAL(tt)

                   !                   write(output_unit_write,*) 'Vden',Vden(F%dim_1D)
                   !                   write(output_unit_write,*) 'post-Vdenave',Vdenave(F%dim_1D)
                   !                   if (pp.eq.9_idef) write(output_unit_write,*) 'Vdenave',Vdenave(F%dim_1D)

                end do !timestep iterator

                !                write(output_unit_write,*) 'Vdenave',Vdenave(F%dim_1D)

                VdenOMP=VdenOMP+Vdenave

                !                write(output_unit_write,*) 'VdenOMP',VdenOMP(F%dim_1D)

                !$OMP SIMD
                do cc=1_idef,pchunk
                   spp(ii)%vars%Y(pp-1+cc,1)=Y_R(cc)
                   spp(ii)%vars%Y(pp-1+cc,2)=Y_PHI(cc)
                   spp(ii)%vars%Y(pp-1+cc,3)=Y_Z(cc)
                   spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                   spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                   spp(ii)%vars%flagCon(pp-1+cc)=flagCon(cc)
                   spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                   spp(ii)%vars%B(pp-1+cc,1) = B_R(cc)
                   spp(ii)%vars%B(pp-1+cc,2) = B_PHI(cc)
                   spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                   spp(ii)%vars%gradB(pp-1+cc,1) = gradB_R(cc)
                   spp(ii)%vars%gradB(pp-1+cc,2) = gradB_PHI(cc)
                   spp(ii)%vars%gradB(pp-1+cc,3) = gradB_Z(cc)

                   spp(ii)%vars%curlb(pp-1+cc,1) = curlb_R(cc)
                   spp(ii)%vars%curlb(pp-1+cc,2) = curlb_PHI(cc)
                   spp(ii)%vars%curlb(pp-1+cc,3) = curlb_Z(cc)

                   spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
                   spp(ii)%vars%PSI_P(pp-1+cc) = PSIp(cc)                
                end do
                !$OMP END SIMD

             else

                call advance_FPinterp_vars(params,Y_R,Y_PHI, &
                     Y_Z,V_PLL,V_MU,m_cache,flagCon,flagCol,F,P,E_PHI,ne,PSIp)

                !$OMP SIMD
                do cc=1_idef,pchunk
                   spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                   spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                   spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                   spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
                end do
                !$OMP END SIMD

             end if


             !$OMP SIMD
             do cc=1_idef,pchunk
                B_R(cc)=spp(ii)%vars%B(pp-1+cc,1)
                B_PHI(cc)=spp(ii)%vars%B(pp-1+cc,2)
                B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)

                Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

                spp(ii)%vars%g(pp-1+cc)=sqrt(1+V_PLL(cc)**2+ &
                     2*V_MU(cc)*Bmag(cc))

                spp(ii)%vars%eta(pp-1+cc) = atan2(sqrt(2*m_cache*Bmag(cc)* &
                     spp(ii)%vars%V(pp-1+cc,2)),spp(ii)%vars%V(pp-1+cc,1))* &
                     180.0_rp/C_PI
             end do
             !$OMP END SIMD

          end do !particle chunk iterator
          !$OMP END PARALLEL DO

          !write(output_unit_write,*) 'VdenOMP',VdenOMP(F%dim_1D)


          call calculate_SC_E1D_FS(params,F,VdenOMP)             


       end do

    end do !species iterator

  end subroutine adv_GCinterp_psi_top_FS

  subroutine adv_GCinterp_psi_top(params,spp,P,F)

    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(INOUT)                                   :: F
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(params%pchunk)               :: Bmag
    REAL(rp),DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_R,E_PHI,E_Z
    REAL(rp),DIMENSION(params%pchunk) :: ne,Te,Zeff    
    REAL(rp),DIMENSION(params%pchunk) :: V_PLL,V_MU
    REAL(rp),DIMENSION(params%pchunk) :: PSIp
    REAL(rp),DIMENSION(params%pchunk) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp),DIMENSION(params%pchunk) :: gradB_R,gradB_PHI,gradB_Z
    INTEGER(is),DIMENSION(params%pchunk) :: flagCon,flagCol
    REAL(rp) :: m_cache,q_cache,B0,EF0,R0,q0,lam,ar


    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                                                    :: cc,pchunk,achunk
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    INTEGER(ip)                                                    :: ttt
    !! time iterator.



    do ii = 1_idef,params%num_species      

       pchunk=params%pchunk
       q_cache=spp(ii)%q
       m_cache=spp(ii)%m

       if (.not.params%LargeCollisions) then

          !$OMP PARALLEL DO default(none) &
          !$OMP& FIRSTPRIVATE(q_cache,m_cache,pchunk) &
          !$OMP& SHARED(params,ii,spp,P,F) &
          !$OMP& PRIVATE(pp,tt,Bmag,cc,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,B_R,B_PHI,B_Z, &
          !$OMP& flagCon,flagCol,E_PHI,PSIp,curlb_R,curlb_PHI,curlb_Z, &
          !$OMP& gradB_R,gradB_PHI,gradB_Z,ne,E_R,E_Z)

          do pp=1_idef,spp(ii)%ppp,pchunk

             !          write(output_unit_write,'("pp: ",I16)') pp

             !$OMP SIMD
             do cc=1_idef,pchunk
                Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)
                Y_PHI(cc)=spp(ii)%vars%Y(pp-1+cc,2)
                Y_Z(cc)=spp(ii)%vars%Y(pp-1+cc,3)

                V_PLL(cc)=spp(ii)%vars%V(pp-1+cc,1)
                V_MU(cc)=spp(ii)%vars%V(pp-1+cc,2)

                PSIp(cc)=spp(ii)%vars%PSI_P(pp-1+cc)

                flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
                flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)           
             end do
             !$OMP END SIMD

             if (.not.params%FokPlan) then
                do tt=1_ip,params%t_skip
                   call advance_GCinterp_psi_vars(pchunk,spp(ii)%vars,pp,tt, &
                        params,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,q_cache,m_cache, &
                        flagCon,flagCol,F,P,B_R,B_PHI,B_Z,E_PHI,PSIp, &
                        curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,ne)

                   if (params%collisions) then       
                      call include_CoulombCollisions_GC_p(tt,params, &
                           Y_R,Y_PHI,Y_Z,V_PLL,V_MU,m_cache,flagCon,flagCol, &
                           F,P,E_PHI,ne,PSIp)                      
                   end if

                end do !timestep iterator


                !$OMP SIMD
                do cc=1_idef,pchunk
                   spp(ii)%vars%Y(pp-1+cc,1)=Y_R(cc)
                   spp(ii)%vars%Y(pp-1+cc,2)=Y_PHI(cc)
                   spp(ii)%vars%Y(pp-1+cc,3)=Y_Z(cc)
                   spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                   spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                   spp(ii)%vars%flagCon(pp-1+cc)=flagCon(cc)
                   spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                   spp(ii)%vars%B(pp-1+cc,1) = B_R(cc)
                   spp(ii)%vars%B(pp-1+cc,2) = B_PHI(cc)
                   spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                   spp(ii)%vars%gradB(pp-1+cc,1) = gradB_R(cc)
                   spp(ii)%vars%gradB(pp-1+cc,2) = gradB_PHI(cc)
                   spp(ii)%vars%gradB(pp-1+cc,3) = gradB_Z(cc)

                   spp(ii)%vars%curlb(pp-1+cc,1) = curlb_R(cc)
                   spp(ii)%vars%curlb(pp-1+cc,2) = curlb_PHI(cc)
                   spp(ii)%vars%curlb(pp-1+cc,3) = curlb_Z(cc)

                   spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
                   spp(ii)%vars%PSI_P(pp-1+cc) = PSIp(cc)

                   spp(ii)%vars%ne(pp-1+cc) = ne(cc)
                end do
                !$OMP END SIMD

             else if (params%FokPlan.and.params%collisions) then

                call advance_FPinterp_vars(params,Y_R,Y_PHI, &
                     Y_Z,V_PLL,V_MU,m_cache,flagCon,flagCol,F,P,E_PHI,ne,PSIp)

                !$OMP SIMD
                do cc=1_idef,pchunk
                   spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                   spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                   spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                   spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)

                   spp(ii)%vars%ne(pp-1+cc) = ne(cc)                
                end do
                !$OMP END SIMD

             else
                do tt=1_ip,params%t_skip
                   call calculate_GCfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
                        E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                        gradB_R,gradB_PHI,gradB_Z,flagCon,PSIp)
                end do

                !$OMP SIMD
                do cc=1_idef,pchunk
                   spp(ii)%vars%B(pp-1+cc,1) = B_R(cc)
                   spp(ii)%vars%B(pp-1+cc,2) = B_PHI(cc)
                   spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                   spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
                   spp(ii)%vars%PSI_P(pp-1+cc) = PSIp(cc)
                end do
                !$OMP END SIMD

             end if

             !$OMP SIMD
             do cc=1_idef,pchunk
                B_R(cc)=spp(ii)%vars%B(pp-1+cc,1)
                B_PHI(cc)=spp(ii)%vars%B(pp-1+cc,2)
                B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)

                Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+ &
                     B_Z(cc)*B_Z(cc))

                spp(ii)%vars%g(pp-1+cc)=sqrt(1+V_PLL(cc)**2+ &
                     2*V_MU(cc)*Bmag(cc))

                spp(ii)%vars%eta(pp-1+cc) = atan2(sqrt(2*m_cache*Bmag(cc)* &
                     spp(ii)%vars%V(pp-1+cc,2)),spp(ii)%vars%V(pp-1+cc,1))* &
                     180.0_rp/C_PI
             end do
             !$OMP END SIMD

          end do !particle chunk iterator
       !$OMP END PARALLEL DO
             
       else

          do tt=1_ip,params%coll_per_dump
             !$OMP PARALLEL DO default(none) &
             !$OMP& FIRSTPRIVATE(q_cache,m_cache,pchunk,achunk) &
             !$OMP& SHARED(params,ii,spp,P,F) &
             !$OMP& PRIVATE(pp,tt,Bmag,cc,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,B_R,B_PHI,B_Z, &
             !$OMP& flagCon,flagCol,E_PHI,PSIp,curlb_R,curlb_PHI,curlb_Z, &
             !$OMP& gradB_R,gradB_PHI,gradB_Z,ne,E_R,E_Z,Te)
             
             do pp=1_idef,spp(ii)%pRE,pchunk

                if ((spp(ii)%pRE-pp).lt.pchunk) then
                   achunk=spp(ii)%pRE-pp+1
                else
                   achunk=pchunk
                end if

                !$OMP SIMD
                do cc=1_idef,achunk
                   Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)
                   Y_PHI(cc)=spp(ii)%vars%Y(pp-1+cc,2)
                   Y_Z(cc)=spp(ii)%vars%Y(pp-1+cc,3)

                   V_PLL(cc)=spp(ii)%vars%V(pp-1+cc,1)
                   V_MU(cc)=spp(ii)%vars%V(pp-1+cc,2)

                   PSIp(cc)=spp(ii)%vars%PSI_P(pp-1+cc)

                   flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
                   flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)           
                end do
                !$OMP END SIMD

                !write(6,*) ''
                !write(6,*) 'pp',pp

                do ttt=1_ip,params%orbits_per_coll
                   call advance_GCinterp_psi_vars(achunk,spp(ii)%vars,pp,tt, &
                        params,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,q_cache,m_cache, &
                        flagCon,flagCol, &
                        F,P,B_R,B_PHI,B_Z,E_PHI,PSIp,curlb_R,curlb_PHI, &
                        curlb_Z,gradB_R,gradB_PHI,gradB_Z,ne)
                end do

                call include_CoulombCollisionsLA_GC_p(spp(ii),achunk, &
                     tt,params,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,m_cache, &
                     flagCon,flagCol,F,P,E_PHI,ne,Te,PSIp)

                !$OMP SIMD
                do cc=1_idef,achunk
                   spp(ii)%vars%Y(pp-1+cc,1)=Y_R(cc)
                   spp(ii)%vars%Y(pp-1+cc,2)=Y_PHI(cc)
                   spp(ii)%vars%Y(pp-1+cc,3)=Y_Z(cc)
                   spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                   spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                   spp(ii)%vars%flagCon(pp-1+cc)=flagCon(cc)
                   spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                   spp(ii)%vars%B(pp-1+cc,1) = B_R(cc)
                   spp(ii)%vars%B(pp-1+cc,2) = B_PHI(cc)
                   spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                   spp(ii)%vars%gradB(pp-1+cc,1) = gradB_R(cc)
                   spp(ii)%vars%gradB(pp-1+cc,2) = gradB_PHI(cc)
                   spp(ii)%vars%gradB(pp-1+cc,3) = gradB_Z(cc)

                   spp(ii)%vars%curlb(pp-1+cc,1) = curlb_R(cc)
                   spp(ii)%vars%curlb(pp-1+cc,2) = curlb_PHI(cc)
                   spp(ii)%vars%curlb(pp-1+cc,3) = curlb_Z(cc)

                   spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
                   spp(ii)%vars%PSI_P(pp-1+cc) = PSIp(cc)

                   spp(ii)%vars%ne(pp-1+cc) = ne(cc)
                   spp(ii)%vars%Te(pp-1+cc) = Te(cc)                
                end do
                !$OMP END SIMD

                !$OMP SIMD
                do cc=1_idef,achunk
                   B_R(cc)=spp(ii)%vars%B(pp-1+cc,1)
                   B_PHI(cc)=spp(ii)%vars%B(pp-1+cc,2)
                   B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)

                   Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+ &
                        B_Z(cc)*B_Z(cc))

                   spp(ii)%vars%g(pp-1+cc)=sqrt(1+V_PLL(cc)**2+ &
                        2*V_MU(cc)*Bmag(cc))

                   spp(ii)%vars%eta(pp-1+cc) = atan2(sqrt(2*m_cache*Bmag(cc)* &
                        spp(ii)%vars%V(pp-1+cc,2)),spp(ii)%vars%V(pp-1+cc,1))* &
                        180.0_rp/C_PI
                end do
                !$OMP END SIMD
                
             end do !particle chunk iterator
             !$OMP END PARALLEL DO

          end do                
          
       end if

    end do !species iterator

  end subroutine adv_GCinterp_psi_top
#endif

#ifdef FIO
  subroutine adv_GCinterp_fio_top(params,spp,P,F)

    USE omp_lib
    IMPLICIT NONE
    
    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(INOUT)                                   :: F
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(params%pchunk)               :: Bmag
    REAL(rp),DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_R,E_PHI,E_Z
    REAL(rp),DIMENSION(params%pchunk) :: ne,Te,Zeff,ni
    REAL(rp),DIMENSION(params%pchunk,params%num_impurity_species) :: nimp    
    REAL(rp),DIMENSION(params%pchunk) :: V_PLL,V_MU
    REAL(rp),DIMENSION(params%pchunk) :: PSIp
    REAL(rp),DIMENSION(params%pchunk) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp),DIMENSION(params%pchunk) :: gradB_R,gradB_PHI,gradB_Z
    INTEGER(is),DIMENSION(params%pchunk) :: flagCon,flagCol
    REAL(rp) :: m_cache,q_cache,B0,EF0,R0,q0,lam,ar
    TYPE(C_PTR), DIMENSION(params%pchunk)  :: hint


    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    INTEGER(ip)                                                    :: ttt
    !! time iterator.
    INTEGER             :: thread_num


    !write(6,*) '2Y_R',spp(1)%vars%Y(1:4,1)*params%cpp%length
    !write(6,*) '2(p/mc)',spp(1)%vars%V(1:4,1)

    do ii = 1_idef,params%num_species      

       pchunk=params%pchunk
       q_cache=spp(ii)%q
       m_cache=spp(ii)%m


       !$OMP PARALLEL DO default(none) &
       !$OMP& FIRSTPRIVATE(q_cache,m_cache,pchunk) &
       !$OMP& SHARED(params,ii,spp,P,F) &
       !$OMP& PRIVATE(pp,tt,Bmag,cc,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,B_R,B_PHI,B_Z, &
       !$OMP& flagCon,flagCol,E_PHI,PSIp,curlb_R,curlb_PHI,curlb_Z, &
       !$OMP& gradB_R,gradB_PHI,gradB_Z,ne,nimp,Te,Zeff,ni,E_R,E_Z,hint, &
       !$OMP& thread_num)

       do pp=1_idef,spp(ii)%ppp,pchunk

          thread_num = OMP_GET_THREAD_NUM()

          !write(6,*) thread_num,'3Y_R',spp(ii)%vars%Y(pp,1)*params%cpp%length
          !write(6,*) thread_num,'3(p/mc)',spp(ii)%vars%V(pp,1)
          
          !          write(output_unit_write,'("pp: ",I16)') pp

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)
             Y_PHI(cc)=spp(ii)%vars%Y(pp-1+cc,2)
             Y_Z(cc)=spp(ii)%vars%Y(pp-1+cc,3)

             B_R(cc)=spp(ii)%vars%B(pp-1+cc,1)
             B_PHI(cc)=spp(ii)%vars%B(pp-1+cc,2)
             B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)

             E_R(cc)=spp(ii)%vars%E(pp-1+cc,1)
             E_PHI(cc)=spp(ii)%vars%E(pp-1+cc,2)
             E_Z(cc)=spp(ii)%vars%E(pp-1+cc,3)

             gradB_R(cc)=spp(ii)%vars%gradB(pp-1+cc,1)
             gradB_PHI(cc)=spp(ii)%vars%gradB(pp-1+cc,2)
             gradB_Z(cc)=spp(ii)%vars%gradB(pp-1+cc,3)

             curlb_R(cc)=spp(ii)%vars%curlb(pp-1+cc,1)
             curlb_PHI(cc)=spp(ii)%vars%curlb(pp-1+cc,2)
             curlb_Z(cc)=spp(ii)%vars%curlb(pp-1+cc,3)

             V_PLL(cc)=spp(ii)%vars%V(pp-1+cc,1)
             V_MU(cc)=spp(ii)%vars%V(pp-1+cc,2)

             PSIp(cc)=spp(ii)%vars%PSI_P(pp-1+cc)

             flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
             flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)

             hint(cc)=spp(ii)%vars%hint(pp-1+cc)

             ne(cc)=spp(ii)%vars%ne(pp-1+cc)
             ni(cc)=spp(ii)%vars%ni(pp-1+cc)
             nimp(cc,:)=spp(ii)%vars%nimp(pp-1+cc,:)
             Te(cc)=spp(ii)%vars%Te(pp-1+cc)
          end do
          !$OMP END SIMD


          do tt=1_ip,params%t_skip

             !write(6,*) thread_num,'4Y_R',Y_R*params%cpp%length
             !write(6,*) thread_num,'4(p/mc)',V_PLL
             
             
             !if (mod(tt,params%t_skip/100).eq.0) then
             !   write(output_unit_write,*) 'iteration',tt
             !   flush(output_unit_write)
             !endif
             
             call advance_GCinterp_fio_vars(spp(ii)%vars,pp,tt, &
                  params,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,q_cache,m_cache, &
                  flagCon,flagCol, &
                  F,P,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,PSIp,curlb_R,curlb_PHI, &
                  curlb_Z,gradB_R,gradB_PHI,gradB_Z,ne,ni,Te,Zeff,nimp,hint)
          end do !timestep iterator


          !$OMP SIMD
          do cc=1_idef,pchunk
             spp(ii)%vars%Y(pp-1+cc,1)=Y_R(cc)
             spp(ii)%vars%Y(pp-1+cc,2)=Y_PHI(cc)
             spp(ii)%vars%Y(pp-1+cc,3)=Y_Z(cc)
             spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
             spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

             spp(ii)%vars%flagCon(pp-1+cc)=flagCon(cc)
             spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

             spp(ii)%vars%B(pp-1+cc,1) = B_R(cc)
             spp(ii)%vars%B(pp-1+cc,2) = B_PHI(cc)
             spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

             spp(ii)%vars%gradB(pp-1+cc,1) = gradB_R(cc)
             spp(ii)%vars%gradB(pp-1+cc,2) = gradB_PHI(cc)
             spp(ii)%vars%gradB(pp-1+cc,3) = gradB_Z(cc)

             spp(ii)%vars%curlb(pp-1+cc,1) = curlb_R(cc)
             spp(ii)%vars%curlb(pp-1+cc,2) = curlb_PHI(cc)
             spp(ii)%vars%curlb(pp-1+cc,3) = curlb_Z(cc)

             spp(ii)%vars%E(pp-1+cc,1) = E_R(cc)
             spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
             spp(ii)%vars%E(pp-1+cc,3) = E_Z(cc)
             spp(ii)%vars%PSI_P(pp-1+cc) = PSIp(cc)

             spp(ii)%vars%ne(pp-1+cc) = ne(cc)
             spp(ii)%vars%ni(pp-1+cc) = ni(cc)
             spp(ii)%vars%nimp(pp-1+cc,:) = nimp(cc,:)
             spp(ii)%vars%Te(pp-1+cc) = Te(cc)
             spp(ii)%vars%Zeff(pp-1+cc) = Zeff(cc)

             spp(ii)%vars%hint(pp-1+cc) = hint(cc)
          end do
          !$OMP END SIMD



          !$OMP SIMD
          do cc=1_idef,pchunk

             Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+ &
                  B_Z(cc)*B_Z(cc))

             spp(ii)%vars%g(pp-1+cc)=sqrt(1+V_PLL(cc)**2+ &
                  2*V_MU(cc)*Bmag(cc))

             spp(ii)%vars%eta(pp-1+cc) = atan2(sqrt(2*m_cache*Bmag(cc)* &
                  spp(ii)%vars%V(pp-1+cc,2)),spp(ii)%vars%V(pp-1+cc,1))* &
                  180.0_rp/C_PI
          end do
          !$OMP END SIMD

       end do !particle chunk iterator
       !$OMP END PARALLEL DO




    end do !species iterator

  end subroutine adv_GCinterp_fio_top
#endif


#ifdef PSPLINE
  subroutine adv_GCinterp_psiwE_top(params,spp,P,F)

    USE omp_lib
    IMPLICIT NONE
    
    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(INOUT)                                   :: F
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(params%pchunk)               :: Bmag
    REAL(rp),DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_R,E_PHI,E_Z
    REAL(rp),DIMENSION(params%pchunk) :: ne,Te,Zeff    
    REAL(rp),DIMENSION(params%pchunk) :: V_PLL,V_MU
    REAL(rp),DIMENSION(params%pchunk) :: PSIp
    REAL(rp),DIMENSION(params%pchunk) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp),DIMENSION(params%pchunk) :: gradB_R,gradB_PHI,gradB_Z
    INTEGER(is),DIMENSION(params%pchunk) :: flagCon,flagCol
    REAL(rp) :: m_cache,q_cache,B0,EF0,R0,q0,lam,ar


    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                 :: cc,pchunk,achunk
    !! Chunk iterator.
    INTEGER(ip)             :: tt,ttt
    !! time iterator.
    INTEGER             :: thread_num

    do ii = 1_idef,params%num_species      

       pchunk=params%pchunk
       q_cache=spp(ii)%q
       m_cache=spp(ii)%m

       if (.not.params%LargeCollisions) then

          !$OMP PARALLEL DO default(none) &
          !$OMP& FIRSTPRIVATE(q_cache,m_cache,pchunk) &
          !$OMP& SHARED(params,ii,spp,P,F) &
          !$OMP& PRIVATE(pp,tt,Bmag,cc,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,B_R,B_PHI,B_Z, &
          !$OMP& flagCon,flagCol,E_PHI,PSIp,curlb_R,curlb_PHI,curlb_Z, &
          !$OMP& gradB_R,gradB_PHI,gradB_Z,ne,E_R,E_Z,thread_num)

          do pp=1_idef,spp(ii)%ppp,pchunk

             thread_num = OMP_GET_THREAD_NUM()

             !          write(output_unit_write,'("pp: ",I16)') pp

             !$OMP SIMD
             do cc=1_idef,pchunk
                Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)
                Y_PHI(cc)=spp(ii)%vars%Y(pp-1+cc,2)
                Y_Z(cc)=spp(ii)%vars%Y(pp-1+cc,3)

                V_PLL(cc)=spp(ii)%vars%V(pp-1+cc,1)
                V_MU(cc)=spp(ii)%vars%V(pp-1+cc,2)

                PSIp(cc)=spp(ii)%vars%PSI_P(pp-1+cc)

                flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
                flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)           
             end do
             !$OMP END SIMD

             if (.not.params%FokPlan) then
                do tt=1_ip,params%t_skip

                   !if (params%t_skip.ge.10) then
                   !   if(mod(tt,params%t_skip/10).eq.0) then
                   !      if((params%mpi_params%rank.eq.0).and. &
                   !           thread_num.eq.0) then
                   !         write(6,'("tt iteration ",I8)') tt
                   !      endif
                   !   end if
                   !end if

                   call advance_GCinterp_psiwE_vars(spp(ii)%vars,pchunk,pp,tt, &
                        params, &
                        Y_R,Y_PHI,Y_Z,V_PLL,V_MU,q_cache,m_cache,flagCon,flagCol, &
                        F,P,B_R,B_PHI,B_Z,E_PHI,PSIp,curlb_R,curlb_PHI, &
                        curlb_Z,gradB_R,gradB_PHI,gradB_Z,ne)

                   if (params%collisions) then       

                      call include_CoulombCollisions_GC_p(tt,params, &
                           Y_R,Y_PHI,Y_Z, V_PLL,V_MU,m_cache, &
                           flagCon,flagCol,F,P,E_PHI,ne,PSIp)

                   end if

                   
                end do !timestep iterator


                !$OMP SIMD
                do cc=1_idef,pchunk
                   spp(ii)%vars%Y(pp-1+cc,1)=Y_R(cc)
                   spp(ii)%vars%Y(pp-1+cc,2)=Y_PHI(cc)
                   spp(ii)%vars%Y(pp-1+cc,3)=Y_Z(cc)
                   spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                   spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                   spp(ii)%vars%flagCon(pp-1+cc)=flagCon(cc)
                   spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                   spp(ii)%vars%B(pp-1+cc,1) = B_R(cc)
                   spp(ii)%vars%B(pp-1+cc,2) = B_PHI(cc)
                   spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                   spp(ii)%vars%gradB(pp-1+cc,1) = gradB_R(cc)
                   spp(ii)%vars%gradB(pp-1+cc,2) = gradB_PHI(cc)
                   spp(ii)%vars%gradB(pp-1+cc,3) = gradB_Z(cc)

                   spp(ii)%vars%curlb(pp-1+cc,1) = curlb_R(cc)
                   spp(ii)%vars%curlb(pp-1+cc,2) = curlb_PHI(cc)
                   spp(ii)%vars%curlb(pp-1+cc,3) = curlb_Z(cc)

                   spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
                   spp(ii)%vars%PSI_P(pp-1+cc) = PSIp(cc)

                   spp(ii)%vars%ne(pp-1+cc) = ne(cc)                
                end do
                !$OMP END SIMD

             else if (params%FokPlan.and.params%collisions) then

                do tt=1_ip,params%t_skip
                   call include_CoulombCollisions_GC_p(tt,params,Y_R,Y_PHI,Y_Z, &
                        V_PLL,V_MU,m_cache,flagCon,flagCol,F,P,E_PHI,ne,PSIp)

                end do

                !$OMP SIMD
                do cc=1_idef,pchunk
                   spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                   spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                   spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                   spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)

                   spp(ii)%vars%ne(pp-1+cc) = ne(cc)                
                end do
                !$OMP END SIMD

             else
                do tt=1_ip,params%t_skip
                   call calculate_GCfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,&
                        E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                        gradB_R,gradB_PHI,gradB_Z,flagCon,PSIp)
                end do

                !$OMP SIMD
                do cc=1_idef,pchunk
                   spp(ii)%vars%B(pp-1+cc,1) = B_R(cc)
                   spp(ii)%vars%B(pp-1+cc,2) = B_PHI(cc)
                   spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                   spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
                   spp(ii)%vars%PSI_P(pp-1+cc) = PSIp(cc)
                end do
                !$OMP END SIMD

             end if


             !$OMP SIMD
             do cc=1_idef,pchunk
                B_R(cc)=spp(ii)%vars%B(pp-1+cc,1)
                B_PHI(cc)=spp(ii)%vars%B(pp-1+cc,2)
                B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)

                Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+ &
                     B_Z(cc)*B_Z(cc))

                spp(ii)%vars%g(pp-1+cc)=sqrt(1+V_PLL(cc)**2+ &
                     2*V_MU(cc)*Bmag(cc))

                spp(ii)%vars%eta(pp-1+cc) = atan2(sqrt(2*m_cache*Bmag(cc)* &
                     spp(ii)%vars%V(pp-1+cc,2)),spp(ii)%vars%V(pp-1+cc,1))* &
                     180.0_rp/C_PI
             end do
             !$OMP END SIMD
             
          end do !particle chunk iterator
          !$OMP END PARALLEL DO

       else

          do tt=1_ip,params%coll_per_dump

#if DBG_CHECK    
    !         if(params%mpi_params%rank.eq.6) then
    !            write(6,*) 'before loop load:ppll',spp(ii)%vars%V(1:8,1),'mu',spp(ii)%vars%V(1:8,2)
    !            write(6,*) 'before loop load:R',spp(ii)%vars%Y(1:8,1),'PHI',spp(ii)%vars%Y(1:8,2),'Z',spp(ii)%vars%Y(1:8,3)
    !         end if
# endif
             
             !if (modulo(tt,params%coll_per_dump/10).eq.0) &
             !     write(6,*) 'mpi',params%mpi_params%rank,', Coll step',tt
             
             !$OMP PARALLEL DO default(none) &
             !$OMP& FIRSTPRIVATE(q_cache,m_cache,pchunk,tt) &
             !$OMP& SHARED(params,ii,spp,P,F) &
             !$OMP& PRIVATE(pp,ttt,Bmag,cc,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
             !$OMP& B_R,B_PHI,B_Z,achunk, &
             !$OMP& flagCon,flagCol,E_PHI,PSIp,curlb_R,curlb_PHI,curlb_Z, &
             !$OMP& gradB_R,gradB_PHI,gradB_Z,ne,Te,E_R,E_Z,thread_num)

             do pp=1_idef,spp(ii)%pRE,pchunk
                
                thread_num = OMP_GET_THREAD_NUM()

                if ((spp(ii)%pRE-pp).lt.pchunk) then
                   achunk=spp(ii)%pRE-pp+1
                else
                   achunk=pchunk
                end if

                
                !          write(output_unit_write,'("pp: ",I16)') pp
                
                !$OMP SIMD
                do cc=1_idef,achunk
                   Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)
                   Y_PHI(cc)=spp(ii)%vars%Y(pp-1+cc,2)
                   Y_Z(cc)=spp(ii)%vars%Y(pp-1+cc,3)

                   V_PLL(cc)=spp(ii)%vars%V(pp-1+cc,1)
                   V_MU(cc)=spp(ii)%vars%V(pp-1+cc,2)

                   PSIp(cc)=spp(ii)%vars%PSI_P(pp-1+cc)

                   flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
                   flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)           
                end do
                !$OMP END SIMD

#if DBG_CHECK    
                do cc=1_idef,achunk
                   if (V_MU(cc).le.0._rp) then
                      write(6,*) 'mu is negative before orbit loop'
                      write(6,*) 'coll_it',tt
                      write(6,*) 'V_PLL,V_MU',V_PLL(cc),V_MU(cc)
                      write(6,*) 'mpi,omp',params%mpi_params%rank,thread_num
                      write(6,*) 'pp,cc,part,pRE',pp,cc,pp-1+cc,spp%pRE
                      call korc_abort(25)
                   endif
                end do
#endif

                !if (params%t_skip.ge.10) then
                !   if(mod(tt,params%t_skip/10).eq.0) then
                !      if((params%mpi_params%rank.eq.0).and. &
                !           thread_num.eq.0) then
                !         write(6,'("tt iteration ",I8)') tt
                !      endif
                !   end if
                !end if

                do ttt=1_ip,params%orbits_per_coll

                   !if (modulo(ttt,params%orbits_per_coll/10).eq.0) &
                   !     write(6,*) 'mpi',params%mpi_params%rank,' , OMP',thread_num,', Orbit step',ttt

#if DBG_CHECK    
                   do cc=1_idef,achunk
                      if (V_MU(cc).le.0._rp) then
                         write(6,*) 'mu is negative in orbit loop'
                         write(6,*) 'coll_it,orb_it',tt,ttt
                         write(6,*) 'V_PLL,V_MU',V_PLL(cc),V_MU(cc)
                         write(6,*) 'mpi,omp',params%mpi_params%rank,thread_num
                         write(6,*) 'pp,cc,part,pRE',pp,cc,pp-1+cc,spp%pRE
                         call korc_abort(25)
                      endif
                   end do
# endif
                   
                   call advance_GCinterp_psiwE_vars(spp(ii)%vars,achunk, &
                        pp,tt,params,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
                        q_cache,m_cache,flagCon,flagCol, &
                        F,P,B_R,B_PHI,B_Z,E_PHI,PSIp,curlb_R,curlb_PHI, &
                        curlb_Z,gradB_R,gradB_PHI,gradB_Z,ne)
                   
                end do

#if DBG_CHECK    
                do cc=1_idef,achunk
                   if (V_MU(cc).le.0._rp) then
                      write(6,*) 'mu is negative after orbit loop'
                      write(6,*) 'coll_it,orb_it',tt,ttt
                      write(6,*) 'V_PLL,V_MU',V_PLL(cc),V_MU(cc)
                      write(6,*) 'mpi,omp',params%mpi_params%rank,thread_num
                      write(6,*) 'part,pRE',pp-1+cc,spp%pRE
                      call korc_abort(25)
                   endif
                end do
# endif
                
                call include_CoulombCollisionsLA_GC_p(spp(ii),achunk, &
                     tt,params,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,m_cache, &
                     flagCon,flagCol,F,P,E_PHI,ne,Te,PSIp)

#if DBG_CHECK    
                do cc=1_idef,achunk
                   if (V_MU(cc).le.0._rp) then
                      write(6,*) 'mu is negative after coll'
                      write(6,*) 'coll_it,orb_it',tt,ttt
                      write(6,*) 'V_PLL,V_MU',V_PLL(cc),V_MU(cc)
                      write(6,*) 'mpi,omp',params%mpi_params%rank,thread_num
                      write(6,*) 'part,pRE',pp-1+cc,spp%pRE
                      call korc_abort(25)
                   endif
                end do                              

                if((params%mpi_params%rank.eq.6).and.(pp.eq.1)) then
                   write(6,*) 'before loop save:ppll',V_PLL,'mu',V_MU
                end if
#endif
                
                !$OMP SIMD
                do cc=1_idef,achunk
                   spp(ii)%vars%Y(pp-1+cc,1)=Y_R(cc)
                   spp(ii)%vars%Y(pp-1+cc,2)=Y_PHI(cc)
                   spp(ii)%vars%Y(pp-1+cc,3)=Y_Z(cc)
                   spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                   spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                   spp(ii)%vars%flagCon(pp-1+cc)=flagCon(cc)
                   spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                   spp(ii)%vars%B(pp-1+cc,1) = B_R(cc)
                   spp(ii)%vars%B(pp-1+cc,2) = B_PHI(cc)
                   spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                   spp(ii)%vars%gradB(pp-1+cc,1) = gradB_R(cc)
                   spp(ii)%vars%gradB(pp-1+cc,2) = gradB_PHI(cc)
                   spp(ii)%vars%gradB(pp-1+cc,3) = gradB_Z(cc)

                   spp(ii)%vars%curlb(pp-1+cc,1) = curlb_R(cc)
                   spp(ii)%vars%curlb(pp-1+cc,2) = curlb_PHI(cc)
                   spp(ii)%vars%curlb(pp-1+cc,3) = curlb_Z(cc)

                   spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
                   spp(ii)%vars%PSI_P(pp-1+cc) = PSIp(cc)

                   spp(ii)%vars%ne(pp-1+cc) = ne(cc)                
                end do
                !$OMP END SIMD

#if DBG_CHECK    
          !      if(params%mpi_params%rank.eq.6.and.(pp.eq.1)) then
          !         write(6,*) 'after loop save 1:ppll',spp(ii)%vars%V(1:8,1),'mu',spp(ii)%vars%V(1:8,2)
          !         write(6,*) 'after loop save 1:R',spp(ii)%vars%Y(1:8,1),'PHI',spp(ii)%vars%Y(1:8,2),'Z',spp(ii)%vars%Y(1:8,3)
          !      end if
#endif
                
             end do !particle chunk iterator
             !$OMP END PARALLEL DO

#if DBG_CHECK    
          !   if(params%mpi_params%rank.eq.6) then
          !      write(6,*) 'after loop save 2:ppll',spp(ii)%vars%V(1:8,1),'mu',spp(ii)%vars%V(1:8,2)
          !      write(6,*) 'after loop save 2:R',spp(ii)%vars%Y(1:8,1),'PHI',spp(ii)%vars%Y(1:8,2),'Z',spp(ii)%vars%Y(1:8,3)
          !   end if
#endif
             
          end do !timestep iterator

          
             !$OMP PARALLEL DO default(none) &
             !$OMP& FIRSTPRIVATE(m_cache,pchunk) &
             !$OMP& SHARED(ii,spp) &
             !$OMP& PRIVATE(pp,Bmag,cc, &
             !$OMP& B_R,B_PHI,B_Z,achunk)

             do pp=1_idef,spp(ii)%pRE,pchunk

                if ((spp(ii)%pRE-pp).lt.pchunk) then
                   achunk=spp(ii)%pRE-pp+1
                else
                   achunk=pchunk
                end if
                
                !$OMP SIMD
                do cc=1_idef,achunk
                   B_R(cc)=spp(ii)%vars%B(pp-1+cc,1)
                   B_PHI(cc)=spp(ii)%vars%B(pp-1+cc,2)
                   B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)

                   Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+ &
                        B_Z(cc)*B_Z(cc))

                   spp(ii)%vars%g(pp-1+cc)=sqrt(1+spp(ii)%vars%V(pp-1+cc,1)**2 &
                        +2*spp(ii)%vars%V(pp-1+cc,2)*Bmag(cc))

                   spp(ii)%vars%eta(pp-1+cc) = atan2(sqrt(2*m_cache*Bmag(cc)* &
                        spp(ii)%vars%V(pp-1+cc,2)),spp(ii)%vars%V(pp-1+cc,1))* &
                        180.0_rp/C_PI
                end do
                !$OMP END SIMD

             end do !particle chunk iterator
             !$OMP END PARALLEL DO
          
       endif


    end do !species iterator

  end subroutine adv_GCinterp_psiwE_top

  subroutine adv_GCinterp_psi2x1t_top(params,spp,P,F)

    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(INOUT)                                   :: F
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(params%pchunk)               :: Bmag
    REAL(rp),DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_R,E_PHI,E_Z
    REAL(rp),DIMENSION(params%pchunk) :: ne,Te,Zeff    
    REAL(rp),DIMENSION(params%pchunk) :: V_PLL,V_MU
    REAL(rp),DIMENSION(params%pchunk) :: PSIp
    REAL(rp),DIMENSION(params%pchunk) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp),DIMENSION(params%pchunk) :: gradB_R,gradB_PHI,gradB_Z
    INTEGER(is),DIMENSION(params%pchunk) :: flagCon,flagCol
    REAL(rp) :: m_cache,q_cache,B0,EF0,R0,q0,lam,ar,time


    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    INTEGER(ip)                                                    :: ttt
    !! time iterator.



    do ii = 1_idef,params%num_species      

       pchunk=params%pchunk
       q_cache=spp(ii)%q
       m_cache=spp(ii)%m


       !$OMP PARALLEL DO default(none) &
       !$OMP& FIRSTPRIVATE(q_cache,m_cache,pchunk) &
       !$OMP& SHARED(params,ii,spp,P,F) &
       !$OMP& PRIVATE(pp,tt,Bmag,cc,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,B_R,B_PHI,B_Z, &
       !$OMP& flagCon,flagCol,E_PHI,PSIp,curlb_R,curlb_PHI,curlb_Z, &
       !$OMP& gradB_R,gradB_PHI,gradB_Z,ne,time,E_R,E_Z)

       do pp=1_idef,spp(ii)%ppp,pchunk

          !          write(output_unit_write,'("pp: ",I16)') pp

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)
             Y_PHI(cc)=spp(ii)%vars%Y(pp-1+cc,2)
             Y_Z(cc)=spp(ii)%vars%Y(pp-1+cc,3)

             V_PLL(cc)=spp(ii)%vars%V(pp-1+cc,1)
             V_MU(cc)=spp(ii)%vars%V(pp-1+cc,2)

             PSIp(cc)=spp(ii)%vars%PSI_P(pp-1+cc)

             flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
             flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)           
          end do
          !$OMP END SIMD

          if (.not.params%FokPlan) then
             do tt=1_ip,params%t_skip
                call advance_GCinterp_psi2x1t_vars(spp(ii)%vars,pp,tt, &
                     params, &
                     Y_R,Y_PHI,Y_Z,V_PLL,V_MU,q_cache,m_cache,flagCon,flagCol, &
                     F,P,B_R,B_PHI,B_Z,E_PHI,PSIp,curlb_R,curlb_PHI, &
                     curlb_Z,gradB_R,gradB_PHI,gradB_Z,ne)


             end do !timestep iterator


             !$OMP SIMD
             do cc=1_idef,pchunk
                spp(ii)%vars%Y(pp-1+cc,1)=Y_R(cc)
                spp(ii)%vars%Y(pp-1+cc,2)=Y_PHI(cc)
                spp(ii)%vars%Y(pp-1+cc,3)=Y_Z(cc)
                spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                spp(ii)%vars%flagCon(pp-1+cc)=flagCon(cc)
                spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                spp(ii)%vars%B(pp-1+cc,1) = B_R(cc)
                spp(ii)%vars%B(pp-1+cc,2) = B_PHI(cc)
                spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                spp(ii)%vars%gradB(pp-1+cc,1) = gradB_R(cc)
                spp(ii)%vars%gradB(pp-1+cc,2) = gradB_PHI(cc)
                spp(ii)%vars%gradB(pp-1+cc,3) = gradB_Z(cc)

                spp(ii)%vars%curlb(pp-1+cc,1) = curlb_R(cc)
                spp(ii)%vars%curlb(pp-1+cc,2) = curlb_PHI(cc)
                spp(ii)%vars%curlb(pp-1+cc,3) = curlb_Z(cc)

                spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
                spp(ii)%vars%PSI_P(pp-1+cc) = PSIp(cc)

                spp(ii)%vars%ne(pp-1+cc) = ne(cc)                
             end do
             !$OMP END SIMD

          else if (params%FokPlan.and.params%collisions) then

             call advance_FPinterp_vars(params,Y_R,Y_PHI, &
                  Y_Z,V_PLL,V_MU,m_cache,flagCon,flagCol,F,P,E_PHI,ne,PSIp)

             !$OMP SIMD
             do cc=1_idef,pchunk
                spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)

                spp(ii)%vars%ne(pp-1+cc) = ne(cc)                
             end do
             !$OMP END SIMD

          else
             do tt=1_ip,params%t_skip
                time=params%init_time+(params%it-1+tt)* &
                     params%dt

                call calculate_GCfields_2x1t_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
                     E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
                     gradB_R,gradB_PHI,gradB_Z,flagCon,PSIp,time)
             end do

             !$OMP SIMD
             do cc=1_idef,pchunk
                spp(ii)%vars%B(pp-1+cc,1) = B_R(cc)
                spp(ii)%vars%B(pp-1+cc,2) = B_PHI(cc)
                spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
                spp(ii)%vars%PSI_P(pp-1+cc) = PSIp(cc)
             end do
             !$OMP END SIMD

          end if


          !$OMP SIMD
          do cc=1_idef,pchunk
             B_R(cc)=spp(ii)%vars%B(pp-1+cc,1)
             B_PHI(cc)=spp(ii)%vars%B(pp-1+cc,2)
             B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)

             Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+ &
                  B_Z(cc)*B_Z(cc))

             spp(ii)%vars%g(pp-1+cc)=sqrt(1+V_PLL(cc)**2+ &
                  2*V_MU(cc)*Bmag(cc))

             spp(ii)%vars%eta(pp-1+cc) = atan2(sqrt(2*m_cache*Bmag(cc)* &
                  spp(ii)%vars%V(pp-1+cc,2)),spp(ii)%vars%V(pp-1+cc,1))* &
                  180.0_rp/C_PI
          end do
          !$OMP END SIMD

       end do !particle chunk iterator
       !$OMP END PARALLEL DO




    end do !species iterator

  end subroutine adv_GCinterp_psi2x1t_top

  subroutine adv_GCinterp_B_top(params,spp,P,F)

    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(params%pchunk)               :: Bmag
    REAL(rp),DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_R,E_PHI,E_Z
    REAL(rp),DIMENSION(params%pchunk) :: ne,Te,Zeff    
    REAL(rp),DIMENSION(params%pchunk) :: V_PLL,V_MU
    REAL(rp),DIMENSION(params%pchunk) :: PSIp
    REAL(rp),DIMENSION(params%pchunk) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp),DIMENSION(params%pchunk) :: gradB_R,gradB_PHI,gradB_Z
    INTEGER(is),DIMENSION(params%pchunk) :: flagCon,flagCol
    REAL(rp) :: m_cache,q_cache,B0,EF0,R0,q0,lam,ar


    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    !! time iterator.


    do ii = 1_idef,params%num_species      

       pchunk=params%pchunk
       q_cache=spp(ii)%q
       m_cache=spp(ii)%m

       !$OMP PARALLEL DO default(none) &
       !$OMP& FIRSTPRIVATE(q_cache,m_cache,pchunk) &
       !$OMP& SHARED(params,ii,spp,P,F) &
       !$OMP& PRIVATE(pp,tt,Bmag,cc,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,B_R,B_PHI,B_Z, &
       !$OMP& flagCon,flagCol,E_PHI,PSIp,curlb_R,curlb_PHI,curlb_Z, &
       !$OMP& gradB_R,gradB_PHI,gradB_Z,ne)
       do pp=1_idef,spp(ii)%ppp,pchunk

          !          write(output_unit_write,'("pp: ",I16)') pp

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)
             Y_PHI(cc)=spp(ii)%vars%Y(pp-1+cc,2)
             Y_Z(cc)=spp(ii)%vars%Y(pp-1+cc,3)

             V_PLL(cc)=spp(ii)%vars%V(pp-1+cc,1)
             V_MU(cc)=spp(ii)%vars%V(pp-1+cc,2)

             PSIp(cc)=spp(ii)%vars%PSI_P(pp-1+cc)

             flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
             flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)           
          end do
          !$OMP END SIMD

          if (.not.params%FokPlan) then       
             do tt=1_ip,params%t_skip
                call advance_GCinterp_B_vars(spp(ii)%vars,pp,tt,params, &
                     Y_R,Y_PHI,Y_Z,V_PLL,V_MU,q_cache,m_cache,flagCon,flagCol, &
                     F,P,B_R,B_PHI,B_Z,E_PHI,curlb_R,curlb_PHI,curlb_Z, &
                     gradB_R,gradB_PHI,gradB_Z,PSIp)
             end do !timestep iterator


             !$OMP SIMD
             do cc=1_idef,pchunk
                spp(ii)%vars%Y(pp-1+cc,1)=Y_R(cc)
                spp(ii)%vars%Y(pp-1+cc,2)=Y_PHI(cc)
                spp(ii)%vars%Y(pp-1+cc,3)=Y_Z(cc)
                spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                spp(ii)%vars%flagCon(pp-1+cc)=flagCon(cc)
                spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                spp(ii)%vars%B(pp-1+cc,1) = B_R(cc)
                spp(ii)%vars%B(pp-1+cc,2) = B_PHI(cc)
                spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                spp(ii)%vars%gradB(pp-1+cc,1) = gradB_R(cc)
                spp(ii)%vars%gradB(pp-1+cc,2) = gradB_PHI(cc)
                spp(ii)%vars%gradB(pp-1+cc,3) = gradB_Z(cc)

                spp(ii)%vars%curlb(pp-1+cc,1) = curlb_R(cc)
                spp(ii)%vars%curlb(pp-1+cc,2) = curlb_PHI(cc)
                spp(ii)%vars%curlb(pp-1+cc,3) = curlb_Z(cc)

                spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
             end do
             !$OMP END SIMD

          else

             call advance_FPinterp_vars(params,Y_R,Y_PHI, &
                  Y_Z,V_PLL,V_MU,m_cache,flagCon,flagCol,F,P,E_PHI,ne,PSIp)

             !$OMP SIMD
             do cc=1_idef,pchunk
                spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
             end do
             !$OMP END SIMD

          end if


          !$OMP SIMD
          do cc=1_idef,pchunk
             B_R(cc)=spp(ii)%vars%B(pp-1+cc,1)
             B_PHI(cc)=spp(ii)%vars%B(pp-1+cc,2)
             B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)

             Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

             spp(ii)%vars%g(pp-1+cc)=sqrt(1+V_PLL(cc)**2+ &
                  2*V_MU(cc)*Bmag(cc))

             spp(ii)%vars%eta(pp-1+cc) = atan2(sqrt(2*m_cache*Bmag(cc)* &
                  spp(ii)%vars%V(pp-1+cc,2)),spp(ii)%vars%V(pp-1+cc,1))* &
                  180.0_rp/C_PI
          end do
          !$OMP END SIMD

       end do !particle chunk iterator
       !$OMP END PARALLEL DO

    end do !species iterator

  end subroutine adv_GCinterp_B_top

  subroutine adv_GCinterp_B2D_top(params,spp,P,F)

    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(params%pchunk)               :: Bmag
    REAL(rp),DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_R,E_PHI,E_Z
    REAL(rp),DIMENSION(params%pchunk) :: ne,Te,Zeff    
    REAL(rp),DIMENSION(params%pchunk) :: V_PLL,V_MU
    REAL(rp),DIMENSION(params%pchunk) :: PSIp
    REAL(rp),DIMENSION(params%pchunk) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp),DIMENSION(params%pchunk) :: gradB_R,gradB_PHI,gradB_Z
    INTEGER(is),DIMENSION(params%pchunk) :: flagCon,flagCol
    REAL(rp) :: m_cache,q_cache,B0,EF0,R0,q0,lam,ar


    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    !! time iterator.


    do ii = 1_idef,params%num_species      

       pchunk=params%pchunk
       q_cache=spp(ii)%q
       m_cache=spp(ii)%m

       !$OMP PARALLEL DO default(none) &
       !$OMP& FIRSTPRIVATE(q_cache,m_cache,pchunk) &
       !$OMP& SHARED(params,ii,spp,P,F) &
       !$OMP& PRIVATE(pp,tt,Bmag,cc,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,B_R,B_PHI,B_Z, &
       !$OMP& flagCon,flagCol,E_PHI,PSIp,curlb_R,curlb_PHI,curlb_Z, &
       !$OMP& gradB_R,gradB_PHI,gradB_Z,ne)
       do pp=1_idef,spp(ii)%ppp,pchunk

          !          write(output_unit_write,'("pp: ",I16)') pp

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)
             Y_PHI(cc)=spp(ii)%vars%Y(pp-1+cc,2)
             Y_Z(cc)=spp(ii)%vars%Y(pp-1+cc,3)

             V_PLL(cc)=spp(ii)%vars%V(pp-1+cc,1)
             V_MU(cc)=spp(ii)%vars%V(pp-1+cc,2)

             PSIp(cc)=spp(ii)%vars%PSI_P(pp-1+cc)

             flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
             flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)           
          end do
          !$OMP END SIMD

          if (.not.params%FokPlan) then       
             do tt=1_ip,params%t_skip
                call advance_GCinterp_B2D_vars(spp(ii)%vars,pp,tt,params, &
                     Y_R,Y_PHI,Y_Z,V_PLL,V_MU,q_cache,m_cache,flagCon,flagCol, &
                     F,P,B_R,B_PHI,B_Z,E_PHI,curlb_R,curlb_PHI,curlb_Z, &
                     gradB_R,gradB_PHI,gradB_Z,PSIp)
             end do !timestep iterator


             !$OMP SIMD
             do cc=1_idef,pchunk
                spp(ii)%vars%Y(pp-1+cc,1)=Y_R(cc)
                spp(ii)%vars%Y(pp-1+cc,2)=Y_PHI(cc)
                spp(ii)%vars%Y(pp-1+cc,3)=Y_Z(cc)
                spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                spp(ii)%vars%flagCon(pp-1+cc)=flagCon(cc)
                spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                spp(ii)%vars%B(pp-1+cc,1) = B_R(cc)
                spp(ii)%vars%B(pp-1+cc,2) = B_PHI(cc)
                spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                spp(ii)%vars%gradB(pp-1+cc,1) = gradB_R(cc)
                spp(ii)%vars%gradB(pp-1+cc,2) = gradB_PHI(cc)
                spp(ii)%vars%gradB(pp-1+cc,3) = gradB_Z(cc)

                spp(ii)%vars%curlb(pp-1+cc,1) = curlb_R(cc)
                spp(ii)%vars%curlb(pp-1+cc,2) = curlb_PHI(cc)
                spp(ii)%vars%curlb(pp-1+cc,3) = curlb_Z(cc)

                spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
             end do
             !$OMP END SIMD

          else

             call advance_FPinterp_vars(params,Y_R,Y_PHI, &
                  Y_Z,V_PLL,V_MU,m_cache,flagCon,flagCol,F,P,E_PHI,ne,PSIp)

             !$OMP SIMD
             do cc=1_idef,pchunk
                spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
             end do
             !$OMP END SIMD

          end if


          !$OMP SIMD
          do cc=1_idef,pchunk
             B_R(cc)=spp(ii)%vars%B(pp-1+cc,1)
             B_PHI(cc)=spp(ii)%vars%B(pp-1+cc,2)
             B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)

             Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

             spp(ii)%vars%g(pp-1+cc)=sqrt(1+V_PLL(cc)**2+ &
                  2*V_MU(cc)*Bmag(cc))

             spp(ii)%vars%eta(pp-1+cc) = atan2(sqrt(2*m_cache*Bmag(cc)* &
                  spp(ii)%vars%V(pp-1+cc,2)),spp(ii)%vars%V(pp-1+cc,1))* &
                  180.0_rp/C_PI
          end do
          !$OMP END SIMD

       end do !particle chunk iterator
       !$OMP END PARALLEL DO

    end do !species iterator

  end subroutine adv_GCinterp_B2D_top

  subroutine adv_GCinterp_2DBdB_top(params,spp,P,F)

    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(params%pchunk)               :: Bmag
    REAL(rp),DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_R,E_PHI,E_Z
    REAL(rp),DIMENSION(params%pchunk) :: ne,Te,Zeff    
    REAL(rp),DIMENSION(params%pchunk) :: V_PLL,V_MU
    REAL(rp),DIMENSION(params%pchunk) :: PSIp
    REAL(rp),DIMENSION(params%pchunk) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp),DIMENSION(params%pchunk) :: gradB_R,gradB_PHI,gradB_Z
    INTEGER(is),DIMENSION(params%pchunk) :: flagCon,flagCol
    REAL(rp) :: m_cache,q_cache,B0,EF0,R0,q0,lam,ar


    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    !! time iterator.


    do ii = 1_idef,params%num_species      

       pchunk=params%pchunk
       q_cache=spp(ii)%q
       m_cache=spp(ii)%m

       !$OMP PARALLEL DO default(none) &
       !$OMP& FIRSTPRIVATE(q_cache,m_cache,pchunk) &
       !$OMP& SHARED(params,ii,spp,P,F) &
       !$OMP& PRIVATE(pp,tt,Bmag,cc,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,B_R,B_PHI,B_Z, &
       !$OMP& flagCon,flagCol,E_PHI,PSIp,curlb_R,curlb_PHI,curlb_Z, &
       !$OMP& gradB_R,gradB_PHI,gradB_Z,ne)
       do pp=1_idef,spp(ii)%ppp,pchunk

          !          write(output_unit_write,'("pp: ",I16)') pp

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)
             Y_PHI(cc)=spp(ii)%vars%Y(pp-1+cc,2)
             Y_Z(cc)=spp(ii)%vars%Y(pp-1+cc,3)

             V_PLL(cc)=spp(ii)%vars%V(pp-1+cc,1)
             V_MU(cc)=spp(ii)%vars%V(pp-1+cc,2)

             PSIp(cc)=spp(ii)%vars%PSI_P(pp-1+cc)

             flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
             flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)           
          end do
          !$OMP END SIMD

          if (.not.params%FokPlan) then       
             do tt=1_ip,params%t_skip
                call advance_GCinterp_2DBdB_vars(spp(ii)%vars,pp,tt,params, &
                     Y_R,Y_PHI,Y_Z,V_PLL,V_MU,q_cache,m_cache,flagCon,flagCol, &
                     F,P,B_R,B_PHI,B_Z,E_PHI,curlb_R,curlb_PHI,curlb_Z, &
                     gradB_R,gradB_PHI,gradB_Z,PSIp)
             end do !timestep iterator


             !$OMP SIMD
             do cc=1_idef,pchunk
                spp(ii)%vars%Y(pp-1+cc,1)=Y_R(cc)
                spp(ii)%vars%Y(pp-1+cc,2)=Y_PHI(cc)
                spp(ii)%vars%Y(pp-1+cc,3)=Y_Z(cc)
                spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                spp(ii)%vars%flagCon(pp-1+cc)=flagCon(cc)
                spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                spp(ii)%vars%B(pp-1+cc,1) = B_R(cc)
                spp(ii)%vars%B(pp-1+cc,2) = B_PHI(cc)
                spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                spp(ii)%vars%gradB(pp-1+cc,1) = gradB_R(cc)
                spp(ii)%vars%gradB(pp-1+cc,2) = gradB_PHI(cc)
                spp(ii)%vars%gradB(pp-1+cc,3) = gradB_Z(cc)

                spp(ii)%vars%curlb(pp-1+cc,1) = curlb_R(cc)
                spp(ii)%vars%curlb(pp-1+cc,2) = curlb_PHI(cc)
                spp(ii)%vars%curlb(pp-1+cc,3) = curlb_Z(cc)

                spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
                spp(ii)%vars%PSI_P(pp-1+cc) = PSIp(cc)
             end do
             !$OMP END SIMD

          else

             call advance_FPinterp_vars(params,Y_R,Y_PHI, &
                  Y_Z,V_PLL,V_MU,m_cache,flagCon,flagCol,F,P,E_PHI,ne,PSIp)

             !$OMP SIMD
             do cc=1_idef,pchunk
                spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                spp(ii)%vars%flagCon(pp-1+cc)=flagCon(cc)
                spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
             end do
             !$OMP END SIMD

          end if


          !$OMP SIMD
          do cc=1_idef,pchunk
             B_R(cc)=spp(ii)%vars%B(pp-1+cc,1)
             B_PHI(cc)=spp(ii)%vars%B(pp-1+cc,2)
             B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)

             Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

             spp(ii)%vars%g(pp-1+cc)=sqrt(1+V_PLL(cc)**2+ &
                  2*V_MU(cc)*Bmag(cc))

             spp(ii)%vars%eta(pp-1+cc) = atan2(sqrt(2*m_cache*Bmag(cc)* &
                  spp(ii)%vars%V(pp-1+cc,2)),spp(ii)%vars%V(pp-1+cc,1))* &
                  180.0_rp/C_PI
          end do
          !$OMP END SIMD

       end do !particle chunk iterator
       !$OMP END PARALLEL DO

    end do !species iterator

  end subroutine adv_GCinterp_2DBdB_top

  subroutine adv_GCinterp_3DBdB1_top(params,spp,P,F)

    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(params%pchunk)               :: Bmag
    REAL(rp),DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_R,E_PHI,E_Z
    REAL(rp),DIMENSION(params%pchunk) :: ne,Te,Zeff    
    REAL(rp),DIMENSION(params%pchunk) :: V_PLL,V_MU
    REAL(rp),DIMENSION(params%pchunk) :: PSIp
    REAL(rp),DIMENSION(params%pchunk) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp),DIMENSION(params%pchunk) :: gradB_R,gradB_PHI,gradB_Z
    INTEGER(is),DIMENSION(params%pchunk) :: flagCon,flagCol
    REAL(rp) :: m_cache,q_cache,B0,EF0,R0,q0,lam,ar


    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    !! time iterator.


    do ii = 1_idef,params%num_species      

       q_cache=spp(ii)%q
       m_cache=spp(ii)%m
       pchunk=params%pchunk

       !$OMP PARALLEL DO default(none) &
       !$OMP& FIRSTPRIVATE(q_cache,m_cache,pchunk) &
       !$OMP& SHARED(params,ii,spp,P,F) &
       !$OMP& PRIVATE(pp,tt,Bmag,cc,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,B_R,B_PHI,B_Z, &
       !$OMP& flagCon,flagCol,E_PHI,PSIp,curlb_R,curlb_PHI,curlb_Z, &
       !$OMP& gradB_R,gradB_PHI,gradB_Z,ne)
       do pp=1_idef,spp(ii)%ppp,pchunk

          !          write(output_unit_write,'("pp: ",I16)') pp

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)
             Y_PHI(cc)=spp(ii)%vars%Y(pp-1+cc,2)
             Y_Z(cc)=spp(ii)%vars%Y(pp-1+cc,3)

             V_PLL(cc)=spp(ii)%vars%V(pp-1+cc,1)
             V_MU(cc)=spp(ii)%vars%V(pp-1+cc,2)

             PSIp(cc)=spp(ii)%vars%PSI_P(pp-1+cc)

             flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
             flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)           
          end do
          !$OMP END SIMD

          if (.not.params%FokPlan) then       
             do tt=1_ip,params%t_skip
                call advance_GCinterp_3DBdB1_vars(spp(ii)%vars,pp,tt,params, &
                     Y_R,Y_PHI,Y_Z,V_PLL,V_MU,q_cache,m_cache,flagCon,flagCol, &
                     F,P,B_R,B_PHI,B_Z,E_PHI,curlb_R,curlb_PHI,curlb_Z, &
                     gradB_R,gradB_PHI,gradB_Z,PSIp)
             end do !timestep iterator


             !$OMP SIMD
             do cc=1_idef,pchunk
                spp(ii)%vars%Y(pp-1+cc,1)=Y_R(cc)
                spp(ii)%vars%Y(pp-1+cc,2)=Y_PHI(cc)
                spp(ii)%vars%Y(pp-1+cc,3)=Y_Z(cc)
                spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                spp(ii)%vars%flagCon(pp-1+cc)=flagCon(cc)
                spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                spp(ii)%vars%B(pp-1+cc,1) = B_R(cc)
                spp(ii)%vars%B(pp-1+cc,2) = B_PHI(cc)
                spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                spp(ii)%vars%gradB(pp-1+cc,1) = gradB_R(cc)
                spp(ii)%vars%gradB(pp-1+cc,2) = gradB_PHI(cc)
                spp(ii)%vars%gradB(pp-1+cc,3) = gradB_Z(cc)

                spp(ii)%vars%curlb(pp-1+cc,1) = curlb_R(cc)
                spp(ii)%vars%curlb(pp-1+cc,2) = curlb_PHI(cc)
                spp(ii)%vars%curlb(pp-1+cc,3) = curlb_Z(cc)

                spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
                spp(ii)%vars%PSI_P(pp-1+cc) = PSIp(cc)
             end do
             !$OMP END SIMD

          else

             call advance_FPinterp_vars(params,Y_R,Y_PHI, &
                  Y_Z,V_PLL,V_MU,m_cache,flagCon,flagCol,F,P,E_PHI,ne,PSIp)

             !$OMP SIMD
             do cc=1_idef,pchunk
                spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
             end do
             !$OMP END SIMD

          end if


          !$OMP SIMD
          do cc=1_idef,pchunk
             B_R(cc)=spp(ii)%vars%B(pp-1+cc,1)
             B_PHI(cc)=spp(ii)%vars%B(pp-1+cc,2)
             B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)

             Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

             spp(ii)%vars%g(pp-1+cc)=sqrt(1+V_PLL(cc)**2+ &
                  2*V_MU(cc)*Bmag(cc))

             spp(ii)%vars%eta(pp-1+cc) = atan2(sqrt(2*m_cache*Bmag(cc)* &
                  spp(ii)%vars%V(pp-1+cc,2)),spp(ii)%vars%V(pp-1+cc,1))* &
                  180.0_rp/C_PI
          end do
          !$OMP END SIMD

       end do !particle chunk iterator
       !$OMP END PARALLEL DO

    end do !species iterator

  end subroutine adv_GCinterp_3DBdB1_top

  subroutine adv_GCinterp_3DBdB_top(params,spp,P,F)

    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp), DIMENSION(params%pchunk)               :: Bmag
    REAL(rp),DIMENSION(params%pchunk) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_R,E_PHI,E_Z
    REAL(rp),DIMENSION(params%pchunk) :: ne,Te,Zeff    
    REAL(rp),DIMENSION(params%pchunk) :: V_PLL,V_MU,PSIp
    REAL(rp),DIMENSION(params%pchunk) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp),DIMENSION(params%pchunk) :: gradB_R,gradB_PHI,gradB_Z
    INTEGER(is),DIMENSION(params%pchunk) :: flagCon,flagCol
    REAL(rp) :: m_cache,q_cache,B0,EF0,R0,q0,lam,ar


    INTEGER                                                    :: ii
    !! Species iterator.
    INTEGER                                                    :: pp
    !! Particles iterator.
    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    INTEGER(ip)                                                    :: tt
    !! time iterator.


    do ii = 1_idef,params%num_species      

       q_cache=spp(ii)%q
       m_cache=spp(ii)%m
       pchunk=params%pchunk

       !$OMP PARALLEL DO default(none) &
       !$OMP& FIRSTPRIVATE(q_cache,m_cache,pchunk) &
       !$OMP& SHARED(params,ii,spp,P,F) &
       !$OMP& PRIVATE(pp,tt,Bmag,cc,Y_R,Y_PHI,Y_Z,V_PLL,V_MU,B_R,B_PHI,B_Z, &
       !$OMP& flagCon,flagCol,E_PHI,curlb_R,curlb_PHI,curlb_Z, &
       !$OMP& gradB_R,gradB_PHI,gradB_Z,ne,PSIp)
       do pp=1_idef,spp(ii)%ppp,pchunk

          !          write(output_unit_write,'("pp: ",I16)') pp

          !$OMP SIMD
          do cc=1_idef,pchunk
             Y_R(cc)=spp(ii)%vars%Y(pp-1+cc,1)
             Y_PHI(cc)=spp(ii)%vars%Y(pp-1+cc,2)
             Y_Z(cc)=spp(ii)%vars%Y(pp-1+cc,3)

             V_PLL(cc)=spp(ii)%vars%V(pp-1+cc,1)
             V_MU(cc)=spp(ii)%vars%V(pp-1+cc,2)

             PSIp(cc)=spp(ii)%vars%PSI_P(pp-1+cc)

             flagCon(cc)=spp(ii)%vars%flagCon(pp-1+cc)
             flagCol(cc)=spp(ii)%vars%flagCol(pp-1+cc)           
          end do
          !$OMP END SIMD

          if (.not.params%FokPlan) then       
             do tt=1_ip,params%t_skip
                call advance_GCinterp_3DBdB_vars(spp(ii)%vars,pp,tt,params, &
                     Y_R,Y_PHI,Y_Z,V_PLL,V_MU,q_cache,m_cache,flagCon,flagCol, &
                     F,P,B_R,B_PHI,B_Z,E_PHI,curlb_R,curlb_PHI,curlb_Z, &
                     gradB_R,gradB_PHI,gradB_Z,PSIp)
             end do !timestep iterator


             !$OMP SIMD
             do cc=1_idef,pchunk
                spp(ii)%vars%Y(pp-1+cc,1)=Y_R(cc)
                spp(ii)%vars%Y(pp-1+cc,2)=Y_PHI(cc)
                spp(ii)%vars%Y(pp-1+cc,3)=Y_Z(cc)
                spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                spp(ii)%vars%flagCon(pp-1+cc)=flagCon(cc)
                spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                spp(ii)%vars%B(pp-1+cc,1) = B_R(cc)
                spp(ii)%vars%B(pp-1+cc,2) = B_PHI(cc)
                spp(ii)%vars%B(pp-1+cc,3) = B_Z(cc)

                spp(ii)%vars%gradB(pp-1+cc,1) = gradB_R(cc)
                spp(ii)%vars%gradB(pp-1+cc,2) = gradB_PHI(cc)
                spp(ii)%vars%gradB(pp-1+cc,3) = gradB_Z(cc)

                spp(ii)%vars%curlb(pp-1+cc,1) = curlb_R(cc)
                spp(ii)%vars%curlb(pp-1+cc,2) = curlb_PHI(cc)
                spp(ii)%vars%curlb(pp-1+cc,3) = curlb_Z(cc)

                spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
             end do
             !$OMP END SIMD

          else

             call advance_FPinterp_vars(params,Y_R,Y_PHI, &
                  Y_Z,V_PLL,V_MU,m_cache,flagCon,flagCol,F,P,E_PHI,ne,PSIp)

             !$OMP SIMD
             do cc=1_idef,pchunk
                spp(ii)%vars%V(pp-1+cc,1)=V_PLL(cc)
                spp(ii)%vars%V(pp-1+cc,2)=V_MU(cc)

                spp(ii)%vars%flagCol(pp-1+cc)=flagCol(cc)

                spp(ii)%vars%E(pp-1+cc,2) = E_PHI(cc)
             end do
             !$OMP END SIMD

          end if


          !$OMP SIMD
          do cc=1_idef,pchunk
             B_R(cc)=spp(ii)%vars%B(pp-1+cc,1)
             B_PHI(cc)=spp(ii)%vars%B(pp-1+cc,2)
             B_Z(cc)=spp(ii)%vars%B(pp-1+cc,3)

             Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

             spp(ii)%vars%g(pp-1+cc)=sqrt(1+V_PLL(cc)**2+ &
                  2*V_MU(cc)*Bmag(cc))

             spp(ii)%vars%eta(pp-1+cc) = atan2(sqrt(2*m_cache*Bmag(cc)* &
                  spp(ii)%vars%V(pp-1+cc,2)),spp(ii)%vars%V(pp-1+cc,1))* &
                  180.0_rp/C_PI
          end do
          !$OMP END SIMD

       end do !particle chunk iterator
       !$OMP END PARALLEL DO

    end do !species iterator

  end subroutine adv_GCinterp_3DBdB_top


  subroutine advance_GCinterp_psi_vars_FS(vars,pp,tt,params,Y_R,Y_PHI,Y_Z, &
       V_PLL,V_MU,q_cache,m_cache,flagCon,flagCol,F,P,B_R,B_PHI,B_Z,E_PHI,PSIp, &
       curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z)
    !! @note Subroutine to advance GC variables \(({\bf X},p_\parallel)\)
    !! @endnote
    !! Comment this section further with evolution equations, numerical
    !! methods, and descriptions of both.
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(PARTICLES), INTENT(INOUT)     :: vars
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    REAL(rp)                                      :: dt
    !! Time step used in the leapfrog step (\(\Delta t\)).

    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    INTEGER(ip),intent(in)                                      :: tt
    !! time iterator.
    INTEGER,intent(in)                                  :: pp


    REAL(rp),DIMENSION(params%pchunk)               :: Bmag
    REAL(rp)              :: a1 = 1./5._rp
    REAL(rp) :: a21 = 3./40._rp,a22=9./40._rp
    REAL(rp) :: a31 = 3./10._rp,a32=-9./10._rp,a33=6./5._rp
    REAL(rp) :: a41 = -11./54._rp,a42=5./2._rp,a43=-70./27._rp,a44=35./27._rp
    REAL(rp) :: a51 = 1631./55296._rp,a52=175./512._rp,a53=575./13824._rp,a54=44275./110592._rp,a55=253./4096._rp
    REAL(rp) :: b1=37./378._rp,b2=0._rp,b3=250./621._rp,b4=125./594._rp,b5=0._rp,b6=512./1771._rp

    REAL(rp),DIMENSION(params%pchunk) :: k1_R,k1_PHI,k1_Z,k1_PLL,k1_MU
    REAL(rp),DIMENSION(params%pchunk) :: k2_R,k2_PHI,k2_Z,k2_PLL,k2_MU
    REAL(rp),DIMENSION(params%pchunk) :: k3_R,k3_PHI,k3_Z,k3_PLL,k3_MU
    REAL(rp),DIMENSION(params%pchunk) :: k4_R,k4_PHI,k4_Z,k4_PLL,k4_MU
    REAL(rp),DIMENSION(params%pchunk) :: k5_R,k5_PHI,k5_Z,k5_PLL,k5_MU
    REAL(rp),DIMENSION(params%pchunk) :: k6_R,k6_PHI,k6_Z,k6_PLL,k6_MU
    REAL(rp),DIMENSION(params%pchunk) :: Y0_R,Y0_PHI,Y0_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_R,E_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: E_PHI
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: PSIp
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: gradB_R,gradB_PHI,gradB_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: V_PLL,V_MU
    REAL(rp),DIMENSION(params%pchunk) :: RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU
    REAL(rp),DIMENSION(params%pchunk) :: V0_PLL,V0_MU
    REAL(rp),DIMENSION(params%pchunk) :: ne,Te,Zeff

    INTEGER(is),DIMENSION(params%pchunk),intent(INOUT) :: flagCon,flagCol
    REAL(rp),intent(IN)  :: q_cache,m_cache

    dt=params%dt
    pchunk=params%pchunk

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU)
    do cc=1_idef,pchunk

       Y0_R(cc)=Y_R(cc)
       Y0_PHI(cc)=Y_PHI(cc)
       Y0_Z(cc)=Y_Z(cc)
       V0_PLL(cc)=V_PLL(cc)
       V0_MU(cc)=V_MU(cc)
    end do
    !$OMP END SIMD

    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfields_p_FS(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    if (params%SC_E_add) then
       call add_interp_SCE_p_FS(params,F,PSIp,E_PHI)
    end if


    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !    write(output_unit_write,*) 'R',Y_R(1)
    !    write(output_unit_write,*) 'PHI',Y_PHI(1)
    !    write(output_unit_write,*) 'Z',Y_Z(1)
    !    write(output_unit_write,*) 'PPLL',V_PLL(1)
    !    write(output_unit_write,*) 'MU',V_MU(1)

    !    write(output_unit_write,*) 'BR',B_R(1)
    !    write(output_unit_write,*) 'BPHI',B_PHI(1)
    !    write(output_unit_write,*) 'BZ',B_Z(1)

    !    write(output_unit_write,*) 'gradBR',gradB_R(1)
    !    write(output_unit_write,*) 'gradBPHI',gradB_PHI(1)
    !    write(output_unit_write,*) 'gradBZ',gradB_Z(1)

    !    write(output_unit_write,*) 'curlBR',curlB_R(1)
    !    write(output_unit_write,*) 'curlBPHI',curlB_PHI(1)
    !    write(output_unit_write,*) 'curlBZ',curlB_Z(1)

    !    write(output_unit_write,*) 'RHS_R',RHS_R(1)
    !    write(output_unit_write,*) 'RHS_PHI',RHS_PHI(1)
    !    write(output_unit_write,*) 'RHS_Z',RHS_Z(1)
    !    write(output_unit_write,*) 'RHS_PLL',RHS_PLL(1)

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k1_R,k1_PHI,k1_Z,k1_PLL,k1_MU)
    do cc=1_idef,pchunk
       k1_R(cc)=dt*RHS_R(cc)              
       k1_PHI(cc)=dt*RHS_PHI(cc)    
       k1_Z(cc)=dt*RHS_Z(cc)    
       k1_PLL(cc)=dt*RHS_PLL(cc)
       k1_MU(cc)=dt*RHS_MU(cc)    

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a1*k1_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a1*k1_MU(cc)
    end do
    !$OMP END SIMD

    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfields_p_FS(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)


    if (params%SC_E_add) then
       call add_interp_SCE_p_FS(params,F,PSIp,E_PHI)
    end if


    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k2_R,k2_PHI,k2_Z,k2_PLL,k2_MU)
    do cc=1_idef,pchunk
       k2_R(cc)=dt*RHS_R(cc)    
       k2_PHI(cc)=dt*RHS_PHI (cc)   
       k2_Z(cc)=dt*RHS_Z(cc)   
       k2_PLL(cc)=dt*RHS_PLL(cc)
       k2_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a21*k1_PLL(cc)+a22*k2_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a21*k1_MU(cc)+a22*k2_MU(cc)
    end do
    !$OMP END SIMD

    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfields_p_FS(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)


    if (params%SC_E_add) then
       call add_interp_SCE_p_FS(params,F,PSIp,E_PHI)
    end if


    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k3_R,k3_PHI,k3_Z,k3_PLL,k3_MU)
    do cc=1_idef,pchunk
       k3_R(cc)=dt*RHS_R(cc)   
       k3_PHI(cc)=dt*RHS_PHI(cc)    
       k3_Z(cc)=dt*RHS_Z(cc)    
       k3_PLL(cc)=dt*RHS_PLL(cc)
       k3_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a31*k1_PLL(cc)+a32*k2_PLL(cc)+a33*k3_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a31*k1_MU(cc)+a32*k2_MU(cc)+a33*k3_MU(cc)
    end do
    !$OMP END SIMD

    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfields_p_FS(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)


    if (params%SC_E_add) then
       call add_interp_SCE_p_FS(params,F,PSIp,E_PHI)
    end if


    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)     

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k4_R,k4_PHI,k4_Z,k4_PLL,k4_MU)
    do cc=1_idef,pchunk
       k4_R(cc)=dt*RHS_R(cc)   
       k4_PHI(cc)=dt*RHS_PHI(cc)    
       k4_Z(cc)=dt*RHS_Z(cc)    
       k4_PLL(cc)=dt*RHS_PLL(cc)
       k4_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a41*k1_R(cc)+a42*k2_R(cc)+a43*k3_R(cc)+ &
            a44*k4_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a41*k1_PHI(cc)+a42*k2_PHI(cc)+ &
            a43*k3_PHI(cc)+a44*k4_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a41*k1_Z(cc)+a42*k2_Z(cc)+a43*k3_Z(cc)+ &
            a44*k4_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a41*k1_PLL(cc)+a42*k2_PLL(cc)+ &
            a43*k3_PLL(cc)+a44*k4_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a41*k1_MU(cc)+a42*k2_MU(cc)+ &
            a43*k3_MU(cc)+a44*k4_MU(cc)
    end do
    !$OMP END SIMD


    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfields_p_FS(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)


    if (params%SC_E_add) then
       call add_interp_SCE_p_FS(params,F,PSIp,E_PHI)
    end if


    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)   

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k5_R,k5_PHI,k5_Z,k5_PLL,k5_MU)
    do cc=1_idef,pchunk
       k5_R(cc)=dt*RHS_R(cc)    
       k5_PHI(cc)=dt*RHS_PHI(cc)    
       k5_Z(cc)=dt*RHS_Z(cc)    
       k5_PLL(cc)=dt*RHS_PLL(cc)
       k5_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a51*k1_R(cc)+a52*k2_R(cc)+a53*k3_R(cc)+ &
            a54*k4_R(cc)+a55*k5_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a51*k1_PHI(cc)+a52*k2_PHI(cc)+ &
            a53*k3_PHI(cc)+a54*k4_PHI(cc)+a55*k5_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a51*k1_Z(cc)+a52*k2_Z(cc)+a53*k3_Z(cc)+ &
            a54*k4_Z(cc)+a55*k5_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a51*k1_PLL(cc)+a52*k2_PLL(cc)+ &
            a53*k3_PLL(cc)+a54*k4_PLL(cc)+a55*k5_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a51*k1_MU(cc)+a52*k2_MU(cc)+ &
            a53*k3_MU(cc)+a54*k4_MU(cc)+a55*k5_MU(cc)
    end do
    !$OMP END SIMD

    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfields_p_FS(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)


    if (params%SC_E_add) then
       call add_interp_SCE_p_FS(params,F,PSIp,E_PHI)
    end if


    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)         

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k6_R,k6_PHI,k6_Z,k6_PLL,k6_MU)
    do cc=1_idef,pchunk
       k6_R(cc)=dt*RHS_R(cc)    
       k6_PHI(cc)=dt*RHS_PHI(cc)    
       k6_Z(cc)=dt*RHS_Z(cc)    
       k6_PLL(cc)=dt*RHS_PLL(cc)
       k6_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+b1*k1_R(cc)+b2*k2_R(cc)+ &
            b3*k3_R(cc)+b4*k4_R(cc)+b5*k5_R(cc)+b6*k6_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+b1*k1_PHI(cc)+b2*k2_PHI(cc)+ &
            b3*k3_PHI(cc)+b4*k4_PHI(cc)+b5*k5_PHI(cc)+b6*k6_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+b1*k1_Z(cc)+b2*k2_Z(cc)+ &
            b3*k3_Z(cc)+b4*k4_Z(cc)+b5*k5_Z(cc)+b6*k6_Z(cc)
       V_PLL(cc)=V0_PLL(cc)+b1*k1_PLL(cc)+b2*k2_PLL(cc)+ &
            b3*k3_PLL(cc)+b4*k4_PLL(cc)+b5*k5_PLL(cc)+b6*k6_PLL(cc)
       V_MU(cc)=V0_MU(cc)+b1*k1_MU(cc)+b2*k2_MU(cc)+ &
            b3*k3_MU(cc)+b4*k4_MU(cc)+b5*k5_MU(cc)+b6*k6_MU(cc)
    end do
    !$OMP END SIMD

    !$OMP SIMD
    !    !$OMP& aligned(Y_R,Y_PHI,Y_Z,V_PLL,V_MU,Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU)
    do cc=1_idef,pchunk

       if ((flagCon(cc).eq.0_is).or.(flagCol(cc).eq.0_is)) then
          Y_R(cc)=Y0_R(cc)
          Y_PHI(cc)=Y0_PHI(cc)
          Y_Z(cc)=Y0_Z(cc)
          V_PLL(cc)=V0_PLL(cc)
          V_MU(cc)=V0_MU(cc)
       end if

    end do
    !$OMP END SIMD

    call calculate_GCfields_p_FS(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    if (params%SC_E_add) then
       call add_interp_SCE_p_FS(params,F,PSIp,E_PHI)
    end if

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !$OMP SIMD
    do cc=1_idef,pchunk
       vars%RHS(pp-1+cc,1)=RHS_R(cc)
       vars%RHS(pp-1+cc,2)=RHS_PHI(cc)
       vars%RHS(pp-1+cc,3)=RHS_Z(cc)
       vars%RHS(pp-1+cc,4)=RHS_PLL(cc)
       vars%RHS(pp-1+cc,5)=RHS_MU(cc)
    end do
    !$OMP END SIMD       

    if (params%collisions) then       

       call include_CoulombCollisions_GC_p(tt,params,Y_R,Y_PHI,Y_Z, &
            V_PLL,V_MU,m_cache,flagCon,flagCol,F,P,E_PHI,ne,PSIp)

    end if


  end subroutine advance_GCinterp_psi_vars_FS

  subroutine advance_GCinterp_psi_vars(pchunk,vars,pp,tt,params,Y_R,Y_PHI,Y_Z, &
       V_PLL,V_MU,q_cache,m_cache,flagCon,flagCol,F,P,B_R,B_PHI,B_Z,E_PHI,PSIp, &
       curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,ne)
    !! @note Subroutine to advance GC variables \(({\bf X},p_\parallel)\)
    !! @endnote
    !! Comment this section further with evolution equations, numerical
    !! methods, and descriptions of both.
    INTEGER,intent(in)                                      :: pchunk
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(PARTICLES), INTENT(INOUT)     :: vars
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    REAL(rp)                                      :: dt
    !! Time step used in the leapfrog step (\(\Delta t\)).

    INTEGER                                                    :: cc
    !! Chunk iterator.
    INTEGER(ip),intent(in)                                      :: tt
    !! time iterator.
    INTEGER,intent(in)                                  :: pp


    REAL(rp),DIMENSION(pchunk)               :: Bmag
    REAL(rp)              :: a1 = 1./5._rp
    REAL(rp) :: a21 = 3./40._rp,a22=9./40._rp
    REAL(rp) :: a31 = 3./10._rp,a32=-9./10._rp,a33=6./5._rp
    REAL(rp) :: a41 = -11./54._rp,a42=5./2._rp,a43=-70./27._rp,a44=35./27._rp
    REAL(rp) :: a51 = 1631./55296._rp,a52=175./512._rp,a53=575./13824._rp,a54=44275./110592._rp,a55=253./4096._rp
    REAL(rp) :: b1=37./378._rp,b2=0._rp,b3=250./621._rp,b4=125./594._rp,b5=0._rp,b6=512./1771._rp

    REAL(rp),DIMENSION(pchunk) :: k1_R,k1_PHI,k1_Z,k1_PLL,k1_MU
    REAL(rp),DIMENSION(pchunk) :: k2_R,k2_PHI,k2_Z,k2_PLL,k2_MU
    REAL(rp),DIMENSION(pchunk) :: k3_R,k3_PHI,k3_Z,k3_PLL,k3_MU
    REAL(rp),DIMENSION(pchunk) :: k4_R,k4_PHI,k4_Z,k4_PLL,k4_MU
    REAL(rp),DIMENSION(pchunk) :: k5_R,k5_PHI,k5_Z,k5_PLL,k5_MU
    REAL(rp),DIMENSION(pchunk) :: k6_R,k6_PHI,k6_Z,k6_PLL,k6_MU
    REAL(rp),DIMENSION(pchunk) :: Y0_R,Y0_PHI,Y0_Z
    REAL(rp),DIMENSION(pchunk),INTENT(INOUT) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(pchunk),INTENT(OUT) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(pchunk) :: E_R,E_Z
    REAL(rp),DIMENSION(pchunk),INTENT(OUT) :: E_PHI
    REAL(rp),DIMENSION(pchunk),INTENT(OUT) :: PSIp
    REAL(rp),DIMENSION(pchunk),INTENT(OUT) :: ne
    REAL(rp),DIMENSION(pchunk),INTENT(OUT) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp),DIMENSION(pchunk),INTENT(OUT) :: gradB_R,gradB_PHI,gradB_Z
    REAL(rp),DIMENSION(pchunk),INTENT(INOUT) :: V_PLL,V_MU
    REAL(rp),DIMENSION(pchunk) :: RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU
    REAL(rp),DIMENSION(pchunk) :: V0_PLL,V0_MU
    REAL(rp),DIMENSION(pchunk) :: Te,Zeff

    INTEGER(is),DIMENSION(pchunk),intent(INOUT) :: flagCon,flagCol
    REAL(rp),intent(IN)  :: q_cache,m_cache

    dt=params%dt

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU)
    do cc=1_idef,pchunk
       Y0_R(cc)=Y_R(cc)
       Y0_PHI(cc)=Y_PHI(cc)
       Y0_Z(cc)=Y_Z(cc)
       V0_PLL(cc)=V_PLL(cc)
       V0_MU(cc)=V_MU(cc)
    end do
    !$OMP END SIMD

        !write(output_unit_write,*) 'R0',Y_R(1)
        !write(output_unit_write,*) 'PHI0',Y_PHI(1)
        !write(output_unit_write,*) 'Z0',Y_Z(1)
        !write(output_unit_write,*) 'PPLL0',V_PLL(1)
        !write(output_unit_write,*) 'MU0',V_MU(1)

    !write(output_unit_write,*) 'ER',E_R(1)
    !write(output_unit_write,*) 'EPHI',E_PHI(1)
    !write(output_unit_write,*) 'EZ',E_Z(1)

    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    !write(output_unit_write,*) 'ER',E_R(1)
    !write(output_unit_write,*) 'EPHI',E_PHI(1)
    !write(output_unit_write,*) 'EZ',E_Z(1)
    
    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    !write(output_unit_write,*) 'ER',E_R(1)
    !write(output_unit_write,*) 'EPHI',E_PHI(1)
    !write(output_unit_write,*) 'EZ',E_Z(1)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

        !write(output_unit_write,*) 'R0',Y_R(1)
        !write(output_unit_write,*) 'PHI0',Y_PHI(1)
        !write(output_unit_write,*) 'Z0',Y_Z(1)
        !write(output_unit_write,*) 'PPLL0',V_PLL(1)
        !write(output_unit_write,*) 'MU0',V_MU(1)

    !    write(output_unit_write,*) 'BR',B_R(1)
    !    write(output_unit_write,*) 'BPHI',B_PHI(1)
    !    write(output_unit_write,*) 'BZ',B_Z(1)

    !    write(output_unit_write,*) 'ER',E_R(1)
    !    write(output_unit_write,*) 'EPHI',E_PHI(1)
    !    write(output_unit_write,*) 'EZ',E_Z(1)

    !    write(output_unit_write,*) 'gradBR',gradB_R(1)
    !    write(output_unit_write,*) 'gradBPHI',gradB_PHI(1)
    !    write(output_unit_write,*) 'gradBZ',gradB_Z(1)

    !    write(output_unit_write,*) 'curlBR',curlB_R(1)
    !    write(output_unit_write,*) 'curlBPHI',curlB_PHI(1)
    !    write(output_unit_write,*) 'curlBZ',curlB_Z(1)

    !    write(output_unit_write,*) 'dt',params%dt
    !    write(output_unit_write,*) 'RHS_R',RHS_R(1)
    !    write(output_unit_write,*) 'RHS_PHI',RHS_PHI(1)
    !    write(output_unit_write,*) 'RHS_Z',RHS_Z(1)
    !    write(output_unit_write,*) 'RHS_PLL',RHS_PLL(1)
    !    write(output_unit_write,*) 'RHS_MU',RHS_MU(1)

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k1_R,k1_PHI,k1_Z,k1_PLL,k1_MU)
    do cc=1_idef,pchunk
       k1_R(cc)=dt*RHS_R(cc)              
       k1_PHI(cc)=dt*RHS_PHI(cc)    
       k1_Z(cc)=dt*RHS_Z(cc)    
       k1_PLL(cc)=dt*RHS_PLL(cc)
       k1_MU(cc)=dt*RHS_MU(cc)    

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a1*k1_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a1*k1_MU(cc)

    end do
    !$OMP END SIMD

       ! write(output_unit_write,*) 'R1',Y_R(1)
       ! write(output_unit_write,*) 'PHI1',Y_PHI(1)
       ! write(output_unit_write,*) 'Z1',Y_Z(1)
       ! write(output_unit_write,*) 'PPLL1',V_PLL(1)
       ! write(output_unit_write,*) 'MU1',V_MU(1)

    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)





    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k2_R,k2_PHI,k2_Z,k2_PLL,k2_MU)
    do cc=1_idef,pchunk
       k2_R(cc)=dt*RHS_R(cc)    
       k2_PHI(cc)=dt*RHS_PHI (cc)   
       k2_Z(cc)=dt*RHS_Z(cc)   
       k2_PLL(cc)=dt*RHS_PLL(cc)
       k2_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a21*k1_PLL(cc)+a22*k2_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a21*k1_MU(cc)+a22*k2_MU(cc)


    end do
    !$OMP END SIMD

    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)





    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k3_R,k3_PHI,k3_Z,k3_PLL,k3_MU)
    do cc=1_idef,pchunk
       k3_R(cc)=dt*RHS_R(cc)   
       k3_PHI(cc)=dt*RHS_PHI(cc)    
       k3_Z(cc)=dt*RHS_Z(cc)    
       k3_PLL(cc)=dt*RHS_PLL(cc)
       k3_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)

       V_PLL(cc)=V0_PLL(cc)   +a31*k1_PLL(cc)+a32*k2_PLL(cc)+a33*k3_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a31*k1_MU(cc)+a32*k2_MU(cc)+a33*k3_MU(cc)


    end do
    !$OMP END SIMD

    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)





    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)     

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k4_R,k4_PHI,k4_Z,k4_PLL,k4_MU)
    do cc=1_idef,pchunk
       k4_R(cc)=dt*RHS_R(cc)   
       k4_PHI(cc)=dt*RHS_PHI(cc)    
       k4_Z(cc)=dt*RHS_Z(cc)    
       k4_PLL(cc)=dt*RHS_PLL(cc)
       k4_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a41*k1_R(cc)+a42*k2_R(cc)+a43*k3_R(cc)+ &
            a44*k4_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a41*k1_PHI(cc)+a42*k2_PHI(cc)+ &
            a43*k3_PHI(cc)+a44*k4_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a41*k1_Z(cc)+a42*k2_Z(cc)+a43*k3_Z(cc)+ &
            a44*k4_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a41*k1_PLL(cc)+a42*k2_PLL(cc)+ &
            a43*k3_PLL(cc)+a44*k4_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a41*k1_MU(cc)+a42*k2_MU(cc)+ &
            a43*k3_MU(cc)+a44*k4_MU(cc)


    end do
    !$OMP END SIMD


    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)





    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)   

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k5_R,k5_PHI,k5_Z,k5_PLL,k5_MU)
    do cc=1_idef,pchunk
       k5_R(cc)=dt*RHS_R(cc)    
       k5_PHI(cc)=dt*RHS_PHI(cc)    
       k5_Z(cc)=dt*RHS_Z(cc)    
       k5_PLL(cc)=dt*RHS_PLL(cc)
       k5_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a51*k1_R(cc)+a52*k2_R(cc)+a53*k3_R(cc)+ &
            a54*k4_R(cc)+a55*k5_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a51*k1_PHI(cc)+a52*k2_PHI(cc)+ &
            a53*k3_PHI(cc)+a54*k4_PHI(cc)+a55*k5_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a51*k1_Z(cc)+a52*k2_Z(cc)+a53*k3_Z(cc)+ &
            a54*k4_Z(cc)+a55*k5_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a51*k1_PLL(cc)+a52*k2_PLL(cc)+ &
            a53*k3_PLL(cc)+a54*k4_PLL(cc)+a55*k5_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a51*k1_MU(cc)+a52*k2_MU(cc)+ &
            a53*k3_MU(cc)+a54*k4_MU(cc)+a55*k5_MU(cc)


    end do
    !$OMP END SIMD

    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)


    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)         

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k6_R,k6_PHI,k6_Z,k6_PLL,k6_MU)
    do cc=1_idef,pchunk
       k6_R(cc)=dt*RHS_R(cc)    
       k6_PHI(cc)=dt*RHS_PHI(cc)    
       k6_Z(cc)=dt*RHS_Z(cc)    
       k6_PLL(cc)=dt*RHS_PLL(cc)
       k6_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+b1*k1_R(cc)+b2*k2_R(cc)+ &
            b3*k3_R(cc)+b4*k4_R(cc)+b5*k5_R(cc)+b6*k6_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+b1*k1_PHI(cc)+b2*k2_PHI(cc)+ &
            b3*k3_PHI(cc)+b4*k4_PHI(cc)+b5*k5_PHI(cc)+b6*k6_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+b1*k1_Z(cc)+b2*k2_Z(cc)+ &
            b3*k3_Z(cc)+b4*k4_Z(cc)+b5*k5_Z(cc)+b6*k6_Z(cc)
       V_PLL(cc)=V0_PLL(cc)+b1*k1_PLL(cc)+b2*k2_PLL(cc)+ &
            b3*k3_PLL(cc)+b4*k4_PLL(cc)+b5*k5_PLL(cc)+b6*k6_PLL(cc)
       V_MU(cc)=V0_MU(cc)+b1*k1_MU(cc)+b2*k2_MU(cc)+ &
            b3*k3_MU(cc)+b4*k4_MU(cc)+b5*k5_MU(cc)+b6*k6_MU(cc)


    end do
    !$OMP END SIMD

    !$OMP SIMD
    !    !$OMP& aligned(Y_R,Y_PHI,Y_Z,V_PLL,V_MU,Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU)
    do cc=1_idef,pchunk

       if ((flagCon(cc).eq.0_is).or.(flagCol(cc).eq.0_is)) then
          Y_R(cc)=Y0_R(cc)
          Y_PHI(cc)=Y0_PHI(cc)
          Y_Z(cc)=Y0_Z(cc)
          V_PLL(cc)=V0_PLL(cc)
          V_MU(cc)=V0_MU(cc)
       end if

    end do
    !$OMP END SIMD

    call calculate_GCfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)



    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !$OMP SIMD
    do cc=1_idef,pchunk
       vars%RHS(pp-1+cc,1)=RHS_R(cc)
       vars%RHS(pp-1+cc,2)=RHS_PHI(cc)
       vars%RHS(pp-1+cc,3)=RHS_Z(cc)
       vars%RHS(pp-1+cc,4)=RHS_PLL(cc)
       vars%RHS(pp-1+cc,5)=RHS_MU(cc)
    end do
    !$OMP END SIMD       




  end subroutine advance_GCinterp_psi_vars
#endif

#ifdef FIO
  subroutine advance_GCinterp_fio_vars(vars,pp,tt,params,Y_R,Y_PHI,Y_Z, &
       V_PLL,V_MU,q_cache,m_cache,flagCon,flagCol,F,P,B_R,B_PHI,B_Z, &
       E_R,E_PHI,E_Z,PSIp,curlb_R,curlb_PHI,curlb_Z, &
       gradB_R,gradB_PHI,gradB_Z,ne,ni,Te,Zeff,nimp,hint)

    USE omp_lib
    IMPLICIT NONE
    
    !! @note Subroutine to advance GC variables \(({\bf X},p_\parallel)\)
    !! @endnote
    !! Comment this section further with evolution equations, numerical
    !! methods, and descriptions of both.
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(PARTICLES), INTENT(INOUT)     :: vars
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    REAL(rp)                                      :: dt
    !! Time step used in the leapfrog step (\(\Delta t\)).

    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    INTEGER(ip),intent(in)                                      :: tt
    !! time iterator.
    INTEGER,intent(in)                                  :: pp


    REAL(rp),DIMENSION(params%pchunk)               :: Bmag
    REAL(rp)              :: a1 = 1./5._rp
    REAL(rp) :: a21 = 3./40._rp,a22=9./40._rp
    REAL(rp) :: a31 = 3./10._rp,a32=-9./10._rp,a33=6./5._rp
    REAL(rp) :: a41 = -11./54._rp,a42=5./2._rp,a43=-70./27._rp,a44=35./27._rp
    REAL(rp) :: a51 = 1631./55296._rp,a52=175./512._rp,a53=575./13824._rp,a54=44275./110592._rp,a55=253./4096._rp
    REAL(rp) :: b1=37./378._rp,b2=0._rp,b3=250./621._rp,b4=125./594._rp,b5=0._rp,b6=512./1771._rp

    REAL(rp),DIMENSION(params%pchunk) :: k1_R,k1_PHI,k1_Z,k1_PLL,k1_MU
    REAL(rp),DIMENSION(params%pchunk) :: k2_R,k2_PHI,k2_Z,k2_PLL,k2_MU
    REAL(rp),DIMENSION(params%pchunk) :: k3_R,k3_PHI,k3_Z,k3_PLL,k3_MU
    REAL(rp),DIMENSION(params%pchunk) :: k4_R,k4_PHI,k4_Z,k4_PLL,k4_MU
    REAL(rp),DIMENSION(params%pchunk) :: k5_R,k5_PHI,k5_Z,k5_PLL,k5_MU
    REAL(rp),DIMENSION(params%pchunk) :: k6_R,k6_PHI,k6_Z,k6_PLL,k6_MU
    REAL(rp),DIMENSION(params%pchunk) :: Y0_R,Y0_PHI,Y0_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: E_PHI,E_R,E_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: PSIp
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: ne,Te,ni
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: Zeff
    REAL(rp),DIMENSION(params%pchunk,params%num_impurity_species)&
         &,INTENT(INOUT) :: nimp
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: curlb_R&
         &,curlb_PHI,curlb_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: gradB_R,gradB_PHI,gradB_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: V_PLL,V_MU
    REAL(rp),DIMENSION(params%pchunk) :: RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU
    REAL(rp),DIMENSION(params%pchunk) :: V0_PLL,V0_MU

    INTEGER(is),DIMENSION(params%pchunk),intent(INOUT) :: flagCon,flagCol
    REAL(rp),intent(IN)  :: q_cache,m_cache
    TYPE(C_PTR), DIMENSION(params%pchunk), INTENT(INOUT)  :: hint
    INTEGER             :: thread_num


    thread_num = OMP_GET_THREAD_NUM()

    
    dt=params%dt
    pchunk=params%pchunk

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU)
    do cc=1_idef,pchunk

       Y0_R(cc)=Y_R(cc)
       Y0_PHI(cc)=Y_PHI(cc)
       Y0_Z(cc)=Y_Z(cc)
       V0_PLL(cc)=V_PLL(cc)
       V0_MU(cc)=V_MU(cc)
    end do
    !$OMP END SIMD

    !write(output_unit_write,*) 'MPI',params%mpi_params%rank,'OMP',thread_num,'R0',Y_R(1)
    !write(output_unit_write,*) 'MPI',params%mpi_params%rank,'OMP',thread_num,'PHI0',Y_PHI(1)
    !write(output_unit_write,*) 'MPI',params%mpi_params%rank,'OMP',thread_num,'Z0',Y_Z(1)
    !write(output_unit_write,*) 'MPI',params%mpi_params%rank,'OMP',thread_num,'PPLL0',V_PLL(1)
    !write(output_unit_write,*) 'MPI',params%mpi_params%rank,'OMP',thread_num,'MU0',V_MU(1)


    call get_fio_GCmagnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       B_R,B_PHI,B_Z,gradB_R,gradB_PHI,gradB_Z, &
       curlb_R,curlb_PHI,curlb_Z,flagCon,hint)
    if (F%FIO_E .ge. 0) then
       call get_fio_GCelectric_fields_p(params,F, &
            Y_R,Y_PHI,Y_Z,E_R,E_PHI,E_Z,flagCon,hint)
    end if
    call get_fio_vector_potential_p(params,F,Y_R,Y_PHI,Y_Z, &
         PSIp,flagCon,hint)

    !write(output_unit_write,*) 'MPI',params%mpi_params%rank,'OMP',thread_num,'B',B_R(1),B_PHI(1),B_Z(1)
    !write(output_unit_write,*) 'MPI',params%mpi_params%rank,'OMP',thread_num,'gradB',gradB_R(1),gradB_PHI(1)
    !write(output_unit_write,*) 'MPI',params%mpi_params%rank,'OMP',thread_num,'curlB',curlB_R(1),curlB_PHI(1),curlB_Z(1)

    call GCEoM1_fio_p(tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
         B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp, &
         ne,ni,nimp,Te,Zeff,flagCon,flagCol,hint)

    !write(output_unit_write,*) 'R',Y_R(1)
    !write(output_unit_write,*) 'PHI',Y_PHI(1)
    !write(output_unit_write,*) 'Z',Y_Z(1)
    !write(output_unit_write,*) 'PPLL',V_PLL(1)
    !write(output_unit_write,*) 'MU',V_MU(1)

    !write(output_unit_write,*) 'BR',B_R(1)
    !write(output_unit_write,*) 'BPHI',B_PHI(1)
    !write(output_unit_write,*) 'BZ',B_Z(1)

    !    write(output_unit_write,*) 'gradBR',gradB_R(1)
    !    write(output_unit_write,*) 'gradBPHI',gradB_PHI(1)
    !    write(output_unit_write,*) 'gradBZ',gradB_Z(1)

    !    write(output_unit_write,*) 'curlBR',curlB_R(1)
    !    write(output_unit_write,*) 'curlBPHI',curlB_PHI(1)
    !    write(output_unit_write,*) 'curlBZ',curlB_Z(1)

    !    write(output_unit_write,*) 'dt',params%dt
    !write(output_unit_write,*) 'RHS_R',RHS_R(1)
    !write(output_unit_write,*) 'RHS_PHI',RHS_PHI(1)
    !write(output_unit_write,*) 'RHS_Z',RHS_Z(1)
    !write(output_unit_write,*) 'RHS_PLL',RHS_PLL(1)
    !write(output_unit_write,*) 'RHS_MU',RHS_MU(1)


    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k1_R,k1_PHI,k1_Z,k1_PLL,k1_MU)
    do cc=1_idef,pchunk
       k1_R(cc)=dt*RHS_R(cc)              
       k1_PHI(cc)=dt*RHS_PHI(cc)    
       k1_Z(cc)=dt*RHS_Z(cc)    
       k1_PLL(cc)=dt*RHS_PLL(cc)
       k1_MU(cc)=dt*RHS_MU(cc)    

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a1*k1_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a1*k1_MU(cc)
    end do
    !$OMP END SIMD
    

    !write(output_unit_write,*) 'R1',Y_R(1)
    !write(output_unit_write,*) 'PHI1',Y_PHI(1)
    !write(output_unit_write,*) 'Z1',Y_Z(1)
    !write(output_unit_write,*) 'PPLL1',V_PLL(1)
    !write(output_unit_write,*) 'MU1',V_MU(1)

    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call get_fio_GCmagnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       B_R,B_PHI,B_Z,gradB_R,gradB_PHI,gradB_Z, &
       curlb_R,curlb_PHI,curlb_Z,flagCon,hint)
    if (F%FIO_E .ge. 0) then
       call get_fio_GCelectric_fields_p(params,F, &
            Y_R,Y_PHI,Y_Z,E_R,E_PHI,E_Z,flagCon,hint)
    end if
    call get_fio_vector_potential_p(params,F,Y_R,Y_PHI,Y_Z, &
         PSIp,flagCon,hint)

    call GCEoM1_fio_p(tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
         B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp, &
         ne,ni,nimp,Te,Zeff,flagCon,flagCol,hint)

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k2_R,k2_PHI,k2_Z,k2_PLL,k2_MU)
    do cc=1_idef,pchunk
       k2_R(cc)=dt*RHS_R(cc)    
       k2_PHI(cc)=dt*RHS_PHI (cc)   
       k2_Z(cc)=dt*RHS_Z(cc)   
       k2_PLL(cc)=dt*RHS_PLL(cc)
       k2_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a21*k1_PLL(cc)+a22*k2_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a21*k1_MU(cc)+a22*k2_MU(cc)
    end do
    !$OMP END SIMD


    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call get_fio_GCmagnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       B_R,B_PHI,B_Z,gradB_R,gradB_PHI,gradB_Z, &
       curlb_R,curlb_PHI,curlb_Z,flagCon,hint)
    if (F%FIO_E .ge. 0) then
       call get_fio_GCelectric_fields_p(params,F, &
            Y_R,Y_PHI,Y_Z,E_R,E_PHI,E_Z,flagCon,hint)
    end if
    call get_fio_vector_potential_p(params,F,Y_R,Y_PHI,Y_Z, &
         PSIp,flagCon,hint)

    call GCEoM1_fio_p(tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
         B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp, &
         ne,ni,nimp,Te,Zeff,flagCon,flagCol,hint) 

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k3_R,k3_PHI,k3_Z,k3_PLL,k3_MU)
    do cc=1_idef,pchunk
       k3_R(cc)=dt*RHS_R(cc)   
       k3_PHI(cc)=dt*RHS_PHI(cc)    
       k3_Z(cc)=dt*RHS_Z(cc)    
       k3_PLL(cc)=dt*RHS_PLL(cc)
       k3_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)

       V_PLL(cc)=V0_PLL(cc)   +a31*k1_PLL(cc)+a32*k2_PLL(cc)+a33*k3_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a31*k1_MU(cc)+a32*k2_MU(cc)+a33*k3_MU(cc)
    end do
    !$OMP END SIMD


    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call get_fio_GCmagnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       B_R,B_PHI,B_Z,gradB_R,gradB_PHI,gradB_Z, &
       curlb_R,curlb_PHI,curlb_Z,flagCon,hint)
    if (F%FIO_E .ge. 0) then
       call get_fio_GCelectric_fields_p(params,F, &
            Y_R,Y_PHI,Y_Z,E_R,E_PHI,E_Z,flagCon,hint)
    end if
    call get_fio_vector_potential_p(params,F,Y_R,Y_PHI,Y_Z, &
         PSIp,flagCon,hint)

    call GCEoM1_fio_p(tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
         B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp, &
         ne,ni,nimp,Te,Zeff,flagCon,flagCol,hint)   

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k4_R,k4_PHI,k4_Z,k4_PLL,k4_MU)
    do cc=1_idef,pchunk
       k4_R(cc)=dt*RHS_R(cc)   
       k4_PHI(cc)=dt*RHS_PHI(cc)    
       k4_Z(cc)=dt*RHS_Z(cc)    
       k4_PLL(cc)=dt*RHS_PLL(cc)
       k4_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a41*k1_R(cc)+a42*k2_R(cc)+a43*k3_R(cc)+ &
            a44*k4_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a41*k1_PHI(cc)+a42*k2_PHI(cc)+ &
            a43*k3_PHI(cc)+a44*k4_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a41*k1_Z(cc)+a42*k2_Z(cc)+a43*k3_Z(cc)+ &
            a44*k4_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a41*k1_PLL(cc)+a42*k2_PLL(cc)+ &
            a43*k3_PLL(cc)+a44*k4_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a41*k1_MU(cc)+a42*k2_MU(cc)+ &
            a43*k3_MU(cc)+a44*k4_MU(cc)
    end do
    !$OMP END SIMD


    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call get_fio_GCmagnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       B_R,B_PHI,B_Z,gradB_R,gradB_PHI,gradB_Z, &
       curlb_R,curlb_PHI,curlb_Z,flagCon,hint)
    if (F%FIO_E .ge. 0) then
       call get_fio_GCelectric_fields_p(params,F, &
            Y_R,Y_PHI,Y_Z,E_R,E_PHI,E_Z,flagCon,hint)
    end if
    call get_fio_vector_potential_p(params,F,Y_R,Y_PHI,Y_Z, &
         PSIp,flagCon,hint)

     call GCEoM1_fio_p(tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
         B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp, &
         ne,ni,nimp,Te,Zeff,flagCon,flagCol,hint)

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k5_R,k5_PHI,k5_Z,k5_PLL,k5_MU)
    do cc=1_idef,pchunk
       k5_R(cc)=dt*RHS_R(cc)    
       k5_PHI(cc)=dt*RHS_PHI(cc)    
       k5_Z(cc)=dt*RHS_Z(cc)    
       k5_PLL(cc)=dt*RHS_PLL(cc)
       k5_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a51*k1_R(cc)+a52*k2_R(cc)+a53*k3_R(cc)+ &
            a54*k4_R(cc)+a55*k5_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a51*k1_PHI(cc)+a52*k2_PHI(cc)+ &
            a53*k3_PHI(cc)+a54*k4_PHI(cc)+a55*k5_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a51*k1_Z(cc)+a52*k2_Z(cc)+a53*k3_Z(cc)+ &
            a54*k4_Z(cc)+a55*k5_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a51*k1_PLL(cc)+a52*k2_PLL(cc)+ &
            a53*k3_PLL(cc)+a54*k4_PLL(cc)+a55*k5_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a51*k1_MU(cc)+a52*k2_MU(cc)+ &
            a53*k3_MU(cc)+a54*k4_MU(cc)+a55*k5_MU(cc)
    end do
    !$OMP END SIMD


    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call get_fio_GCmagnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       B_R,B_PHI,B_Z,gradB_R,gradB_PHI,gradB_Z, &
       curlb_R,curlb_PHI,curlb_Z,flagCon,hint)
    if (F%FIO_E .ge. 0) then
       call get_fio_GCelectric_fields_p(params,F, &
            Y_R,Y_PHI,Y_Z,E_R,E_PHI,E_Z,flagCon,hint)
    end if
    call get_fio_vector_potential_p(params,F,Y_R,Y_PHI,Y_Z, &
         PSIp,flagCon,hint)
    call GCEoM1_fio_p(tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
         B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp, &
         ne,ni,nimp,Te,Zeff,flagCon,flagCol,hint)         

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k6_R,k6_PHI,k6_Z,k6_PLL,k6_MU)
    do cc=1_idef,pchunk
       k6_R(cc)=dt*RHS_R(cc)    
       k6_PHI(cc)=dt*RHS_PHI(cc)    
       k6_Z(cc)=dt*RHS_Z(cc)    
       k6_PLL(cc)=dt*RHS_PLL(cc)
       k6_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+b1*k1_R(cc)+b2*k2_R(cc)+ &
            b3*k3_R(cc)+b4*k4_R(cc)+b5*k5_R(cc)+b6*k6_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+b1*k1_PHI(cc)+b2*k2_PHI(cc)+ &
            b3*k3_PHI(cc)+b4*k4_PHI(cc)+b5*k5_PHI(cc)+b6*k6_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+b1*k1_Z(cc)+b2*k2_Z(cc)+ &
            b3*k3_Z(cc)+b4*k4_Z(cc)+b5*k5_Z(cc)+b6*k6_Z(cc)
       V_PLL(cc)=V0_PLL(cc)+b1*k1_PLL(cc)+b2*k2_PLL(cc)+ &
            b3*k3_PLL(cc)+b4*k4_PLL(cc)+b5*k5_PLL(cc)+b6*k6_PLL(cc)
       V_MU(cc)=V0_MU(cc)+b1*k1_MU(cc)+b2*k2_MU(cc)+ &
            b3*k3_MU(cc)+b4*k4_MU(cc)+b5*k5_MU(cc)+b6*k6_MU(cc)
    end do
    !$OMP END SIMD
    

    !$OMP SIMD
    !    !$OMP& aligned(Y_R,Y_PHI,Y_Z,V_PLL,V_MU,Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU)
    do cc=1_idef,pchunk

       if ((flagCon(cc).eq.0_is).or.(flagCol(cc).eq.0_is)) then
          Y_R(cc)=Y0_R(cc)
          Y_PHI(cc)=Y0_PHI(cc)
          Y_Z(cc)=Y0_Z(cc)
          V_PLL(cc)=V0_PLL(cc)
          V_MU(cc)=V0_MU(cc)
       end if

    end do
    !$OMP END SIMD

    call get_fio_GCmagnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       B_R,B_PHI,B_Z,gradB_R,gradB_PHI,gradB_Z, &
       curlb_R,curlb_PHI,curlb_Z,flagCon,hint)
    if (F%FIO_E .ge. 0) then
       call get_fio_GCelectric_fields_p(params,F, &
            Y_R,Y_PHI,Y_Z,E_R,E_PHI,E_Z,flagCon,hint)
    end if
    call get_fio_vector_potential_p(params,F,Y_R,Y_PHI,Y_Z, &
         PSIp,flagCon,hint)

    call GCEoM1_fio_p(tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
         B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp, &
         ne,ni,nimp,Te,Zeff,flagCon,flagCol,hint) 

    !$OMP SIMD
    do cc=1_idef,pchunk
       vars%RHS(pp-1+cc,1)=RHS_R(cc)
       vars%RHS(pp-1+cc,2)=RHS_PHI(cc)
       vars%RHS(pp-1+cc,3)=RHS_Z(cc)
       vars%RHS(pp-1+cc,4)=RHS_PLL(cc)
       vars%RHS(pp-1+cc,5)=RHS_MU(cc)
    end do
    !$OMP END SIMD

    !write(6,*) E_PHI
    
    if (params%collisions) then       

       call include_CoulombCollisions_GCfio_p(tt,params,Y_R,Y_PHI,Y_Z, &
            V_PLL,V_MU,m_cache,flagCon,flagCol,F,P,E_PHI,ne,ni,Te,Zeff&
            &,nimp,PSIp,hint)

    end if


  end subroutine advance_GCinterp_fio_vars
#endif

#ifdef PSPLINE
  subroutine advance_GCinterp_psiwE_vars(vars,pchunk,pp,tt,params,Y_R,Y_PHI,Y_Z, &
       V_PLL,V_MU,q_cache,m_cache,flagCon,flagCol,F,P,B_R,B_PHI,B_Z,E_PHI,PSIp, &
       curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,ne)
    !! @note Subroutine to advance GC variables \(({\bf X},p_\parallel)\)
    !! @endnote
    !! Comment this section further with evolution equations, numerical
    !! methods, and descriptions of both.
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(PARTICLES), INTENT(INOUT)     :: vars
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    REAL(rp)                                      :: dt
    !! Time step used in the leapfrog step (\(\Delta t\)).

    INTEGER                                                    :: cc
    INTEGER, INTENT(IN)                 :: pchunk
    !! Chunk iterator.
    INTEGER(ip),intent(in)                                      :: tt
    !! time iterator.
    INTEGER,intent(in)                                  :: pp


    REAL(rp),DIMENSION(pchunk)               :: Bmag
    REAL(rp)              :: a1 = 1./5._rp
    REAL(rp) :: a21 = 3./40._rp,a22=9./40._rp
    REAL(rp) :: a31 = 3./10._rp,a32=-9./10._rp,a33=6./5._rp
    REAL(rp) :: a41 = -11./54._rp,a42=5./2._rp,a43=-70./27._rp,a44=35./27._rp
    REAL(rp) :: a51 = 1631./55296._rp,a52=175./512._rp,a53=575./13824._rp,a54=44275./110592._rp,a55=253./4096._rp
    REAL(rp) :: b1=37./378._rp,b2=0._rp,b3=250./621._rp,b4=125./594._rp,b5=0._rp,b6=512./1771._rp

    REAL(rp),DIMENSION(pchunk) :: k1_R,k1_PHI,k1_Z,k1_PLL,k1_MU
    REAL(rp),DIMENSION(pchunk) :: k2_R,k2_PHI,k2_Z,k2_PLL,k2_MU
    REAL(rp),DIMENSION(pchunk) :: k3_R,k3_PHI,k3_Z,k3_PLL,k3_MU
    REAL(rp),DIMENSION(pchunk) :: k4_R,k4_PHI,k4_Z,k4_PLL,k4_MU
    REAL(rp),DIMENSION(pchunk) :: k5_R,k5_PHI,k5_Z,k5_PLL,k5_MU
    REAL(rp),DIMENSION(pchunk) :: k6_R,k6_PHI,k6_Z,k6_PLL,k6_MU
    REAL(rp),DIMENSION(pchunk) :: Y0_R,Y0_PHI,Y0_Z
    REAL(rp),DIMENSION(pchunk),INTENT(INOUT) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(pchunk),INTENT(OUT) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(pchunk) :: E_R,E_Z
    REAL(rp),DIMENSION(pchunk),INTENT(OUT) :: E_PHI
    REAL(rp),DIMENSION(pchunk),INTENT(OUT) :: PSIp
    REAL(rp),DIMENSION(pchunk),INTENT(OUT) :: ne
    REAL(rp),DIMENSION(pchunk),INTENT(OUT) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp),DIMENSION(pchunk),INTENT(OUT) :: gradB_R,gradB_PHI,gradB_Z
    REAL(rp),DIMENSION(pchunk),INTENT(INOUT) :: V_PLL,V_MU
    REAL(rp),DIMENSION(pchunk) :: RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU
    REAL(rp),DIMENSION(pchunk) :: V0_PLL,V0_MU
    REAL(rp),DIMENSION(pchunk) :: Te,Zeff

    INTEGER(is),DIMENSION(pchunk),intent(INOUT) :: flagCon,flagCol
    REAL(rp),intent(IN)  :: q_cache,m_cache

    dt=params%dt

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU)
    do cc=1_idef,pchunk

       Y0_R(cc)=Y_R(cc)
       Y0_PHI(cc)=Y_PHI(cc)
       Y0_Z(cc)=Y_Z(cc)
       V0_PLL(cc)=V_PLL(cc)
       V0_MU(cc)=V_MU(cc)



    end do
    !$OMP END SIMD

    !    write(output_unit_write,*) 'R0',Y_R(1)
    !    write(output_unit_write,*) 'PHI0',Y_PHI(1)
    !    write(output_unit_write,*) 'Z0',Y_Z(1)
    !    write(output_unit_write,*) 'PPLL0',V_PLL(1)
    !    write(output_unit_write,*) 'MU0',V_MU(1)


    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfieldswE_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)


    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !    write(output_unit_write,*) 'R0',Y_R(1)
    !    write(output_unit_write,*) 'PHI0',Y_PHI(1)
    !    write(output_unit_write,*) 'Z0',Y_Z(1)
    !    write(output_unit_write,*) 'PPLL0',V_PLL(1)
    !    write(output_unit_write,*) 'MU0',V_MU(1)

    !    write(output_unit_write,*) 'BR',B_R(1)
    !    write(output_unit_write,*) 'BPHI',B_PHI(1)
    !    write(output_unit_write,*) 'BZ',B_Z(1)

    !    write(output_unit_write,*) 'gradBR',gradB_R(1)
    !    write(output_unit_write,*) 'gradBPHI',gradB_PHI(1)
    !    write(output_unit_write,*) 'gradBZ',gradB_Z(1)

    !    write(output_unit_write,*) 'curlBR',curlB_R(1)
    !    write(output_unit_write,*) 'curlBPHI',curlB_PHI(1)
    !    write(output_unit_write,*) 'curlBZ',curlB_Z(1)

    !    write(output_unit_write,*) 'dt',params%dt
    !    write(output_unit_write,*) 'RHS_R',RHS_R(1)
    !    write(output_unit_write,*) 'RHS_PHI',RHS_PHI(1)
    !    write(output_unit_write,*) 'RHS_Z',RHS_Z(1)
    !    write(output_unit_write,*) 'RHS_PLL',RHS_PLL(1)
    !    write(output_unit_write,*) 'RHS_MU',RHS_MU(1)

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k1_R,k1_PHI,k1_Z,k1_PLL,k1_MU)
    do cc=1_idef,pchunk
       k1_R(cc)=dt*RHS_R(cc)              
       k1_PHI(cc)=dt*RHS_PHI(cc)    
       k1_Z(cc)=dt*RHS_Z(cc)    
       k1_PLL(cc)=dt*RHS_PLL(cc)
       k1_MU(cc)=dt*RHS_MU(cc)    

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a1*k1_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a1*k1_MU(cc)

    end do
    !$OMP END SIMD

    !    write(output_unit_write,*) 'R1',Y_R(1)
    !    write(output_unit_write,*) 'PHI1',Y_PHI(1)
    !!    write(output_unit_write,*) 'Z1',Y_Z(1)
    !   write(output_unit_write,*) 'PPLL1',V_PLL(1)
    !   write(output_unit_write,*) 'MU1',V_MU(1)

    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfieldswE_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)
   

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k2_R,k2_PHI,k2_Z,k2_PLL,k2_MU)
    do cc=1_idef,pchunk
       k2_R(cc)=dt*RHS_R(cc)    
       k2_PHI(cc)=dt*RHS_PHI (cc)   
       k2_Z(cc)=dt*RHS_Z(cc)   
       k2_PLL(cc)=dt*RHS_PLL(cc)
       k2_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a21*k1_PLL(cc)+a22*k2_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a21*k1_MU(cc)+a22*k2_MU(cc)


    end do
    !$OMP END SIMD

    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfieldswE_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)


    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k3_R,k3_PHI,k3_Z,k3_PLL,k3_MU)
    do cc=1_idef,pchunk
       k3_R(cc)=dt*RHS_R(cc)   
       k3_PHI(cc)=dt*RHS_PHI(cc)    
       k3_Z(cc)=dt*RHS_Z(cc)    
       k3_PLL(cc)=dt*RHS_PLL(cc)
       k3_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)

       V_PLL(cc)=V0_PLL(cc)   +a31*k1_PLL(cc)+a32*k2_PLL(cc)+a33*k3_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a31*k1_MU(cc)+a32*k2_MU(cc)+a33*k3_MU(cc)


    end do
    !$OMP END SIMD

    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfieldswE_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)     

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k4_R,k4_PHI,k4_Z,k4_PLL,k4_MU)
    do cc=1_idef,pchunk
       k4_R(cc)=dt*RHS_R(cc)   
       k4_PHI(cc)=dt*RHS_PHI(cc)    
       k4_Z(cc)=dt*RHS_Z(cc)    
       k4_PLL(cc)=dt*RHS_PLL(cc)
       k4_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a41*k1_R(cc)+a42*k2_R(cc)+a43*k3_R(cc)+ &
            a44*k4_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a41*k1_PHI(cc)+a42*k2_PHI(cc)+ &
            a43*k3_PHI(cc)+a44*k4_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a41*k1_Z(cc)+a42*k2_Z(cc)+a43*k3_Z(cc)+ &
            a44*k4_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a41*k1_PLL(cc)+a42*k2_PLL(cc)+ &
            a43*k3_PLL(cc)+a44*k4_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a41*k1_MU(cc)+a42*k2_MU(cc)+ &
            a43*k3_MU(cc)+a44*k4_MU(cc)


    end do
    !$OMP END SIMD


    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfieldswE_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)   

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k5_R,k5_PHI,k5_Z,k5_PLL,k5_MU)
    do cc=1_idef,pchunk
       k5_R(cc)=dt*RHS_R(cc)    
       k5_PHI(cc)=dt*RHS_PHI(cc)    
       k5_Z(cc)=dt*RHS_Z(cc)    
       k5_PLL(cc)=dt*RHS_PLL(cc)
       k5_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a51*k1_R(cc)+a52*k2_R(cc)+a53*k3_R(cc)+ &
            a54*k4_R(cc)+a55*k5_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a51*k1_PHI(cc)+a52*k2_PHI(cc)+ &
            a53*k3_PHI(cc)+a54*k4_PHI(cc)+a55*k5_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a51*k1_Z(cc)+a52*k2_Z(cc)+a53*k3_Z(cc)+ &
            a54*k4_Z(cc)+a55*k5_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a51*k1_PLL(cc)+a52*k2_PLL(cc)+ &
            a53*k3_PLL(cc)+a54*k4_PLL(cc)+a55*k5_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a51*k1_MU(cc)+a52*k2_MU(cc)+ &
            a53*k3_MU(cc)+a54*k4_MU(cc)+a55*k5_MU(cc)


    end do
    !$OMP END SIMD

    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfieldswE_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)         

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k6_R,k6_PHI,k6_Z,k6_PLL,k6_MU)
    do cc=1_idef,pchunk
       k6_R(cc)=dt*RHS_R(cc)    
       k6_PHI(cc)=dt*RHS_PHI(cc)    
       k6_Z(cc)=dt*RHS_Z(cc)    
       k6_PLL(cc)=dt*RHS_PLL(cc)
       k6_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+b1*k1_R(cc)+b2*k2_R(cc)+ &
            b3*k3_R(cc)+b4*k4_R(cc)+b5*k5_R(cc)+b6*k6_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+b1*k1_PHI(cc)+b2*k2_PHI(cc)+ &
            b3*k3_PHI(cc)+b4*k4_PHI(cc)+b5*k5_PHI(cc)+b6*k6_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+b1*k1_Z(cc)+b2*k2_Z(cc)+ &
            b3*k3_Z(cc)+b4*k4_Z(cc)+b5*k5_Z(cc)+b6*k6_Z(cc)
       V_PLL(cc)=V0_PLL(cc)+b1*k1_PLL(cc)+b2*k2_PLL(cc)+ &
            b3*k3_PLL(cc)+b4*k4_PLL(cc)+b5*k5_PLL(cc)+b6*k6_PLL(cc)
       V_MU(cc)=V0_MU(cc)+b1*k1_MU(cc)+b2*k2_MU(cc)+ &
            b3*k3_MU(cc)+b4*k4_MU(cc)+b5*k5_MU(cc)+b6*k6_MU(cc)


    end do
    !$OMP END SIMD

    !$OMP SIMD
    !    !$OMP& aligned(Y_R,Y_PHI,Y_Z,V_PLL,V_MU,Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU)
    do cc=1_idef,pchunk

       if ((flagCon(cc).eq.0_is).or.(flagCol(cc).eq.0_is)) then
          Y_R(cc)=Y0_R(cc)
          Y_PHI(cc)=Y0_PHI(cc)
          Y_Z(cc)=Y0_Z(cc)
          V_PLL(cc)=V0_PLL(cc)
          V_MU(cc)=V0_MU(cc)
       end if

    end do
    !$OMP END SIMD

    call calculate_GCfieldswE_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !$OMP SIMD
    do cc=1_idef,pchunk
       vars%RHS(pp-1+cc,1)=RHS_R(cc)
       vars%RHS(pp-1+cc,2)=RHS_PHI(cc)
       vars%RHS(pp-1+cc,3)=RHS_Z(cc)
       vars%RHS(pp-1+cc,4)=RHS_PLL(cc)
       vars%RHS(pp-1+cc,5)=RHS_MU(cc)
    end do
    !$OMP END SIMD       


  end subroutine advance_GCinterp_psiwE_vars

  subroutine advance_GCinterp_psi2x1t_vars(vars,pp,tt,params,Y_R,Y_PHI,Y_Z, &
       V_PLL,V_MU,q_cache,m_cache,flagCon,flagCol,F,P,B_R,B_PHI,B_Z,E_PHI,PSIp, &
       curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,ne)
    !! @note Subroutine to advance GC variables \(({\bf X},p_\parallel)\)
    !! @endnote
    !! Comment this section further with evolution equations, numerical
    !! methods, and descriptions of both.
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(PARTICLES), INTENT(INOUT)     :: vars
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    REAL(rp)                                      :: dt,time
    !! Time step used in the leapfrog step (\(\Delta t\)).

    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    INTEGER(ip),intent(in)                                      :: tt
    !! time iterator.
    INTEGER,intent(in)                                  :: pp


    REAL(rp),DIMENSION(params%pchunk)               :: Bmag
    REAL(rp)              :: a1 = 1./5._rp
    REAL(rp) :: a21 = 3./40._rp,a22=9./40._rp
    REAL(rp) :: a31 = 3./10._rp,a32=-9./10._rp,a33=6./5._rp
    REAL(rp) :: a41 = -11./54._rp,a42=5./2._rp,a43=-70./27._rp,a44=35./27._rp
    REAL(rp) :: a51 = 1631./55296._rp,a52=175./512._rp,a53=575./13824._rp,a54=44275./110592._rp,a55=253./4096._rp
    REAL(rp) :: b1=37./378._rp,b2=0._rp,b3=250./621._rp,b4=125./594._rp,b5=0._rp,b6=512./1771._rp

    REAL(rp),DIMENSION(params%pchunk) :: k1_R,k1_PHI,k1_Z,k1_PLL,k1_MU
    REAL(rp),DIMENSION(params%pchunk) :: k2_R,k2_PHI,k2_Z,k2_PLL,k2_MU
    REAL(rp),DIMENSION(params%pchunk) :: k3_R,k3_PHI,k3_Z,k3_PLL,k3_MU
    REAL(rp),DIMENSION(params%pchunk) :: k4_R,k4_PHI,k4_Z,k4_PLL,k4_MU
    REAL(rp),DIMENSION(params%pchunk) :: k5_R,k5_PHI,k5_Z,k5_PLL,k5_MU
    REAL(rp),DIMENSION(params%pchunk) :: k6_R,k6_PHI,k6_Z,k6_PLL,k6_MU
    REAL(rp),DIMENSION(params%pchunk) :: Y0_R,Y0_PHI,Y0_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_R,E_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: E_PHI
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: PSIp
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: ne
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: gradB_R,gradB_PHI,gradB_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: V_PLL,V_MU
    REAL(rp),DIMENSION(params%pchunk) :: RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU
    REAL(rp),DIMENSION(params%pchunk) :: V0_PLL,V0_MU
    REAL(rp),DIMENSION(params%pchunk) :: Te,Zeff

    INTEGER(is),DIMENSION(params%pchunk),intent(INOUT) :: flagCon,flagCol
    REAL(rp),intent(IN)  :: q_cache,m_cache

    dt=params%dt
    pchunk=params%pchunk
    time=params%init_time+(params%it-1+tt)*params%dt

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU)
    do cc=1_idef,pchunk

       Y0_R(cc)=Y_R(cc)
       Y0_PHI(cc)=Y_PHI(cc)
       Y0_Z(cc)=Y_Z(cc)
       V0_PLL(cc)=V_PLL(cc)
       V0_MU(cc)=V_MU(cc)



    end do
    !$OMP END SIMD

    !    write(output_unit_write,*) 'R0',Y_R(1)
    !    write(output_unit_write,*) 'PHI0',Y_PHI(1)
    !    write(output_unit_write,*) 'Z0',Y_Z(1)
    !    write(output_unit_write,*) 'PPLL0',V_PLL(1)
    !    write(output_unit_write,*) 'MU0',V_MU(1)


    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfields_2x1t_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp,time)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)


    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !    write(output_unit_write,*) 'R0',Y_R(1)
    !    write(output_unit_write,*) 'PHI0',Y_PHI(1)
    !    write(output_unit_write,*) 'Z0',Y_Z(1)
    !    write(output_unit_write,*) 'PPLL0',V_PLL(1)
    !    write(output_unit_write,*) 'MU0',V_MU(1)

    !    write(output_unit_write,*) 'BR',B_R(1)
    !    write(output_unit_write,*) 'BPHI',B_PHI(1)
    !    write(output_unit_write,*) 'BZ',B_Z(1)

    !    write(output_unit_write,*) 'gradBR',gradB_R(1)
    !    write(output_unit_write,*) 'gradBPHI',gradB_PHI(1)
    !    write(output_unit_write,*) 'gradBZ',gradB_Z(1)

    !    write(output_unit_write,*) 'curlBR',curlB_R(1)
    !    write(output_unit_write,*) 'curlBPHI',curlB_PHI(1)
    !    write(output_unit_write,*) 'curlBZ',curlB_Z(1)

    !    write(output_unit_write,*) 'dt',params%dt
    !    write(output_unit_write,*) 'RHS_R',RHS_R(1)
    !    write(output_unit_write,*) 'RHS_PHI',RHS_PHI(1)
    !    write(output_unit_write,*) 'RHS_Z',RHS_Z(1)
    !    write(output_unit_write,*) 'RHS_PLL',RHS_PLL(1)
    !    write(output_unit_write,*) 'RHS_MU',RHS_MU(1)

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k1_R,k1_PHI,k1_Z,k1_PLL,k1_MU)
    do cc=1_idef,pchunk
       k1_R(cc)=dt*RHS_R(cc)              
       k1_PHI(cc)=dt*RHS_PHI(cc)    
       k1_Z(cc)=dt*RHS_Z(cc)    
       k1_PLL(cc)=dt*RHS_PLL(cc)
       k1_MU(cc)=dt*RHS_MU(cc)    

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a1*k1_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a1*k1_MU(cc)

    end do
    !$OMP END SIMD

    !    write(output_unit_write,*) 'R1',Y_R(1)
    !    write(output_unit_write,*) 'PHI1',Y_PHI(1)
    !!    write(output_unit_write,*) 'Z1',Y_Z(1)
    !   write(output_unit_write,*) 'PPLL1',V_PLL(1)
    !   write(output_unit_write,*) 'MU1',V_MU(1)

    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfields_2x1t_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp,time)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)





    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k2_R,k2_PHI,k2_Z,k2_PLL,k2_MU)
    do cc=1_idef,pchunk
       k2_R(cc)=dt*RHS_R(cc)    
       k2_PHI(cc)=dt*RHS_PHI (cc)   
       k2_Z(cc)=dt*RHS_Z(cc)   
       k2_PLL(cc)=dt*RHS_PLL(cc)
       k2_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a21*k1_PLL(cc)+a22*k2_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a21*k1_MU(cc)+a22*k2_MU(cc)


    end do
    !$OMP END SIMD

    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfields_2x1t_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp,time)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)





    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k3_R,k3_PHI,k3_Z,k3_PLL,k3_MU)
    do cc=1_idef,pchunk
       k3_R(cc)=dt*RHS_R(cc)   
       k3_PHI(cc)=dt*RHS_PHI(cc)    
       k3_Z(cc)=dt*RHS_Z(cc)    
       k3_PLL(cc)=dt*RHS_PLL(cc)
       k3_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)

       V_PLL(cc)=V0_PLL(cc)   +a31*k1_PLL(cc)+a32*k2_PLL(cc)+a33*k3_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a31*k1_MU(cc)+a32*k2_MU(cc)+a33*k3_MU(cc)


    end do
    !$OMP END SIMD

    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfields_2x1t_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp,time)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)





    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)     

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k4_R,k4_PHI,k4_Z,k4_PLL,k4_MU)
    do cc=1_idef,pchunk
       k4_R(cc)=dt*RHS_R(cc)   
       k4_PHI(cc)=dt*RHS_PHI(cc)    
       k4_Z(cc)=dt*RHS_Z(cc)    
       k4_PLL(cc)=dt*RHS_PLL(cc)
       k4_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a41*k1_R(cc)+a42*k2_R(cc)+a43*k3_R(cc)+ &
            a44*k4_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a41*k1_PHI(cc)+a42*k2_PHI(cc)+ &
            a43*k3_PHI(cc)+a44*k4_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a41*k1_Z(cc)+a42*k2_Z(cc)+a43*k3_Z(cc)+ &
            a44*k4_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a41*k1_PLL(cc)+a42*k2_PLL(cc)+ &
            a43*k3_PLL(cc)+a44*k4_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a41*k1_MU(cc)+a42*k2_MU(cc)+ &
            a43*k3_MU(cc)+a44*k4_MU(cc)


    end do
    !$OMP END SIMD


    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfields_2x1t_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp,time)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)





    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)   

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k5_R,k5_PHI,k5_Z,k5_PLL,k5_MU)
    do cc=1_idef,pchunk
       k5_R(cc)=dt*RHS_R(cc)    
       k5_PHI(cc)=dt*RHS_PHI(cc)    
       k5_Z(cc)=dt*RHS_Z(cc)    
       k5_PLL(cc)=dt*RHS_PLL(cc)
       k5_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a51*k1_R(cc)+a52*k2_R(cc)+a53*k3_R(cc)+ &
            a54*k4_R(cc)+a55*k5_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a51*k1_PHI(cc)+a52*k2_PHI(cc)+ &
            a53*k3_PHI(cc)+a54*k4_PHI(cc)+a55*k5_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a51*k1_Z(cc)+a52*k2_Z(cc)+a53*k3_Z(cc)+ &
            a54*k4_Z(cc)+a55*k5_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a51*k1_PLL(cc)+a52*k2_PLL(cc)+ &
            a53*k3_PLL(cc)+a54*k4_PLL(cc)+a55*k5_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a51*k1_MU(cc)+a52*k2_MU(cc)+ &
            a53*k3_MU(cc)+a54*k4_MU(cc)+a55*k5_MU(cc)


    end do
    !$OMP END SIMD

    !    call interp_fields_p(F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
    call calculate_GCfields_2x1t_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp,time)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)


    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)         

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k6_R,k6_PHI,k6_Z,k6_PLL,k6_MU)
    do cc=1_idef,pchunk
       k6_R(cc)=dt*RHS_R(cc)    
       k6_PHI(cc)=dt*RHS_PHI(cc)    
       k6_Z(cc)=dt*RHS_Z(cc)    
       k6_PLL(cc)=dt*RHS_PLL(cc)
       k6_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+b1*k1_R(cc)+b2*k2_R(cc)+ &
            b3*k3_R(cc)+b4*k4_R(cc)+b5*k5_R(cc)+b6*k6_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+b1*k1_PHI(cc)+b2*k2_PHI(cc)+ &
            b3*k3_PHI(cc)+b4*k4_PHI(cc)+b5*k5_PHI(cc)+b6*k6_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+b1*k1_Z(cc)+b2*k2_Z(cc)+ &
            b3*k3_Z(cc)+b4*k4_Z(cc)+b5*k5_Z(cc)+b6*k6_Z(cc)
       V_PLL(cc)=V0_PLL(cc)+b1*k1_PLL(cc)+b2*k2_PLL(cc)+ &
            b3*k3_PLL(cc)+b4*k4_PLL(cc)+b5*k5_PLL(cc)+b6*k6_PLL(cc)
       V_MU(cc)=V0_MU(cc)+b1*k1_MU(cc)+b2*k2_MU(cc)+ &
            b3*k3_MU(cc)+b4*k4_MU(cc)+b5*k5_MU(cc)+b6*k6_MU(cc)


    end do
    !$OMP END SIMD

    !$OMP SIMD
    !    !$OMP& aligned(Y_R,Y_PHI,Y_Z,V_PLL,V_MU,Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU)
    do cc=1_idef,pchunk

       if ((flagCon(cc).eq.0_is).or.(flagCol(cc).eq.0_is)) then
          Y_R(cc)=Y0_R(cc)
          Y_PHI(cc)=Y0_PHI(cc)
          Y_Z(cc)=Y0_Z(cc)
          V_PLL(cc)=V0_PLL(cc)
          V_MU(cc)=V0_MU(cc)
       end if

    end do
    !$OMP END SIMD

    call calculate_GCfields_2x1t_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp,time)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)



    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !$OMP SIMD
    do cc=1_idef,pchunk
       vars%RHS(pp-1+cc,1)=RHS_R(cc)
       vars%RHS(pp-1+cc,2)=RHS_PHI(cc)
       vars%RHS(pp-1+cc,3)=RHS_Z(cc)
       vars%RHS(pp-1+cc,4)=RHS_PLL(cc)
       vars%RHS(pp-1+cc,5)=RHS_MU(cc)
    end do
    !$OMP END SIMD       

    if (params%collisions) then       

       call include_CoulombCollisions_GC_p(tt,params,Y_R,Y_PHI,Y_Z, &
            V_PLL,V_MU,m_cache,flagCon,flagCol,F,P,E_PHI,ne,PSIp)

    end if


  end subroutine advance_GCinterp_psi2x1t_vars

  subroutine advance_GCinterp_B2D_vars(vars,pp,tt,params,Y_R,Y_PHI,Y_Z, &
       V_PLL,V_MU,q_cache,m_cache,flagCon,flagCol,F,P,B_R,B_PHI,B_Z,E_PHI, &
       curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,PSIp)
    !! @note Subroutine to advance GC variables \(({\bf X},p_\parallel)\)
    !! @endnote
    !! Comment this section further with evolution equations, numerical
    !! methods, and descriptions of both.
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(PARTICLES), INTENT(INOUT)     :: vars
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    REAL(rp)                                      :: dt
    !! Time step used in the leapfrog step (\(\Delta t\)).

    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    INTEGER(ip),intent(in)                                      :: tt
    !! time iterator.
    INTEGER,intent(in)                                  :: pp


    REAL(rp),DIMENSION(params%pchunk)               :: Bmag
    REAL(rp)              :: a1 = 1./5._rp
    REAL(rp) :: a21 = 3./40._rp,a22=9./40._rp
    REAL(rp) :: a31 = 3./10._rp,a32=-9./10._rp,a33=6./5._rp
    REAL(rp) :: a41 = -11./54._rp,a42=5./2._rp,a43=-70./27._rp,a44=35./27._rp
    REAL(rp) :: a51 = 1631./55296._rp,a52=175./512._rp,a53=575./13824._rp,a54=44275./110592._rp,a55=253./4096._rp
    REAL(rp) :: b1=37./378._rp,b2=0._rp,b3=250./621._rp,b4=125./594._rp,b5=0._rp,b6=512./1771._rp

    REAL(rp),DIMENSION(params%pchunk) :: k1_R,k1_PHI,k1_Z,k1_PLL,k1_MU
    REAL(rp),DIMENSION(params%pchunk) :: k2_R,k2_PHI,k2_Z,k2_PLL,k2_MU
    REAL(rp),DIMENSION(params%pchunk) :: k3_R,k3_PHI,k3_Z,k3_PLL,k3_MU
    REAL(rp),DIMENSION(params%pchunk) :: k4_R,k4_PHI,k4_Z,k4_PLL,k4_MU
    REAL(rp),DIMENSION(params%pchunk) :: k5_R,k5_PHI,k5_Z,k5_PLL,k5_MU
    REAL(rp),DIMENSION(params%pchunk) :: k6_R,k6_PHI,k6_Z,k6_PLL,k6_MU
    REAL(rp),DIMENSION(params%pchunk) :: Y0_R,Y0_PHI,Y0_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: Y_R,Y_PHI,Y_Z,PSIp
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_R,E_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: E_PHI
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: gradB_R,gradB_PHI,gradB_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: V_PLL,V_MU
    REAL(rp),DIMENSION(params%pchunk) :: RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU
    REAL(rp),DIMENSION(params%pchunk) :: V0_PLL,V0_MU
    REAL(rp),DIMENSION(params%pchunk) :: ne,Te,Zeff

    INTEGER(is),DIMENSION(params%pchunk),intent(INOUT) :: flagCon,flagCol
    REAL(rp),intent(IN)  :: q_cache,m_cache

    dt=params%dt
    pchunk=params%pchunk

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU)
    do cc=1_idef,pchunk

       Y0_R(cc)=Y_R(cc)
       Y0_PHI(cc)=Y_PHI(cc)
       Y0_Z(cc)=Y_Z(cc)
       V0_PLL(cc)=V_PLL(cc)
       V0_MU(cc)=V_MU(cc)
    end do
    !$OMP END SIMD

    call interp_fields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k1_R,k1_PHI,k1_Z,k1_PLL,k1_MU)
    do cc=1_idef,pchunk
       k1_R(cc)=dt*RHS_R(cc)              
       k1_PHI(cc)=dt*RHS_PHI(cc)    
       k1_Z(cc)=dt*RHS_Z(cc)    
       k1_PLL(cc)=dt*RHS_PLL(cc)
       k1_MU(cc)=dt*RHS_MU(cc)    

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a1*k1_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a1*k1_MU(cc)
    end do
    !$OMP END SIMD

    call interp_fields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k2_R,k2_PHI,k2_Z,k2_PLL,k2_MU)
    do cc=1_idef,pchunk
       k2_R(cc)=dt*RHS_R(cc)    
       k2_PHI(cc)=dt*RHS_PHI (cc)   
       k2_Z(cc)=dt*RHS_Z(cc)   
       k2_PLL(cc)=dt*RHS_PLL(cc)
       k2_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a21*k1_PLL(cc)+a22*k2_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a21*k1_MU(cc)+a22*k2_MU(cc)
    end do
    !$OMP END SIMD

    call interp_fields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k3_R,k3_PHI,k3_Z,k3_PLL,k3_MU)
    do cc=1_idef,pchunk
       k3_R(cc)=dt*RHS_R(cc)   
       k3_PHI(cc)=dt*RHS_PHI(cc)    
       k3_Z(cc)=dt*RHS_Z(cc)    
       k3_PLL(cc)=dt*RHS_PLL(cc)
       k3_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a31*k1_PLL(cc)+a32*k2_PLL(cc)+a33*k3_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a31*k1_MU(cc)+a32*k2_MU(cc)+a33*k3_MU(cc)
    end do
    !$OMP END SIMD

    call interp_fields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)     

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k4_R,k4_PHI,k4_Z,k4_PLL,k4_MU)
    do cc=1_idef,pchunk
       k4_R(cc)=dt*RHS_R(cc)   
       k4_PHI(cc)=dt*RHS_PHI(cc)    
       k4_Z(cc)=dt*RHS_Z(cc)    
       k4_PLL(cc)=dt*RHS_PLL(cc)
       k4_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a41*k1_R(cc)+a42*k2_R(cc)+a43*k3_R(cc)+ &
            a44*k4_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a41*k1_PHI(cc)+a42*k2_PHI(cc)+ &
            a43*k3_PHI(cc)+a44*k4_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a41*k1_Z(cc)+a42*k2_Z(cc)+a43*k3_Z(cc)+ &
            a44*k4_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a41*k1_PLL(cc)+a42*k2_PLL(cc)+ &
            a43*k3_PLL(cc)+a44*k4_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a41*k1_MU(cc)+a42*k2_MU(cc)+ &
            a43*k3_MU(cc)+a44*k4_MU(cc)
    end do
    !$OMP END SIMD


    call interp_fields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)   

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k5_R,k5_PHI,k5_Z,k5_PLL,k5_MU)
    do cc=1_idef,pchunk
       k5_R(cc)=dt*RHS_R(cc)    
       k5_PHI(cc)=dt*RHS_PHI(cc)    
       k5_Z(cc)=dt*RHS_Z(cc)    
       k5_PLL(cc)=dt*RHS_PLL(cc)
       k5_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a51*k1_R(cc)+a52*k2_R(cc)+a53*k3_R(cc)+ &
            a54*k4_R(cc)+a55*k5_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a51*k1_PHI(cc)+a52*k2_PHI(cc)+ &
            a53*k3_PHI(cc)+a54*k4_PHI(cc)+a55*k5_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a51*k1_Z(cc)+a52*k2_Z(cc)+a53*k3_Z(cc)+ &
            a54*k4_Z(cc)+a55*k5_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a51*k1_PLL(cc)+a52*k2_PLL(cc)+ &
            a53*k3_PLL(cc)+a54*k4_PLL(cc)+a55*k5_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a51*k1_MU(cc)+a52*k2_MU(cc)+ &
            a53*k3_MU(cc)+a54*k4_MU(cc)+a55*k5_MU(cc)
    end do
    !$OMP END SIMD

    call interp_fields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)         

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k6_R,k6_PHI,k6_Z,k6_PLL,k6_MU)
    do cc=1_idef,pchunk
       k6_R(cc)=dt*RHS_R(cc)    
       k6_PHI(cc)=dt*RHS_PHI(cc)    
       k6_Z(cc)=dt*RHS_Z(cc)    
       k6_PLL(cc)=dt*RHS_PLL(cc)
       k6_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+b1*k1_R(cc)+b2*k2_R(cc)+ &
            b3*k3_R(cc)+b4*k4_R(cc)+b5*k5_R(cc)+b6*k6_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+b1*k1_PHI(cc)+b2*k2_PHI(cc)+ &
            b3*k3_PHI(cc)+b4*k4_PHI(cc)+b5*k5_PHI(cc)+b6*k6_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+b1*k1_Z(cc)+b2*k2_Z(cc)+ &
            b3*k3_Z(cc)+b4*k4_Z(cc)+b5*k5_Z(cc)+b6*k6_Z(cc)
       V_PLL(cc)=V0_PLL(cc)+b1*k1_PLL(cc)+b2*k2_PLL(cc)+ &
            b3*k3_PLL(cc)+b4*k4_PLL(cc)+b5*k5_PLL(cc)+b6*k6_PLL(cc)
       V_MU(cc)=V0_MU(cc)+b1*k1_MU(cc)+b2*k2_MU(cc)+ &
            b3*k3_MU(cc)+b4*k4_MU(cc)+b5*k5_MU(cc)+b6*k6_MU(cc)
    end do
    !$OMP END SIMD

    !$OMP SIMD
    !    !$OMP& aligned(Y_R,Y_PHI,Y_Z,V_PLL,V_MU,Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU)
    do cc=1_idef,pchunk

       if ((flagCon(cc).eq.0_is).or.(flagCol(cc).eq.0_is)) then
          Y_R(cc)=Y0_R(cc)
          Y_PHI(cc)=Y0_PHI(cc)
          Y_Z(cc)=Y0_Z(cc)
          V_PLL(cc)=V0_PLL(cc)
          V_MU(cc)=V0_MU(cc)
       end if

    end do
    !$OMP END SIMD

    call interp_fields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !$OMP SIMD
    do cc=1_idef,pchunk
       vars%RHS(pp-1+cc,1)=RHS_R(cc)
       vars%RHS(pp-1+cc,2)=RHS_PHI(cc)
       vars%RHS(pp-1+cc,3)=RHS_Z(cc)
       vars%RHS(pp-1+cc,4)=RHS_PLL(cc)
       vars%RHS(pp-1+cc,5)=RHS_MU(cc)
    end do
    !$OMP END SIMD



    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    if (params%collisions) then       

       call include_CoulombCollisions_GC_p(tt,params,Y_R,Y_PHI,Y_Z, &
            V_PLL,V_MU,m_cache,flagCon,flagCol,F,P,E_PHI,ne,PSIp)

    end if


  end subroutine advance_GCinterp_B2D_vars

  subroutine advance_GCinterp_2DBdB_vars(vars,pp,tt,params,Y_R,Y_PHI,Y_Z, &
       V_PLL,V_MU,q_cache,m_cache,flagCon,flagCol,F,P,B_R,B_PHI,B_Z,E_PHI, &
       curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,PSIp)
    !! @note Subroutine to advance GC variables \(({\bf X},p_\parallel)\)
    !! @endnote
    !! Comment this section further with evolution equations, numerical
    !! methods, and descriptions of both.
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(PARTICLES), INTENT(INOUT)     :: vars
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    REAL(rp)                                      :: dt
    !! Time step used in the leapfrog step (\(\Delta t\)).

    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    INTEGER(ip),intent(in)                                      :: tt
    !! time iterator.
    INTEGER,intent(in)                                  :: pp


    REAL(rp),DIMENSION(params%pchunk)               :: Bmag
    REAL(rp)              :: a1 = 1./5._rp
    REAL(rp) :: a21 = 3./40._rp,a22=9./40._rp
    REAL(rp) :: a31 = 3./10._rp,a32=-9./10._rp,a33=6./5._rp
    REAL(rp) :: a41 = -11./54._rp,a42=5./2._rp,a43=-70./27._rp,a44=35./27._rp
    REAL(rp) :: a51 = 1631./55296._rp,a52=175./512._rp,a53=575./13824._rp,a54=44275./110592._rp,a55=253./4096._rp
    REAL(rp) :: b1=37./378._rp,b2=0._rp,b3=250./621._rp,b4=125./594._rp,b5=0._rp,b6=512./1771._rp

    REAL(rp),DIMENSION(params%pchunk) :: k1_R,k1_PHI,k1_Z,k1_PLL,k1_MU
    REAL(rp),DIMENSION(params%pchunk) :: k2_R,k2_PHI,k2_Z,k2_PLL,k2_MU
    REAL(rp),DIMENSION(params%pchunk) :: k3_R,k3_PHI,k3_Z,k3_PLL,k3_MU
    REAL(rp),DIMENSION(params%pchunk) :: k4_R,k4_PHI,k4_Z,k4_PLL,k4_MU
    REAL(rp),DIMENSION(params%pchunk) :: k5_R,k5_PHI,k5_Z,k5_PLL,k5_MU
    REAL(rp),DIMENSION(params%pchunk) :: k6_R,k6_PHI,k6_Z,k6_PLL,k6_MU
    REAL(rp),DIMENSION(params%pchunk) :: Y0_R,Y0_PHI,Y0_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_R,E_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: E_PHI,PSIp
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: gradB_R,gradB_PHI,gradB_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: V_PLL,V_MU
    REAL(rp),DIMENSION(params%pchunk) :: RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU
    REAL(rp),DIMENSION(params%pchunk) :: V0_PLL,V0_MU
    REAL(rp),DIMENSION(params%pchunk) :: ne,Te,Zeff

    INTEGER(is),DIMENSION(params%pchunk),intent(INOUT) :: flagCon,flagCol
    REAL(rp),intent(IN)  :: q_cache,m_cache

    dt=params%dt
    pchunk=params%pchunk

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU)
    do cc=1_idef,pchunk

       Y0_R(cc)=Y_R(cc)
       Y0_PHI(cc)=Y_PHI(cc)
       Y0_Z(cc)=Y_Z(cc)
       V0_PLL(cc)=V_PLL(cc)
       V0_MU(cc)=V_MU(cc)
    end do
    !$OMP END SIMD

    call calculate_2DBdBfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !write(output_unit_write,*) 'R',Y_R(1)
    !write(output_unit_write,*) 'PHI',Y_PHI(1)
    !write(output_unit_write,*) 'Z',Y_Z(1)
    !write(output_unit_write,*) 'PPLL',V_PLL(1)
    !write(output_unit_write,*) 'MU',V_MU(1)

    !write(output_unit_write,*) 'BR',B_R(1)
    !write(output_unit_write,*) 'BPHI',B_PHI(1)
    !write(output_unit_write,*) 'BZ',B_Z(1)

    !write(output_unit_write,*) 'gradBR',gradB_R(1)
    !write(output_unit_write,*) 'gradBPHI',gradB_PHI(1)
    !write(output_unit_write,*) 'gradBZ',gradB_Z(1)

    !write(output_unit_write,*) 'curlBR',curlB_R(1)
    !write(output_unit_write,*) 'curlBPHI',curlB_PHI(1)
    !write(output_unit_write,*) 'curlBZ',curlB_Z(1)

    !write(output_unit_write,*) 'RHS_R',RHS_R(1)
    !write(output_unit_write,*) 'RHS_PHI',RHS_PHI(1)
    !write(output_unit_write,*) 'RHS_Z',RHS_Z(1)
    !write(output_unit_write,*) 'RHS_PLL',RHS_PLL(1)


    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k1_R,k1_PHI,k1_Z,k1_PLL,k1_MU)
    do cc=1_idef,pchunk
       k1_R(cc)=dt*RHS_R(cc)              
       k1_PHI(cc)=dt*RHS_PHI(cc)    
       k1_Z(cc)=dt*RHS_Z(cc)    
       k1_PLL(cc)=dt*RHS_PLL(cc)
       k1_MU(cc)=dt*RHS_MU(cc)    

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a1*k1_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a1*k1_MU(cc)
    end do
    !$OMP END SIMD

    call calculate_2DBdBfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k2_R,k2_PHI,k2_Z,k2_PLL,k2_MU)
    do cc=1_idef,pchunk
       k2_R(cc)=dt*RHS_R(cc)    
       k2_PHI(cc)=dt*RHS_PHI (cc)   
       k2_Z(cc)=dt*RHS_Z(cc)   
       k2_PLL(cc)=dt*RHS_PLL(cc)
       k2_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a21*k1_PLL(cc)+a22*k2_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a21*k1_MU(cc)+a22*k2_MU(cc)
    end do
    !$OMP END SIMD

    call calculate_2DBdBfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k3_R,k3_PHI,k3_Z,k3_PLL,k3_MU)
    do cc=1_idef,pchunk
       k3_R(cc)=dt*RHS_R(cc)   
       k3_PHI(cc)=dt*RHS_PHI(cc)    
       k3_Z(cc)=dt*RHS_Z(cc)    
       k3_PLL(cc)=dt*RHS_PLL(cc)
       k3_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a31*k1_PLL(cc)+a32*k2_PLL(cc)+a33*k3_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a31*k1_MU(cc)+a32*k2_MU(cc)+a33*k3_MU(cc)
    end do
    !$OMP END SIMD

    call calculate_2DBdBfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)     

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k4_R,k4_PHI,k4_Z,k4_PLL,k4_MU)
    do cc=1_idef,pchunk
       k4_R(cc)=dt*RHS_R(cc)   
       k4_PHI(cc)=dt*RHS_PHI(cc)    
       k4_Z(cc)=dt*RHS_Z(cc)    
       k4_PLL(cc)=dt*RHS_PLL(cc)
       k4_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a41*k1_R(cc)+a42*k2_R(cc)+a43*k3_R(cc)+ &
            a44*k4_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a41*k1_PHI(cc)+a42*k2_PHI(cc)+ &
            a43*k3_PHI(cc)+a44*k4_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a41*k1_Z(cc)+a42*k2_Z(cc)+a43*k3_Z(cc)+ &
            a44*k4_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a41*k1_PLL(cc)+a42*k2_PLL(cc)+ &
            a43*k3_PLL(cc)+a44*k4_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a41*k1_MU(cc)+a42*k2_MU(cc)+ &
            a43*k3_MU(cc)+a44*k4_MU(cc)
    end do
    !$OMP END SIMD


    call calculate_2DBdBfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)   

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k5_R,k5_PHI,k5_Z,k5_PLL,k5_MU)
    do cc=1_idef,pchunk
       k5_R(cc)=dt*RHS_R(cc)    
       k5_PHI(cc)=dt*RHS_PHI(cc)    
       k5_Z(cc)=dt*RHS_Z(cc)    
       k5_PLL(cc)=dt*RHS_PLL(cc)
       k5_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a51*k1_R(cc)+a52*k2_R(cc)+a53*k3_R(cc)+ &
            a54*k4_R(cc)+a55*k5_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a51*k1_PHI(cc)+a52*k2_PHI(cc)+ &
            a53*k3_PHI(cc)+a54*k4_PHI(cc)+a55*k5_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a51*k1_Z(cc)+a52*k2_Z(cc)+a53*k3_Z(cc)+ &
            a54*k4_Z(cc)+a55*k5_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a51*k1_PLL(cc)+a52*k2_PLL(cc)+ &
            a53*k3_PLL(cc)+a54*k4_PLL(cc)+a55*k5_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a51*k1_MU(cc)+a52*k2_MU(cc)+ &
            a53*k3_MU(cc)+a54*k4_MU(cc)+a55*k5_MU(cc)
    end do
    !$OMP END SIMD

    call calculate_2DBdBfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)         

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k6_R,k6_PHI,k6_Z,k6_PLL,k6_MU)
    do cc=1_idef,pchunk
       k6_R(cc)=dt*RHS_R(cc)    
       k6_PHI(cc)=dt*RHS_PHI(cc)    
       k6_Z(cc)=dt*RHS_Z(cc)    
       k6_PLL(cc)=dt*RHS_PLL(cc)
       k6_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+b1*k1_R(cc)+b2*k2_R(cc)+ &
            b3*k3_R(cc)+b4*k4_R(cc)+b5*k5_R(cc)+b6*k6_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+b1*k1_PHI(cc)+b2*k2_PHI(cc)+ &
            b3*k3_PHI(cc)+b4*k4_PHI(cc)+b5*k5_PHI(cc)+b6*k6_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+b1*k1_Z(cc)+b2*k2_Z(cc)+ &
            b3*k3_Z(cc)+b4*k4_Z(cc)+b5*k5_Z(cc)+b6*k6_Z(cc)
       V_PLL(cc)=V0_PLL(cc)+b1*k1_PLL(cc)+b2*k2_PLL(cc)+ &
            b3*k3_PLL(cc)+b4*k4_PLL(cc)+b5*k5_PLL(cc)+b6*k6_PLL(cc)
       V_MU(cc)=V0_MU(cc)+b1*k1_MU(cc)+b2*k2_MU(cc)+ &
            b3*k3_MU(cc)+b4*k4_MU(cc)+b5*k5_MU(cc)+b6*k6_MU(cc)
    end do
    !$OMP END SIMD

    !$OMP SIMD
    !    !$OMP& aligned(Y_R,Y_PHI,Y_Z,V_PLL,V_MU,Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU)
    do cc=1_idef,pchunk

       if ((flagCon(cc).eq.0_is).or.(flagCol(cc).eq.0_is)) then
          Y_R(cc)=Y0_R(cc)
          Y_PHI(cc)=Y0_PHI(cc)
          Y_Z(cc)=Y0_Z(cc)
          V_PLL(cc)=V0_PLL(cc)
          V_MU(cc)=V0_MU(cc)
       end if

    end do
    !$OMP END SIMD

    call calculate_2DBdBfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !$OMP SIMD
    do cc=1_idef,pchunk
       vars%RHS(pp-1+cc,1)=RHS_R(cc)
       vars%RHS(pp-1+cc,2)=RHS_PHI(cc)
       vars%RHS(pp-1+cc,3)=RHS_Z(cc)
       vars%RHS(pp-1+cc,4)=RHS_PLL(cc)
       vars%RHS(pp-1+cc,5)=RHS_MU(cc)
    end do
    !$OMP END SIMD



    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    if (params%collisions) then       

       call include_CoulombCollisions_GC_p(tt,params,Y_R,Y_PHI,Y_Z, &
            V_PLL,V_MU,m_cache,flagCon,flagCol,F,P,E_PHI,ne,PSIp)

    end if


  end subroutine advance_GCinterp_2DBdB_vars

  subroutine advance_GCinterp_3DBdB_vars(vars,pp,tt,params,Y_R,Y_PHI,Y_Z, &
       V_PLL,V_MU,q_cache,m_cache,flagCon,flagCol,F,P,B_R,B_PHI,B_Z,E_PHI, &
       curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,PSIp)
    !! @note Subroutine to advance GC variables \(({\bf X},p_\parallel)\)
    !! @endnote
    !! Comment this section further with evolution equations, numerical
    !! methods, and descriptions of both.
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(PARTICLES), INTENT(INOUT)     :: vars
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    REAL(rp)                                      :: dt
    !! Time step used in the leapfrog step (\(\Delta t\)).

    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    INTEGER(ip),intent(in)                                      :: tt
    !! time iterator.
    INTEGER,intent(in)                                  :: pp


    REAL(rp),DIMENSION(params%pchunk)               :: Bmag
    REAL(rp)              :: a1 = 1./5._rp
    REAL(rp) :: a21 = 3./40._rp,a22=9./40._rp
    REAL(rp) :: a31 = 3./10._rp,a32=-9./10._rp,a33=6./5._rp
    REAL(rp) :: a41 = -11./54._rp,a42=5./2._rp,a43=-70./27._rp,a44=35./27._rp
    REAL(rp) :: a51 = 1631./55296._rp,a52=175./512._rp,a53=575./13824._rp,a54=44275./110592._rp,a55=253./4096._rp
    REAL(rp) :: b1=37./378._rp,b2=0._rp,b3=250./621._rp,b4=125./594._rp,b5=0._rp,b6=512./1771._rp

    REAL(rp),DIMENSION(params%pchunk) :: k1_R,k1_PHI,k1_Z,k1_PLL,k1_MU
    REAL(rp),DIMENSION(params%pchunk) :: k2_R,k2_PHI,k2_Z,k2_PLL,k2_MU
    REAL(rp),DIMENSION(params%pchunk) :: k3_R,k3_PHI,k3_Z,k3_PLL,k3_MU
    REAL(rp),DIMENSION(params%pchunk) :: k4_R,k4_PHI,k4_Z,k4_PLL,k4_MU
    REAL(rp),DIMENSION(params%pchunk) :: k5_R,k5_PHI,k5_Z,k5_PLL,k5_MU
    REAL(rp),DIMENSION(params%pchunk) :: k6_R,k6_PHI,k6_Z,k6_PLL,k6_MU
    REAL(rp),DIMENSION(params%pchunk) :: Y0_R,Y0_PHI,Y0_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: Y_R,Y_PHI,Y_Z,PSIp
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_R,E_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: E_PHI
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: gradB_R,gradB_PHI,gradB_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: V_PLL,V_MU
    REAL(rp),DIMENSION(params%pchunk) :: RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU
    REAL(rp),DIMENSION(params%pchunk) :: V0_PLL,V0_MU
    REAL(rp),DIMENSION(params%pchunk) :: ne,Te,Zeff

    INTEGER(is),DIMENSION(params%pchunk),intent(INOUT) :: flagCon,flagCol
    REAL(rp),intent(IN)  :: q_cache,m_cache

    dt=params%dt
    pchunk=params%pchunk

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU)
    do cc=1_idef,pchunk

       Y0_R(cc)=Y_R(cc)
       Y0_PHI(cc)=Y_PHI(cc)
       Y0_Z(cc)=Y_Z(cc)
       V0_PLL(cc)=V_PLL(cc)
       V0_MU(cc)=V_MU(cc)
    end do
    !$OMP END SIMD

    call calculate_3DBdBfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !write(output_unit_write,*) 'R',Y_R(1)
    !write(output_unit_write,*) 'PHI',Y_PHI(1)
    !write(output_unit_write,*) 'Z',Y_Z(1)
    !write(output_unit_write,*) 'PPLL',V_PLL(1)
    !write(output_unit_write,*) 'MU',V_MU(1)

    !write(output_unit_write,*) 'BR',B_R(1)
    !write(output_unit_write,*) 'BPHI',B_PHI(1)
    !write(output_unit_write,*) 'BZ',B_Z(1)

    !write(output_unit_write,*) 'gradBR',gradB_R(1)
    !write(output_unit_write,*) 'gradBPHI',gradB_PHI(1)
    !write(output_unit_write,*) 'gradBZ',gradB_Z(1)

    !write(output_unit_write,*) 'curlBR',curlB_R(1)
    !write(output_unit_write,*) 'curlBPHI',curlB_PHI(1)
    !write(output_unit_write,*) 'curlBZ',curlB_Z(1)

    !write(output_unit_write,*) 'RHS_R',RHS_R(1)
    !write(output_unit_write,*) 'RHS_PHI',RHS_PHI(1)
    !write(output_unit_write,*) 'RHS_Z',RHS_Z(1)
    !write(output_unit_write,*) 'RHS_PLL',RHS_PLL(1)


    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k1_R,k1_PHI,k1_Z,k1_PLL,k1_MU)
    do cc=1_idef,pchunk
       k1_R(cc)=dt*RHS_R(cc)              
       k1_PHI(cc)=dt*RHS_PHI(cc)    
       k1_Z(cc)=dt*RHS_Z(cc)    
       k1_PLL(cc)=dt*RHS_PLL(cc)
       k1_MU(cc)=dt*RHS_MU(cc)    

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a1*k1_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a1*k1_MU(cc)
    end do
    !$OMP END SIMD

    call calculate_3DBdBfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k2_R,k2_PHI,k2_Z,k2_PLL,k2_MU)
    do cc=1_idef,pchunk
       k2_R(cc)=dt*RHS_R(cc)    
       k2_PHI(cc)=dt*RHS_PHI (cc)   
       k2_Z(cc)=dt*RHS_Z(cc)   
       k2_PLL(cc)=dt*RHS_PLL(cc)
       k2_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a21*k1_PLL(cc)+a22*k2_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a21*k1_MU(cc)+a22*k2_MU(cc)
    end do
    !$OMP END SIMD

    call calculate_3DBdBfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k3_R,k3_PHI,k3_Z,k3_PLL,k3_MU)
    do cc=1_idef,pchunk
       k3_R(cc)=dt*RHS_R(cc)   
       k3_PHI(cc)=dt*RHS_PHI(cc)    
       k3_Z(cc)=dt*RHS_Z(cc)    
       k3_PLL(cc)=dt*RHS_PLL(cc)
       k3_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a31*k1_PLL(cc)+a32*k2_PLL(cc)+a33*k3_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a31*k1_MU(cc)+a32*k2_MU(cc)+a33*k3_MU(cc)
    end do
    !$OMP END SIMD

    call calculate_3DBdBfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)     

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k4_R,k4_PHI,k4_Z,k4_PLL,k4_MU)
    do cc=1_idef,pchunk
       k4_R(cc)=dt*RHS_R(cc)   
       k4_PHI(cc)=dt*RHS_PHI(cc)    
       k4_Z(cc)=dt*RHS_Z(cc)    
       k4_PLL(cc)=dt*RHS_PLL(cc)
       k4_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a41*k1_R(cc)+a42*k2_R(cc)+a43*k3_R(cc)+ &
            a44*k4_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a41*k1_PHI(cc)+a42*k2_PHI(cc)+ &
            a43*k3_PHI(cc)+a44*k4_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a41*k1_Z(cc)+a42*k2_Z(cc)+a43*k3_Z(cc)+ &
            a44*k4_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a41*k1_PLL(cc)+a42*k2_PLL(cc)+ &
            a43*k3_PLL(cc)+a44*k4_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a41*k1_MU(cc)+a42*k2_MU(cc)+ &
            a43*k3_MU(cc)+a44*k4_MU(cc)
    end do
    !$OMP END SIMD


    call calculate_3DBdBfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)   

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k5_R,k5_PHI,k5_Z,k5_PLL,k5_MU)
    do cc=1_idef,pchunk
       k5_R(cc)=dt*RHS_R(cc)    
       k5_PHI(cc)=dt*RHS_PHI(cc)    
       k5_Z(cc)=dt*RHS_Z(cc)    
       k5_PLL(cc)=dt*RHS_PLL(cc)
       k5_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a51*k1_R(cc)+a52*k2_R(cc)+a53*k3_R(cc)+ &
            a54*k4_R(cc)+a55*k5_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a51*k1_PHI(cc)+a52*k2_PHI(cc)+ &
            a53*k3_PHI(cc)+a54*k4_PHI(cc)+a55*k5_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a51*k1_Z(cc)+a52*k2_Z(cc)+a53*k3_Z(cc)+ &
            a54*k4_Z(cc)+a55*k5_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a51*k1_PLL(cc)+a52*k2_PLL(cc)+ &
            a53*k3_PLL(cc)+a54*k4_PLL(cc)+a55*k5_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a51*k1_MU(cc)+a52*k2_MU(cc)+ &
            a53*k3_MU(cc)+a54*k4_MU(cc)+a55*k5_MU(cc)
    end do
    !$OMP END SIMD

    call calculate_3DBdBfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)         

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k6_R,k6_PHI,k6_Z,k6_PLL,k6_MU)
    do cc=1_idef,pchunk
       k6_R(cc)=dt*RHS_R(cc)    
       k6_PHI(cc)=dt*RHS_PHI(cc)    
       k6_Z(cc)=dt*RHS_Z(cc)    
       k6_PLL(cc)=dt*RHS_PLL(cc)
       k6_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+b1*k1_R(cc)+b2*k2_R(cc)+ &
            b3*k3_R(cc)+b4*k4_R(cc)+b5*k5_R(cc)+b6*k6_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+b1*k1_PHI(cc)+b2*k2_PHI(cc)+ &
            b3*k3_PHI(cc)+b4*k4_PHI(cc)+b5*k5_PHI(cc)+b6*k6_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+b1*k1_Z(cc)+b2*k2_Z(cc)+ &
            b3*k3_Z(cc)+b4*k4_Z(cc)+b5*k5_Z(cc)+b6*k6_Z(cc)
       V_PLL(cc)=V0_PLL(cc)+b1*k1_PLL(cc)+b2*k2_PLL(cc)+ &
            b3*k3_PLL(cc)+b4*k4_PLL(cc)+b5*k5_PLL(cc)+b6*k6_PLL(cc)
       V_MU(cc)=V0_MU(cc)+b1*k1_MU(cc)+b2*k2_MU(cc)+ &
            b3*k3_MU(cc)+b4*k4_MU(cc)+b5*k5_MU(cc)+b6*k6_MU(cc)
    end do
    !$OMP END SIMD

    !$OMP SIMD
    !    !$OMP& aligned(Y_R,Y_PHI,Y_Z,V_PLL,V_MU,Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU)
    do cc=1_idef,pchunk

       if ((flagCon(cc).eq.0_is).or.(flagCol(cc).eq.0_is)) then
          Y_R(cc)=Y0_R(cc)
          Y_PHI(cc)=Y0_PHI(cc)
          Y_Z(cc)=Y0_Z(cc)
          V_PLL(cc)=V0_PLL(cc)
          V_MU(cc)=V0_MU(cc)
       end if

    end do
    !$OMP END SIMD

    call calculate_3DBdBfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !$OMP SIMD
    do cc=1_idef,pchunk
       vars%RHS(pp-1+cc,1)=RHS_R(cc)
       vars%RHS(pp-1+cc,2)=RHS_PHI(cc)
       vars%RHS(pp-1+cc,3)=RHS_Z(cc)
       vars%RHS(pp-1+cc,4)=RHS_PLL(cc)
       vars%RHS(pp-1+cc,5)=RHS_MU(cc)
    end do
    !$OMP END SIMD



    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    if (params%collisions) then       

       call include_CoulombCollisions_GC_p(tt,params,Y_R,Y_PHI,Y_Z, &
            V_PLL,V_MU,m_cache,flagCon,flagCol,F,P,E_PHI,ne,PSIp)

    end if


  end subroutine advance_GCinterp_3DBdB_vars


  subroutine advance_GCinterp_3DBdB1_vars(vars,pp,tt,params,Y_R,Y_PHI,Y_Z, &
       V_PLL,V_MU,q_cache,m_cache,flagCon,flagCol,F,P,B_R,B_PHI,B_Z,E_PHI, &
       curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,PSIp)
    !! @note Subroutine to advance GC variables \(({\bf X},p_\parallel)\)
    !! @endnote
    !! Comment this section further with evolution equations, numerical
    !! methods, and descriptions of both.
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(PARTICLES), INTENT(INOUT)     :: vars
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    REAL(rp)                                      :: dt
    !! Time step used in the leapfrog step (\(\Delta t\)).

    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    INTEGER(ip),intent(in)                                      :: tt
    !! time iterator.
    INTEGER,intent(in)                                  :: pp


    REAL(rp),DIMENSION(params%pchunk)               :: Bmag
    REAL(rp)              :: a1 = 1./5._rp
    REAL(rp) :: a21 = 3./40._rp,a22=9./40._rp
    REAL(rp) :: a31 = 3./10._rp,a32=-9./10._rp,a33=6./5._rp
    REAL(rp) :: a41 = -11./54._rp,a42=5./2._rp,a43=-70./27._rp,a44=35./27._rp
    REAL(rp) :: a51 = 1631./55296._rp,a52=175./512._rp,a53=575./13824._rp,a54=44275./110592._rp,a55=253./4096._rp
    REAL(rp) :: b1=37./378._rp,b2=0._rp,b3=250./621._rp,b4=125./594._rp,b5=0._rp,b6=512./1771._rp

    REAL(rp),DIMENSION(params%pchunk) :: k1_R,k1_PHI,k1_Z,k1_PLL,k1_MU
    REAL(rp),DIMENSION(params%pchunk) :: k2_R,k2_PHI,k2_Z,k2_PLL,k2_MU
    REAL(rp),DIMENSION(params%pchunk) :: k3_R,k3_PHI,k3_Z,k3_PLL,k3_MU
    REAL(rp),DIMENSION(params%pchunk) :: k4_R,k4_PHI,k4_Z,k4_PLL,k4_MU
    REAL(rp),DIMENSION(params%pchunk) :: k5_R,k5_PHI,k5_Z,k5_PLL,k5_MU
    REAL(rp),DIMENSION(params%pchunk) :: k6_R,k6_PHI,k6_Z,k6_PLL,k6_MU
    REAL(rp),DIMENSION(params%pchunk) :: Y0_R,Y0_PHI,Y0_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_R,E_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: E_PHI,PSIp
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: gradB_R,gradB_PHI,gradB_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: V_PLL,V_MU
    REAL(rp),DIMENSION(params%pchunk) :: RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU
    REAL(rp),DIMENSION(params%pchunk) :: V0_PLL,V0_MU
    REAL(rp),DIMENSION(params%pchunk) :: ne,Te,Zeff

    INTEGER(is),DIMENSION(params%pchunk),intent(INOUT) :: flagCon,flagCol
    REAL(rp),intent(IN)  :: q_cache,m_cache

    dt=params%dt
    pchunk=params%pchunk

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU)
    do cc=1_idef,pchunk

       Y0_R(cc)=Y_R(cc)
       Y0_PHI(cc)=Y_PHI(cc)
       Y0_Z(cc)=Y_Z(cc)
       V0_PLL(cc)=V_PLL(cc)
       V0_MU(cc)=V_MU(cc)
    end do
    !$OMP END SIMD

    call calculate_3DBdBfields1_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !write(output_unit_write,*) 'R',Y_R(1)
    !write(output_unit_write,*) 'PHI',Y_PHI(1)
    !write(output_unit_write,*) 'Z',Y_Z(1)
    !write(output_unit_write,*) 'PPLL',V_PLL(1)
    !write(output_unit_write,*) 'MU',V_MU(1)

    !write(output_unit_write,*) 'BR',B_R(1)
    !write(output_unit_write,*) 'BPHI',B_PHI(1)
    !write(output_unit_write,*) 'BZ',B_Z(1)

    !write(output_unit_write,*) 'gradBR',gradB_R(1)
    !write(output_unit_write,*) 'gradBPHI',gradB_PHI(1)
    !write(output_unit_write,*) 'gradBZ',gradB_Z(1)

    !write(output_unit_write,*) 'curlBR',curlB_R(1)
    !write(output_unit_write,*) 'curlBPHI',curlB_PHI(1)
    !write(output_unit_write,*) 'curlBZ',curlB_Z(1)

    !write(output_unit_write,*) 'RHS_R',RHS_R(1)
    !write(output_unit_write,*) 'RHS_PHI',RHS_PHI(1)
    !write(output_unit_write,*) 'RHS_Z',RHS_Z(1)
    !write(output_unit_write,*) 'RHS_PLL',RHS_PLL(1)


    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k1_R,k1_PHI,k1_Z,k1_PLL,k1_MU)
    do cc=1_idef,pchunk
       k1_R(cc)=dt*RHS_R(cc)              
       k1_PHI(cc)=dt*RHS_PHI(cc)    
       k1_Z(cc)=dt*RHS_Z(cc)    
       k1_PLL(cc)=dt*RHS_PLL(cc)
       k1_MU(cc)=dt*RHS_MU(cc)    

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a1*k1_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a1*k1_MU(cc)
    end do
    !$OMP END SIMD

    call calculate_3DBdBfields1_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k2_R,k2_PHI,k2_Z,k2_PLL,k2_MU)
    do cc=1_idef,pchunk
       k2_R(cc)=dt*RHS_R(cc)    
       k2_PHI(cc)=dt*RHS_PHI (cc)   
       k2_Z(cc)=dt*RHS_Z(cc)   
       k2_PLL(cc)=dt*RHS_PLL(cc)
       k2_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a21*k1_PLL(cc)+a22*k2_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a21*k1_MU(cc)+a22*k2_MU(cc)
    end do
    !$OMP END SIMD

    call calculate_3DBdBfields1_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k3_R,k3_PHI,k3_Z,k3_PLL,k3_MU)
    do cc=1_idef,pchunk
       k3_R(cc)=dt*RHS_R(cc)   
       k3_PHI(cc)=dt*RHS_PHI(cc)    
       k3_Z(cc)=dt*RHS_Z(cc)    
       k3_PLL(cc)=dt*RHS_PLL(cc)
       k3_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a31*k1_PLL(cc)+a32*k2_PLL(cc)+a33*k3_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a31*k1_MU(cc)+a32*k2_MU(cc)+a33*k3_MU(cc)
    end do
    !$OMP END SIMD

    call calculate_3DBdBfields1_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)     

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k4_R,k4_PHI,k4_Z,k4_PLL,k4_MU)
    do cc=1_idef,pchunk
       k4_R(cc)=dt*RHS_R(cc)   
       k4_PHI(cc)=dt*RHS_PHI(cc)    
       k4_Z(cc)=dt*RHS_Z(cc)    
       k4_PLL(cc)=dt*RHS_PLL(cc)
       k4_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a41*k1_R(cc)+a42*k2_R(cc)+a43*k3_R(cc)+ &
            a44*k4_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a41*k1_PHI(cc)+a42*k2_PHI(cc)+ &
            a43*k3_PHI(cc)+a44*k4_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a41*k1_Z(cc)+a42*k2_Z(cc)+a43*k3_Z(cc)+ &
            a44*k4_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a41*k1_PLL(cc)+a42*k2_PLL(cc)+ &
            a43*k3_PLL(cc)+a44*k4_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a41*k1_MU(cc)+a42*k2_MU(cc)+ &
            a43*k3_MU(cc)+a44*k4_MU(cc)
    end do
    !$OMP END SIMD


    call calculate_3DBdBfields1_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)   

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k5_R,k5_PHI,k5_Z,k5_PLL,k5_MU)
    do cc=1_idef,pchunk
       k5_R(cc)=dt*RHS_R(cc)    
       k5_PHI(cc)=dt*RHS_PHI(cc)    
       k5_Z(cc)=dt*RHS_Z(cc)    
       k5_PLL(cc)=dt*RHS_PLL(cc)
       k5_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a51*k1_R(cc)+a52*k2_R(cc)+a53*k3_R(cc)+ &
            a54*k4_R(cc)+a55*k5_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a51*k1_PHI(cc)+a52*k2_PHI(cc)+ &
            a53*k3_PHI(cc)+a54*k4_PHI(cc)+a55*k5_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a51*k1_Z(cc)+a52*k2_Z(cc)+a53*k3_Z(cc)+ &
            a54*k4_Z(cc)+a55*k5_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a51*k1_PLL(cc)+a52*k2_PLL(cc)+ &
            a53*k3_PLL(cc)+a54*k4_PLL(cc)+a55*k5_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a51*k1_MU(cc)+a52*k2_MU(cc)+ &
            a53*k3_MU(cc)+a54*k4_MU(cc)+a55*k5_MU(cc)
    end do
    !$OMP END SIMD

    call calculate_3DBdBfields1_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)         

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k6_R,k6_PHI,k6_Z,k6_PLL,k6_MU)
    do cc=1_idef,pchunk
       k6_R(cc)=dt*RHS_R(cc)    
       k6_PHI(cc)=dt*RHS_PHI(cc)    
       k6_Z(cc)=dt*RHS_Z(cc)    
       k6_PLL(cc)=dt*RHS_PLL(cc)
       k6_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+b1*k1_R(cc)+b2*k2_R(cc)+ &
            b3*k3_R(cc)+b4*k4_R(cc)+b5*k5_R(cc)+b6*k6_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+b1*k1_PHI(cc)+b2*k2_PHI(cc)+ &
            b3*k3_PHI(cc)+b4*k4_PHI(cc)+b5*k5_PHI(cc)+b6*k6_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+b1*k1_Z(cc)+b2*k2_Z(cc)+ &
            b3*k3_Z(cc)+b4*k4_Z(cc)+b5*k5_Z(cc)+b6*k6_Z(cc)
       V_PLL(cc)=V0_PLL(cc)+b1*k1_PLL(cc)+b2*k2_PLL(cc)+ &
            b3*k3_PLL(cc)+b4*k4_PLL(cc)+b5*k5_PLL(cc)+b6*k6_PLL(cc)
       V_MU(cc)=V0_MU(cc)+b1*k1_MU(cc)+b2*k2_MU(cc)+ &
            b3*k3_MU(cc)+b4*k4_MU(cc)+b5*k5_MU(cc)+b6*k6_MU(cc)
    end do
    !$OMP END SIMD

    !$OMP SIMD
    !    !$OMP& aligned(Y_R,Y_PHI,Y_Z,V_PLL,V_MU,Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU)
    do cc=1_idef,pchunk

       if ((flagCon(cc).eq.0_is).or.(flagCol(cc).eq.0_is)) then
          Y_R(cc)=Y0_R(cc)
          Y_PHI(cc)=Y0_PHI(cc)
          Y_Z(cc)=Y0_Z(cc)
          V_PLL(cc)=V0_PLL(cc)
          V_MU(cc)=V0_MU(cc)
       end if

    end do
    !$OMP END SIMD

    call calculate_3DBdBfields1_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon,PSIp)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !$OMP SIMD
    do cc=1_idef,pchunk
       vars%RHS(pp-1+cc,1)=RHS_R(cc)
       vars%RHS(pp-1+cc,2)=RHS_PHI(cc)
       vars%RHS(pp-1+cc,3)=RHS_Z(cc)
       vars%RHS(pp-1+cc,4)=RHS_PLL(cc)
       vars%RHS(pp-1+cc,5)=RHS_MU(cc)
    end do
    !$OMP END SIMD



    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    if (params%collisions) then       

       call include_CoulombCollisions_GC_p(tt,params,Y_R,Y_PHI,Y_Z, &
            V_PLL,V_MU,m_cache,flagCon,flagCol,F,P,E_PHI,ne,PSIp)

    end if


  end subroutine advance_GCinterp_3DBdB1_vars

  subroutine advance_GCinterp_B_vars(vars,pp,tt,params,Y_R,Y_PHI,Y_Z, &
       V_PLL,V_MU,q_cache,m_cache,flagCon,flagCol,F,P,B_R,B_PHI,B_Z,E_PHI, &
       curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,PSIp)
    !! @note Subroutine to advance GC variables \(({\bf X},p_\parallel)\)
    !! @endnote
    !! Comment this section further with evolution equations, numerical
    !! methods, and descriptions of both.
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(PARTICLES), INTENT(INOUT)     :: vars
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    REAL(rp)                                      :: dt
    !! Time step used in the leapfrog step (\(\Delta t\)).

    INTEGER                                                    :: cc,pchunk
    !! Chunk iterator.
    INTEGER(ip),intent(in)                                      :: tt
    !! time iterator.
    INTEGER,intent(in)                                  :: pp


    REAL(rp),DIMENSION(params%pchunk)               :: Bmag
    REAL(rp)              :: a1 = 1./5._rp
    REAL(rp) :: a21 = 3./40._rp,a22=9./40._rp
    REAL(rp) :: a31 = 3./10._rp,a32=-9./10._rp,a33=6./5._rp
    REAL(rp) :: a41 = -11./54._rp,a42=5./2._rp,a43=-70./27._rp,a44=35./27._rp
    REAL(rp) :: a51 = 1631./55296._rp,a52=175./512._rp,a53=575./13824._rp,a54=44275./110592._rp,a55=253./4096._rp
    REAL(rp) :: b1=37./378._rp,b2=0._rp,b3=250./621._rp,b4=125./594._rp,b5=0._rp,b6=512./1771._rp

    REAL(rp),DIMENSION(params%pchunk) :: k1_R,k1_PHI,k1_Z,k1_PLL,k1_MU
    REAL(rp),DIMENSION(params%pchunk) :: k2_R,k2_PHI,k2_Z,k2_PLL,k2_MU
    REAL(rp),DIMENSION(params%pchunk) :: k3_R,k3_PHI,k3_Z,k3_PLL,k3_MU
    REAL(rp),DIMENSION(params%pchunk) :: k4_R,k4_PHI,k4_Z,k4_PLL,k4_MU
    REAL(rp),DIMENSION(params%pchunk) :: k5_R,k5_PHI,k5_Z,k5_PLL,k5_MU
    REAL(rp),DIMENSION(params%pchunk) :: k6_R,k6_PHI,k6_Z,k6_PLL,k6_MU
    REAL(rp),DIMENSION(params%pchunk) :: Y0_R,Y0_PHI,Y0_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: Y_R,Y_PHI,Y_Z,PSIp
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: B_R,B_PHI,B_Z
    REAL(rp),DIMENSION(params%pchunk) :: E_R,E_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: E_PHI
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: curlb_R,curlb_PHI,curlb_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: gradB_R,gradB_PHI,gradB_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(INOUT) :: V_PLL,V_MU
    REAL(rp),DIMENSION(params%pchunk) :: RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU
    REAL(rp),DIMENSION(params%pchunk) :: V0_PLL,V0_MU
    REAL(rp),DIMENSION(params%pchunk) :: ne,Te,Zeff

    INTEGER(is),DIMENSION(params%pchunk),intent(INOUT) :: flagCon,flagCol
    REAL(rp),intent(IN)  :: q_cache,m_cache

    dt=params%dt
    pchunk=params%pchunk

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU)
    do cc=1_idef,pchunk

       Y0_R(cc)=Y_R(cc)
       Y0_PHI(cc)=Y_PHI(cc)
       Y0_Z(cc)=Y_Z(cc)
       V0_PLL(cc)=V_PLL(cc)
       V0_MU(cc)=V_MU(cc)
    end do
    !$OMP END SIMD

    call interp_fields_3D_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k1_R,k1_PHI,k1_Z,k1_PLL,k1_MU)
    do cc=1_idef,pchunk
       k1_R(cc)=dt*RHS_R(cc)              
       k1_PHI(cc)=dt*RHS_PHI(cc)    
       k1_Z(cc)=dt*RHS_Z(cc)    
       k1_PLL(cc)=dt*RHS_PLL(cc)
       k1_MU(cc)=dt*RHS_MU(cc)    

       Y_R(cc)=Y0_R(cc)+a1*k1_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a1*k1_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a1*k1_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a1*k1_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a1*k1_MU(cc)
    end do
    !$OMP END SIMD

    call interp_fields_3D_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k2_R,k2_PHI,k2_Z,k2_PLL,k2_MU)
    do cc=1_idef,pchunk
       k2_R(cc)=dt*RHS_R(cc)    
       k2_PHI(cc)=dt*RHS_PHI (cc)   
       k2_Z(cc)=dt*RHS_Z(cc)   
       k2_PLL(cc)=dt*RHS_PLL(cc)
       k2_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a21*k1_R(cc)+a22*k2_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a21*k1_PHI(cc)+a22*k2_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a21*k1_Z(cc)+a22*k2_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a21*k1_PLL(cc)+a22*k2_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a21*k1_MU(cc)+a22*k2_MU(cc)
    end do
    !$OMP END SIMD

    call interp_fields_3D_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k3_R,k3_PHI,k3_Z,k3_PLL,k3_MU)
    do cc=1_idef,pchunk
       k3_R(cc)=dt*RHS_R(cc)   
       k3_PHI(cc)=dt*RHS_PHI(cc)    
       k3_Z(cc)=dt*RHS_Z(cc)    
       k3_PLL(cc)=dt*RHS_PLL(cc)
       k3_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a31*k1_R(cc)+a32*k2_R(cc)+a33*k3_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a31*k1_PHI(cc)+a32*k2_PHI(cc)+ &
            a33*k3_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a31*k1_Z(cc)+a32*k2_Z(cc)+a33*k3_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a31*k1_PLL(cc)+a32*k2_PLL(cc)+a33*k3_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a31*k1_MU(cc)+a32*k2_MU(cc)+a33*k3_MU(cc)
    end do
    !$OMP END SIMD

    call interp_fields_3D_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)     

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k4_R,k4_PHI,k4_Z,k4_PLL,k4_MU)
    do cc=1_idef,pchunk
       k4_R(cc)=dt*RHS_R(cc)   
       k4_PHI(cc)=dt*RHS_PHI(cc)    
       k4_Z(cc)=dt*RHS_Z(cc)    
       k4_PLL(cc)=dt*RHS_PLL(cc)
       k4_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a41*k1_R(cc)+a42*k2_R(cc)+a43*k3_R(cc)+ &
            a44*k4_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a41*k1_PHI(cc)+a42*k2_PHI(cc)+ &
            a43*k3_PHI(cc)+a44*k4_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a41*k1_Z(cc)+a42*k2_Z(cc)+a43*k3_Z(cc)+ &
            a44*k4_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a41*k1_PLL(cc)+a42*k2_PLL(cc)+ &
            a43*k3_PLL(cc)+a44*k4_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a41*k1_MU(cc)+a42*k2_MU(cc)+ &
            a43*k3_MU(cc)+a44*k4_MU(cc)
    end do
    !$OMP END SIMD


    call interp_fields_3D_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)   

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k5_R,k5_PHI,k5_Z,k5_PLL,k5_MU)
    do cc=1_idef,pchunk
       k5_R(cc)=dt*RHS_R(cc)    
       k5_PHI(cc)=dt*RHS_PHI(cc)    
       k5_Z(cc)=dt*RHS_Z(cc)    
       k5_PLL(cc)=dt*RHS_PLL(cc)
       k5_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+a51*k1_R(cc)+a52*k2_R(cc)+a53*k3_R(cc)+ &
            a54*k4_R(cc)+a55*k5_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+a51*k1_PHI(cc)+a52*k2_PHI(cc)+ &
            a53*k3_PHI(cc)+a54*k4_PHI(cc)+a55*k5_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+a51*k1_Z(cc)+a52*k2_Z(cc)+a53*k3_Z(cc)+ &
            a54*k4_Z(cc)+a55*k5_Z(cc)
       V_PLL(cc)=V0_PLL(cc)   +a51*k1_PLL(cc)+a52*k2_PLL(cc)+ &
            a53*k3_PLL(cc)+a54*k4_PLL(cc)+a55*k5_PLL(cc)
       V_MU(cc)=V0_MU(cc)   +a51*k1_MU(cc)+a52*k2_MU(cc)+ &
            a53*k3_MU(cc)+a54*k4_MU(cc)+a55*k5_MU(cc)
    end do
    !$OMP END SIMD

    call interp_fields_3D_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon)         

    !$OMP SIMD
    !    !$OMP& aligned(Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
    !    !$OMP& RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& k6_R,k6_PHI,k6_Z,k6_PLL,k6_MU)
    do cc=1_idef,pchunk
       k6_R(cc)=dt*RHS_R(cc)    
       k6_PHI(cc)=dt*RHS_PHI(cc)    
       k6_Z(cc)=dt*RHS_Z(cc)    
       k6_PLL(cc)=dt*RHS_PLL(cc)
       k6_MU(cc)=dt*RHS_MU(cc)

       Y_R(cc)=Y0_R(cc)+b1*k1_R(cc)+b2*k2_R(cc)+ &
            b3*k3_R(cc)+b4*k4_R(cc)+b5*k5_R(cc)+b6*k6_R(cc)
       Y_PHI(cc)=Y0_PHI(cc)+b1*k1_PHI(cc)+b2*k2_PHI(cc)+ &
            b3*k3_PHI(cc)+b4*k4_PHI(cc)+b5*k5_PHI(cc)+b6*k6_PHI(cc)
       Y_Z(cc)=Y0_Z(cc)+b1*k1_Z(cc)+b2*k2_Z(cc)+ &
            b3*k3_Z(cc)+b4*k4_Z(cc)+b5*k5_Z(cc)+b6*k6_Z(cc)
       V_PLL(cc)=V0_PLL(cc)+b1*k1_PLL(cc)+b2*k2_PLL(cc)+ &
            b3*k3_PLL(cc)+b4*k4_PLL(cc)+b5*k5_PLL(cc)+b6*k6_PLL(cc)
       V_MU(cc)=V0_MU(cc)+b1*k1_MU(cc)+b2*k2_MU(cc)+ &
            b3*k3_MU(cc)+b4*k4_MU(cc)+b5*k5_MU(cc)+b6*k6_MU(cc)
    end do
    !$OMP END SIMD

    !$OMP SIMD
    !    !$OMP& aligned(Y_R,Y_PHI,Y_Z,V_PLL,V_MU,Y0_R,Y0_PHI,Y0_Z,V0_PLL,V0_MU)
    do cc=1_idef,pchunk

       if ((flagCon(cc).eq.0_is).or.(flagCol(cc).eq.0_is)) then
          Y_R(cc)=Y0_R(cc)
          Y_PHI(cc)=Y0_PHI(cc)
          Y_Z(cc)=Y0_Z(cc)
          V_PLL(cc)=V0_PLL(cc)
          V_MU(cc)=V0_MU(cc)
       end if

    end do
    !$OMP END SIMD

    call interp_fields_3D_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI, &
         E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
         flagCon)

    call GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU,B_R,B_PHI, &
         B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R, &
         gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flagCon) 

    !$OMP SIMD
    do cc=1_idef,pchunk
       vars%RHS(pp-1+cc,1)=RHS_R(cc)
       vars%RHS(pp-1+cc,2)=RHS_PHI(cc)
       vars%RHS(pp-1+cc,3)=RHS_Z(cc)
       vars%RHS(pp-1+cc,4)=RHS_PLL(cc)
       vars%RHS(pp-1+cc,5)=RHS_MU(cc)
    end do
    !$OMP END SIMD



    call add_analytical_E_p(params,tt,F,E_PHI,Y_R,Y_Z)

    if (params%collisions) then       

       call include_CoulombCollisions_GC_p(tt,params,Y_R,Y_PHI,Y_Z, &
            V_PLL,V_MU,m_cache,flagCon,flagCol,F,P,E_PHI,ne,PSIp)

    end if


  end subroutine advance_GCinterp_B_vars
#endif

  subroutine advance_FPinterp_vars(params,Y_R,Y_PHI,Y_Z,V_PLL,V_MU, &
       m_cache,flagCon,flagCol,F,P,E_PHI,ne,PSIp)    
    TYPE(KORC_PARAMS), INTENT(INOUT)                              :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(IN)                                 :: P
    TYPE(FIELDS), INTENT(IN)                                   :: F
    INTEGER(ip)                                                    :: tt
    !! time iterator.
    REAL(rp),DIMENSION(params%pchunk), INTENT(IN)  :: Y_R,Y_PHI,Y_Z
    REAL(rp),DIMENSION(params%pchunk), INTENT(INOUT)  :: V_PLL,V_MU,PSIp
    REAL(rp),DIMENSION(params%pchunk), INTENT(OUT)  :: E_PHI
    REAL(rp),intent(in) :: m_cache
    INTEGER(is),DIMENSION(params%pchunk),intent(INOUT) :: flagCon,flagCol
    REAL(rp),DIMENSION(params%pchunk), INTENT(OUT) :: ne

    !    write(output_unit_write,'("E_PHI_FP: ",E17.10)') E_PHI

    do tt=1_ip,params%t_skip
       call include_CoulombCollisions_GC_p(tt,params,Y_R,Y_PHI,Y_Z, &
            V_PLL,V_MU,m_cache,flagCon,flagCol,F,P,E_PHI,ne,PSIp)

       !       write(output_unit_write,'("Collision Loop in FP")')

    end do

    !    write(output_unit_write,'("V_PLL: ",E17.10)') V_PLL
    !    write(output_unit_write,'("V_MU: ",E17.10)') V_MU

  end subroutine advance_FPinterp_vars



  subroutine GCEoM_p(params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,B_R,B_PHI, &
       B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI, &
       gradB_Z,V_PLL,V_MU,Y_R,q_cache,m_cache)
    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    REAL(rp),DIMENSION(params%pchunk)  :: Bmag,bhat_R,bhat_PHI,bhat_Z,Bst_R,Bst_PHI
    REAL(rp),DIMENSION(params%pchunk)  :: BstdotE,BstdotgradB,EcrossB_R,EcrossB_PHI,bdotBst
    REAL(rp),DIMENSION(params%pchunk)  :: bcrossgradB_R,bcrossgradB_PHI,bcrossgradB_Z,gamgc
    REAL(rp),DIMENSION(params%pchunk)  :: EcrossB_Z,Bst_Z
    REAL(rp),DIMENSION(params%pchunk)  :: pm,xi,tau_R
    REAL(rp),DIMENSION(params%pchunk),INTENT(in) :: gradB_R,gradB_PHI,gradB_Z,curlb_R
    REAL(rp),DIMENSION(params%pchunk),INTENT(in) :: curlb_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: RHS_R,RHS_PHI,RHS_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: RHS_PLL
    REAL(rp),DIMENSION(params%pchunk),INTENT(IN) :: V_PLL,V_MU,Y_R,curlb_PHI
    REAL(rp),INTENT(in) :: q_cache,m_cache
    INTEGER(ip)  :: cc,pchunk

    pchunk=params%pchunk

    !$OMP SIMD
    !    !$OMP& aligned(gradB_R,gradB_PHI,gradB_Z,curlb_R,curlb_Z, &
    !    !$OMP& B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,RHS_R,RHS_PHI,RHS_Z,RHS_PLL, &
    !    !$OMP& V_PLL,V_MU,Y_R,curlb_PHI)
    do cc=1_idef,pchunk
       Bmag(cc) = SQRT(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

       bhat_R(cc) = B_R(cc)/Bmag(cc)
       bhat_PHI(cc) = B_PHI(cc)/Bmag(cc)
       bhat_Z(cc) = B_Z(cc)/Bmag(cc)

       Bst_R(cc)=q_cache*B_R(cc)+V_PLL(cc)*curlb_R(cc)
       Bst_PHI(cc)=q_cache*B_PHI(cc)+V_PLL(cc)*curlb_PHI(cc)
       Bst_Z(cc)=q_cache*B_Z(cc)+V_PLL(cc)*curlb_Z(cc)

       bdotBst(cc)=bhat_R(cc)*Bst_R(cc)+bhat_PHI(cc)*Bst_PHI(cc)+ &
            bhat_Z(cc)*Bst_Z(cc)
       BstdotE(cc)=Bst_R(cc)*E_R(cc)+Bst_PHI(cc)*E_PHI(cc)+Bst_Z(cc)*E_Z(cc)   
       BstdotgradB(cc)=Bst_R(cc)*gradB_R(cc)+Bst_PHI(cc)*gradB_PHI(cc)+ &
            Bst_Z(cc)*gradB_Z(cc)

       Ecrossb_R(cc)=E_PHI(cc)*bhat_Z(cc)-E_Z(cc)*bhat_PHI(cc)
       Ecrossb_PHI(cc)=E_Z(cc)*bhat_R(cc)-E_R(cc)*bhat_Z(cc)
       Ecrossb_Z(cc)=E_R(cc)*bhat_PHI(cc)-E_PHI(cc)*bhat_R(cc)


       bcrossgradB_R(cc)=bhat_PHI(cc)*gradB_Z(cc)-bhat_Z(cc)*gradB_PHI(cc)
       bcrossgradB_PHI(cc)=bhat_Z(cc)*gradB_R(cc)-bhat_R(cc)*gradB_Z(cc)
       bcrossgradB_Z(cc)=bhat_R(cc)*gradB_PHI(cc)-bhat_PHI(cc)*gradB_R(cc)

       gamgc(cc)=sqrt(1+V_PLL(cc)*V_PLL(cc)+2*V_MU(cc)*Bmag(cc))

       pm(cc)=sqrt(gamgc(cc)**2-1)
       xi(cc)=V_PLL(cc)/pm(cc)

       RHS_R(cc)=(q_cache*Ecrossb_R(cc)+(m_cache*V_MU(cc)* &
            bcrossgradB_R(cc)+V_PLL(cc)*Bst_R(cc))/(m_cache*gamgc(cc)))/ &
            bdotBst(cc)
       RHS_PHI(cc)=(q_cache*Ecrossb_PHI(cc)+(m_cache*V_MU(cc)* &
            bcrossgradB_PHI(cc)+V_PLL(cc)*Bst_PHI(cc))/(m_cache*gamgc(cc)))/ &
            (Y_R(cc)*bdotBst(cc))
       RHS_Z(cc)=(q_cache*Ecrossb_Z(cc)+(m_cache*V_MU(cc)* &
            bcrossgradB_Z(cc)+V_PLL(cc)*Bst_Z(cc))/(m_cache*gamgc(cc)))/ &
            bdotBst(cc)
       RHS_PLL(cc)=(q_cache*BstdotE(cc)-V_MU(cc)*BstdotgradB(cc)/gamgc(cc))/ &
            bdotBst(cc)

    end do
    !$OMP END SIMD

    !    write(output_unit_write,*) 'RHS_R: ',RHS_R(1)
    !    write(output_unit_write,*) 'RHS_PHI: ',RHS_PHI(1)
    !    write(output_unit_write,*) 'RHS_Z: ',RHS_Z(1)
    !    write(output_unit_write,*) 'RHS_PLL: ',RHS_PLL(1)

  end subroutine GCEoM_p

  subroutine GCEoM1_p(pchunk,tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
       B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
       gradB_R,gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp,ne,flag_cache)

    USE omp_lib
    IMPLICIT NONE
    
    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)      :: F
    TYPE(PROFILES), INTENT(IN)                                 :: P
    INTEGER, INTENT(IN) :: pchunk
    REAL(rp),DIMENSION(pchunk)  :: Bmag,bhat_R,bhat_PHI,bhat_Z,Bst_R,Bst_PHI
    REAL(rp),DIMENSION(pchunk)  :: BstdotE,BstdotgradB,EcrossB_R,EcrossB_PHI,bdotBst
    REAL(rp),DIMENSION(pchunk)  :: bcrossgradB_R,bcrossgradB_PHI,bcrossgradB_Z,gamgc
    REAL(rp),DIMENSION(pchunk)  :: EcrossB_Z,Bst_Z
    REAL(rp),DIMENSION(pchunk)  :: pm,xi,tau_R
    REAL(rp),DIMENSION(pchunk)  :: SR_PLL,SR_MU,BREM_PLL,BREM_MU,BREM_P
    REAL(rp),DIMENSION(pchunk),INTENT(in) :: gradB_R,gradB_PHI,gradB_Z,curlb_R
    REAL(rp),DIMENSION(pchunk),INTENT(in) :: curlb_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z
    REAL(rp),DIMENSION(pchunk),INTENT(OUT) :: RHS_R,RHS_PHI,RHS_Z
    REAL(rp),DIMENSION(pchunk),INTENT(OUT) :: RHS_PLL,RHS_MU
    REAL(rp),DIMENSION(pchunk),INTENT(IN) :: V_PLL,V_MU,Y_R,Y_PHI,Y_Z,curlb_PHI
    REAL(rp),DIMENSION(pchunk),INTENT(IN) :: PSIp
    REAL(rp),INTENT(in) :: q_cache,m_cache
    INTEGER  :: cc
    INTEGER(ip),INTENT(IN)  :: tt
    INTEGER(is),DIMENSION(pchunk),INTENT(OUT)  :: flag_cache
    REAL(rp)  :: time,re_cache,alpha_cache
    REAL(rp), DIMENSION(pchunk) 			:: Zeff,Te
    REAL(rp), DIMENSION(pchunk),INTENT(OUT) 		:: ne
    INTEGER             :: thread_num


    thread_num = OMP_GET_THREAD_NUM()

    !$OMP SIMD
    !    !$OMP& aligned(gradB_R,gradB_PHI,gradB_Z,curlb_R,curlb_Z, &
    !    !$OMP& B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& V_PLL,V_MU,Y_R,curlb_PHI,tau_R)
    do cc=1_idef,pchunk

       ne(cc)=-1._rp
       Te(cc)=-1._rp
       Zeff(cc)=-1._rp
       
       Bmag(cc) = SQRT(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))
       
       bhat_R(cc) = B_R(cc)/Bmag(cc)
       bhat_PHI(cc) = B_PHI(cc)/Bmag(cc)
       bhat_Z(cc) = B_Z(cc)/Bmag(cc)

       Bst_R(cc)=q_cache*B_R(cc)+V_PLL(cc)*curlb_R(cc)
       Bst_PHI(cc)=q_cache*B_PHI(cc)+V_PLL(cc)*curlb_PHI(cc)
       Bst_Z(cc)=q_cache*B_Z(cc)+V_PLL(cc)*curlb_Z(cc)

      ! write(output_unit_write,*) 'bmag',Bmag(cc),'bhat',bhat_R(cc),bhat_PHI(cc),bhat_Z(cc),'Bst',Bst_R(cc),Bst_PHI(cc),Bst_Z(cc)
       
       bdotBst(cc)=bhat_R(cc)*Bst_R(cc)+bhat_PHI(cc)*Bst_PHI(cc)+ &
            bhat_Z(cc)*Bst_Z(cc)
       BstdotE(cc)=Bst_R(cc)*E_R(cc)+Bst_PHI(cc)*E_PHI(cc)+Bst_Z(cc)*E_Z(cc)   
       BstdotgradB(cc)=Bst_R(cc)*gradB_R(cc)+Bst_PHI(cc)*gradB_PHI(cc)+ &
            Bst_Z(cc)*gradB_Z(cc)

       !write(output_unit_write,*) 'bdotBst',bdotBst(cc),BstdotE(cc),BstdotgradB(cc)
       
       Ecrossb_R(cc)=E_PHI(cc)*bhat_Z(cc)-E_Z(cc)*bhat_PHI(cc)
       Ecrossb_PHI(cc)=E_Z(cc)*bhat_R(cc)-E_R(cc)*bhat_Z(cc)
       Ecrossb_Z(cc)=E_R(cc)*bhat_PHI(cc)-E_PHI(cc)*bhat_R(cc)

       !write(output_unit_write,*) 'Ecrossb',Ecrossb_R(cc),Ecrossb_PHI(cc),Ecrossb_Z(cc)
       
       bcrossgradB_R(cc)=bhat_PHI(cc)*gradB_Z(cc)-bhat_Z(cc)*gradB_PHI(cc)
       bcrossgradB_PHI(cc)=bhat_Z(cc)*gradB_R(cc)-bhat_R(cc)*gradB_Z(cc)
       bcrossgradB_Z(cc)=bhat_R(cc)*gradB_PHI(cc)-bhat_PHI(cc)*gradB_R(cc)

      ! write(output_unit_write,*) 'bcrossgradB',bcrossgradB_R(cc),bcrossgradB_PHI(cc),bcrossgradB_Z(cc)
       
       gamgc(cc)=sqrt(1+V_PLL(cc)*V_PLL(cc)+2*V_MU(cc)*Bmag(cc))
       
       pm(cc)=sqrt(gamgc(cc)**2-1)
       xi(cc)=V_PLL(cc)/pm(cc)

       RHS_R(cc)=(q_cache*Ecrossb_R(cc)+(m_cache*V_MU(cc)* &
            bcrossgradB_R(cc)+V_PLL(cc)*Bst_R(cc))/(m_cache*gamgc(cc)))/ &
            bdotBst(cc)
       RHS_PHI(cc)=(q_cache*Ecrossb_PHI(cc)+(m_cache*V_MU(cc)* &
            bcrossgradB_PHI(cc)+V_PLL(cc)*Bst_PHI(cc))/(m_cache*gamgc(cc)))/ &
            (Y_R(cc)*bdotBst(cc))
       RHS_Z(cc)=(q_cache*Ecrossb_Z(cc)+(m_cache*V_MU(cc)* &
            bcrossgradB_Z(cc)+V_PLL(cc)*Bst_Z(cc))/(m_cache*gamgc(cc)))/ &
            bdotBst(cc)
       RHS_PLL(cc)=(q_cache*BstdotE(cc)-V_MU(cc)*BstdotgradB(cc)/gamgc(cc))/ &
            bdotBst(cc)
       RHS_MU(cc)=0._rp

    end do
    !$OMP END SIMD

    !write(output_unit_write,*) 'bmag',Bmag(1),'bhat',bhat_R(1),bhat_PHI(1),bhat_Z(1),'Bst',Bst_R(1),Bst_PHI(1),Bst_Z(1)
    !write(output_unit_write,*) 'bdotBst',bdotBst(1),BstdotE(1),BstdotgradB(1)
    !write(output_unit_write,*) 'Ecrossb',Ecrossb_R(1),Ecrossb_PHI(1),Ecrossb_Z(1)
    !write(output_unit_write,*) 'bcrossgradB',bcrossgradB_R(1),bcrossgradB_PHI(1),bcrossgradB_Z(1)

    if (params%radiation.and.(params%GC_rad_model.eq.'SDE')) then

       !       write(output_unit_write,*) 'RHS_PLL',RHS_PLL(1)

       re_cache=C_RE/params%cpp%length
       alpha_cache=C_a

       if (.not.params%LargeCollisions) then
          time=params%init_time+(params%it-1+tt)*params%dt
       else
          time=params%init_time+params%it*params%dt+ &
               tt*params%coll_per_dump_dt
       end if

       if (params%profile_model(1:10).eq.'ANALYTICAL') then
          call analytical_profiles_p(pchunk,time,params,Y_R,Y_Z,P,F,ne,Te,Zeff,PSIp)
       else
#ifdef PSPLINE
          call interp_FOcollision_p(pchunk,Y_R,Y_PHI,Y_Z,ne,Te,Zeff,flag_cache)
#endif
       endif

       !write(6,*) 'Y: ',Y_R,Y_PHI,Y_Z
       !write(6,*) 'Profs: ',ne,Te,Zeff
       
       
       !$OMP SIMD
       !       !$OMP& aligned(tau_R,Bmag,RHS_PLL,V_PLL,xi,gamgc,RHS_MU,V_MU)
       do cc=1_idef,pchunk

          tau_R(cc)=6*C_PI*E0/(Bmag(cc)*Bmag(cc))

          SR_PLL(cc)=V_PLL(cc)*(1._rp-xi(cc)*xi(cc))/tau_R(cc)* &
               (1._rp/gamgc(cc)-gamgc(cc))
          SR_MU(cc)=-2._rp*V_MU(cc)/tau_R(cc)* &
               (gamgc(cc)*(1-xi(cc)*xi(cc))+xi(cc)*xi(cc)/gamgc(cc))

          !Normalizations done here
          BREM_P(cc)=-4._rp*re_cache**2*ne(cc)* &
               Zeff(cc)*(Zeff(cc)+1._rp)*alpha_cache* &
               (gamgc(cc)-1._rp)*(log(2._rp*gamgc(cc))-1._rp/3._rp)
          BREM_PLL(cc)=xi(cc)*BREM_P(cc)
          BREM_MU(cc)=(1._rp-xi(cc)*xi(cc))*V_PLL(cc)/ &
               (Bmag(cc)*xi(cc))*BREM_P(cc)

          RHS_PLL(cc)=RHS_PLL(cc)+SR_PLL(cc)+BREM_PLL(cc)
          RHS_MU(cc)=SR_MU(cc)+BREM_MU(cc)

       end do
       !$OMP END SIMD

    end if   

    
#if DBG_CHECK    
    !$OMP SIMD
    do cc=1_idef,pchunk 
       if(isnan(gamgc(cc)).and.flag_cache(cc)==1) then           
          write(6,*) 'gamgc is NaN'
          write(6,*) 'V_PLL',V_PLL(cc)
          write(6,*) 'V_MU',V_MU(cc)
          write(6,*) 'Bmag',Bmag(cc)
          
          stop 'gamgc is a NaN'
       endif
       if(isnan(RHS_R(cc)).and.flag_cache(cc)==1) then  

          write(6,*) thread_num,'flags',flag_cache(cc)
          write(6,*) thread_num,'Y',Y_R(cc)*params%cpp%length,Y_PHI(cc),Y_Z(cc)*params%cpp%length
          write(6,*) thread_num,'B',B_R(cc),B_PHI(cc),B_Z(cc)
          write(6,*) thread_num,'E',E_R(cc),E_PHI(cc),E_Z(cc)
          write(6,*) thread_num,'gradB',gradB_R(cc),gradB_PHI(cc),gradB_Z(cc)
          write(6,*) thread_num,'curlb',curlb_R(cc),curlb_PHI(cc),curlb_Z(cc)
          write(6,*) thread_num,'V',V_PLL(cc),V_MU(cc)
          write(6,*) thread_num,'Exb',Ecrossb_R(cc)
          write(6,*) thread_num,'bxgradB',bcrossgradB_R(cc)
          write(6,*) thread_num,'Bst',Bst_R(cc)
          write(6,*) thread_num,'bdotBst',bdotBst(cc)
          write(6,*) thread_num,'gamma',gamgc(cc)
          
          stop 'RHS_R1 is a NaN'
       endif
       if(isnan(RHS_PHI(cc)).and.flag_cache(cc)==1) then
          stop 'RHS_PHI1 is a NaN'
       end if
       if(isnan(RHS_Z(cc)).and. flag_cache(cc)==1) then
          stop 'RHS_Z1 is a NaN'
       endif
       if(isnan(RHS_PLL(cc)).and. flag_cache(cc)==1) then

          write(6,*) params%mpi_params%rank,thread_num,'flags',flag_cache(cc)
          write(6,*) params%mpi_params%rank,thread_num,'Y',Y_R(cc)*params%cpp%length,Y_PHI(cc),Y_Z(cc)*params%cpp%length
          write(6,*) params%mpi_params%rank,thread_num,'V',V_PLL(cc),V_MU(cc)
          write(6,*) params%mpi_params%rank,thread_num,'gamma',gamgc(cc)
          write(6,*) params%mpi_params%rank,thread_num,'xi',xi(cc)
          stop 'RHS_PLL1 is a NaN'
       endif
       if(isnan(RHS_MU(cc)).and. flag_cache(cc)==1) then
          stop 'RHS_MU1 is a NaN'
       endif
    end do
    !$OMP END SIMD

    if (params%radiation.and.(params%GC_rad_model.eq.'SDE')) then
       !$OMP SIMD
       do cc=1_idef,pchunk
          if(isnan(ne(cc)).and. flag_cache(cc)==1) then
             stop 'ne is a NaN'
          endif
          if(isnan(Zeff(cc)).and. flag_cache(cc)==1) then
             stop 'Zeff is a NaN'
          endif
          if(isnan(BREM_P(cc)).and. flag_cache(cc)==1) then
             stop 'BREM_P is a NaN'
          endif
          if(isnan(BREM_PLL(cc)).and. flag_cache(cc)==1) then
             stop 'BREM_PLL is a NaN'
          endif
          if(isnan(BREM_MU(cc)).and. flag_cache(cc)==1) then
             stop 'BREM_MU is a NaN'
          endif
          if(isnan(SR_PLL(cc)).and. flag_cache(cc)==1) then
             stop 'SR_PLL is a NaN'
          endif
          if(isnan(SR_MU(cc)).and. flag_cache(cc)==1) then
             stop 'SR_MU is a NaN'
          endif
       end do
       !$OMP END SIMD
    end if
#endif
    
    !    write(output_unit_write,*) 'RHS_R: ',RHS_R(1)
    !    write(output_unit_write,*) 'RHS_PHI: ',RHS_PHI(1)
    !    write(output_unit_write,*) 'RHS_Z: ',RHS_Z(1)
    !    write(output_unit_write,*) 'RHS_PLL: ',RHS_PLL(1)
    !    write(output_unit_write,*) 'RHS_MU: ',RHS_MU(1)

  end subroutine GCEoM1_p

#ifdef FIO
  subroutine GCEoM1_fio_p(tt,P,F,params,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
       B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z, &
       gradB_R,gradB_PHI,gradB_Z,V_PLL,V_MU,Y_R,Y_PHI,Y_Z,q_cache,m_cache,PSIp, &
       ne,ni,nimp,Te,Zeff,flagCon,flagCol,hint)

    USE omp_lib
    IMPLICIT NONE
    
    TYPE(KORC_PARAMS), INTENT(INOUT)                           :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)      :: F
    TYPE(PROFILES), INTENT(IN)                                 :: P
    REAL(rp),DIMENSION(params%pchunk)  :: Bmag,bhat_R,bhat_PHI,bhat_Z,Bst_R,Bst_PHI
    REAL(rp),DIMENSION(params%pchunk)  :: BstdotE,BstdotgradB,EcrossB_R,EcrossB_PHI,bdotBst
    REAL(rp),DIMENSION(params%pchunk)  :: bcrossgradB_R,bcrossgradB_PHI,bcrossgradB_Z,gamgc
    REAL(rp),DIMENSION(params%pchunk)  :: EcrossB_Z,Bst_Z
    REAL(rp),DIMENSION(params%pchunk)  :: pm,xi,tau_R
    REAL(rp),DIMENSION(params%pchunk)  :: SR_PLL,SR_MU,BREM_PLL,BREM_MU,BREM_P
    REAL(rp),DIMENSION(params%pchunk),INTENT(in) :: gradB_R,gradB_PHI,gradB_Z,curlb_R
    REAL(rp),DIMENSION(params%pchunk),INTENT(in) :: curlb_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: RHS_R,RHS_PHI,RHS_Z
    REAL(rp),DIMENSION(params%pchunk),INTENT(OUT) :: RHS_PLL,RHS_MU
    REAL(rp),DIMENSION(params%pchunk),INTENT(IN) :: V_PLL,V_MU,Y_R,Y_PHI,Y_Z,curlb_PHI
    REAL(rp),DIMENSION(params%pchunk),INTENT(IN) :: PSIp
    REAL(rp),INTENT(in) :: q_cache,m_cache
    INTEGER(ip)  :: cc,pchunk,ii
    INTEGER(ip),INTENT(IN)  :: tt
    REAL(rp)  :: time,re_cache,alpha_cache
    REAL(rp), DIMENSION(params%pchunk),INTENT(OUT) 		:: ne,Te,ni,nimp,Zeff
    TYPE(C_PTR), DIMENSION(params%pchunk), INTENT(INOUT)  :: hint
    INTEGER(is),DIMENSION(params%pchunk),intent(INOUT) :: flagCon,flagCol
    INTEGER             :: thread_num

    pchunk=params%pchunk

    thread_num = OMP_GET_THREAD_NUM()
    
!    !$OMP SIMD
!    do cc=1_idef,pchunk 
!       if(isnan(B_R(cc))) then
!          write(6,*) thread_num,'0Y',Y_R(cc)*params%cpp%length,Y_PHI(cc),Y_Z(cc)*params%cpp%length
!          write(6,*) thread_num,'0B',B_R(cc),B_PHI(cc),B_Z(cc)
!          write(6,*) thread_num,'0E',E_R(cc),E_PHI(cc),E_Z(cc)
!          write(6,*) thread_num,'0gradB',gradB_R(cc),gradB_PHI(cc),gradB_Z(cc)
!          write(6,*) thread_num,'0curlb',curlb_R(cc),curlb_PHI(cc),curlb_Z(cc)
!          write(6,*) thread_num,'0V',V_PLL(cc),V_MU(cc)
!          write(6,*) 'Exb',Ecrossb_R(cc)
!          write(6,*) 'bxgradB',bcrossgradB_R(cc)
!          write(6,*) 'Bst',Bst_R(cc)
!          write(6,*) 'bdotBst',bdotBst(cc)
!          write(6,*) thread_num,'0gamma',gamgc(cc)

!          stop 'B_R is a NaN'
!       endif
!    end do
!    !$OMP END SIMD

    
    !$OMP SIMD
    !    !$OMP& aligned(gradB_R,gradB_PHI,gradB_Z,curlb_R,curlb_Z, &
    !    !$OMP& B_R,B_PHI,B_Z,E_R,E_PHI,E_Z,RHS_R,RHS_PHI,RHS_Z,RHS_PLL,RHS_MU, &
    !    !$OMP& V_PLL,V_MU,Y_R,curlb_PHI,tau_R)
    do cc=1_idef,pchunk

       ne(cc)=1._rp
       Te(cc)=1._rp
       Zeff(cc)=1._rp
       
       Bmag(cc) = SQRT(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))
       
       bhat_R(cc) = B_R(cc)/Bmag(cc)
       bhat_PHI(cc) = B_PHI(cc)/Bmag(cc)
       bhat_Z(cc) = B_Z(cc)/Bmag(cc)

       Bst_R(cc)=q_cache*B_R(cc)+V_PLL(cc)*curlb_R(cc)
       Bst_PHI(cc)=q_cache*B_PHI(cc)+V_PLL(cc)*curlb_PHI(cc)
       Bst_Z(cc)=q_cache*B_Z(cc)+V_PLL(cc)*curlb_Z(cc)

      ! write(output_unit_write,*) 'bmag',Bmag(cc),'bhat',bhat_R(cc),bhat_PHI(cc),bhat_Z(cc),'Bst',Bst_R(cc),Bst_PHI(cc),Bst_Z(cc)
       
       bdotBst(cc)=bhat_R(cc)*Bst_R(cc)+bhat_PHI(cc)*Bst_PHI(cc)+ &
            bhat_Z(cc)*Bst_Z(cc)
       BstdotE(cc)=Bst_R(cc)*E_R(cc)+Bst_PHI(cc)*E_PHI(cc)+Bst_Z(cc)*E_Z(cc)   
       BstdotgradB(cc)=Bst_R(cc)*gradB_R(cc)+Bst_PHI(cc)*gradB_PHI(cc)+ &
            Bst_Z(cc)*gradB_Z(cc)

       !write(output_unit_write,*) 'bdotBst',bdotBst(cc),BstdotE(cc),BstdotgradB(cc)
       
       Ecrossb_R(cc)=E_PHI(cc)*bhat_Z(cc)-E_Z(cc)*bhat_PHI(cc)
       Ecrossb_PHI(cc)=E_Z(cc)*bhat_R(cc)-E_R(cc)*bhat_Z(cc)
       Ecrossb_Z(cc)=E_R(cc)*bhat_PHI(cc)-E_PHI(cc)*bhat_R(cc)

       !write(output_unit_write,*) 'Ecrossb',Ecrossb_R(cc),Ecrossb_PHI(cc),Ecrossb_Z(cc)
       
       bcrossgradB_R(cc)=bhat_PHI(cc)*gradB_Z(cc)-bhat_Z(cc)*gradB_PHI(cc)
       bcrossgradB_PHI(cc)=bhat_Z(cc)*gradB_R(cc)-bhat_R(cc)*gradB_Z(cc)
       bcrossgradB_Z(cc)=bhat_R(cc)*gradB_PHI(cc)-bhat_PHI(cc)*gradB_R(cc)

      ! write(output_unit_write,*) 'bcrossgradB',bcrossgradB_R(cc),bcrossgradB_PHI(cc),bcrossgradB_Z(cc)
       
       gamgc(cc)=sqrt(1+V_PLL(cc)*V_PLL(cc)+2*V_MU(cc)*Bmag(cc))
       
       pm(cc)=sqrt(gamgc(cc)**2-1)
       xi(cc)=V_PLL(cc)/pm(cc)

       RHS_R(cc)=(q_cache*Ecrossb_R(cc)+(m_cache*V_MU(cc)* &
            bcrossgradB_R(cc)+V_PLL(cc)*Bst_R(cc))/(m_cache*gamgc(cc)))/ &
            bdotBst(cc)
       RHS_PHI(cc)=(q_cache*Ecrossb_PHI(cc)+(m_cache*V_MU(cc)* &
            bcrossgradB_PHI(cc)+V_PLL(cc)*Bst_PHI(cc))/(m_cache*gamgc(cc)))/ &
            (Y_R(cc)*bdotBst(cc))
       RHS_Z(cc)=(q_cache*Ecrossb_Z(cc)+(m_cache*V_MU(cc)* &
            bcrossgradB_Z(cc)+V_PLL(cc)*Bst_Z(cc))/(m_cache*gamgc(cc)))/ &
            bdotBst(cc)
       RHS_PLL(cc)=(q_cache*BstdotE(cc)-V_MU(cc)*BstdotgradB(cc)/gamgc(cc))/ &
            bdotBst(cc)
       RHS_MU(cc)=0._rp

    end do
    !$OMP END SIMD

    !write(output_unit_write,*) 'bmag',Bmag(1),'bhat',bhat_R(1),bhat_PHI(1),bhat_Z(1),'Bst',Bst_R(1),Bst_PHI(1),Bst_Z(1)
    !write(output_unit_write,*) 'bdotBst',bdotBst(1),BstdotE(1),BstdotgradB(1)
    !write(output_unit_write,*) 'Ecrossb',Ecrossb_R(1),Ecrossb_PHI(1),Ecrossb_Z(1)
    !write(output_unit_write,*) 'bcrossgradB',bcrossgradB_R(1),bcrossgradB_PHI(1),bcrossgradB_Z(1)

    !    !$OMP SIMD
    !    do cc=1_idef,8
    !       if(isnan(RHS_R(cc))) stop 'RHS_R0 is a NaN'
    !       if(isnan(RHS_PHI(cc))) stop 'RHS_PHI0 is a NaN'
    !       if(isnan(RHS_Z(cc))) stop 'RHS_Z0 is a NaN'
    !       if(isnan(RHS_PLL(cc))) stop 'RHS_PLL0 is a NaN'
    !       if(isnan(RHS_MU(cc))) stop 'RHS_MU0 is a NaN'
    !    end do
    !    !$OMP END SIMD

    if (params%radiation.and.(params%GC_rad_model.eq.'SDE')) then

       !       write(output_unit_write,*) 'RHS_PLL',RHS_PLL(1)

       re_cache=C_RE/params%cpp%length
       alpha_cache=C_a

       call get_fio_profile_p(params,P,Y_R,Y_PHI,Y_Z, &
            ne,Te,flagCon,hint)

       call get_fio_ion_p(params,P,Y_R,Y_PHI,Y_Z, &
            ne,ni,nimp,Zeff,flagCon,hint)
          
       !$OMP SIMD
       !       !$OMP& aligned(tau_R,Bmag,RHS_PLL,V_PLL,xi,gamgc,RHS_MU,V_MU)
       do cc=1_idef,pchunk

          tau_R(cc)=6*C_PI*E0/(Bmag(cc)*Bmag(cc))

          SR_PLL(cc)=V_PLL(cc)*(1._rp-xi(cc)*xi(cc))/tau_R(cc)* &
               (1._rp/gamgc(cc)-gamgc(cc))
          SR_MU(cc)=-2._rp*V_MU(cc)/tau_R(cc)* &
               (gamgc(cc)*(1-xi(cc)*xi(cc))+xi(cc)*xi(cc)/gamgc(cc))
          
          !Normalizations done here
          BREM_P(cc)=-4._rp*re_cache**2*ne(cc)* &
               Zeff(cc)*(Zeff(cc)+1._rp)*alpha_cache* &
               (gamgc(cc)-1._rp)*(log(2._rp*gamgc(cc))-1._rp/3._rp)
          BREM_PLL(cc)=xi(cc)*BREM_P(cc)
          BREM_MU(cc)=(1._rp-xi(cc)*xi(cc))*V_PLL(cc)/ &
               (Bmag(cc)*xi(cc))*BREM_P(cc)

          RHS_PLL(cc)=RHS_PLL(cc)+SR_PLL(cc)+BREM_PLL(cc)
          RHS_MU(cc)=SR_MU(cc)+BREM_MU(cc)

       end do
       !$OMP END SIMD


    end if

#if DBG_CHECK    
    !$OMP SIMD
    do cc=1_idef,pchunk 
       if(isnan(gamgc(cc)).and. &
            ((flagCon(cc)==1).and.(flagCol(cc)==1))) then           
          write(6,*) 'gamgc is NaN'
          write(6,*) 'V_PLL',V_PLL(cc)
          write(6,*) 'V_MU',V_MU(cc)
          write(6,*) 'Bmag',Bmag(cc)
          
          stop 'gamgc is a NaN'
       endif
       if(isnan(RHS_R(cc)).and. &
            ((flagCon(cc)==1).and.(flagCol(cc)==1))) then  

          write(6,*) thread_num,'flags',flagCon(cc),flagCol(cc)
          write(6,*) thread_num,'Y',Y_R(cc)*params%cpp%length,Y_PHI(cc),Y_Z(cc)*params%cpp%length
          write(6,*) thread_num,'B',B_R(cc),B_PHI(cc),B_Z(cc)
          write(6,*) thread_num,'E',E_R(cc),E_PHI(cc),E_Z(cc)
          write(6,*) thread_num,'gradB',gradB_R(cc),gradB_PHI(cc),gradB_Z(cc)
          write(6,*) thread_num,'curlb',curlb_R(cc),curlb_PHI(cc),curlb_Z(cc)
          write(6,*) thread_num,'V',V_PLL(cc),V_MU(cc)
          write(6,*) thread_num,'Exb',Ecrossb_R(cc)
          write(6,*) thread_num,'bxgradB',bcrossgradB_R(cc)
          write(6,*) thread_num,'Bst',Bst_R(cc)
          write(6,*) thread_num,'bdotBst',bdotBst(cc)
          write(6,*) thread_num,'gamma',gamgc(cc)
          
          stop 'RHS_R1 is a NaN'
       endif
       if(isnan(RHS_PHI(cc)).and. &
            ((flagCon(cc)==1).and.(flagCol(cc)==1))) then
          stop 'RHS_PHI1 is a NaN'
       end if
       if(isnan(RHS_Z(cc)).and. &
            ((flagCon(cc)==1).and.(flagCol(cc)==1))) then
          stop 'RHS_Z1 is a NaN'
       endif
       if(isnan(RHS_PLL(cc)).and. &
            ((flagCon(cc)==1).and.(flagCol(cc)==1))) then
          stop 'RHS_PLL1 is a NaN'
       endif
       if(isnan(RHS_MU(cc)).and. &
            ((flagCon(cc)==1).and.(flagCol(cc)==1))) then
          stop 'RHS_MU1 is a NaN'
       endif
    end do
    !$OMP END SIMD

    if (params%radiation.and.(params%GC_rad_model.eq.'SDE')) then
       !$OMP SIMD
       do cc=1_idef,pchunk
          if(isnan(ne(cc)).and. &
               ((flagCon(cc)==1).and.(flagCol(cc)==1))) then
             stop 'ne is a NaN'
          endif
          if(isnan(Zeff(cc)).and. &
               ((flagCon(cc)==1).and.(flagCol(cc)==1))) then
             stop 'Zeff is a NaN'
          endif
          if(isnan(BREM_P(cc)).and. &
               ((flagCon(cc)==1).and.(flagCol(cc)==1))) then
             stop 'BREM_P is a NaN'
          endif
          if(isnan(BREM_PLL(cc)).and. &
               ((flagCon(cc)==1).and.(flagCol(cc)==1))) then
             stop 'BREM_PLL is a NaN'
          endif
          if(isnan(BREM_MU(cc)).and. &
               ((flagCon(cc)==1).and.(flagCol(cc)==1))) then
             stop 'BREM_MU is a NaN'
          endif
          if(isnan(SR_PLL(cc)).and. &
               ((flagCon(cc)==1).and.(flagCol(cc)==1))) then
             stop 'SR_PLL is a NaN'
          endif
          if(isnan(SR_MU(cc)).and. &
               ((flagCon(cc)==1).and.(flagCol(cc)==1))) then
             stop 'SR_MU is a NaN'
          endif
       end do
       !$OMP END SIMD
    end if
#endif
    !    write(output_unit_write,*) 'RHS_R: ',RHS_R(1)
    !    write(output_unit_write,*) 'RHS_PHI: ',RHS_PHI(1)
    !    write(output_unit_write,*) 'RHS_Z: ',RHS_Z(1)
    !    write(output_unit_write,*) 'RHS_PLL: ',RHS_PLL(1)
    !    write(output_unit_write,*) 'RHS_MU: ',RHS_MU(1)

  end subroutine GCEoM1_fio_p
#endif
  
  subroutine aux_fields(pp,spp,gradB,curlb,Bmag)
    TYPE(SPECIES), INTENT(IN)    :: spp
    !! An instance of the derived type SPECIES containing all the parameters
    !! and simulation variables of the different species in the simulation.
    REAL(rp),DIMENSION(3),INTENT(INOUT) :: gradB
    REAL(rp),DIMENSION(3),INTENT(INOUT) :: curlb
    REAL(rp),INTENT(IN) :: Bmag
    REAL(rp) :: dRB
    REAL(rp) :: dPHIB
    REAL(rp) :: dZB
    INTEGER  :: pp

    dRB=(spp%vars%B(pp,1)*spp%vars%BR(pp,1)+ &
         spp%vars%B(pp,2)*spp%vars%BPHI(pp,1)+ &
         spp%vars%B(pp,3)*spp%vars%BZ(pp,1))/Bmag
    dPHIB=(spp%vars%B(pp,1)*spp%vars%BR(pp,2)+ &
         spp%vars%B(pp,2)*spp%vars%BPHI(pp,2)+ &
         spp%vars%B(pp,3)*spp%vars%BZ(pp,2))/Bmag
    dZB=(spp%vars%B(pp,1)*spp%vars%BR(pp,3)+ &
         spp%vars%B(pp,2)*spp%vars%BPHI(pp,3)+ &
         spp%vars%B(pp,3)*spp%vars%BZ(pp,3))/Bmag

    gradB(1)=dRB
    gradB(2)=dPHIB/spp%vars%Y(pp,1)
    gradB(3)=dZB

    curlb(1)=((Bmag*spp%vars%BZ(pp,2)-spp%vars%B(pp,3)*dPHIB)/spp%vars%Y(pp,1)- &
         (Bmag*spp%vars%BPHI(pp,3)-spp%vars%B(pp,2)*dZB))/Bmag**2
    curlb(2)=((Bmag*spp%vars%BR(pp,3)-spp%vars%B(pp,1)*dZB)- &
         (Bmag*spp%vars%BZ(pp,1)-spp%vars%B(pp,3)*dRB))/Bmag**2
    curlb(3)=((Bmag*spp%vars%BPHI(pp,1)-spp%vars%B(pp,2)*dRB) - &
         (Bmag*spp%vars%BPHI(pp,1)-spp%vars%B(pp,1)*dPHIB)/ &
         spp%vars%Y(pp,1))/Bmag**2+ &
         spp%vars%B(pp,2)/(Bmag*spp%vars%Y(pp,1))  

  end subroutine aux_fields


end module korc_ppusher
