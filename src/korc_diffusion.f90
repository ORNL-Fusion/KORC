module korc_diffusion
  use korc_types
  use korc_constants
  use korc_HDF5
  use korc_interp
  use korc_profiles
  use korc_fields
  use korc_input
  use korc_random

#ifdef __NVCOMPILER
  use ieee_arithmetic
#endif

  IMPLICIT NONE

  PUBLIC :: include_CoulombCollisions_FO_p,&
       define_collisions_time_step

  contains

  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !
  ! * SUBROUTINES FOR INITIALIZING COLLISIONS PARAMS * !
  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !

  subroutine define_collisions_time_step(params,params_ACC,F,init)
    TYPE(KORC_PARAMS), INTENT(INOUT) 	:: params
    TYPE(KORC_PARAMS_ACC) :: params_ACC
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

       if (params%snapshot_frequency.gt.0._rp) then

          !write(6,*) 'params%snapshot_frequency',params%snapshot_frequency*params%cpp%time
          !write(6,*) 'cparams_ss%dTau*Tau',cparams_ss%dTau*Tau*params%cpp%time
          !write(6,*) 'FLOOR(params%snapshot_frequency/cparams_ss%dTau*Tau)', &
          !     FLOOR(params%snapshot_frequency/ &
          !     (cparams_ss%dTau*Tau),ip)

          params%coll_per_dump=ceiling(params%snapshot_frequency/ &
               (cparams_ss%dTau*Tau))

          cparams_ss%coll_per_dump_dt=params%snapshot_frequency/params%coll_per_dump

          params%coll_per_dump_dt=cparams_ss%coll_per_dump_dt
          params_ACC%coll_per_dump_dt=cparams_ss%coll_per_dump_dt

          params%orbits_per_coll=ceiling(cparams_ss%coll_per_dump_dt/ &
               params%dt)


#ifdef __NVCOMPILER
          params%dt=cparams_ss%coll_per_dump_dt/real(params%orbits_per_coll)
#else
          params%dt=cparams_ss%coll_per_dump_dt/float(params%orbits_per_coll)
#endif     
          

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

          !if (params%LargeCollisions) then

             write(output_unit_write,'("Number of collision steps per dump step: ",I16)') params%coll_per_dump

             write(output_unit_write,'("Collision time step: ",E17.10)') params%coll_per_dump_dt*params%cpp%time

             write(output_unit_write,'("Number of orbit steps per collision step: ",I16)') params%orbits_per_coll

             write(output_unit_write,'("Orbit time step: ",E17.10)') params%dt*params%cpp%time

          !end if


          write(output_unit_write,'("* * * * * * * * * * * * * * * * * * * * &
               * * * * * * * * * * * * * * *",/)')
       end if

        if (params%mpi_params%rank .EQ. 0) then
          flush(output_unit_write)
        end if

        if (params%coll_per_dump.gt.params%t_skip) then
          write(6,*) 'more collisional iterations than orbit iterations, decrease orbit timestep!'
          call korc_abort(26)
        endif

    else if (params%orbit_model(1:2).eq.'GC'.and.params%field_eval.eq.'eqn' &
      .and..not.params%field_model.eq.'M3D_C1') then

      params%coll_per_dump=1   
      params_ACC%coll_per_dump=1 

      params%orbits_per_coll=params%t_skip
      params_ACC%orbits_per_coll=params%t_skip

    end if
  end subroutine define_collisions_time_step

  subroutine include_CoulombCollisions_FO_p(tt,params,random,X_X,X_Y,X_Z, &
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
    CLASS(random_context), POINTER, INTENT(INOUT) :: random
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

       CALL random%uniform%set(0.0_rp, 1.0_rp)

       !$OMP SIMD
   !       !$OMP& aligned(rnd1,dW,CAL,dCAL,CFL,CBL,vm,ne,Te,Zeff,dpm, &
   !       !$OMP& flagCon,flagCol,dxi,xi,pm,dphi,um,Ub_X,Ub_Y,Ub_Z,U_X,U_Y,U_Z, &
   !       !$OMP& b1_X,b1_Y,b1_Z,b2_X,b2_Y,b2_Z,b3_X,b3_Y,b3_Z)
       do cc=1_idef,pchunk

          ! uses C library to generate normal_distribution random variables,
          ! preserving parallelization where Fortran random number generator
          ! does not
          rnd1(cc,1) = random%uniform%get()
          rnd1(cc,2) = random%uniform%get()
          rnd1(cc,3) = random%uniform%get()

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

end module korc_collisions
