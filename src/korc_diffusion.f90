module korc_diffusion
  use korc_types
  use korc_constants
  use korc_input
  use korc_random
  use korc_coords

#ifdef __NVCOMPILER
  use ieee_arithmetic
#endif

  IMPLICIT NONE

  PUBLIC :: include_diffusion_p,&
       define_diffusion_time_step

  contains

  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !
  ! * SUBROUTINES FOR INITIALIZING COLLISIONS PARAMS * !
  ! * * * * * * * * * * * *  * * * * * * * * * * * * * !

  subroutine define_diffusion_time_step(spp,params,params_ACC,F,init)
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)       :: spp
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


    if (params%diffusion) then
      E = 1. + spp(1)%Eo
      v = SQRT(1.0_rp - (1./E)**2)
      vpll=v*cos(C_PI*spp(1)%etao/180.)

      q=F%AB%qo*(1+(F%AB%mu_turb/F%AB%lambda)**2)
      Lc=C_PI*F%Ro*q

      DRR=abs(vpll)*Lc*F%AB%A_turb
      dDRR=abs(vpll)**(Lc+C_PI*F%Ro*F%AB%qo*2*F%AB%mu_turb/F%AB%lambda**2*F%AB%A_turb)

      mu=dDRR+DRR/F%AB%mu_turb

      nu=(/mu/F%AB%mu_turb,2*DRR/F%AB%mu_turb**2/)

      Tau = MINVAL( 1.0_rp/nu )

      params%diff_subcycling_iterations = ceiling(0.05*Tau/params%dt,ip)

      if (params%snapshot_frequency.gt.0._rp) then

        params%diff_per_dump=ceiling(params%snapshot_frequency/ &
            (0.05*Tau))

        params%diff_per_dump_dt=params%snapshot_frequency/params%diff_per_dump
        params_ACC%diff_per_dump_dt=params%diff_per_dump_dt

        params%orbits_per_diff=ceiling(params%diff_per_dump_dt/ &
              params%dt)

#ifdef __NVCOMPILER
        params%dt=params%diff_per_dump_dt/real(params%orbits_per_diff)
#else
        params%dt=params%diff_per_dump_dt/float(params%orbits_per_diff)
#endif         

      end if

      if (init) num_diffusion_in_simulation = params%simulation_time/Tau

      if (params%mpi_params%rank .EQ. 0) then
        write(output_unit_write,'("* * * * * * * * * * * SUBCYCLING FOR  &
          DIFFUSION * * * * * * * * * * *")')

        write(output_unit_write,'("Drift frequency: ",E17.10)') &
          nu(1)/params%cpp%time
        write(output_unit_write,'("RR diffusion: ",E17.10)') &
          nu(2)/params%cpp%time

        write(output_unit_write,'("The shorter diffusion time in the simulations  &
          is: ",E17.10," s")') Tau*params%cpp%time
        write(output_unit_write,'("Number of KORC iterations per collision: ",I16)')  &
          params%diff_subcycling_iterations
        if (init) then
          write(output_unit_write,'("Number of collisions in simulated time: ",E17.10)')  &
            num_diffusion_in_simulation
        end if

        write(output_unit_write,'("Number of diffusion steps per dump step: ",I16)') params%diff_per_dump

        write(output_unit_write,'("Diffusion time step: ",E17.10)') params%diff_per_dump_dt*params%cpp%time

        write(output_unit_write,'("Number of orbit steps per diffusion step: ",I16)') params%orbits_per_diff

        write(output_unit_write,'("Orbit time step: ",E17.10)') params%dt*params%cpp%time

        write(output_unit_write,'("* * * * * * * * * * * * * * * * * * * * &
          * * * * * * * * * * * * * * *",/)')
      end if

      if (params%mpi_params%rank .EQ. 0) then
        flush(output_unit_write)
      end if

      if (params%diff_per_dump.gt.params%t_skip) then
        write(6,*) 'more diffusion iterations than orbit iterations, decrease orbit timestep!'
        call korc_abort(26)
      endif

    endif
  
  end subroutine define_diffusion_time_step

  subroutine include_diffusion_p(tt,params,random,X_X,X_Y,X_Z, &
    U_X,U_Y,U_Z,B_X,B_Y,B_Z,me,F)
    !! This subroutine performs a Stochastic collision process consistent
    !! with the Fokker-Planck model for relativitic electron colliding with
    !! a thermal (Maxwellian) plasma. The collision operator is in spherical
    !! coordinates of the form found in Papp et al., NF (2011). CA
    !! corresponds to the parallel (speed diffusion) process, CF corresponds
    !! to a slowing down (momentum loss) process, and CB corresponds to a
    !! perpendicular diffusion process. Ordering of the processes are
    !! $$ \sqrt{CB}\gg CB \gg CF \sim \sqrt{CA} \gg CA,$$
    !! and only the dominant terms are kept.
    TYPE(FIELDS), INTENT(IN)      :: F
    TYPE(KORC_PARAMS), INTENT(IN) 		:: params
    CLASS(random_context), POINTER, INTENT(INOUT) :: random
    REAL(rp), DIMENSION(params%pchunk), INTENT(INOUT) 	:: X_X,X_Y,X_Z
    REAL(rp), DIMENSION(params%pchunk)  	:: T_R,T_T,T_Z
    REAL(rp), DIMENSION(params%pchunk), INTENT(IN) 	:: U_X,U_Y,U_Z
    REAL(rp), INTENT(IN)  :: me
    INTEGER(ip), INTENT(IN) 			:: tt
    REAL(rp), DIMENSION(params%pchunk), INTENT(IN) 		:: B_X,B_Y,B_Z

    REAL(rp), DIMENSION(params%pchunk) 		:: b_unit_X,b_unit_Y,b_unit_Z
    REAL(rp), DIMENSION(params%pchunk) 		:: Bmag
    REAL(rp), DIMENSION(params%pchunk) 			:: dW
    !! 3D Weiner process
    REAL(rp), DIMENSION(params%pchunk) 			:: rnd1
    REAL(rp) 					:: dt,time
    REAL(rp), DIMENSION(params%pchunk) 					:: um
    REAL(rp), DIMENSION(params%pchunk) 					:: dr
    REAL(rp), DIMENSION(params%pchunk) 					:: vm
    REAL(rp), DIMENSION(params%pchunk) 					:: pm
    REAL(rp), DIMENSION(params%pchunk) 			:: xi
    integer :: cc,pchunk
    INTEGER(ip) ::flagCon

    pchunk=params%pchunk

    if (MODULO(params%it+tt,params%diff_subcycling_iterations) .EQ. 0_ip) then
      dt = REAL(params%diff_subcycling_iterations,rp)*params%dt

      call cart_to_tor_check_if_confined_p(pchunk,F%AB%a,F%AB%Ro,F%AB%kappa, &
        X_X,X_Y,X_Z,T_R,T_T,T_Z,flagCon)

      CALL random%uniform%set(0.0_rp, 1.0_rp)

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

        rnd1(cc) = random%uniform%get()
        dW(cc) = SQRT(3*dt)*(-1+2*rnd1(cc))

        q(cc) = F%AB%qo*(1+(T_R(cc)/F%AB%lambda)**2)
        Lc(cc) = C_PI*F%AB%Ro*q(cc)
        dLc(cc) = C_PI*F%AB%Ro*F%AB%qo*2*T_R(cc)/F%AB%lambda**2

        g_r(cc) = F%AB%A_turb * exp(-0.5 * ((T_R(cc) - F%AB%mu_turb)/F%AB%sigma_turb)**2)
        ballooning(cc) = 0.25 * (1.0 + cos(T_T(cc)))**2
        dBr_norm_squared(cc)=g_r(cc)*ballooning(cc)
        ddBr_norm_squared(cc)=dBr_norm_squared(cc)*(-(T_R(cc) - F%AB%mu_turb)/F%AB%sigma_turb**2)

        DRR(cc)=abs(vm(cc)*xi(cc))*Lc(cc)*dBr_norm_squared(cc)
        dDRR(cc)=abs(vm(cc)*xi(cc))*(Lc(cc)*ddBr_norm_squared(cc) + &
          dLc(cc)*dBr_norm_squared(cc))

        dr(cc)=REAL(flagCon(cc))* &
              ((dDRR(cc)+DRR(cc)/T_R(cc))*dt+ &
              sqrt(2.0_rp*DRR(cc))*dW(cc))

        T_R(cc)=T_R(cc)+dr(cc)

      end do

      call tor_to_cart_p(pchunk,F%AB%Ro,F%AB%kappa, &
        X_X,X_Y,X_Z,T_R,T_T,T_Z)

    end if

  end subroutine include_diffusion_p

end module korc_diffusion
