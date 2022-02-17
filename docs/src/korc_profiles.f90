module korc_profiles
  !! @note Module that contain subroutines for calculating analytical plasma 
  !! profiles and calls to subroutines for interpolating external plasma 
  !! profiles to the particles positions. @endnote
  use korc_types
  use korc_hpc
  use korc_coords
  use korc_interp
  use korc_HDF5
  use korc_input

  IMPLICIT NONE

  PUBLIC :: get_profiles,&
       initialize_profiles,&
       analytical_profiles_p,&
       DEALLOCATE_PROFILES_ARRAYS
  PRIVATE :: get_analytical_profiles,&
       uniform_profiles,&
       load_profiles_data_from_hdf5,&
       ALLOCATE_2D_PROFILES_ARRAYS,&
       ALLOCATE_3D_PROFILES_ARRAYS

CONTAINS


  subroutine initialize_profiles(params,P,F)
    !! @note Subroutine that initializes the parameters of analytical 
    !! or pre-computed plasma profiles for being used in the
    !! simulation. @endnote
    !! KORC can run using either analytical and pre-computed plasma
    !! profiles. Pre-computed plasma profiles, as in the case of
    !! pre-computed electric or magnetic fields, are interpolated to
    !! electrons' position in [[korc_profiles]].
    !!
    !! There are two types of analytical plasma profiles that can be
    !! used in KORC: 3rd degree polynomial radial plasma profiles,
    !!
    !! $$f(r) = a_3r^3 + a_2r^2 +a_1r + a_0,$$
    !!
    !! and radial plasma profiles with a \(\tanh(r)\) dependency:
    !!
    !! $$f(r) = f_0\left[1 - \tanh^n\left(\frac{2r}{a}\right)\right]$$,
    !!
    !! where \(r\) is the radial coordinate in toroidal coordinates,
    !! \(f_0\) is a given plasma parameter at the magnetic axis,
    !! and \(a\) is the plasma radius as measured from the magnetic
    !! axis to the last closed flux surface. Notice that the larger
    !! \(n\) is, the more uniform the radial profiles are.

    TYPE(KORC_PARAMS), INTENT(IN)   :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(OUT)     :: P
    !! An instance of KORC's derived type PROFILES containing all 
    !! the information about the plasma profiles used in the
    !! simulation. See [[korc_types]] and [[korc_profiles]].
    TYPE(FIELDS), INTENT(IN)     :: F
    !CHARACTER(MAX_STRING_LENGTH)    :: ne_profile
    !! String containing the type of electron density profile 
    !! to be used in the simulation.
    !CHARACTER(MAX_STRING_LENGTH)    :: Te_profile
    !! String containing the type of electron temperature profile 
    !! to be used in the simulation.
    !CHARACTER(MAX_STRING_LENGTH)    :: Zeff_profile  
    !! String containing the type of \(Z_{eff}\) profile to be used
    !! in the simulation.
    !CHARACTER(MAX_STRING_LENGTH)    :: filename
    !! Full path to the HDF5 file containing the pre-computed
    !! plasma profiles.
    !REAL(rp)                        :: radius_profile
    !! Plasma radius \(a\) as measured from the magnetic axis.
    !REAL(rp)                        :: neo
    !! Electron density at the magnetic axis \(f_0 = n_{e,0}\).
    !REAL(rp)                        :: Teo
    !! Electron temperature at the magnetic axis \(f_0 = T_{e,0}\).
    !REAL(rp)                        :: Zeffo
    !! \(Z_{eff}\) at the magnetic axis \(f_0 = Z_{eff,0}\).
    !REAL(rp)                        :: n_ne
    !! Exponent \(n\) used in \(\tanh^n(r)\) of the electron
    !! density profile.
    !REAL(rp)                        :: n_Te
    !! Exponent \(n\) used in \(\tanh^n(r)\) of the electron
    !! temperature profile.
    !REAL(rp)                        :: n_Zeff
    !! Exponent \(n\) used in \(\tanh^n(r)\) of the \(Z_{eff}\) profile.
    !REAL(rp), DIMENSION(4)          :: a_ne
    !! Coefficients of the polynomial electron density profile. 
    !! See detailed description above,
    !! a_ne=(\(a_{0}\),\(a_{2}\),\(a_{3}\),\(a_{4}\)).
    !REAL(rp), DIMENSION(4)          :: a_Te
    !! Coefficients of the polynomial electron temperature profile. 
    !! See detailed description above,
    !! a_Te=(\(a_{0}\),\(a_{2}\),\(a_{3}\),\(a_{4}\)).
    !REAL(rp), DIMENSION(4)          :: a_Zeff
    !! Coefficients of the \(Z_{eff}\) profile. See detailed
    !! description above, a_Zeff=(\(a_{0}\),\(a_{2}\),\(a_{3}\),\(a_{4}\)).
    !LOGICAL                         :: axisymmetric
    !! Flag to indicate if the plasma profiles are axisymmetric.
    INTEGER  :: ii,kk
    !REAL(rp)  ::  n_REr0
    !REAL(rp)  ::  n_tauion
    !REAL(rp)  ::  n_lamfront
    !REAL(rp)  ::  n_lamback,n_lamshelf,n_shelfdelay,n_tauin,n_tauout,n_shelf
    REAL(rp)  ::  rm,r_a!,psiN_0

    !NAMELIST /plasmaProfiles/ radius_profile,ne_profile,neo,n_ne,a_ne, &
    !     Te_profile,Teo,n_Te,a_Te,n_REr0,n_tauion,n_lamfront,n_lamback, &
    !     Zeff_profile,Zeffo,n_Zeff,a_Zeff,filename,axisymmetric, &
    !     n_lamshelf,n_shelfdelay,n_tauin,n_tauout,n_shelf,psiN_0

    if (params%mpi_params%rank .EQ. 0) then
       write(output_unit_write,'("* * * * * * * * INITIALIZING PROFILES * * * * * * * *")')
    end if
    
    if (params%profile_model(1:10).eq.'ANALYTICAL') then
       !open(unit=default_unit_open,file=TRIM(params%path_to_inputs), &
       !     status='OLD',form='formatted')
       !read(default_unit_open,nml=plasmaProfiles)
       !close(default_unit_open)

       P%a = radius_profile
       P%R0=F%Ro
       P%Z0=F%Zo 

       
       P%ne_profile = TRIM(ne_profile)
       P%neo = neo
       P%n_ne = n_ne
       P%a_ne = a_ne

       P%n_REr0=n_REr0
       P%n_tauion=n_tauion
       P%n_tauin=n_tauin
       P%n_tauout=n_tauout
       P%n_shelfdelay=n_shelfdelay
       P%n_lamfront=n_lamfront
       P%n_lamback=n_lamback
       P%n_lamshelf=n_lamshelf
       P%n_shelf=n_shelf
       P%psiN_0=psiN_0
        
       P%Te_profile = TRIM(Te_profile)
       P%Teo = Teo*C_E ! Converted to Joules
       P%n_Te = n_Te
       P%a_Te = a_Te

       P%Zeff_profile = TRIM(Zeff_profile)
       P%Zeffo = Zeffo
       P%n_Zeff = n_Zeff
       P%a_Zeff = a_Zeff


       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'("ANALYTICAL")')
          write(output_unit_write,'("ne profile: ",A20)') P%ne_profile
          write(output_unit_write,'("Te profile: ",A20)') P%Te_profile
          write(output_unit_write,'("Zeff profile: ",A20)') P%Zeff_profile
       end if
       
       if (params%field_eval.eq.'interp') then

          P%axisymmetric = axisymmetric

          P%dims(1) = F%dims(1)
          P%dims(3) = F%dims(3)

          call ALLOCATE_2D_PROFILES_ARRAYS(params,P)

          P%X%R=F%X%R
          P%X%Z=F%X%Z

          do ii=1_idef,P%dims(1)
             do kk=1_idef,P%dims(3)

                rm=sqrt((P%X%R(ii)-P%R0)**2+(P%X%Z(kk)-P%Z0)**2)
                r_a=rm/P%a

                SELECT CASE (TRIM(P%ne_profile))
                CASE('FLAT')
                   P%ne_2D(ii,kk) = P%neo
                CASE('SPONG')
                   P%ne_2D(ii,kk) = P%neo*(1._rp-0.2*r_a**8)+P%n_ne
                CASE('RE-EVO')                   
                   !flat profile placeholder, updates every timestep
                   P%ne_2D(ii,kk) = P%neo
                CASE('RE-EVO1')                   
                   !flat profile placeholder, updates every timestep
                   P%ne_2D(ii,kk) = P%neo
                CASE('RE-EVO-PSI')                   
                   !flat profile placeholder, updates every timestep
                   P%ne_2D(ii,kk) = P%neo
                CASE('RE-EVO-PSIN-SG')                   
                   !flat profile placeholder, updates every timestep
                   P%ne_2D(ii,kk) = P%neo
                CASE('RE-EVO-PSIP-G')                   
                   !flat profile placeholder, updates every timestep
                   P%ne_2D(ii,kk) = P%neo
                CASE('RE-EVO-PSIP-G1')                   
                   !flat profile placeholder, updates every timestep
                   P%ne_2D(ii,kk) = P%neo
                CASE('MST_FSA')
                   P%ne_2D(ii,kk) = P%neo*(1._rp-r_a**4._rp)**4._rp
                CASE DEFAULT
                   P%ne_2D(ii,kk) = P%neo
                END SELECT

                SELECT CASE (TRIM(P%Te_profile))
                CASE('FLAT')       
                   P%Te_2D(ii,kk) = P%Teo
                CASE('SPONG')
                   P%Te_2D(ii,kk) = (P%Teo-P%n_Te)*(1._rp-0.6*r_a**2)**2+P%n_Te
                CASE('MST_FSA')
                   P%Te_2D(ii,kk) = P%Teo*(1._rp-r_a**8._rp)**4._rp
                CASE DEFAULT
                   P%Te_2D(ii,kk) = P%Teo
                END SELECT

                SELECT CASE (TRIM(P%Zeff_profile))
                CASE('FLAT') 
                   P%Zeff_2D(ii,kk) = P%Zeffo
                CASE('SPONG')
                   P%Zeff_2D(ii,kk) = P%Zeffo
                CASE DEFAULT
                   P%Zeff_2D(ii,kk) = P%Zeffo
                END SELECT

                
             end do
          end do


          
       end if

    else if (params%profile_model(1:8).eq.'EXTERNAL') then
       !open(unit=default_unit_open,file=TRIM(params%path_to_inputs), &
       !     status='OLD',form='formatted')
       !read(default_unit_open,nml=plasmaProfiles)
       !close(default_unit_open)

       P%a = radius_profile
       P%R0=F%Ro
       P%Z0=F%Zo 
       
       P%ne_profile = TRIM(ne_profile)
       P%neo = neo
       P%Te_profile = TRIM(Te_profile)
       P%Teo = Teo*C_E ! Converted to Joules
       P%Zeff_profile = TRIM(Zeff_profile)
       P%Zeffo = Zeffo

       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'("EXTERNAL")')
       end if
       
       P%filename = TRIM(filename)
       P%axisymmetric = axisymmetric

       call load_profiles_data_from_hdf5(params,P)
    else if (params%profile_model.eq.'UNIFORM') then
       !open(unit=default_unit_open,file=TRIM(params%path_to_inputs), &
       !     status='OLD',form='formatted')
       !read(default_unit_open,nml=plasmaProfiles)
       !close(default_unit_open)

       P%a = radius_profile
       P%R0=F%Ro
       P%Z0=F%Zo 
       
       P%ne_profile = TRIM(ne_profile)
       P%neo = neo
       P%n_ne = 0.0_rp
       P%a_ne = (/0.0_rp,0.0_rp,0.0_rp,0.0_rp/)

       P%Te_profile = TRIM(Te_profile)
       P%Teo = Teo*C_E ! Converted to Joules
       P%n_Te = 0.0_rp
       P%a_Te = (/0.0_rp,0.0_rp,0.0_rp,0.0_rp/)

       P%Zeff_profile = TRIM(Zeff_profile)
       P%Zeffo = Zeffo
       P%n_Zeff = 0.0_rp
       P%a_Zeff = (/0.0_rp,0.0_rp,0.0_rp,0.0_rp/)

       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'("UNIFORM")')
       end if
       
    endif

    if (params%mpi_params%rank .EQ. 0) then
       write(output_unit_write,'("* * * * * * * * * * * * * * * * * * * * * * * * *")')
    end if
    
  end subroutine initialize_profiles


  subroutine uniform_profiles(vars,P)
    !! @note Subroutine that returns the value of uniform plasma p
    !! arameters. @endnote
    !! This subroutie is used only when the simulation is ran for a 'UNIFORM'
    !! plasma. As a convention, in a uniform plasma we set
    !! \(n_e = n_{e,0}\), \(T_e = T_{e,0}\), and \(Z_{eff} = Z_{eff,0}\).
    TYPE(PROFILES), INTENT(IN)     :: P
    !! An instance of KORC's derived type PROFILES containing all the information
    !! about the plasma profiles used in the simulation. See [[korc_types]]
    !! and [[korc_profiles]].
    TYPE(PARTICLES), INTENT(INOUT) :: vars
    !! An instance of PARTICLES containing the variables of a given species.

    vars%ne = P%neo
    vars%Te = P%Teo
    vars%Zeff = P%Zeffo
  end subroutine uniform_profiles

  subroutine analytical_profiles_p(pchunk,time,params,Y_R,Y_Z,P,F,ne,Te,Zeff,PSIp)
    !! @note Subroutine that calculates the analytical plasma profiles at
    !! the particles' position. @endnote
    TYPE(KORC_PARAMS), INTENT(IN)                           :: params
    INTEGER, INTENT(IN)                             :: pchunk
    REAL(rp), DIMENSION(params%pchunk), INTENT(IN)  :: Y_R,Y_Z,PSIp
    REAL(rp), INTENT(IN)  :: time
    TYPE(PROFILES), INTENT(IN)                         :: P
    !! An instance of KORC's derived type PROFILES containing all the
    !! information about the plasma profiles used in the simulation.
    !! See [[korc_types]] and [[korc_profiles]].
    TYPE(FIELDS), INTENT(IN)      :: F
    REAL(rp), DIMENSION(params%pchunk),INTENT(OUT) :: ne
    !! Background electron density seen by simulated particles.
    REAL(rp), DIMENSION(params%pchunk),INTENT(OUT) :: Te
    !! Backgroun temperature density seen by simulated particles.
    REAL(rp), DIMENSION(params%pchunk),INTENT(OUT) :: Zeff
    !! Effective atomic charge seen by simulated particles.
    INTEGER(ip)                                        :: cc
    !! Particle iterator.
    REAL(rp) :: R0,Z0,a,ne0,n_ne,Te0,n_Te,Zeff0,R0a
    REAL(rp) :: R0_RE,Z0_RE,sigmaR_RE,sigmaZ_RE,psimax_RE
    REAL(rp) :: n_REr0,n_tauion,n_lamfront,n_lamback,n_lamshelf
    REAL(rp) :: n_psifront,n_psiback,n_psishelf
    REAL(rp) :: n_tauin,n_tauout,n_shelfdelay,n_shelf
    REAL(rp) :: n0t,n_taut
    REAL(rp) :: PSIp0,PSIp_lim,psiN_0
    REAL(rp), DIMENSION(params%pchunk) :: r_a,rm,rm_RE,PSIpN,PSIp_temp
    
    R0=P%R0
    Z0=P%Z0
    a=P%a
    R0a=F%AB%Ro
    
    ne0=P%neo
    n_ne=P%n_ne

    Te0=P%Teo
    n_Te=P%n_Te

    Zeff0=P%Zeffo

    R0_RE=P%R0_RE
    Z0_RE=P%Z0_RE
    n_REr0=P%n_REr0
    n_tauion=P%n_tauion
    n_tauin=P%n_tauin
    n_tauout=P%n_tauout
    n_shelfdelay=P%n_shelfdelay
    n_lamfront=P%n_lamfront
    n_lamback=P%n_lamback
    n_lamshelf=P%n_lamshelf
    n_psifront=P%n_lamfront*params%cpp%length
    n_psiback=P%n_lamback*params%cpp%length
    n_psishelf=P%n_lamshelf*params%cpp%length
    n_shelf=P%n_shelf
    
    PSIp_lim=F%PSIp_lim
    PSIp0=F%PSIP_min
    psiN_0=P%psiN_0
    

!    write(output_unit_write,*) 'PSIp',PSIp(1)*(params%cpp%Bo*params%cpp%length**2)
!    write(output_unit_write,*) 'PSIp_lim',PSIp_lim*(params%cpp%Bo*params%cpp%length**2)   
!    write(output_unit_write,*) 'PSIp0',PSIp0*(params%cpp%Bo*params%cpp%length**2)
    
!    write(output_unit_write,'("R0_RE: "E17.10)') R0_RE
!    write(output_unit_write,'("Z0_RE: "E17.10)') Z0_RE
!    write(output_unit_write,'("n_REr0: "E17.10)') n_REr0

    
    SELECT CASE (TRIM(P%ne_profile))
    CASE('FLAT')
       !$OMP SIMD
       do cc=1_idef,pchunk
          ne(cc) = ne0
       end do
       !$OMP END SIMD
    CASE('FLAT-RAMP')
       !$OMP SIMD
       do cc=1_idef,pchunk
          ne(cc) = n_ne+(ne0-n_ne)*time/n_tauion
       end do
       !$OMP END SIMD

    CASE('TANH-RAMP')
       !$OMP SIMD
       do cc=1_idef,pchunk
          ne(cc) = n_ne+(ne0-n_ne)/2*(tanh((time-n_shelfdelay)/n_tauion)+1._rp)
       end do
       !$OMP END SIMD
    CASE('SINE')
       !$OMP SIMD
       do cc=1_idef,pchunk
          ne(cc) = n_ne+(ne0-n_ne)*sin(time/n_tauion)
       end do
       !$OMP END SIMD
    CASE('SPONG')
       !$OMP SIMD
       do cc=1_idef,pchunk
          rm(cc)=sqrt((Y_R(cc)-R0)**2+(Y_Z(cc)-Z0)**2)
          r_a(cc)=rm(cc)/a
          ne(cc) = ne0*(1._rp-0.2*r_a(cc)**8)+n_ne
       end do
       !$OMP END SIMD
    CASE('MST_FSA')       
       !$OMP SIMD
       do cc=1_idef,pchunk
          rm(cc)=sqrt((Y_R(cc)-R0a)**2+(Y_Z(cc)-Z0)**2)
          r_a(cc)=rm(cc)/a
          ne(cc) = (ne0-n_ne)*(1._rp-r_a(cc)**4._rp)**4._rp+n_ne

          !write(6,*) 'R',Y_R*params%cpp%length,'R0',R0*params%cpp%length, &
          !     'Z',Y_Z*params%cpp%length,'Z0',Z0*params%cpp%length, &
          !     'a',a*params%cpp%length
          !write(6,*) 'r_a',r_a,'ne',ne(cc)*params%cpp%density
          
       end do
       !$OMP END SIMD
    CASE('RE-EVO')
       !$OMP SIMD
       do cc=1_idef,pchunk
          rm_RE(cc)=sqrt((Y_R(cc)-R0_RE)**2+(Y_Z(cc)-Z0_RE)**2)
          ne(cc) = (ne0-n_ne)/4._rp*(1+tanh((rm_RE(cc)+ &
               n_REr0*(time/n_tauion-1))/n_lamfront))* &
               (1+tanh(-(rm_RE(cc)-n_REr0)/n_lamback))+n_ne
       end do
       !$OMP END SIMD
    CASE('RE-EVO1')
       !$OMP SIMD
       do cc=1_idef,pchunk
          rm_RE(cc)=sqrt((Y_R(cc)-R0_RE)**2+(Y_Z(cc)-Z0_RE)**2)
          ne(cc) = (ne0-n_ne)/8._rp*(1+tanh((rm_RE(cc)+ &
               n_REr0*(time/n_tauion-1))/n_lamfront))* &
               (1+tanh(-(rm_RE(cc)-n_REr0)/n_lamback))* &
               (2*(n_shelf-n_ne)/(ne0-n_ne)+(ne0-n_shelf)/(ne0-n_ne)* &
               (1-tanh((rm_RE(cc)+n_REr0*((time-n_shelfdelay)/n_tauin-1))/ &
               n_lamshelf)))+n_ne
       end do
       !$OMP END SIMD

       
    CASE('RE-EVO-PSI')
       !$OMP SIMD
       do cc=1_idef,pchunk
          PSIpN(cc)=(PSIp(cc)-PSIp0)/(PSIp_lim-PSIp0)
          ne(cc) = (ne0-n_ne)/8._rp*(1+tanh((sqrt(abs(PSIpN(cc)))+ &
               sqrt(abs(psiN_0))*(time/n_tauion-1))/n_psifront))* &
               (1+tanh(-(sqrt(abs(PSIpN(cc)))-sqrt(abs(psiN_0)))/n_psiback))* &
               (2*(n_shelf-n_ne)/(ne0-n_ne)+(ne0-n_shelf)/(ne0-n_ne)* &
               (1-tanh((sqrt(abs(PSIpN(cc)))+ sqrt(abs(psiN_0))* &
               ((time-n_shelfdelay)/n_tauin-1))/n_psishelf)))+n_ne
       end do
       !$OMP END SIMD

       !       write(output_unit_write,*) 'at time ',time*params%cpp%time, &
       !' ne: ',ne(1)/params%cpp%length**3
      
!       !$OMP SIMD
!       do cc=1_idef,8
!          if(isnan(ne(cc))) then
!             write(output_unit_write,*) 'PSIp: ',PSIp(cc)
!             write(output_unit_write,*) 'PSIp0: ',PSIp0
!             write(output_unit_write,*) 'PSIp_lim: ',PSIp_lim
!             write(output_unit_write,*) 'PSIpN: ',PSIpN(cc)

!             stop 'ne_eval is a NaN'
!          end if
!       end do
       !       !$OMP END SIMD

    CASE('RE-EVO-PSIN-SG')

       n0t=(ne0-n_ne)/2._rp*(tanh(time/n_tauin)- &
            tanh((time-n_shelfdelay)/n_tauin))
       n_taut=n_psishelf*erf((time+params%dt/100._rp)/n_tauion)
       
       !$OMP SIMD
       do cc=1_idef,pchunk
          PSIpN(cc)=(PSIp(cc)-PSIp0)/(PSIp_lim-PSIp0)
          ne(cc) = n0t*exp(-(sqrt(abs(PSIpN(cc)))-sqrt(abs(psiN_0)))**2._rp/ &
               (2._rp*n_taut**2._rp))*(1._rp+erf(-10._rp* &
               (sqrt(abs(PSIpN(cc)))-sqrt(abs(psiN_0)))/ &
               (sqrt(2._rp)*n_taut)))/2._rp+n_ne
       end do
       !$OMP END SIMD
       
    CASE('RE-EVO-PSIP-G')

!       write(output_unit_write,*) 'time: ',time*params%cpp%time
       
       n0t=(ne0-n_ne)/2._rp*(tanh((time-n_tauin)/n_tauin)- &
            tanh((time-n_shelfdelay)/n_tauout))
       n_taut=n_psishelf*erf((time+params%dt/100._rp)/n_tauion)
       
       !$OMP SIMD
       do cc=1_idef,pchunk
          PSIp_temp(cc)=PSIp(cc)*(params%cpp%Bo*params%cpp%length**2)
          ne(cc) = n0t*exp(-(sqrt(abs(PSIp_temp(cc)))-sqrt(abs(psiN_0)))**2._rp/ &
               (2._rp*n_taut**2._rp))+n_ne
       end do
       !$OMP END SIMD

    CASE('RE-EVO-PSIP-G1')

!       write(output_unit_write,*) 'time: ',time*params%cpp%time
       
       n0t=(ne0-n_ne)/2._rp*(tanh((time-1.75*n_tauin)/n_tauin)- &
            tanh((time-n_shelfdelay)/n_tauout))
       n_taut=n_psishelf*erf((time+params%dt/100._rp)/n_tauion)
       
       !$OMP SIMD
       do cc=1_idef,pchunk
          PSIp_temp(cc)=PSIp(cc)*(params%cpp%Bo*params%cpp%length**2)
          ne(cc) = n0t*exp(-(sqrt(abs(PSIp_temp(cc)))-sqrt(abs(psiN_0)))**2._rp/ &
               (2._rp*n_taut**2._rp))+n_ne
       end do
       !$OMP END SIMD
       
    CASE DEFAULT
       !$OMP SIMD
       do cc=1_idef,pchunk
          ne(cc) = ne0
       end do
       !$OMP END SIMD
    END SELECT

    SELECT CASE (TRIM(P%Te_profile))
    CASE('FLAT')
       !$OMP SIMD
       do cc=1_idef,pchunk
          Te(cc) = Te0
       end do
       !$OMP END SIMD
    CASE('SPONG')
       !$OMP SIMD
       do cc=1_idef,pchunk
          rm(cc)=sqrt((Y_R(cc)-R0)**2+(Y_Z(cc)-Z0)**2)
          r_a(cc)=rm(cc)/a
          Te(cc) = Te0*(1._rp-0.6*r_a(cc)**2)**2+Te0*n_Te
       end do
       !$OMP END SIMD
    CASE('MST_FSA')
       !$OMP SIMD
       do cc=1_idef,pchunk
          rm(cc)=sqrt((Y_R(cc)-R0a)**2+(Y_Z(cc)-Z0)**2)
          r_a(cc)=rm(cc)/a
          Te(cc) = (Te0-n_Te)*(1._rp-r_a(cc)**8._rp)**4._rp+n_Te

          !write(6,*) 'T_e',Te(cc)*params%cpp%temperature/C_E         
          
       end do
       !$OMP END SIMD
    CASE DEFAULT
       !$OMP SIMD
       do cc=1_idef,pchunk
          Te(cc) = P%Teo
       end do
       !$OMP END SIMD
    END SELECT

    SELECT CASE (TRIM(P%Zeff_profile))
    CASE('FLAT')
       !$OMP SIMD
       do cc=1_idef,pchunk
          Zeff(cc) = P%Zeffo
       end do
       !$OMP END SIMD
    CASE('SPONG')
       !$OMP SIMD
       do cc=1_idef,pchunk
          Zeff(cc) = P%Zeffo
       end do
       !$OMP END SIMD
    CASE DEFAULT
       !$OMP SIMD
       do cc=1_idef,pchunk
          Zeff(cc) = P%Zeffo
       end do
       !$OMP END SIMD
    END SELECT
          

!    write(output_unit_write,*) PSIpN(1)
    
!    write(output_unit_write,'("ne: "E17.10)') ne(1)/params%cpp%length**3
!    write(output_unit_write,'("rm_RE: "E17.10)') rm_RE(1)
    
  end subroutine analytical_profiles_p

  subroutine get_analytical_profiles(P,Y,ne,Te,Zeff,flag)
    !! @note Subroutine that calculates the analytical plasma profiles at
    !! the particles' position. @endnote
    TYPE(PROFILES), INTENT(IN)                         :: P
    !! An instance of KORC's derived type PROFILES containing all the
    !! information about the plasma profiles used in the simulation.
    !! See [[korc_types]] and [[korc_profiles]].
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)  :: Y
    !! Particles' position in toroidal coordinates; Y(1,:) = \(r\), Y(2,:)
    !! = \(\theta\), Y(3,:) = \(\zeta\).
    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ne
    !! Background electron density seen by simulated particles.
    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Te
    !! Backgroun temperature density seen by simulated particles.
    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Zeff
    !! Effective atomic charge seen by simulated particles.
    INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: flag
    !! Flag for each particle to decide whether it is being followed
    !! (flag=T) or not (flag=F).
    REAL(rp)                                           :: r_a
    !! Normalized toroidal radial position of simulated particles
    !! \(r/a\), where \(a\) is the plasma radius.
    REAL(rp)                                           :: fr
    !! Calculated radial profile.
    INTEGER(ip)                                        :: pp
    !! Particle iterator.
    INTEGER(ip)                                        :: ss
    !! Species iterator.

    if (Y(2,1).eq.0) then
       ss=1_idef
    else
       ss = size(Y,1)
    end if
   

    !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,fr,r_a) &
    !$OMP& SHARED(P,Y,ne,Te,Zeff,flag)
    do pp=1_idef,ss
       if ( flag(pp) .EQ. 1_is ) then

          r_a = Y(pp,1)/P%a

!          write(output_unit_write,'("r: ",E17.10)') r_a
          
          SELECT CASE (TRIM(P%ne_profile))
          CASE('TANH')
             fr = 1_ip - TANH(2.0_rp*r_a)**P%n_ne
             ne(pp) = P%neo*fr
          CASE('FLAT')
             ne(pp) = P%neo
          CASE('POLYNOMIAL')
             fr = P%a_ne(1) + P%a_ne(2)*r_a + P%a_ne(3)*r_a**2 + &
                  P%a_ne(4)*r_a**3
             ne(pp) = P%neo*fr
          CASE('SPONG')
             fr = P%neo*(1._rp-0.2*r_a**8)+P%neo*P%n_ne
             ne(pp) = P%neo*fr
          CASE DEFAULT
             ne(pp) = P%neo
          END SELECT

          SELECT CASE (TRIM(P%Te_profile))
          CASE('TANH')
             fr = 1_ip - TANH(2.0_rp*r_a)**P%n_Te
             Te(pp) = P%Teo*fr
          CASE('FLAT')
             Te(pp) = P%Teo
          CASE('POLYNOMIAL')
             fr = P%a_Te(1) + P%a_Te(2)*r_a + P%a_Te(3)*r_a**2 + &
                  P%a_Te(4)*r_a**3
             Te(pp) = P%Teo*fr
          CASE('SPONG')
             fr = P%Teo*(1._rp-0.6*r_a**2)**2+P%Teo*P%n_Te
             ne(pp) = P%neo*fr
          CASE DEFAULT
             Te(pp) = P%Teo
          END SELECT

          SELECT CASE (TRIM(P%Zeff_profile))
          CASE('TANH')
             fr = 1_ip - TANH(2.0_rp*r_a)**P%n_Zeff
             Zeff(pp) = P%Zeffo*fr
          CASE('FLAT')
             Zeff(pp) = P%Zeffo
          CASE('POLYNOMIAL')
             fr = P%a_Zeff(1) + P%a_Zeff(2)*r_a + P%a_Zeff(3)*r_a**2 + &
                  P%a_Zeff(4)*r_a**3
             Zeff(pp) = P%Zeffo*fr
          CASE('SPONG')
             Zeff(pp) = P%Zeffo
          CASE DEFAULT
             Zeff(pp) = P%Zeffo
          END SELECT
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine get_analytical_profiles


  subroutine get_profiles(params,vars,P,F)
    !! @note Subrotuine that calls the appropriate subroutine for calculating
    !! or interpolating the plasma profiles at the particles' position. @endnote
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    !! Core KORC simulation parameters.
    TYPE(PARTICLES), INTENT(INOUT) :: vars
    !! An instance of PARTICLES containing the variables of a given species.
    TYPE(PROFILES), INTENT(IN)     :: P
    !! An instance of KORC's derived type PROFILES containing all
    !! the information about the plasma profiles used in the
    !! simulation. See [[korc_types] and [[korc_profiles]].
    TYPE(FIELDS), INTENT(IN)                                   :: F
    !! An instance of the KORC derived type FIELDS.

    SELECT CASE (TRIM(params%profile_model))
    CASE('ANALYTICAL')

!       write(output_unit_write,'("Y in: ",E17.10)') vars%Y(1,:)
       
       call cyl_to_cart(vars%Y, vars%X)

!       write(output_unit_write,'("X getprof: ",E17.10)') vars%X(1,:)
       
       call cart_to_tor_check_if_confined(vars%X,F,vars%Y,vars%flagCon)

!       write(output_unit_write,'("flag: ",I15)') vars%flagCon(1) 

       call get_analytical_profiles(P,vars%Y,vars%ne,vars%Te, &
            vars%Zeff,vars%flagCon)

       call cart_to_cyl(vars%X, vars%Y)

!       write(output_unit_write,'("Y out: ",E17.10)') vars%Y(1,:)
    CASE('EXTERNAL')
       call interp_profiles(params,vars,P)
    CASE('UNIFORM')
       call uniform_profiles(vars,P)
    CASE DEFAULT
    END SELECT

  end subroutine get_profiles


  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  ! Subroutines for getting the profiles data from HDF5 files
  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


  subroutine load_profiles_data_from_hdf5(params,P)
    !! @note Subroutine that loads pre-computed plasma profiles' data
    !! from an input HDF5 file. @endnote
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(INOUT)  :: P
    !! An instance of KORC's derived type PROFILES containing all the
    !! information about the plasma profiles used in the
    !! simulation. See [[korc_types]] and [[korc_profiles]].
    CHARACTER(MAX_STRING_LENGTH)   :: filename
    !!String containing the name of the input HDF5 file.
    CHARACTER(MAX_STRING_LENGTH)   :: gname
    !! String containing the group name of a parameter in the HDF5 file.
    CHARACTER(MAX_STRING_LENGTH)   :: subgname
    !! String containing the subgroup name of a parameter in the HDF5 file.
    CHARACTER(MAX_STRING_LENGTH)   :: dset
    !!Name of data set to read from file.
    INTEGER(HID_T)                 :: h5file_id
    !! HDF5 file identifier.
    INTEGER(HID_T)                 :: group_id
    !! HDF5 group identifier.
    INTEGER(HID_T)                 :: subgroup_id
    !!HDF5 subgroup identifier.
    REAL(rp)                       :: rdatum
    !!
    INTEGER                        :: h5error
    !! HDF5 error status.

    filename = TRIM(P%filename)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: load_profiles_data_from_hdf5 --> h5fopen_f")')
    end if


    dset = "/NR"
    call load_from_hdf5(h5file_id,dset,rdatum)
    P%dims(1) = INT(rdatum)

    if (P%axisymmetric) then
       P%dims(2) = 0
    else
       dset = "/NPHI"
       call load_from_hdf5(h5file_id,dset,rdatum)
       P%dims(2) = INT(rdatum)
    end if

    dset = "/NZ"
    call load_from_hdf5(h5file_id,dset,rdatum)
    P%dims(3) = INT(rdatum)

    if (P%axisymmetric) then
       call ALLOCATE_2D_PROFILES_ARRAYS(params,P)
    else
       call ALLOCATE_3D_PROFILES_ARRAYS(P)
    end if

    dset = "/R"
    call load_array_from_hdf5(h5file_id,dset,P%X%R)

    if (.NOT.P%axisymmetric) then
       dset = "/PHI"
       call load_array_from_hdf5(h5file_id,dset,P%X%PHI)
    end if

    dset = "/Z"
    call load_array_from_hdf5(h5file_id,dset,P%X%Z)

    dset = "/FLAG"
    if (P%axisymmetric) then
       call load_array_from_hdf5(h5file_id,dset,P%FLAG2D)
    else
       call load_array_from_hdf5(h5file_id,dset,P%FLAG3D)
    end if

    dset = "/ne"
    if (P%axisymmetric) then
       call load_array_from_hdf5(h5file_id,dset,P%ne_2D)
    else
       call load_array_from_hdf5(h5file_id,dset,P%ne_3D)
    end if

    dset = "/Te"
    if (P%axisymmetric) then
       call load_array_from_hdf5(h5file_id,dset,P%Te_2D)
       P%Te_2D = P%Te_2D*C_E
    else
       call load_array_from_hdf5(h5file_id,dset,P%Te_3D)
       P%Te_3D = P%Te_3D*C_E
    end if

    !write(output_unit_write,'("Te: ",E17.10)') P%Te_2D(1,1)

    dset = "/Zeff"
    if (P%axisymmetric) then
       call load_array_from_hdf5(h5file_id,dset,P%Zeff_2D)
    else
       call load_array_from_hdf5(h5file_id,dset,P%Zeff_3D)
    end if

    if (params%profile_model(10:10).eq.'H') then

       dset = "/RHON"
       call load_array_from_hdf5(h5file_id,dset,P%RHON)
       dset = "/nRE"
       call load_array_from_hdf5(h5file_id,dset,P%nRE_2D)
       dset = "/nAr0"
       call load_array_from_hdf5(h5file_id,dset,P%nAr0_2D)
       dset = "/nAr1"
       call load_array_from_hdf5(h5file_id,dset,P%nAr1_2D)
       dset = "/nAr2"
       call load_array_from_hdf5(h5file_id,dset,P%nAr2_2D)
       dset = "/nAr3"
       call load_array_from_hdf5(h5file_id,dset,P%nAr3_2D)
       dset = "/nD"
       call load_array_from_hdf5(h5file_id,dset,P%nD_2D)
       dset = "/nD1"
       call load_array_from_hdf5(h5file_id,dset,P%nD1_2D)
       
    end if
   
    call h5fclose_f(h5file_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: load_profiles_data_from_hdf5 --> h5fclose_f")')
    end if
  end subroutine load_profiles_data_from_hdf5


  !> @brief Subroutine that allocates the mesh information and 2-D arrays for keeping the data of pre-computed plasma profiles.
  !!
  !! @param[out] P An instance of KORC's derived type PROFILES containing all the information about the plasma profiles used in the
  !! simulation. See korc_types.f90 and korc_profiles.f90.
  subroutine ALLOCATE_2D_PROFILES_ARRAYS(params,P)
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    TYPE(PROFILES), INTENT(INOUT) :: P

    ALLOCATE(P%X%R(P%dims(1)))
    ALLOCATE(P%X%Z(P%dims(3)))
    ALLOCATE(P%FLAG2D(P%dims(1),P%dims(3)))
    ALLOCATE(P%ne_2D(P%dims(1),P%dims(3)))
    ALLOCATE(P%Te_2D(P%dims(1),P%dims(3)))
    ALLOCATE(P%Zeff_2D(P%dims(1),P%dims(3)))

    if (params%profile_model(10:10).eq.'H') then
       ALLOCATE(P%RHON(P%dims(1),P%dims(3)))
       ALLOCATE(P%nRE_2D(P%dims(1),P%dims(3)))
       ALLOCATE(P%nAr0_2D(P%dims(1),P%dims(3)))
       ALLOCATE(P%nAr1_2D(P%dims(1),P%dims(3)))
       ALLOCATE(P%nAr2_2D(P%dims(1),P%dims(3)))
       ALLOCATE(P%nAr3_2D(P%dims(1),P%dims(3)))
       ALLOCATE(P%nD_2D(P%dims(1),P%dims(3)))
       ALLOCATE(P%nD1_2D(P%dims(1),P%dims(3)))
    end if
    
  end subroutine ALLOCATE_2D_PROFILES_ARRAYS

  subroutine ALLOCATE_3D_PROFILES_ARRAYS(P)
    !! @note Subroutine that allocates the mesh information and 3-D arrays
    !! for keeping the data of pre-computed plasma profiles. @endnote
    TYPE(PROFILES), INTENT(INOUT) :: P
    !! @param[out] P An instance of KORC's derived type PROFILES containing
    !! all the information about the plasma profiles used in the
    !! simulation. See [[korc_types]] and [[korc_profiles]].
    ALLOCATE(P%X%R(P%dims(1)))
    ALLOCATE(P%X%PHI(P%dims(2)))
    ALLOCATE(P%X%Z(P%dims(3)))
    ALLOCATE(P%FLAG3D(P%dims(1),P%dims(2),P%dims(3)))
    ALLOCATE(P%ne_3D(P%dims(1),P%dims(2),P%dims(3)))
    ALLOCATE(P%Te_3D(P%dims(1),P%dims(2),P%dims(3)))
    ALLOCATE(P%Zeff_3D(P%dims(1),P%dims(2),P%dims(3)))
  end subroutine ALLOCATE_3D_PROFILES_ARRAYS

  subroutine DEALLOCATE_PROFILES_ARRAYS(P)
    TYPE(PROFILES), INTENT(INOUT)              :: P

    if (ALLOCATED(P%X%R)) DEALLOCATE(P%X%R)
    if (ALLOCATED(P%X%PHI)) DEALLOCATE(P%X%PHI)
    if (ALLOCATED(P%X%Z)) DEALLOCATE(P%X%Z)

    if (ALLOCATED(P%FLAG2D)) DEALLOCATE(P%FLAG2D)
    if (ALLOCATED(P%FLAG3D)) DEALLOCATE(P%FLAG3D)

    if (ALLOCATED(P%ne_2D)) DEALLOCATE(P%ne_2D)
    if (ALLOCATED(P%Te_2D)) DEALLOCATE(P%Te_2D)
    if (ALLOCATED(P%Zeff_2D)) DEALLOCATE(P%Zeff_2D)
    if (ALLOCATED(P%ne_3D)) DEALLOCATE(P%ne_3D)
    if (ALLOCATED(P%Te_3D)) DEALLOCATE(P%Te_3D)
    if (ALLOCATED(P%Zeff_3D)) DEALLOCATE(P%Zeff_3D)

    if (ALLOCATED(P%RHON)) DEALLOCATE(P%RHON)
    if (ALLOCATED(P%nRE_2D)) DEALLOCATE(P%nRE_2D)
    if (ALLOCATED(P%nAr0_2D)) DEALLOCATE(P%nAr0_2D)
    if (ALLOCATED(P%nAr1_2D)) DEALLOCATE(P%nAr1_2D)
    if (ALLOCATED(P%nAr2_2D)) DEALLOCATE(P%nAr2_2D)
    if (ALLOCATED(P%nAr3_2D)) DEALLOCATE(P%nAr3_2D)
    if (ALLOCATED(P%nD_2D)) DEALLOCATE(P%nD_2D)
    if (ALLOCATED(P%nD1_2D)) DEALLOCATE(P%nD1_2D)
    
#ifdef FIO
    if (ALLOCATED(P%FIO_nimp)) DEALLOCATE(P%FIO_nimp)
#endif
    
  end subroutine DEALLOCATE_PROFILES_ARRAYS
  
end module korc_profiles
