MODULE korc_simple_equilibrium_pdf
  USE korc_types
  USE korc_constants
  USE korc_HDF5
  USE korc_hpc
  USE special_functions
  USE korc_input

  IMPLICIT NONE

  TYPE, PRIVATE :: PARAMS
     REAL(rp) :: E ! Parallel electric field normalized using the critical electric field
     REAL(rp) :: Zeff ! Effective atomic number of impurities

     REAL(rp) :: max_pitch_angle ! Maximum pitch angle of sampled PDF in degrees
     REAL(rp) :: min_pitch_angle ! Minimum pitch angle of sampled PDF in degrees
     REAL(rp) :: po ! Momentum of sampled PDF in units of mc

     REAL(rp) :: Bo
     REAL(rp) :: lambda
  END TYPE PARAMS

  TYPE(PARAMS), PRIVATE :: pdf_params
  REAL(rp), PRIVATE, PARAMETER :: xo = (C_ME*C_C**2/C_E)/1.0E6
  REAL(rp), PRIVATE, PARAMETER :: Tol = 1.0E-5_rp
  REAL(rp), PRIVATE, PARAMETER :: minmax_buffer_size = 10.0_rp

  PUBLIC :: get_equilibrium_distribution
  PRIVATE :: initialize_params,&
       save_params,&
       sample_distribution,&
       deg2rad,&
       rad2deg,&
       fRE,&
       random_norm,&
       PR,&
       P_integral,&
       IntK,&
       IntBesselK

CONTAINS

  SUBROUTINE get_equilibrium_distribution(params,eta,go,etao)
    TYPE(KORC_PARAMS), INTENT(IN) :: params
    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: eta
    REAL(rp), INTENT(IN) :: go
    REAL(rp), INTENT(OUT) :: etao

    call initialize_params(params,go)

    call save_params(params)

    call sample_distribution(params,eta,etao)
  END SUBROUTINE get_equilibrium_distribution


  SUBROUTINE initialize_params(params,go)
    TYPE(KORC_PARAMS), INTENT(IN) :: params
    REAL(rp), INTENT(IN) :: go
    !REAL(rp) :: max_pitch_angle
    !REAL(rp) :: min_pitch_angle
    !REAL(rp) :: Zeff
    !REAL(rp) :: E
    !REAL(rp) :: Bo
    !REAL(rp) :: lambda
    !NAMELIST /SimpleEquilibriumPDF/ max_pitch_angle,min_pitch_angle,Zeff,E,Bo,lambda

    !open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
    !read(default_unit_open,nml=SimpleEquilibriumPDF)
    !close(default_unit_open)

    pdf_params%max_pitch_angle = max_pitch_angle_simple
    pdf_params%min_pitch_angle = min_pitch_angle_simple
    pdf_params%Zeff = Zeff_simple
    pdf_params%E = E_simple
    pdf_params%Bo = Bo_simple
    pdf_params%lambda = lambda_simple

    pdf_params%po = sqrt(go**2 - 1.0_rp)
  END SUBROUTINE initialize_params


  FUNCTION deg2rad(x)
    REAL(rp), INTENT(IN) :: x
    REAL(rp) :: deg2rad

    deg2rad = C_PI*x/180.0_rp
  END FUNCTION deg2rad


  FUNCTION rad2deg(x)
    REAL(rp), INTENT(IN) :: x
    REAL(rp) :: rad2deg

    rad2deg = 180.0_rp*x/C_PI
  END FUNCTION rad2deg


  FUNCTION fRE(eta,p)
    REAL(rp), INTENT(IN) :: eta ! pitch angle in degrees
    REAL(rp), INTENT(IN) :: p ! momentum in units of mc
    REAL(rp) :: fRE
    REAL(rp) :: A

    A = (2.0_rp*pdf_params%E/(pdf_params%Zeff + 1.0_rp))*(p**2/SQRT(p**2.0_rp + 1.0_rp))
    fRE = 0.5_rp*A*EXP(A*COS(deg2rad(eta)))/SINH(A)
    !	fRE = fRE*PR(eta,p,pdf_params%Bo,pdf_params%lambda)
  END FUNCTION fRE


  FUNCTION random_norm(mean,sigma)
    REAL(rp), INTENT(IN) :: mean
    REAL(rp), INTENT(IN) :: sigma
    REAL(rp) :: random_norm
    REAL(rp) :: rand1, rand2

    call RANDOM_NUMBER(rand1)
    call RANDOM_NUMBER(rand2)

    random_norm = SQRT(-2.0_rp*LOG(1.0_rp-rand1))*COS(2.0_rp*C_PI*rand2);
  END FUNCTION random_norm


  FUNCTION IntK(v,x)
    IMPLICIT NONE
    REAL(rp) :: IntK
    REAL(rp), INTENT(IN) :: v
    REAL(rp), INTENT(IN) :: x

    IntK = (C_PI/SQRT(2.0_rp))*(1.0_rp - 0.25_rp*(4.0_rp*v**2 - 1.0_rp))*ERFC(SQRT(x))&
         + 0.25_rp*(4.0_rp*v**2 - 1.0_rp)*SQRT(0.5_rp*C_PI/x)*EXP(-x)
  END FUNCTION IntK


  FUNCTION besselk(v,x)
    IMPLICIT NONE
    REAL(rp) :: besselk
    REAL(rp), INTENT(IN) :: x
    REAL(rp), INTENT(IN) :: v
    REAL(4) :: ri,rk,rip,rkp

    call bessik(REAL(x,4),REAL(v,4),ri,rk,rip,rkp)
    besselk = REAL(rk,rp)
  END FUNCTION besselk


  FUNCTION IntBesselK(a,b)
    REAL(rp), INTENT(IN) :: a
    REAL(rp), INTENT(IN) :: b
    REAL(rp) :: IntBesselK
    REAL(rp) :: Iold
    REAL(rp) :: Inew
    REAL(rp) :: rerr
    REAL(rp) :: sum_f
    REAL(rp) :: v,h,z
    INTEGER :: ii,jj,npoints
    LOGICAL :: flag

    v = 5.0_rp/3.0_rp
    h = b - a
    sum_f = 0.5*(besselk(v,a) + besselk(v,b))

    Iold = 0.0_rp
    Inew = sum_f*h

    ii = 1_idef
    flag = .TRUE.
    do while (flag)
       Iold = Inew

       ii = ii + 1_idef
       npoints = 2_idef**(ii-2_idef)
       h = 0.5_rp*(b-a)/REAL(npoints,rp)
       sum_f = 0.0_rp
       do jj=1_idef,npoints
          z = a + h + 2.0_rp*(REAL(jj,rp) - 1.0_rp)*h
          sum_f = sum_f + besselk(v,z)
       end do

       Inew = 0.5_rp*Iold + sum_f*h
       rerr = ABS((Inew - Iold)/Iold)
       flag = .NOT.(rerr.LT.Tol)
    end do
    IntBesselK = Inew
  END FUNCTION IntBesselK

  SUBROUTINE P_integral(z,P)
    REAL(rp), INTENT(OUT) :: P
    REAL(rp), INTENT(IN) :: z
    REAL(rp) :: a

    P = 0.0_rp

    IF (z .LT. 0.5_rp) THEN
       a = (2.16_rp/2.0_rp**(2.0_rp/3.0_rp))*z**(1.0_rp/3.0_rp)
       P = IntBesselK(z,a) + IntK(5.0_rp/3.0_rp,a)
    ELSE IF ((z .GE. 0.5_rp).AND.(z .LT. 2.5_rp)) THEN
       a = 0.72_rp*(z + 1.0_rp)
       P = IntBesselK(z,a) + IntK(5.0_rp/3.0_rp,a)
    ELSE
       P = IntK(5.0_rp/3.0_rp,z)
    END IF
  END SUBROUTINE P_integral


  FUNCTION PR(eta,p,Bo,l)
    REAL(rp), INTENT(IN) :: eta ! in radians
    REAL(rp), INTENT(IN) :: p ! dimensionless (in units of mc)
    REAL(rp), INTENT(IN) :: Bo
    REAL(rp), INTENT(IN) :: l
    REAL(rp) :: PR
    REAL(rp) :: g
    REAL(rp) :: v
    REAL(rp) :: k
    REAL(rp) :: lc
    REAL(rp) :: z
    REAL(rp) :: Pi

    g = SQRT(p**2 + 1.0_rp)
    v = C_C*SQRT(1.0_rp - 1.0_rp/g**2)

    k = C_E*Bo*SIN(deg2rad(eta))/(g*C_ME*v)

    lc = (4.0_rp*C_PI/3.0_rp)/(k*g**3) ! Critical wavelength

    z = lc/l

    call P_integral(z,Pi)

    PR = (C_C*C_E**2)*Pi/(SQRT(3.0_rp)*C_E0*g**2*l**3)
  END FUNCTION PR


  SUBROUTINE sample_distribution(params,eta,etao)
    TYPE(KORC_PARAMS), INTENT(IN) :: params
    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: eta
    REAL(rp), INTENT(OUT) :: etao
    REAL(rp) :: go_root
    REAL(rp) :: etao_root
    REAL(rp) :: eta_buffer, eta_test
    REAL(rp) :: ratio, rand_unif
    REAL(rp), DIMENSION(:), ALLOCATABLE :: eta_samples
    REAL(rp), DIMENSION(:), ALLOCATABLE :: eta_tmp
    REAL(rp) :: minmax, min_pitch_angle, max_pitch_angle
    REAL(rp) :: deta
    LOGICAL :: leta
    INTEGER :: num_accepted
    INTEGER :: ii,jj,ppp,nsamples
    INTEGER :: mpierr

    ppp = SIZE(eta)
    nsamples = ppp*params%mpi_params%nmpi

    deta = (pdf_params%max_pitch_angle - pdf_params%min_pitch_angle)/100.0_rp

    if (pdf_params%min_pitch_angle.GE.korc_zero) then
       do jj=1_idef,INT(minmax_buffer_size,idef)
          minmax = pdf_params%min_pitch_angle -  REAL(jj,rp)*deta
          if (minmax.GT.0.0_rp) then
             min_pitch_angle = minmax
          end if
       end do
    else
       min_pitch_angle = pdf_params%min_pitch_angle
    end if

    do jj=1_idef,INT(minmax_buffer_size,idef)
       minmax = pdf_params%max_pitch_angle + REAL(jj,rp)*deta
       if (minmax.LE.90.0_rp) then
          max_pitch_angle = minmax
       end if
    end do

    !	write(output_unit_write,*) min_pitch_angle,max_pitch_angle

    if (params%mpi_params%rank.EQ.0_idef) then
       ALLOCATE(eta_samples(nsamples))! Number of samples to distribute among all MPI processes
       ALLOCATE(eta_tmp(nsamples))! Number of samples to distribute among all MPI processes

       !* * * Transient * * *!
       call RANDOM_SEED()
       call RANDOM_NUMBER(rand_unif)
       eta_buffer = pdf_params%min_pitch_angle + (pdf_params%max_pitch_angle - pdf_params%min_pitch_angle)*rand_unif

       ii=2_idef
       do while (ii .LE. 1000_idef)
          eta_test = eta_buffer + random_norm(0.0_rp,deta)
          do while ((ABS(eta_test) .GT. pdf_params%max_pitch_angle).OR.(ABS(eta_test) .LT. pdf_params%min_pitch_angle))
             eta_test = eta_buffer + random_norm(0.0_rp,deta)
          end do

          ratio = fRE(eta_test,pdf_params%po)/fRE(eta_buffer,pdf_params%po)

          if (ratio .GE. 1.0_rp) then
             eta_buffer = eta_test
             ii = ii + 1_idef
          else
             call RANDOM_NUMBER(rand_unif)
             if (rand_unif .LT. ratio) then
                eta_buffer = eta_test
                ii = ii + 1_idef
             end if
          end if
       end do
       !* * * Transient * * *!


       eta_tmp(1) = eta_buffer

       call RANDOM_SEED()
       call RANDOM_NUMBER(rand_unif)

       num_accepted = 0_idef
       do while(num_accepted.LT.nsamples)
          ii=2_idef
          do while (ii .LE. nsamples)
             eta_test = eta_tmp(ii-1) + random_norm(0.0_rp,deta)
             do while ((ABS(eta_test) .GT. max_pitch_angle).OR.(ABS(eta_test) .LT. min_pitch_angle))
                eta_test = eta_tmp(ii-1) + random_norm(0.0_rp,deta)
             end do

             ratio = fRE(eta_test,pdf_params%po)/fRE(eta_tmp(ii-1),pdf_params%po)

             if (ratio .GE. 1.0_rp) then
                eta_tmp(ii) = eta_test
                ii = ii + 1_idef
             else
                call RANDOM_NUMBER(rand_unif)
                if (rand_unif .LT. ratio) then
                   eta_tmp(ii) = eta_test
                   ii = ii + 1_idef
                end if
             end if
          end do

          eta_tmp = ABS(eta_tmp)

          ii = 1_idef
          do while ( (ii.LT.nsamples).AND.(num_accepted.LT.nsamples) )
             leta = (eta_tmp(ii).LE.pdf_params%max_pitch_angle).AND.(eta_tmp(ii).GE.pdf_params%min_pitch_angle)
             if (leta) then
                num_accepted = num_accepted + 1_idef
                eta_samples(num_accepted) = eta_tmp(ii)
             end if
             ii = ii + 1_idef
          end do
       end do

       etao = SUM(eta_samples)/nsamples
    end if

    CALL MPI_SCATTER(eta_samples,ppp,MPI_REAL8,eta,ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

    CALL MPI_BCAST(etao,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

    if (params%mpi_params%rank.EQ.0_idef) then
       DEALLOCATE(eta_samples)
       DEALLOCATE(eta_tmp)
    end if

  END SUBROUTINE sample_distribution


  SUBROUTINE save_params(params)
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
       filename = TRIM(params%path_to_outputs) // "simple_equilibrium_pdf.h5"
       call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)

       gname = "pdf_params"
       call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

       dset = TRIM(gname) // "/max_pitch_angle"
       attr = "Maximum pitch angle in avalanche PDF (degrees)"
       call save_to_hdf5(h5file_id,dset,pdf_params%max_pitch_angle,attr)

       dset = TRIM(gname) // "/min_pitch_angle"
       attr = "Minimum pitch angle in avalanche PDF (degrees)"
       call save_to_hdf5(h5file_id,dset,pdf_params%min_pitch_angle,attr)

       dset = TRIM(gname) // "/Zeff"
       attr = "Effective atomic number of ions."
       call save_to_hdf5(h5file_id,dset,pdf_params%Zeff,attr)

       dset = TRIM(gname) // "/E"
       attr = "Parallel electric field in (Ec)"
       call save_to_hdf5(h5file_id,dset,pdf_params%E,attr)

       dset = TRIM(gname) // "/Bo"
       attr = "Characteristic magnetic field in T (in case of using PR)"
       call save_to_hdf5(h5file_id,dset,pdf_params%Bo,attr)

       dset = TRIM(gname) // "/lambda"
       attr = "Characteristic wavelength in m (in case of using PR)"
       call save_to_hdf5(h5file_id,dset,pdf_params%lambda,attr)

       call h5gclose_f(group_id, h5error)

       call h5fclose_f(h5file_id, h5error)
    end if
  END SUBROUTINE save_params

END MODULE korc_simple_equilibrium_pdf
