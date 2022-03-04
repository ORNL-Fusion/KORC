!> @brief Module that contains subroutines and functions to sample various energy distributions.
MODULE korc_energy_pdfs
  USE korc_types
  USE korc_constants
  USE korc_HDF5
  USE korc_hpc
  USE korc_input
  
  IMPLICIT NONE


  !> @brief KORC derived type that contains information about a given Gamma distribution function @f$f_\Gamma(x,\kappa,\theta)@f$.
  !! @details We write a given Gamma distribution function in terms of its shape factor @f$\kappa@f$ and scale factor
  !! @f$\theta@f$, so that:
  !!
  !!
  !! @f$f_\Gamma(x,\kappa,\theta) = \frac{1}{\Gamma(\kappa) \theta^\kappa}x^{\kappa-1}\exp{\left(-x/\theta\right)}@f$.
  TYPE, PRIVATE :: GAMMA_PARAMS
     REAL(rp) :: min_energy !< Minimum energy of sampled @f$f_\Gamma(x,\kappa,\theta)@f$ in MeV.
     REAL(rp) :: max_energy !< Maximum energy of sampled @f$f_\Gamma(x,\kappa,\theta)@f$ in MeV.
     REAL(rp) :: min_p !< Minimum momentum of sampled @f$f_\Gamma(x,\kappa,\theta)@f$.
     REAL(rp) :: max_p !< Maximum momentum of sampled @f$f_\Gamma(x,\kappa,\theta)@f$.
     REAL(rp) :: k !< Shape factor @f$\kappa@f$.
     REAL(rp) :: t !< Scale factor @f$\theta@f$.
  END TYPE GAMMA_PARAMS

  TYPE(GAMMA_PARAMS), PRIVATE 	:: gamma_pdf_params !< An instance of the KORC derived type GAMMA_PARAMS
  REAL(rp), PRIVATE, PARAMETER 	:: co = (C_E*1.0E6)/(C_ME*C_C**2) !< Constant with units to transform @f$\mathcal{E}@f$ in
  !! @f$f_\Gamma(\mathcal{E},\kappa,\theta)@f$ from MeV to Joules.
  REAL(rp), PRIVATE, PARAMETER 	:: minmax_buffer_size = 10.0_rp !< This is the size of the buffer zone in each direction when
  !!using a Metropolis-Hastings method to sample a distribution.

  PUBLIC 	:: get_gamma_distribution
  PRIVATE :: initialize_gamma_params,&
       save_gamma_params,&
       sample_gamma_distribution,&
       deg2rad,&
       fRE,&
       random_norm,&
       fGamma

CONTAINS

  !> @brief Subroutine that contains calls to subroutine to generate a gamma distribution for the energy distribution of a given
  !! species in the simulation.
  !!
  !! @param[in] params Core KORC simulation parameters.
  !! @param[in,out] g Relativistic gamma factor @f$\gamma@f$ of the particles in a given species in the simulation. These are so that, they follow
  !! a Gamma distribution in energy. The parameters of the Gamma distributions are given by the user.
  !! @param[out] go Mean value of @f$\gamma@f$ of the particles in a given species. used to calculate the minimum required time step to
  !! resolve in detail the full-orbit dynamics of the particles.
  SUBROUTINE get_gamma_distribution(params,g,go)
    TYPE(KORC_PARAMS), INTENT(IN) 						:: params
    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) 	:: g
    REAL(rp), INTENT(OUT) 								:: go

    call initialize_gamma_params(params)

    call save_gamma_params(params)

    call sample_gamma_distribution(params,g,go)
  END SUBROUTINE get_gamma_distribution


  !> @brief Subroutine that reads from the input file the parameters of the Gamma distribution
  !! @f$f_\Gamma(x,\kappa,\theta) = \frac{1}{\Gamma(\kappa) \theta^\kappa}x^{\kappa-1}\exp{\left(-x/\theta\right)}@f$.
  !!
  !! @param[in] params Core KORC simulation parameters.
  !! @param max_energy Maximum energy of sampled @f$f_\Gamma(x,\kappa,\theta)@f$ in MeV.
  !! @param min_energy Minimum energy of sampled @f$f_\Gamma(x,\kappa,\theta)@f$ in MeV.
  !! @param k Shape factor @f$\kappa@f$.
  !! @param t Scale factor @f$\theta@f$.
  SUBROUTINE initialize_gamma_params(params)
    TYPE(KORC_PARAMS), INTENT(IN) 	:: params
    REAL(rp) 						:: max_energy
    REAL(rp) 						:: min_energy
    REAL(rp) 						:: k
    REAL(rp) 						:: t

    !NAMELIST /EnergyGammaPDF/ max_energy,min_energy,k,t

    !open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
    !read(default_unit_open,nml=EnergyGammaPDF)
    !close(default_unit_open)

    gamma_pdf_params%min_energy = min_energy_gamma*C_E ! In Joules
    gamma_pdf_params%max_energy = max_energy_gamma*C_E ! In Joules
    gamma_pdf_params%k = k_gamma
    gamma_pdf_params%t = t_gamma

    gamma_pdf_params%max_p = SQRT((gamma_pdf_params%max_energy/(C_ME*C_C**2))**2 - 1.0_rp) ! In units of mc
    gamma_pdf_params%min_p = SQRT((gamma_pdf_params%min_energy/(C_ME*C_C**2))**2 - 1.0_rp) ! In units of mc
  END SUBROUTINE initialize_gamma_params

  !> @brief Function that converts @f$x@f$ from degrees to radians.
  !!
  !! @param x Angle @f$x@f$ in degrees.
  !! @param deg2rad Angle @f$x@f$ converted to radians.
  FUNCTION deg2rad(x)
    REAL(rp), INTENT(IN) 	:: x
    REAL(rp) 				:: deg2rad

    deg2rad = C_PI*x/180.0_rp
  END FUNCTION deg2rad


  !> @brief Function that calculates the value of the Gamma distribution @f$f_\Gamma(x,\kappa,\theta) =
  !! \frac{1}{\Gamma(\kappa) \theta^\kappa}x^{\kappa-1}\exp{\left(-x/\theta\right)}@f$.
  !!
  !! @param x Variable @f$x@f$ of @f$f_\Gamma(x,\kappa,\theta)@f$.
  !! @param k Shape factor @f$\kappa@f$ of @f$f_\Gamma(x,\kappa,\theta)@f$.
  !! @param t Scale factor @f$\theta@f$ of @f$f_\Gamma(x,\kappa,\theta)@f$.
  !! @param fGamma Computed value of @f$f_\Gamma(x,\kappa,\theta)@f$.
  FUNCTION fGamma(x,k,t)
    REAL(rp), INTENT(IN) 	:: x
    REAL(rp), INTENT(IN) 	:: k
    REAL(rp), INTENT(IN) 	:: t
    REAL(rp)				:: fGamma

    fGamma = x**(k - 1.0_rp)*EXP(-x/t)/(GAMMA(k)*t**k)
  END FUNCTION fGamma


  !> Evaluation of the energy distribution function @f$f_{RE}(\mathcal{E})@f$ of runaway electrons as function of the normalized momentum
  !! @f$p' = p/m_ec@f$. Here, @f$p'@f$ is the normalized momentum and @f$m_e@f$ and @f$c@f$ are the electron mass and the speed of light.
  !!
  !! @param p Normalized momentum @f$p' = p/m_ec@f$ of a given electron in the simulation.
  !! @param fRE Computed value of the energy distribution function of runaway electrons.
  !! @param Eo Normalized energy @f$\mathcal{E}' = \sqrt{1 + p'}@f$ of the the electron with normalized momentum @f$p'@f$.
  FUNCTION fRE(p)
    REAL(rp), INTENT(IN) 	:: p ! momentum in units of mc
    REAL(rp) 				:: fRE
    REAL(rp) 				:: Eo ! In units of mc^2

    Eo = SQRT(p**2.0_rp + 1.0_rp)

    fRE = fGamma(Eo,gamma_pdf_params%k,gamma_pdf_params%t*co)
  END FUNCTION fRE


  !> @brief Gaussian random number generator.
  !! @details This function returns a deviate of a Gaussian distribution @f$f_G(x;\mu,\sigma) = \frac{1}{\sigma\sqrt{2\pi}} \exp{\left( -(x-\mu)^2/2\sigma^2 \right)}@f$,
  !! with mean @f$\mu@f$, and standard deviation @f$\sigma@f$.
  !!
  !! We use the Inverse Transform Sampling Method for sampling @f$x@f$. With this method we get @f$x = \sqrt{-2\log{(1-y)}}\cos(2\pi z)@f$,
  !! where @f$y@f$ and @f$z@f$ are uniform random numbers in the interval @f$[0,1]@f$.
  !!
  !! @param[in] mu Mean value @f$\mu@f$ of the Gaussian distribution.
  !! @param[in] mu Standard deviation @f$\sigma@f$ of the Gaussian distribution.
  !! @param random_norm Sampled number @f$x@f$ from the Gaussian distribution @f$f_G(x;\mu,\sigma)@f$.
  !! @param rand1 Uniform random number in the interval @f$[0,1]@f$.
  !! @param rand2 Uniform random number in the interval @f$[0,1]@f$.
  FUNCTION random_norm(mean,sigma)
    REAL(rp), INTENT(IN) 	:: mean
    REAL(rp), INTENT(IN) 	:: sigma
    REAL(rp) 				:: random_norm
    REAL(rp) 				:: rand1
    REAL(rp) 				:: rand2

    call RANDOM_NUMBER(rand1)
    call RANDOM_NUMBER(rand2)

    random_norm = SQRT(-2.0_rp*LOG(1.0_rp-rand1))*COS(2.0_rp*C_PI*rand2);
  END FUNCTION random_norm


  !> @brief Subroutine that samples a Gamma distribution representing the runaways' (marginal) energy distribution function.
  !! @details This subroutine uses the Metropolis-Hastings method for sampling the Gamma distribution representing the runaways'
  !! (marginal) energy distribution function. Unlike the typical Metropolis-Hasting method, after setting the boundaries of the region
  !! we want to sample, we perform a sampling in a larger region that contains the original sampling area plus a buffer region.
  !! After finishing the first sampling, we only keep the particles in the original sampling region, the particles in the p_buffer
  !! are sampled again until all of them lie within the original sampling region. This method ensures that the boundaries are
  !! well sampled.
  !!
  !! @param[in] params Core KORC simulation parameters.
  !! @param[in,out] g Relativistic gamma factor @f$\gamma@f$ of the particles in a given species in the simulation. These are so that,
  !! they follow a Gamma distribution in energy. The parameters of the Gamma distributions are given by the user.
  !! @param[out] go Mean value of @f$\gamma@f$ of the particles in a given species. used to calculate the minimum required time step to
  !! resolve in detail the full-orbit dynamics of the particles.
  !! @param p Sampled normalized momentum @f$p' = p/m_ec@f$ of particles in a given particle species.
  !! @param p_buffer Size along the momentum axis of the buffer used in the Metropolis-Hastings method.
  !! @param p_test The test value of the normalized momentum used in the Metropolis-Hastings method.
  !! @param ratio Ratio of probabilities used to determine when a move in during the sampling is kept as part of the sampled chain.
  !! @param rand_unif A deviate of a uniform random distribution in the interval @f$[0,1]@f$.
  !! @param p_samples Temporary array to keep the sampled normalized momentum.
  !! @param deta Step size along the pitch-angle direction of the random walk used in the Metropolis-Hastings sampling.
  !! @param dp  Step size along the momentum direction of the random walk used in the Metropolis-Hastings sampling.
  !! @param ii Iterator.
  !! @param ppp Number of particles per MPI processes.
  !! @param nsamples Number of total samples in the initial condition of a given simulation.
  !! @param mpierr MPI error status.
  SUBROUTINE sample_gamma_distribution(params,g,go)
    TYPE(KORC_PARAMS), INTENT(IN) 						:: params
    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) 	:: g
    REAL(rp), INTENT(OUT) 								:: go
    REAL(rp), DIMENSION(:), ALLOCATABLE 				:: p
    REAL(rp) 											:: p_buffer
    REAL(rp) 											:: p_test
    REAL(rp) 											:: ratio
    REAL(rp) 											:: rand_unif
    REAL(rp), DIMENSION(:), ALLOCATABLE 				:: p_samples
    REAL(rp) 											:: deta
    REAL(rp) 											:: dp
    INTEGER 											:: ii
    INTEGER 											:: ppp
    INTEGER 											:: nsamples
    INTEGER 											:: mpierr

    ppp = SIZE(g)
    nsamples = ppp*params%mpi_params%nmpi
    ALLOCATE(p(ppp))

    dp = 1.0_rp

    if (params%mpi_params%rank.EQ.0_idef) then
       ALLOCATE(p_samples(nsamples))! Number of samples to distribute among all MPI processes

       call RANDOM_SEED()
       call RANDOM_NUMBER(rand_unif)
       p_buffer = gamma_pdf_params%min_p + (gamma_pdf_params%max_p - gamma_pdf_params%min_p)*rand_unif

       ii=2_idef
       do while (ii .LE. 1000000_idef)
          p_test = p_buffer + random_norm(0.0_rp,dp)
          do while ((p_test.LT.gamma_pdf_params%min_p).OR.(p_test.GT.gamma_pdf_params%max_p))
             p_test = p_buffer + random_norm(0.0_rp,dp)
          end do

          ratio = fRE(p_test)/fRE(p_buffer)

          if (ratio .GE. 1.0_rp) then
             p_buffer = p_test
             ii = ii + 1_idef
          else
             call RANDOM_NUMBER(rand_unif)
             if (rand_unif .LT. ratio) then
                p_buffer = p_test
                ii = ii + 1_idef
             end if
          end if
       end do

       call RANDOM_SEED()
       call RANDOM_NUMBER(rand_unif)
       p_samples(1) = p_buffer

       ii=2_idef
       do while (ii .LE. nsamples)
          p_test = p_samples(ii-1) + random_norm(0.0_rp,dp)
          do while ((p_test.LT.gamma_pdf_params%min_p).OR.(p_test.GT.gamma_pdf_params%max_p))
             p_test = p_samples(ii-1) + random_norm(0.0_rp,dp)
          end do

          ratio = fRE(p_test)/fRE(p_samples(ii-1))

          if (ratio .GE. 1.0_rp) then
             p_samples(ii) = p_test
             ii = ii + 1_idef
          else
             call RANDOM_NUMBER(rand_unif)
             if (rand_unif .LT. ratio) then
                p_samples(ii) = p_test
                ii = ii + 1_idef
             end if
          end if
       end do

       go = SUM(SQRT(1.0_rp + p_samples**2))/nsamples
    end if

    CALL MPI_SCATTER(p_samples,ppp,MPI_REAL8,p,ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

    CALL MPI_BCAST(go,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

    g = SQRT(1.0_rp + p**2)

    DEALLOCATE(p)
    if (params%mpi_params%rank.EQ.0_idef) then
       DEALLOCATE(p_samples)
    end if
  END SUBROUTINE sample_gamma_distribution


  !> @brief Surboutine that saves the Gamma distribution parameters to the HDF5 file <i>gamma_distribution_parameters.h5</i>.
  !!
  !! @param[in] params Core KORC simulation parameters.
  !! @param filename String containing the name of the HDF5 file.
  !! @param gname String containing the group name of a set of KORC parameters.
  !! @param attr_array An 1-D array with attributes of 1-D real or integer arrays that are passed to KORC interfaces of HDF5 I/O subroutines.
  !! @param dset Name of data set to be saved to file.
  !! @param attr A single attributes of real or integer data that is passed to KORC interfaces of HDF5 I/O subroutines.
  !! @param h5file_id HDF5 file identifier.
  !! @param group_id HDF5 group identifier.
  !! @param h5error HDF5 error status.
  !! @param units Temporary variable used to add physical units to KORC parameters.
  SUBROUTINE save_gamma_params(params)
    TYPE(KORC_PARAMS), INTENT(IN) 							:: params
    CHARACTER(MAX_STRING_LENGTH) 							:: filename
    CHARACTER(MAX_STRING_LENGTH) 							:: gname
    CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
    CHARACTER(MAX_STRING_LENGTH) 							:: dset
    CHARACTER(MAX_STRING_LENGTH) 							:: attr
    INTEGER(HID_T) 											:: h5file_id
    INTEGER(HID_T) 											:: group_id
    INTEGER 												:: h5error
    REAL(rp) 												:: units

    if (params%mpi_params%rank .EQ. 0) then
       filename = TRIM(params%path_to_outputs) // "gamma_distribution_parameters.h5"
       call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)

       gname = "params"
       call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

       dset = TRIM(gname) // "/min_energy"
       attr = "Minimum energy in avalanche PDF (eV)"
       units = 1.0_rp/C_E
       call save_to_hdf5(h5file_id,dset,units*gamma_pdf_params%min_energy,attr)

       dset = TRIM(gname) // "/max_energy"
       attr = "Maximum energy in avalanche PDF (eV)"
       units = 1.0_rp/C_E
       call save_to_hdf5(h5file_id,dset,units*gamma_pdf_params%max_energy,attr)

       dset = TRIM(gname) // "/max_p"
       attr = "Maximum momentum in avalanche PDF (me*c^2)"
       call save_to_hdf5(h5file_id,dset,gamma_pdf_params%max_p,attr)

       dset = TRIM(gname) // "/min_p"
       attr = "Maximum momentum in avalanche PDF (me*c^2)"
       call save_to_hdf5(h5file_id,dset,gamma_pdf_params%min_p,attr)

       dset = TRIM(gname) // "/k"
       attr = "Shape factor"
       call save_to_hdf5(h5file_id,dset,gamma_pdf_params%k,attr)

       dset = TRIM(gname) // "/t"
       attr = "Scale factor"
       call save_to_hdf5(h5file_id,dset,gamma_pdf_params%t,attr)

       call h5gclose_f(group_id, h5error)

       call h5fclose_f(h5file_id, h5error)
    end if
  END SUBROUTINE save_gamma_params

END MODULE korc_energy_pdfs
