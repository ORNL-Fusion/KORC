module korc_finalize
  !! @note Module containing subroutines to terminate parallel
  !! communications and free memory.
  use korc_types
  use korc_fields
  use korc_profiles
  use korc_hpc

  IMPLICIT NONE

  PUBLIC :: finalize_communications,&
       deallocate_variables

CONTAINS

  
  subroutine finalize_communications(params)
    !! @note Interface to function that finalizes MPI communications.
    !! See [[korc_hpc]].
    TYPE(KORC_PARAMS), INTENT(IN) :: params
    !! Core KORC simulation parameters.  

    call finalize_mpi(params)
  end subroutine finalize_communications


  subroutine deallocate_variables(params,F,P,spp)
    !! @note Subroutine to free allocatable simulation variables.    
    TYPE(KORC_PARAMS), INTENT(INOUT) 			:: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(INOUT) 			:: F
    !! An instance of KORC's derived type FIELDS containing all the
    !! information about the fields used in the simulation. See
    !! [[korc_types]] and [[korc_fields]].
    TYPE(PROFILES), INTENT(INOUT)              :: P
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
    !! An instance of KORC's derived type SPECIES containing all the
    !! information of different electron species. See [[korc_types]].
    INTEGER 						:: ii
    !! Iterator of the spp array
    
    DEALLOCATE(params%outputs_list)

    do ii=1_idef,params%num_species
       DEALLOCATE(spp(ii)%vars%X)
       DEALLOCATE(spp(ii)%vars%V)
       DEALLOCATE(spp(ii)%vars%Rgc)
       DEALLOCATE(spp(ii)%vars%Y)
       DEALLOCATE(spp(ii)%vars%Yborn)
       DEALLOCATE(spp(ii)%vars%E)
       DEALLOCATE(spp(ii)%vars%B)
       DEALLOCATE(spp(ii)%vars%PSI_P)
       DEALLOCATE(spp(ii)%vars%ne)
       DEALLOCATE(spp(ii)%vars%ni)
       if (ALLOCATED(spp(ii)%vars%nimp)) DEALLOCATE(spp(ii)%vars%nimp)
       DEALLOCATE(spp(ii)%vars%Te)
       DEALLOCATE(spp(ii)%vars%Zeff)
       DEALLOCATE(spp(ii)%vars%g)
       DEALLOCATE(spp(ii)%vars%eta)
       DEALLOCATE(spp(ii)%vars%mu)
       DEALLOCATE(spp(ii)%vars%Prad)
       DEALLOCATE(spp(ii)%vars%flagCon)
       DEALLOCATE(spp(ii)%vars%flagCol)
       DEALLOCATE(spp(ii)%vars%flagRE)
       DEALLOCATE(spp(ii)%vars%initLCFS)
       DEALLOCATE(spp(ii)%vars%wt)

       if (params%orbit_model(1:2).eq.'GC') then
          DEALLOCATE(spp(ii)%vars%Y0)
          DEALLOCATE(spp(ii)%vars%Y1)
          DEALLOCATE(spp(ii)%vars%V0)
          DEALLOCATE(spp(ii)%vars%k1)
          DEALLOCATE(spp(ii)%vars%k2)
          DEALLOCATE(spp(ii)%vars%k3)
          DEALLOCATE(spp(ii)%vars%k4)
          DEALLOCATE(spp(ii)%vars%RHS)
       end if

       if (ALLOCATED(spp(ii)%BMC_ra)) DEALLOCATE(spp(ii)%BMC_ra)
       if (ALLOCATED(spp(ii)%BMC_nRE)) DEALLOCATE(spp(ii)%BMC_nRE)
    end do

    DEALLOCATE(spp)

    call DEALLOCATE_FIELDS_ARRAYS(F)
    call DEALLOCATE_PROFILES_ARRAYS(P)
  end subroutine deallocate_variables

end module korc_finalize
