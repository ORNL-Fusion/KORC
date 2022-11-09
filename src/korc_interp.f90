module korc_interp
  !! @note Module containing functions and subroutines for performing
  !! interpolations using the PSPLINE library. @endnote
  !! For a detailed documentation of the PSPLINE library we refer the
  !! user to "https://w3.pppl.gov/ntcc/PSPLINE/".
  use korc_types
  use korc_coords
  use korc_rnd_numbers
  use korc_hpc

#ifdef PSPLINE
  use EZspline_obj	! psplines module
  use EZspline		! psplines module
#endif

#ifdef FIO
  use korc_fio
#endif

  !$ use OMP_LIB

  IMPLICIT NONE

#ifdef PSPLINE
#ifdef DOUBLE_PRECISION


  TYPE, PRIVATE :: KORC_3D_FIELDS_INTERPOLANT
     !! @note Derived type containing 3-D PSPLINE interpolants for
     !! cylindrical components of vector fields
     !! \(\mathbf{F}(R,\phi,Z) = F_R\hat{e}_R + F_\phi\hat{e}_phi +
     !! F_Z\hat{e}_Z\). Real precision of 8 bytes. @endnote
     TYPE(EZspline3)    :: A     !! Interpolant of a scalar field \(A(R,Z)\).
     TYPE(EZspline3)    :: R
     !! Interpolant of \(F_R(R,\phi,Z)\).
     TYPE(EZspline3)    :: PHI
     !! Interpolant of \(F_\phi(R,\phi,Z)\).
     TYPE(EZspline3)    :: Z
     !! Interpolant of \(F_Z(R,\phi,Z)\).

     INTEGER               :: NR
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NPHI
     !! Size of mesh containing the field data along the \(\phi\)-axis.
     INTEGER               :: NZ
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCSPHI = (/ -1, -1 /)
     !! Periodic boundary condition for the interpolants at both
     !! ends of the \(\phi\) direction.
     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(Z\) direction.
  END TYPE KORC_3D_FIELDS_INTERPOLANT

  TYPE, PRIVATE :: KORC_2X1T_FIELDS_INTERPOLANT
     !! @note Derived type containing 3-D PSPLINE interpolants for
     !! cylindrical components of vector fields
     !! \(\mathbf{F}(R,\phi,Z) = F_R\hat{e}_R + F_\phi\hat{e}_phi +
     !! F_Z\hat{e}_Z\). Real precision of 8 bytes. @endnote
     TYPE(EZspline3)    :: A
     !! Interpolant of a scalar field \(A(R,Z)\).
     TYPE(EZspline3)    :: R
     !! Interpolant of \(F_R(R,\phi,Z)\).
     TYPE(EZspline3)    :: T
     !! Interpolant of \(F_\phi(R,\phi,Z)\).
     TYPE(EZspline3)    :: Z
     !! Interpolant of \(F_Z(R,\phi,Z)\).

     INTEGER               :: NR
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NT
     !! Size of mesh containing the field data along the \(\phi\)-axis.
     INTEGER               :: NZ
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCST = (/ 0, 0 /)
     !! Periodic boundary condition for the interpolants at both
     !! ends of the \(\phi\) direction.
     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(Z\) direction.
  END TYPE KORC_2X1T_FIELDS_INTERPOLANT

  TYPE, PRIVATE :: KORC_2DX_FIELDS_INTERPOLANT
     !! @note Derived type containing 2-D PSPLINE interpolants for
     !! cylindrical components of vector fields \(\mathbf{F}(R,Z) =
     !! F_R\hat{e}_R + F_\phi\hat{e}_phi+ F_Z\hat{e}_Z\).
     !! Real precision of 8 bytes. @endnote
     TYPE(EZspline2)    :: A
     !! Interpolant of a scalar field \(A(R,Z)\).
     TYPE(EZspline2)    :: X
     !! Interpolant of \(F_X(R,Z)\).
     TYPE(EZspline2)    :: Y
     !! Interpolant of \(F_Y(R,Z)\).
     TYPE(EZspline2)    :: Z
     !! Interpolant of \(F_Z(R,Z)\).

     INTEGER               :: NR
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NZ
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(Z\) direction.
  END TYPE KORC_2DX_FIELDS_INTERPOLANT

  TYPE, PRIVATE :: KORC_2D_FIELDS_INTERPOLANT
     !! @note Derived type containing 2-D PSPLINE interpolants for
     !! cylindrical components of vector fields \(\mathbf{F}(R,Z) =
     !! F_R\hat{e}_R + F_\phi\hat{e}_phi+ F_Z\hat{e}_Z\).
     !! Real precision of 8 bytes. @endnote
     TYPE(EZspline2)    :: A
     !! Interpolant of a scalar field \(A(R,Z)\).
     TYPE(EZspline2)    :: R
     !! Interpolant of \(F_R(R,Z)\).
     TYPE(EZspline2)    :: PHI
     !! Interpolant of \(F_\phi(R,Z)\).
     TYPE(EZspline2)    :: Z
     !! Interpolant of \(F_Z(R,Z)\).

     INTEGER               :: NR
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NZ
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(Z\) direction.
  END TYPE KORC_2D_FIELDS_INTERPOLANT

  TYPE, PRIVATE :: KORC_1D_FIELDS_INTERPOLANT
     !! @note Derived type containing 2-D PSPLINE interpolants for
     !! cylindrical components of vector fields \(\mathbf{F}(R,Z) =
     !! F_R\hat{e}_R + F_\phi\hat{e}_phi+ F_Z\hat{e}_Z\).
     !! Real precision of 8 bytes. @endnote
     TYPE(EZspline1)    :: A
     !! Interpolant of a scalar field \(A(R,Z)\).
     TYPE(EZspline1)    :: R
     !! Interpolant of \(F_R(R,Z)\).
     TYPE(EZspline1)    :: PHI
     !! Interpolant of \(F_\phi(R,Z)\).
     TYPE(EZspline1)    :: Z
     !! Interpolant of \(F_Z(R,Z)\).

     INTEGER               :: Nrm
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER, DIMENSION(2) :: BCSrm = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(R\) direction.
     INTEGER               :: NPSIP
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER, DIMENSION(2) :: BCSPSIP = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(R\) direction.
  END TYPE KORC_1D_FIELDS_INTERPOLANT

  TYPE, PRIVATE :: KORC_3D_PROFILES_INTERPOLANT
     !! @note Derived type containing 3-D PSPLINE interpolants for cylindrical
     !! components of the density \(n_e(R,\phi,Z)\),
     !! temperature \(T_e(R,\phi,Z)\), and effective charge number
     !! \(Z_{eff}(R,\phi,Z)\) profiles. Real precision of 8 bytes. @endnote
     TYPE(EZspline3)    :: ne
     !! Interpolant of background electron density \(n_e(R,\phi,Z)\).
     TYPE(EZspline3)    :: Te
     !! Interpolant of background electron temperature \(T_e(R,\phi,Z)\).
     TYPE(EZspline3)    :: Zeff
     !! Interpolant of effective charge number \(Z_{eff}(R,\phi,Z)\).

     INTEGER               :: NR
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NPHI
     !! Size of mesh containing the field data along the \(\phi\)-axis.
     INTEGER               :: NZ
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends
     !! of the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCSPHI = (/ -1, -1 /)
     !! Periodic boundary condition for the interpolants at both ends of
     !! the \(\phi\) direction.
     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends
     !! of the \(Z\) direction.
  END TYPE KORC_3D_PROFILES_INTERPOLANT


  TYPE, PRIVATE :: KORC_2D_PROFILES_INTERPOLANT
     !! @note Derived type containing 2-D PSPLINE interpolants for cylindrical
     !! components of the density \(n_e(R,Z)\), temperature \(T_e(R,Z)\), and
     !! effective charge number \(Z_{eff}(R,Z)\) profiles.
     !! Real precision of 8 bytes. @endnote
     TYPE(EZspline2)    :: ne
     !! Interpolant of background electron density \(n_e(R,Z)\).
     TYPE(EZspline2)    :: Te
     !! Interpolant of background electron temperature \(T_e(R,Z)\).
     TYPE(EZspline2)    :: Zeff
     !! Interpolant of effective charge number \(Z_{eff}(R,Z)\)
     TYPE(EZspline2)    :: RHON
     TYPE(EZspline2)    :: nRE
     TYPE(EZspline2)    :: nAr0
     TYPE(EZspline2)    :: nAr1
     TYPE(EZspline2)    :: nAr2
     TYPE(EZspline2)    :: nAr3
     TYPE(EZspline2)    :: nD
     TYPE(EZspline2)    :: nD1

     INTEGER               :: NR
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NZ
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(Z\) direction.
  END TYPE KORC_2D_PROFILES_INTERPOLANT

  TYPE, PRIVATE :: KORC_2D_HOLLMANN_INTERPOLANT
     !! @note Derived type containing 2-D PSPLINE interpolants for cylindrical
     !! components of the density \(n_e(R,Z)\), temperature \(T_e(R,Z)\), and
     !! effective charge number \(Z_{eff}(R,Z)\) profiles.
     !! Real precision of 8 bytes. @endnote
     TYPE(EZspline2)    :: fRE_E
     TYPE(EZspline2)    :: fRE_pitch

     INTEGER               :: NE
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NRHO
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCRHO = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCE = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(Z\) direction.
  END TYPE KORC_2D_HOLLMANN_INTERPOLANT


#elif SINGLE_PRECISION


  TYPE, PRIVATE :: KORC_3D_FIELDS_INTERPOLANT
     !! @note Derived type containing 3-D PSPLINE interpolants for cylindrical
     !! components of vector fields \(\mathbf{F}(R,\phi,Z) = F_R\hat{e}_R +
     !! F_\phi\hat{e}_phi + F_Z\hat{e}_Z\).
     !! Real precision of 4 bytes. @endnote
     TYPE(EZspline3_r4)    :: R
     !! Interpolant of \(F_R(R,\phi,Z)\).
     TYPE(EZspline3_r4)    :: PHI
     !! Interpolant of \(F_\phi(R,\phi,Z)\).
     TYPE(EZspline3_r4)    :: Z
     !! Interpolant of \(F_Z(R,\phi,Z)\).

     INTEGER               :: NR
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NPHI
     !! Size of mesh containing the field data along the \(\phi\)-axis.
     INTEGER               :: NZ
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends
     !! of the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCSPHI = (/ -1, -1 /)
     !! Periodic boundary condition for the interpolants at both ends of
     !! the \(\phi\) direction.
     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends of
     !! the \(Z\) direction.
  END TYPE KORC_3D_FIELDS_INTERPOLANT

    TYPE, PRIVATE :: KORC_2X1T_FIELDS_INTERPOLANT
     !! @note Derived type containing 3-D PSPLINE interpolants for cylindrical
     !! components of vector fields \(\mathbf{F}(R,\phi,Z) = F_R\hat{e}_R +
     !! F_\phi\hat{e}_phi + F_Z\hat{e}_Z\).
     !! Real precision of 4 bytes. @endnote
     TYPE(EZspline3_r4)    :: R
     !! Interpolant of \(F_R(R,\phi,Z)\).
     TYPE(EZspline3_r4)    :: T
     !! Interpolant of \(F_\phi(R,\phi,Z)\).
     TYPE(EZspline3_r4)    :: Z
     !! Interpolant of \(F_Z(R,\phi,Z)\).

     INTEGER               :: NR
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NT
     !! Size of mesh containing the field data along the \(\phi\)-axis.
     INTEGER               :: NZ
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends
     !! of the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCST = (/ 0, 0 /)
     !! Periodic boundary condition for the interpolants at both ends of
     !! the \(\phi\) direction.
     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends of
     !! the \(Z\) direction.
  END TYPE KORC_2X1T_FIELDS_INTERPOLANT

  TYPE, PRIVATE :: KORC_2D_FIELDS_INTERPOLANT
     !! @note Derived type containing 2-D PSPLINE interpolants for cylindrical
     !! components of vector fields \(\mathbf{F}(R,Z) = F_R\hat{e}_R +
     !! F_\phi\hat{e}_phi+ F_Z\hat{e}_Z\).
     !! Real precision of 4 bytes. @endnote
     TYPE(EZspline2_r4)    :: A
     !! Interpolant of a scalar field \(A(R,Z)\).
     TYPE(EZspline2_r4)    :: R
     !! Interpolant of \(F_R(R,Z)\).
     TYPE(EZspline2_r4)    :: PHI
     !! Interpolant of \(F_\phi(R,Z)\).
     TYPE(EZspline2_r4)    :: Z
     !! Interpolant of \(F_Z(R,Z)\).

     INTEGER               :: NR
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NZ
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends
     !! of the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends
     !! of the \(Z\) direction.
  END TYPE KORC_2D_FIELDS_INTERPOLANT

  TYPE, PRIVATE :: KORC_1D_FIELDS_INTERPOLANT
     !! @note Derived type containing 2-D PSPLINE interpolants for
     !! cylindrical components of vector fields \(\mathbf{F}(R,Z) =
     !! F_R\hat{e}_R + F_\phi\hat{e}_phi+ F_Z\hat{e}_Z\).
     !! Real precision of 8 bytes. @endnote
     TYPE(EZspline1_r4)    :: A
     !! Interpolant of a scalar field \(A(R,Z)\).
     TYPE(EZspline1_r4)    :: R
     !! Interpolant of \(F_R(R,Z)\).
     TYPE(EZspline1_r4)    :: PHI
     !! Interpolant of \(F_\phi(R,Z)\).
     TYPE(EZspline1_r4)    :: Z
     !! Interpolant of \(F_Z(R,Z)\).

     INTEGER               :: Nrm
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER, DIMENSION(2) :: BCSrm = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both
     !! ends of the \(R\) direction.
  END TYPE KORC_1D_FIELDS_INTERPOLANT

  TYPE, PRIVATE :: KORC_3D_PROFILES_INTERPOLANT
     !! @note Derived type containing 3-D PSPLINE interpolants for cylindrical
     !! components of the density \(n_e(R,\phi,Z)\),
     !! temperature \(T_e(R,\phi,Z)\), and effective charge number
     !! \(Z_{eff}(R,\phi,Z)\) profiles.
     !! Real precision of 4 bytes. @endnote
     TYPE(EZspline3_r4)    :: ne
     !! Interpolant of background electron density \(n_e(R,\phi,Z)\).
     TYPE(EZspline3_r4)    :: Te
     !! Interpolant of background electron temperature \(T_e(R,\phi,Z)\).
     TYPE(EZspline3_r4)    :: Zeff
     !! Interpolant of effective charge number \(Z_{eff}(R,\phi,Z)\).

     INTEGER               :: NR
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NPHI
     !! Size of mesh containing the field data along the \(\phi\)-axis.
     INTEGER               :: NZ
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends of
     !! the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCSPHI = (/ -1, -1 /)
     !! Periodic boundary condition for the interpolants at both ends of
     !! the \(\phi\) direction.
     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends
     !! of the \(Z\) direction.
  END TYPE KORC_3D_PROFILES_INTERPOLANT



  TYPE, PRIVATE :: KORC_2D_PROFILES_INTERPOLANT
     !! @note Derived type containing 2-D PSPLINE interpolants for
     !! cylindrical components of the density \(n_e(R,Z)\),
     !! temperature \(T_e(R,Z)\), and effective charge number \(Z_{eff}(R,Z)\) profiles.
     !! Real precision of 8 bytes. @endnote
     TYPE(EZspline2_r4)    :: ne
     !! Interpolant of background electron density \(n_e(R,Z)\).
     TYPE(EZspline2_r4)    :: Te
     !! Interpolant of background electron temperature \(T_e(R,Z)\).
     TYPE(EZspline2_r4)    :: Zeff
     !! Interpolant of effective charge number \(Z_{eff}(R,Z)\).

     INTEGER               :: NR
     !! Size of mesh containing the field data along the \(R\)-axis.
     INTEGER               :: NZ
     !! Size of mesh containing the field data along the \(Z\)-axis.
     INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends
     !! of the \(R\) direction.
     INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
     !! Not-a-knot boundary condition for the interpolants at both ends
     !! of the \(Z\) direction.
  END TYPE KORC_2D_PROFILES_INTERPOLANT


#endif


  TYPE, PRIVATE :: KORC_INTERPOLANT_DOMAIN
     !! @note Derived type containing 2-D and 3-D arrays with the information of
     !! the spatial domain where the fields and profiles are known.
     !! This info is used for detecting when a particle is lost, and therefore not
     !! followed anymore. @endnote
     INTEGER(KIND=1), DIMENSION(:), ALLOCATABLE      :: FLAG1D
     !! 2-D array with info of the spatial domain where the axisymmetric fields
     !! and plasma profiles are known.
     INTEGER(KIND=1), DIMENSION(:,:), ALLOCATABLE      :: FLAG2D
     !! 2-D array with info of the spatial domain where the axisymmetric fields
     !! and plasma profiles are known.
     INTEGER(KIND=1), DIMENSION(:,:,:), ALLOCATABLE    :: FLAG3D
     !! 3-D array with info of the spatial domain where the 3-D fields and plasma
     !! profiles are known.
     INTEGER(KIND=1), DIMENSION(:,:), ALLOCATABLE      :: LCFS2D
     !! 2-D array with info of the spatial domain where the axisymmetric fields
     !! and plasma profiles are known.

     REAL(rp)                                          :: Ro
     !! Smaller radial position of the fields and profiles domain.
     REAL(rp)                                          :: Zo
     !! Smaller vertical position of the fields and profiles domain
     REAL(rp)                                          :: To

     REAL(rp)                                          :: Drm
     REAL(rp)                                          :: DPSIP
     REAL(rp)                                          :: DR
     !! Separation between grid points along the radial direction.
     REAL(rp)                                          :: DPHI  !
     ! Separation between grid points along the azimuthal direction.
     REAL(rp)                                          :: DT  !
     ! Separation between grid points along the azimuthal direction.
     REAL(rp)                                          :: DZ
     !! Separation between grid points along the vertical direction.

  END TYPE KORC_INTERPOLANT_DOMAIN

  TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE      :: bfield_2d
  !! An instance of KORC_2D_FIELDS_INTERPOLANT for interpolating
  !! the magnetic field.
  TYPE(KORC_3D_FIELDS_INTERPOLANT), PRIVATE      :: bfield_3d
  !! An instance of KORC_3D_FIELDS_INTERPOLANT for interpolating
  !! the magnetic field.
  TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE      :: b1Refield_2d
  TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE      :: b1Imfield_2d
  TYPE(KORC_2DX_FIELDS_INTERPOLANT), PRIVATE      :: b1Refield_2dx
  TYPE(KORC_2DX_FIELDS_INTERPOLANT), PRIVATE      :: b1Imfield_2dx
  TYPE(KORC_2DX_FIELDS_INTERPOLANT), PRIVATE      :: e1Refield_2dx
  TYPE(KORC_2DX_FIELDS_INTERPOLANT), PRIVATE      :: e1Imfield_2dx
  TYPE(KORC_2X1T_FIELDS_INTERPOLANT), PRIVATE      :: bfield_2X1T
  TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE      :: dbdR_2d
  TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE      :: dbdPHI_2d
  TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE      :: dbdZ_2d
  !! An instance of KORC_2D_FIELDS_INTERPOLANT for interpolating
  !! the magnetic field.
  TYPE(KORC_3D_FIELDS_INTERPOLANT), PRIVATE      :: dbdR_3d
  TYPE(KORC_3D_FIELDS_INTERPOLANT), PRIVATE      :: dbdPHI_3d
  TYPE(KORC_3D_FIELDS_INTERPOLANT), PRIVATE      :: dbdZ_3d
  !! An instance of KORC_3D_FIELDS_INTERPOLANT for interpolating
  !! the magnetic field.
  TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE      :: efield_2d
  !! An instance of KORC_2D_FIELDS_INTERPOLANT for interpolating
  !! the electric field.
  TYPE(KORC_3D_FIELDS_INTERPOLANT), PRIVATE      :: efield_3d
  !! An instance of KORC_3D_FIELDS_INTERPOLANT for interpolating
  !! the electric field.
  TYPE(KORC_1D_FIELDS_INTERPOLANT), PRIVATE      :: efield_SC1d
  !! An instance of KORC_1D_FIELDS_INTERPOLANT for interpolating
  !! the self-consistent electric field.
  TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE      :: gradB_2d
  !! An instance of KORC_2D_FIELDS_INTERPOLANT for interpolating
  !! the magnetic field.
  TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE      :: curlb_2d
  !! An instance of KORC_2D_FIELDS_INTERPOLANT for interpolating
  !! the magnetic field.
  TYPE(KORC_3D_FIELDS_INTERPOLANT), PRIVATE      :: gradB_3d
  !! An instance of KORC_2D_FIELDS_INTERPOLANT for interpolating
  !! the magnetic field.
  TYPE(KORC_3D_FIELDS_INTERPOLANT), PRIVATE      :: curlb_3d
  !! An instance of KORC_2D_FIELDS_INTERPOLANT for interpolating
  !! the magnetic field.
  TYPE(KORC_INTERPOLANT_DOMAIN), PRIVATE         :: fields_domain
  !! An instance of KORC_INTERPOLANT_DOMAIN used for interpolating fields.
  TYPE(KORC_2D_PROFILES_INTERPOLANT), PRIVATE    :: profiles_2d
  !! An instance of KORC_2D_PROFILES_INTERPOLANT for interpolating plasma
  !! profiles.
  TYPE(KORC_3D_PROFILES_INTERPOLANT), PRIVATE    :: profiles_3d
  !! An instance of KORC_3D_PROFILES_INTERPOLANT for interpolating plasma
  !! profiles.
  TYPE(KORC_2D_HOLLMANN_INTERPOLANT), PRIVATE    :: hollmann_2d
  !! An instance of KORC_2D_PROFILES_INTERPOLANT for interpolating plasma
  !! profiles.
  TYPE(KORC_INTERPOLANT_DOMAIN), PRIVATE         :: profiles_domain
  !! An instance of KORC_INTERPOLANT_DOMAIN used for interpolating plasma
  !! profiles.
  INTEGER                                        :: ezerr
  !! Error status during PSPLINE interpolations.

#endif



#ifdef PSPLINE
  PUBLIC :: interp_fields,&
       interp_fields_p,&
       interp_fields_3D_p,&
       interp_FOfields_p,&
       interp_FOfields1_p,&
       interp_bmag_p,&
       interp_collision_p,&
!       interp_fields_FO_p,&
       interp_profiles,&
       initialize_fields_interpolant,&
       initialize_profiles_interpolant,&
       finalize_interpolants,&
       calculate_initial_magnetic_field,&
       calculate_magnetic_field_p,&
       calculate_GCfields_p,&
       calculate_GCfieldswE_p,&
       calculate_GCfields_2x1t_p,&
       calculate_GCfields_p_FS,&
       calculate_2DBdBfields_p,&
       calculate_3DBdBfields_p,&
       calculate_3DBdBfields1_p,&
       sample_poloidal_flux,&
       initialize_SC1D_field_interpolant,&
       add_interp_SCE_p,&
       initialize_SC1D_field_interpolant_FS,&
       add_interp_SCE_p_FS
  PRIVATE :: interp_3D_bfields,&
       interp_2D_bfields,&
       interp_3D_efields,&
       interp_2D_efields,&
       interp_2D_profiles,&
       interp_3D_profiles,&
       check_if_in_fields_domain,&
       check_if_in_LCFS,&
       check_if_in_profiles_domain,&
       check_if_in_profiles_domain_p,&
       check_if_in_fields_domain_p,&
       interp_2D_gradBfields,&
       interp_2D_curlbfields,&
       gradient_2D_Bfields
#else
  PUBLIC :: interp_fields

#endif

CONTAINS

#ifdef PSPLINE
  subroutine initialize_fields_interpolant(params,F)
    !! @note Subroutine that initializes fields interpolants. @endnote
    !! This subroutine initializes either 2-D or 3-D PSPLINE interpolants
    !! using the data of fields in the KORC-dervied-type variable F.
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(INOUT)       :: F
    !! An instance of KORC's derived type FIELDS containing all the information
    !! about the fields used in the simulation.
    !! See [[korc_types]] and [[korc_fields]].
    integer :: ii,jj

    if (((params%field_model(1:8) .EQ. 'EXTERNAL').or. &
         (params%field_eval.eq.'interp')).and. &
         (.not.TRIM(params%field_model).eq.'M3D_C1'.and. &
         .not.TRIM(params%field_model).eq.'NIMROD')) then

       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'("* * * * INITIALIZING FIELDS INTERPOLANT * * * *")')
       end if

       ! * * * * * * * * MAGNETIC FIELD * * * * * * * * !
       if (F%Bflux.or.F%ReInterp_2x1t.or. &
            ((params%orbit_model(1:2).eq.'FO').and.F%axisymmetric_fields)) then

          write(output_unit_write,*) '2D poloidal flux function'

          if(F%ReInterp_2x1t) then

             if (.not.(EZspline_allocated(bfield_2d%A))) then

                bfield_2d%NR = F%dims(1)
                bfield_2d%NZ = F%dims(3)

                ! Initializing poloidal flux interpolant
                call EZspline_init(bfield_2d%A,bfield_2d%NR,bfield_2d%NZ, &
                     bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)

                call EZspline_error(ezerr)

                bfield_2d%A%x1 = F%X%R
                bfield_2d%A%x2 = F%X%Z
             end if

             !write(output_unit_write,'("R",E17.10)') F%X%R
             !write(output_unit_write,'("Z",E17.10)') F%X%Z

             call EZspline_setup(bfield_2d%A, F%PSIp3D(:,F%ind_2x1t,:), &
                  ezerr, .TRUE.)
             call EZspline_error(ezerr)

             !write(output_unit_write,'("bfield_2d%A: ",E17.10)') bfield_2d%A%fspl(1,:,:)

             if (.not.ALLOCATED(fields_domain%FLAG2D)) &
                  ALLOCATE(fields_domain%FLAG2D(bfield_2d%NR,bfield_2d%NZ))

             fields_domain%FLAG2D = F%FLAG3D(:,F%ind_2x1t,:)

             if (.not.ALLOCATED(fields_domain%LCFS2D)) &
                  ALLOCATE(fields_domain%LCFS2D(bfield_2d%NR,bfield_2d%NZ))

             fields_domain%LCFS2D = F%LCFS3D(:,F%ind_2x1t,:)

             fields_domain%DR = ABS(F%X%R(2) - F%X%R(1))
             fields_domain%DZ = ABS(F%X%Z(2) - F%X%Z(1))

             F%Bflux3D = .FALSE.

          else

             if (EZspline_allocated(bfield_2d%A)) &
                  call Ezspline_free(bfield_2d%A, ezerr)

             bfield_2d%NR = F%dims(1)
             bfield_2d%NZ = F%dims(3)

             ! Initializing poloidal flux interpolant
             call EZspline_init(bfield_2d%A,bfield_2d%NR,bfield_2d%NZ, &
                  bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)

             call EZspline_error(ezerr)

             bfield_2d%A%x1 = F%X%R
             bfield_2d%A%x2 = F%X%Z

#if DBG_CHECK
             !write(6,*) 'R',F%X%R
             !write(6,*) 'Z',F%X%Z
#endif

             call EZspline_setup(bfield_2d%A, F%PSIp, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             !write(output_unit_write,'("bfield_2d%A: ",E17.10)') bfield_2d%A%fspl(1,:,:)

             if (.not.ALLOCATED(fields_domain%FLAG2D)) &
                  ALLOCATE(fields_domain%FLAG2D(bfield_2d%NR,bfield_2d%NZ))

             fields_domain%FLAG2D = F%FLAG2D

             if (.not.ALLOCATED(fields_domain%LCFS2D)) &
                  ALLOCATE(fields_domain%LCFS2D(bfield_2d%NR,bfield_2d%NZ))

             fields_domain%LCFS2D = F%LCFS2D

             fields_domain%DR = ABS(F%X%R(2) - F%X%R(1))
             fields_domain%DZ = ABS(F%X%Z(2) - F%X%Z(1))
          endif
       end if

       if(F%B1field) then

          if (params%field_model(10:13).eq.'MARS') then

             write(output_unit_write,*) '2D n=1 MARS magnetic fields'

             b1Refield_2d%NR = F%dims(1)
             b1Refield_2d%NZ = F%dims(3)

             ! Initializing BR1Re interpolant
             call EZspline_init(b1Refield_2d%R,b1Refield_2d%NR, &
                  b1Refield_2d%NZ,b1Refield_2d%BCSR,b1Refield_2d%BCSZ,ezerr)

             call EZspline_error(ezerr)

             b1Refield_2d%R%x1 = F%X%R
             b1Refield_2d%R%x2 = F%X%Z

             !write(output_unit_write,'("R",E17.10)') F%X%R
             !write(output_unit_write,'("Z",E17.10)') F%X%Z

             call EZspline_setup(b1Refield_2d%R, F%B1Re_2D%R, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             !write(6,*) 'size',size(F%B1Re_2D%R(:,200))
             !write(6,*) 'B1Re_2D%R',F%B1Re_2D%R(:,200)*params%cpp%Bo

             ! Initializing BPHI1Re interpolant
             call EZspline_init(b1Refield_2d%PHI,b1Refield_2d%NR, &
                  b1Refield_2d%NZ,b1Refield_2d%BCSR,b1Refield_2d%BCSZ,ezerr)

             call EZspline_error(ezerr)

             b1Refield_2d%PHI%x1 = F%X%R
             b1Refield_2d%PHI%x2 = F%X%Z

             !write(output_unit_write,'("R",E17.10)') F%X%R
             !write(output_unit_write,'("Z",E17.10)') F%X%Z

             call EZspline_setup(b1Refield_2d%PHI, F%B1Re_2D%PHI, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ! Initializing BZ1_Re interpolant
             call EZspline_init(b1Refield_2d%Z,b1Refield_2d%NR, &
                  b1Refield_2d%NZ,b1Refield_2d%BCSR,b1Refield_2d%BCSZ,ezerr)

             call EZspline_error(ezerr)

             b1Refield_2d%Z%x1 = F%X%R
             b1Refield_2d%Z%x2 = F%X%Z

             !write(output_unit_write,'("R",E17.10)') F%X%R
             !write(output_unit_write,'("Z",E17.10)') F%X%Z

             call EZspline_setup(b1Refield_2d%Z, F%B1Re_2D%Z, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             b1Imfield_2d%NR = F%dims(1)
             b1Imfield_2d%NZ = F%dims(3)

             ! Initializing BR1RIm interpolant
             call EZspline_init(b1Imfield_2d%R,b1Imfield_2d%NR, &
                  b1Imfield_2d%NZ,b1Imfield_2d%BCSR,b1Imfield_2d%BCSZ,ezerr)

             call EZspline_error(ezerr)

             b1Imfield_2d%R%x1 = F%X%R
             b1Imfield_2d%R%x2 = F%X%Z

             !write(output_unit_write,'("R",E17.10)') F%X%R
             !write(output_unit_write,'("Z",E17.10)') F%X%Z

             call EZspline_setup(b1Imfield_2d%R, F%B1Im_2D%R, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ! Initializing BPHI1Im interpolant
             call EZspline_init(b1Imfield_2d%PHI,b1Imfield_2d%NR, &
                  b1Imfield_2d%NZ,b1Imfield_2d%BCSR,b1Imfield_2d%BCSZ,ezerr)

             call EZspline_error(ezerr)

             b1Imfield_2d%PHI%x1 = F%X%R
             b1Imfield_2d%PHI%x2 = F%X%Z

             !write(output_unit_write,'("R",E17.10)') F%X%R
             !write(output_unit_write,'("Z",E17.10)') F%X%Z

             call EZspline_setup(b1Imfield_2d%PHI, F%B1Im_2D%PHI, ezerr, .TRUE.)
             call EZspline_error(ezerr)


             ! Initializing BZ1_Im interpolant
             call EZspline_init(b1Imfield_2d%Z,b1Imfield_2d%NR, &
                  b1Imfield_2d%NZ,b1Imfield_2d%BCSR,b1Imfield_2d%BCSZ,ezerr)

             call EZspline_error(ezerr)

             b1Imfield_2d%Z%x1 = F%X%R
             b1Imfield_2d%Z%x2 = F%X%Z

             !write(output_unit_write,'("R",E17.10)') F%X%R
             !write(output_unit_write,'("Z",E17.10)') F%X%Z

             call EZspline_setup(b1Imfield_2d%Z, F%B1Im_2D%Z, ezerr, .TRUE.)
             call EZspline_error(ezerr)

          else if (params%field_model(10:14).eq.'AORSA') then

             write(output_unit_write,*) '2D AORSA magnetic fields'
             flush(output_unit_write)

             b1Refield_2dx%NR = F%dims(1)
             b1Refield_2dx%NZ = F%dims(3)

             ! Initializing BX1Re interpolant
             call EZspline_init(b1Refield_2dx%X,b1Refield_2dx%NR, &
                  b1Refield_2dx%NZ,b1Refield_2dx%BCSR,b1Refield_2dx%BCSZ,ezerr)

             call EZspline_error(ezerr)

             b1Refield_2dx%X%x1 = F%X%R
             b1Refield_2dx%X%x2 = F%X%Z

             !write(output_unit_write,'("R",E17.10)') F%X%R
             !write(output_unit_write,'("Z",E17.10)') F%X%Z

             call EZspline_setup(b1Refield_2dx%X, F%B1Re_2DX%X, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             !write(6,*) 'size',size(F%B1Re_2D%R(:,200))
             !write(6,*) 'B1Re_2DX%X',F%B1Re_2DX%X(:,200)*params%cpp%Bo


             ! Initializing BY1Re interpolant
             call EZspline_init(b1Refield_2dx%Y,b1Refield_2dx%NR, &
                  b1Refield_2dx%NZ,b1Refield_2dx%BCSR,b1Refield_2dx%BCSZ,ezerr)

             call EZspline_error(ezerr)

             b1Refield_2dx%Y%x1 = F%X%R
             b1Refield_2dx%Y%x2 = F%X%Z

             !write(output_unit_write,'("R",E17.10)') F%X%R
             !write(output_unit_write,'("Z",E17.10)') F%X%Z

             call EZspline_setup(b1Refield_2dx%Y, F%B1Re_2DX%Y, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             !write(6,*) 'size',size(F%B1Re_2D%R(:,200))
             !write(6,*) 'B1Re_2D%R',F%B1Re_2D%R(:,200)*params%cpp%Bo

             ! Initializing BZ1Re interpolant
             call EZspline_init(b1Refield_2dx%Z,b1Refield_2dx%NR, &
                  b1Refield_2dx%NZ,b1Refield_2dx%BCSR,b1Refield_2dx%BCSZ,ezerr)

             call EZspline_error(ezerr)

             b1Refield_2dx%Z%x1 = F%X%R
             b1Refield_2dx%Z%x2 = F%X%Z

             !write(output_unit_write,'("R",E17.10)') F%X%R
             !write(output_unit_write,'("Z",E17.10)') F%X%Z

             call EZspline_setup(b1Refield_2dx%Z, F%B1Re_2DX%Z, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             !write(6,*) 'size',size(F%B1Re_2D%R(:,200))
             !write(6,*) 'B1Re_2D%R',F%B1Re_2D%R(:,200)*params%cpp%Bo

             b1Imfield_2dx%NR = F%dims(1)
             b1Imfield_2dx%NZ = F%dims(3)

             ! Initializing BX1Im interpolant
             call EZspline_init(b1Imfield_2dx%X,b1Imfield_2dx%NR, &
                  b1Imfield_2dx%NZ,b1Imfield_2dx%BCSR,b1Imfield_2dx%BCSZ,ezerr)

             call EZspline_error(ezerr)

             b1Imfield_2dx%X%x1 = F%X%R
             b1Imfield_2dx%X%x2 = F%X%Z

             !write(output_unit_write,'("R",E17.10)') F%X%R
             !write(output_unit_write,'("Z",E17.10)') F%X%Z

             call EZspline_setup(b1Imfield_2dx%X, F%B1Im_2DX%X, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             !write(6,*) 'size',size(F%B1Im_2D%R(:,200))
             !write(6,*) 'B1Im_2D%R',F%B1Im_2D%R(:,200)*params%cpp%Bo


             ! Initializing BY1Im interpolant
             call EZspline_init(b1Imfield_2dx%Y,b1Imfield_2dx%NR, &
                  b1Imfield_2dx%NZ,b1Imfield_2dx%BCSR,b1Imfield_2dx%BCSZ,ezerr)

             call EZspline_error(ezerr)

             b1Imfield_2dx%Y%x1 = F%X%R
             b1Imfield_2dx%Y%x2 = F%X%Z

             !write(output_unit_write,'("R",E17.10)') F%X%R
             !write(output_unit_write,'("Z",E17.10)') F%X%Z

             call EZspline_setup(b1Imfield_2dx%Y, F%B1Im_2DX%Y, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             !write(6,*) 'size',size(F%B1Im_2D%R(:,200))
             !write(6,*) 'B1Im_2D%R',F%B1Im_2D%R(:,200)*params%cpp%Bo

             ! Initializing BZ1Im interpolant
             call EZspline_init(b1Imfield_2dx%Z,b1Imfield_2dx%NR, &
                  b1Imfield_2dx%NZ,b1Imfield_2dx%BCSR,b1Imfield_2dx%BCSZ,ezerr)

             call EZspline_error(ezerr)

             b1Imfield_2dx%Z%x1 = F%X%R
             b1Imfield_2dx%Z%x2 = F%X%Z

             !write(output_unit_write,'("R",E17.10)') F%X%R
             !write(output_unit_write,'("Z",E17.10)') F%X%Z

             call EZspline_setup(b1Imfield_2dx%Z, F%B1Im_2DX%Z, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             !write(6,*) 'size',size(F%B1Im_2D%R(:,200))
             !write(6,*) 'B1Im_2D%R',F%B1Im_2D%R(:,200)*params%cpp%Bo

          endif

       end if

       if (F%Bflux3D) then

          if (params%mpi_params%rank .EQ. 0) then
             write(output_unit_write,'("3D psi potential")')
          end if

          if(F%Dim2x1t) then

             write(output_unit_write,*) '2X1T poloidal flux function'

             bfield_2X1T%NR = F%dims(1)
             bfield_2X1T%NT = F%dims(2)
             bfield_2X1T%NZ = F%dims(3)

             call EZspline_init(bfield_2X1T%A, bfield_2X1T%NR, bfield_2X1T%NT, &
                  bfield_2X1T%NZ,&
                  bfield_2X1T%BCSR, bfield_2X1T%BCST, bfield_2X1T%BCSZ, ezerr)
             call EZspline_error(ezerr)

             bfield_2X1T%A%x1 = F%X%R
             bfield_2X1T%A%x2 = F%X%PHI
             bfield_2X1T%A%x3 = F%X%Z

             !write(output_unit_write,*) F%X%PHI

             call EZspline_setup(bfield_2X1T%A, F%PSIp3D, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             if (.not.ALLOCATED(fields_domain%FLAG3D)) &
                  ALLOCATE(fields_domain%FLAG3D(bfield_2X1T%NR,bfield_2X1T%NT, &
                  bfield_2X1T%NZ))
             fields_domain%FLAG3D = F%FLAG3D

             fields_domain%DR = ABS(F%X%R(2) - F%X%R(1))
             fields_domain%DT = ABS(F%X%PHI(2) - F%X%PHI(1))
             fields_domain%DZ = ABS(F%X%Z(2) - F%X%Z(1))

             fields_domain%To = F%X%PHI(1)

          else

             write(output_unit_write,*) '3D poloidal flux function'

             bfield_3d%NR = F%dims(1)
             bfield_3d%NPHI = F%dims(2)
             bfield_3d%NZ = F%dims(3)

             ! Initializing R component of interpolant
             call EZspline_init(bfield_3d%A, bfield_3d%NR, bfield_3d%NPHI, &
                  bfield_3d%NZ,&
                  bfield_3d%BCSR, bfield_3d%BCSPHI, bfield_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             bfield_3d%A%x1 = F%X%R
             ! bfield_3d%R%x2 = F%X%PHI
             bfield_3d%A%x3 = F%X%Z

             call EZspline_setup(bfield_3d%A, F%PSIp3D, ezerr, .TRUE.)
             call EZspline_error(ezerr)
          end if

       end if


       if (F%Bfield) then
          if (F%axisymmetric_fields) then

             write(output_unit_write,*) '2D magnetic field'

             bfield_2d%NR = F%dims(1)
             bfield_2d%NZ = F%dims(3)

             ! Initializing R component
             call EZspline_init(bfield_2d%R,bfield_2d%NR,bfield_2d%NZ, &
                  bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             bfield_2d%R%x1 = F%X%R
             bfield_2d%R%x2 = F%X%Z

             call EZspline_setup(bfield_2d%R, F%B_2D%R, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ! Initializing PHI component
             call EZspline_init(bfield_2d%PHI,bfield_2d%NR,bfield_2d%NZ, &
                  bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             bfield_2d%PHI%x1 = F%X%R
             bfield_2d%PHI%x2 = F%X%Z

             call EZspline_setup(bfield_2d%PHI, F%B_2D%PHI, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ! Initializing Z component
             call EZspline_init(bfield_2d%Z,bfield_2d%NR,bfield_2d%NZ, &
                  bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             bfield_2d%Z%x1 = F%X%R
             bfield_2d%Z%x2 = F%X%Z

             call EZspline_setup(bfield_2d%Z, F%B_2D%Z, ezerr, .TRUE.)
             call EZspline_error(ezerr)

!             do ii=1_idef,bfield_2d%PHI%n1
!                do jj=1_idef,bfield_2d%PHI%n2
!                   write(output_unit_write,'("BPHI_spline1 at R ",E17.10,", Z ",E17.10,": ",E17.10)') &
!                        bfield_2d%PHI%x1(ii)*params%cpp%length, &
!                        bfield_2d%PHI%x2(jj)*params%cpp%length, &
!                        bfield_2d%PHI%fspl(1,ii,jj)*params%cpp%Bo
!                end do
!             end do


             if (params%orbit_model.eq.'GCpre') then

                write(output_unit_write,*) '2D precalculated auxiliary fields'

                gradB_2d%NR = F%dims(1)
                gradB_2d%NZ = F%dims(3)

                ! Initializing GRADBR component
                call EZspline_init(gradB_2d%R,gradB_2d%NR,gradB_2d%NZ, &
                     gradB_2d%BCSR,gradB_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                gradB_2d%R%x1 = F%X%R
                gradB_2d%R%x2 = F%X%Z

                call EZspline_setup(gradB_2d%R, F%gradB_2D%R, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                ! Initializing GRADBPHI component
                call EZspline_init(gradB_2d%PHI,gradB_2d%NR,gradB_2d%NZ, &
                     gradB_2d%BCSR,gradB_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                gradB_2d%PHI%x1 = F%X%R
                gradB_2d%PHI%x2 = F%X%Z

                call EZspline_setup(gradB_2d%PHI, F%gradB_2D%PHI, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                ! Initializing GRADBZ component
                call EZspline_init(gradB_2d%Z,gradB_2d%NR,gradB_2d%NZ, &
                     gradB_2d%BCSR,gradB_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                gradB_2d%Z%x1 = F%X%R
                gradB_2d%Z%x2 = F%X%Z

                call EZspline_setup(gradB_2d%Z, F%gradB_2D%Z, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                curlb_2d%NR = F%dims(1)
                curlb_2d%NZ = F%dims(3)

                ! Initializing CURLBR component
                call EZspline_init(curlb_2d%R,curlb_2d%NR,curlb_2d%NZ, &
                     curlb_2d%BCSR,curlb_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                curlb_2d%R%x1 = F%X%R
                curlb_2d%R%x2 = F%X%Z

                call EZspline_setup(curlb_2d%R, F%curlb_2D%R, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                ! Initializing CURLBPHI component
                call EZspline_init(curlb_2d%PHI,curlb_2d%NR,curlb_2d%NZ, &
                     curlb_2d%BCSR,curlb_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                curlb_2d%PHI%x1 = F%X%R
                curlb_2d%PHI%x2 = F%X%Z

                call EZspline_setup(curlb_2d%PHI, F%curlb_2D%PHI, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                ! Initializing CURLBZ component
                call EZspline_init(curlb_2d%Z,curlb_2d%NR,curlb_2d%NZ, &
                     curlb_2d%BCSR,curlb_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                curlb_2d%Z%x1 = F%X%R
                curlb_2d%Z%x2 = F%X%Z

                call EZspline_setup(curlb_2d%Z, F%curlb_2D%Z, ezerr, .TRUE.)
                call EZspline_error(ezerr)

             end if

             if (.not.ALLOCATED(fields_domain%FLAG2D)) &
                  ALLOCATE(fields_domain%FLAG2D(bfield_2d%NR,bfield_2d%NZ))

             fields_domain%FLAG2D = F%FLAG2D

             if (.not.ALLOCATED(fields_domain%LCFS2D)) &
                  ALLOCATE(fields_domain%LCFS2D(bfield_2d%NR,bfield_2d%NZ))

             fields_domain%LCFS2D = F%LCFS2D

             fields_domain%DR = ABS(F%X%R(2) - F%X%R(1))
             fields_domain%DZ = ABS(F%X%Z(2) - F%X%Z(1))
          else

             write(output_unit_write,*) '3D magnetic field'

             bfield_3d%NR = F%dims(1)
             bfield_3d%NPHI = F%dims(2)
             bfield_3d%NZ = F%dims(3)

             ! Initializing R component of interpolant
             call EZspline_init(bfield_3d%R, bfield_3d%NR, bfield_3d%NPHI, &
                  bfield_3d%NZ,&
                  bfield_3d%BCSR, bfield_3d%BCSPHI, bfield_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             bfield_3d%R%x1 = F%X%R
             bfield_3d%R%x2 = F%X%PHI
             bfield_3d%R%x3 = F%X%Z

             call EZspline_setup(bfield_3d%R, F%B_3D%R, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             !write(6,'("bfield_3d%R: ",E17.10)') bfield_3d%R%fspl(1,:,:,:)
             !do ii=1_idef,bfield_3d%R%n1
             !   do jj=1_idef,bfield_3d%R%n3
             !      write(6,'("BR_spline1 at R ",E17.10,", Z ",E17.10,": ",E17.10)') &
             !          bfield_3d%R%x1(ii)*params%cpp%length, &
             !          bfield_3d%R%x3(jj)*params%cpp%length, &
             !          bfield_3d%R%fspl(1,ii,1,jj)*params%cpp%Bo
             !  end do
             !end do

             ! Initializing PHI component of interpolant
             call EZspline_init(bfield_3d%PHI, bfield_3d%NR, bfield_3d%NPHI, &
                  bfield_3d%NZ,&
                  bfield_3d%BCSR, bfield_3d%BCSPHI, bfield_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             bfield_3d%PHI%x1 = F%X%R
             bfield_3d%PHI%x2 = F%X%PHI
             bfield_3d%PHI%x3 = F%X%Z

             call EZspline_setup(bfield_3d%PHI, F%B_3D%PHI, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             !write(output_unit_write,*) bfield_3d%PHI%x2


             ! Initializing Z component of interpolant
             call EZspline_init(bfield_3d%Z, bfield_3d%NR, bfield_3d%NPHI, &
                  bfield_3d%NZ,&
                  bfield_3d%BCSR, bfield_3d%BCSPHI, bfield_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             bfield_3d%Z%x1 = F%X%R
             bfield_3d%Z%x2 = F%X%PHI
             bfield_3d%Z%x3 = F%X%Z

             call EZspline_setup(bfield_3d%Z, F%B_3D%Z, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             if (params%orbit_model.eq.'GCpre') then

                write(output_unit_write,*) '3D precalculated auxiliary fields'

                gradB_3d%NR = F%dims(1)
                gradB_3d%NPHI = F%dims(2)
                gradB_3d%NZ = F%dims(3)

                ! Initializing GRADBR component
                call EZspline_init(gradB_3d%R,gradB_3d%NR,gradB_3d%NPHI,&
                     gradB_3d%NZ,gradB_3d%BCSR,gradB_3d%BCSPHI, &
                     gradB_3d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                gradB_3d%R%x1 = F%X%R
                !gradB_3d%R%x2 = F%X%PHI
                gradB_3d%R%x3 = F%X%Z

                call EZspline_setup(gradB_3d%R, F%gradB_3D%R, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                ! Initializing GRADBPHI component
                call EZspline_init(gradB_3d%PHI,gradB_3d%NR,gradB_3d%NPHI,&
                     gradB_3d%NZ,gradB_3d%BCSR,gradB_3d%BCSPHI, &
                     gradB_3d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                gradB_3d%PHI%x1 = F%X%R
                !gradB_3d%PHI%x2 = F%X%PHI
                gradB_3d%PHI%x3 = F%X%Z

                call EZspline_setup(gradB_3d%PHI, F%gradB_3D%PHI, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                ! Initializing GRADBZ component
                call EZspline_init(gradB_3d%Z,gradB_3d%NR,gradB_3d%NPHI,&
                     gradB_3d%NZ,gradB_3d%BCSR,gradB_3d%BCSPHI, &
                     gradB_3d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                gradB_3d%Z%x1 = F%X%R
                !gradB_3d%Z%x2 = F%X%PHI
                gradB_3d%Z%x3 = F%X%Z

                call EZspline_setup(gradB_3d%Z, F%gradB_3D%Z, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                curlb_3d%NR = F%dims(1)
                curlb_3d%NPHI = F%dims(2)
                curlb_3d%NZ = F%dims(3)

                ! Initializing CURLBR component
                call EZspline_init(curlb_3d%R,curlb_3d%NR,curlb_3d%NPHI,&
                     curlb_3d%NZ,curlb_3d%BCSR,curlb_3d%BCSPHI, &
                     curlb_3d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                curlb_3d%R%x1 = F%X%R
                !curlb_3d%R%x2 = F%X%PHI
                curlb_3d%R%x3 = F%X%Z

                call EZspline_setup(curlb_3d%R, F%curlb_3D%R, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                ! Initializing CURLBPHI component
                call EZspline_init(curlb_3d%PHI,curlb_3d%NR,curlb_3d%NPHI,&
                     curlb_3d%NZ,curlb_3d%BCSR,curlb_3d%BCSPHI, &
                     curlb_3d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                curlb_3d%PHI%x1 = F%X%R
                !curlb_3d%PHI%x2 = F%X%PHI
                curlb_3d%PHI%x3 = F%X%Z

                call EZspline_setup(curlb_3d%PHI, F%curlb_3D%PHI, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                ! Initializing CURLBZ component
                call EZspline_init(curlb_3d%Z,curlb_3d%NR,curlb_3d%NPHI,&
                     curlb_3d%NZ,curlb_3d%BCSR,curlb_3d%BCSPHI, &
                     curlb_3d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                curlb_3d%Z%x1 = F%X%R
                !curlb_3d%Z%x2 = F%X%PHI
                curlb_3d%Z%x3 = F%X%Z

                call EZspline_setup(curlb_3d%Z, F%curlb_3D%Z, ezerr, .TRUE.)
                call EZspline_error(ezerr)

             end if

             ALLOCATE(fields_domain%FLAG3D(bfield_3d%NR,bfield_3d%NPHI, &
                  bfield_3d%NZ))
             fields_domain%FLAG3D = F%FLAG3D

             fields_domain%DR = ABS(F%X%R(2) - F%X%R(1))
             fields_domain%DPHI = 2.0_rp*C_PI/bfield_3d%NPHI
             fields_domain%DZ = ABS(F%X%Z(2) - F%X%Z(1))
          end if
       end if

       if (F%dBfield) then
          if (F%axisymmetric_fields) then

             write(output_unit_write,*) '2D precalculated gradient fields'

             ! dBdR
             dbdR_2d%NR = F%dims(1)
             dbdR_2d%NZ = F%dims(3)

             ! Initializing R component
             call EZspline_init(dbdR_2d%R,dbdR_2d%NR,dbdR_2d%NZ, &
                  dbdR_2d%BCSR,dbdR_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             dbdR_2d%R%x1 = F%X%R
             dbdR_2d%R%x2 = F%X%Z

             call EZspline_setup(dbdR_2d%R, F%dBdR_2D%R, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ! Initializing PHI component
             call EZspline_init(dbdR_2d%PHI,dbdR_2d%NR,dbdR_2d%NZ, &
                  dbdR_2d%BCSR,dbdR_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             dbdR_2d%PHI%x1 = F%X%R
             dbdR_2d%PHI%x2 = F%X%Z

             call EZspline_setup(dbdR_2d%PHI, F%dBdR_2D%PHI, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ! Initializing Z component
             call EZspline_init(dbdR_2d%Z,dbdR_2d%NR,dbdR_2d%NZ, &
                  dbdR_2d%BCSR,dbdR_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             dbdR_2d%Z%x1 = F%X%R
             dbdR_2d%Z%x2 = F%X%Z

             call EZspline_setup(dbdR_2d%Z, F%dBdR_2D%Z, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             !dBdPHI
             dbdPHI_2d%NR = F%dims(1)
             dbdPHI_2d%NZ = F%dims(3)
             ! Initializing R component
             call EZspline_init(dbdPHI_2d%R,dbdPHI_2d%NR,dbdPHI_2d%NZ, &
                  dbdPHI_2d%BCSR,dbdPHI_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             dbdPHI_2d%R%x1 = F%X%R
             dbdPHI_2d%R%x2 = F%X%Z

             call EZspline_setup(dbdPHI_2d%R, F%dBdPHI_2D%R, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ! Initializing PHI component
             call EZspline_init(dbdPHI_2d%PHI,dbdPHI_2d%NR,dbdPHI_2d%NZ, &
                  dbdPHI_2d%BCSR,dbdPHI_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             dbdPHI_2d%PHI%x1 = F%X%R
             dbdPHI_2d%PHI%x2 = F%X%Z

             call EZspline_setup(dbdPHI_2d%PHI, F%dBdPHI_2D%PHI, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ! Initializing Z component
             call EZspline_init(dbdPHI_2d%Z,dbdPHI_2d%NR,dbdPHI_2d%NZ, &
                  dbdPHI_2d%BCSR,dbdPHI_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             dbdPHI_2d%Z%x1 = F%X%R
             dbdPHI_2d%Z%x2 = F%X%Z

             call EZspline_setup(dbdPHI_2d%Z, F%dBdPHI_2D%Z, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             !dBdZ
             dbdZ_2d%NR = F%dims(1)
             dbdZ_2d%NZ = F%dims(3)
             ! Initializing R component
             call EZspline_init(dbdZ_2d%R,dbdZ_2d%NR,dbdZ_2d%NZ, &
                  dbdZ_2d%BCSR,dbdZ_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             dbdZ_2d%R%x1 = F%X%R
             dbdZ_2d%R%x2 = F%X%Z

             call EZspline_setup(dbdZ_2d%R, F%dBdZ_2D%R, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ! Initializing PHI component
             call EZspline_init(dbdZ_2d%PHI,dbdZ_2d%NR,dbdZ_2d%NZ, &
                  dbdZ_2d%BCSR,dbdZ_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             dbdZ_2d%PHI%x1 = F%X%R
             dbdZ_2d%PHI%x2 = F%X%Z

             call EZspline_setup(dbdZ_2d%PHI, F%dBdZ_2D%PHI, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ! Initializing Z component
             call EZspline_init(dbdZ_2d%Z,dbdZ_2d%NR,dbdZ_2d%NZ, &
                  dbdZ_2d%BCSR,dbdZ_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             dbdZ_2d%Z%x1 = F%X%R
             dbdZ_2d%Z%x2 = F%X%Z

             call EZspline_setup(dbdZ_2d%Z, F%dBdZ_2D%Z, ezerr, .TRUE.)
             call EZspline_error(ezerr)

          else

             write(output_unit_write,*) '3D precalculated gradient fields'

             ! dBdR
             dbdR_3d%NR = F%dims(1)
             dbdR_3d%NPHI = F%dims(2)
             dbdR_3d%NZ = F%dims(3)

             ! Initializing R component of interpolant
             call EZspline_init(dbdR_3d%R, dbdR_3d%NR, dbdR_3d%NPHI, &
                  dbdR_3d%NZ,&
                  dbdR_3d%BCSR, dbdR_3d%BCSPHI, dbdR_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             dbdR_3d%R%x1 = F%X%R
             ! dbdR_3d%R%x2 = F%X%PHI
             dbdR_3d%R%x3 = F%X%Z

             call EZspline_setup(dbdR_3d%R, F%dBdR_3D%R, ezerr)
             call EZspline_error(ezerr)

             ! Initializing PHI component of interpolant
             call EZspline_init(dbdR_3d%PHI, dbdR_3d%NR, dbdR_3d%NPHI, &
                  dbdR_3d%NZ,&
                  dbdR_3d%BCSR, dbdR_3d%BCSPHI, dbdR_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             dbdR_3d%PHI%x1 = F%X%R
             ! dbdR_3d%PHI%x2 = F%X%PHI
             dbdR_3d%PHI%x3 = F%X%Z

             call EZspline_setup(dbdR_3d%PHI, F%dBdR_3D%PHI, ezerr)
             call EZspline_error(ezerr)

             !write(output_unit_write,*) dbdR_3d%PHI%x2


             ! Initializing Z component of interpolant
             call EZspline_init(dbdR_3d%Z, dbdR_3d%NR, dbdR_3d%NPHI, &
                  dbdR_3d%NZ,&
                  dbdR_3d%BCSR, dbdR_3d%BCSPHI, dbdR_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             dbdR_3d%Z%x1 = F%X%R
             ! dbdR_3d%Z%x2 = F%X%PHI
             dbdR_3d%Z%x3 = F%X%Z

             call EZspline_setup(dbdR_3d%Z, F%dBdR_3D%Z, ezerr)
             call EZspline_error(ezerr)

             !dBdPHI
             dbdPHI_3d%NR = F%dims(1)
             dbdPHI_3d%NPHI = F%dims(2)
             dbdPHI_3d%NZ = F%dims(3)

             ! Initializing R component of interpolant
             call EZspline_init(dbdPHI_3d%R, dbdPHI_3d%NR, dbdPHI_3d%NPHI, &
                  dbdPHI_3d%NZ,&
                  dbdPHI_3d%BCSR, dbdPHI_3d%BCSPHI, dbdPHI_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             dbdPHI_3d%R%x1 = F%X%R
             ! dbdPHI_3d%R%x2 = F%X%PHI
             dbdPHI_3d%R%x3 = F%X%Z

             call EZspline_setup(dbdPHI_3d%R, F%dBdPHI_3D%R, ezerr)
             call EZspline_error(ezerr)

             ! Initializing PHI component of interpolant
             call EZspline_init(dbdPHI_3d%PHI, dbdPHI_3d%NR, dbdPHI_3d%NPHI, &
                  dbdPHI_3d%NZ,&
                  dbdPHI_3d%BCSR, dbdPHI_3d%BCSPHI, dbdPHI_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             dbdPHI_3d%PHI%x1 = F%X%R
             ! dbdPHI_3d%PHI%x2 = F%X%PHI
             dbdPHI_3d%PHI%x3 = F%X%Z

             call EZspline_setup(dbdPHI_3d%PHI, F%dBdPHI_3D%PHI, ezerr)
             call EZspline_error(ezerr)

             !write(output_unit_write,*) dbdPHI_3d%PHI%x2


             ! Initializing Z component of interpolant
             call EZspline_init(dbdPHI_3d%Z, dbdPHI_3d%NR, dbdPHI_3d%NPHI, &
                  dbdPHI_3d%NZ,&
                  dbdPHI_3d%BCSR, dbdPHI_3d%BCSPHI, dbdPHI_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             dbdPHI_3d%Z%x1 = F%X%R
             ! dbdPHI_3d%Z%x2 = F%X%PHI
             dbdPHI_3d%Z%x3 = F%X%Z

             call EZspline_setup(dbdPHI_3d%Z, F%dBdPHI_3D%Z, ezerr)
             call EZspline_error(ezerr)

             !dBdZ
             dbdZ_3d%NR = F%dims(1)
             dbdZ_3d%NPHI = F%dims(2)
             dbdZ_3d%NZ = F%dims(3)

             ! Initializing R component of interpolant
             call EZspline_init(dbdZ_3d%R, dbdZ_3d%NR, dbdZ_3d%NPHI, &
                  dbdZ_3d%NZ,&
                  dbdZ_3d%BCSR, dbdZ_3d%BCSPHI, dbdZ_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             dbdZ_3d%R%x1 = F%X%R
             ! dbdZ_3d%R%x2 = F%X%PHI
             dbdZ_3d%R%x3 = F%X%Z

             call EZspline_setup(dbdZ_3d%R, F%dBdZ_3D%R, ezerr)
             call EZspline_error(ezerr)

             ! Initializing PHI component of interpolant
             call EZspline_init(dbdZ_3d%PHI, dbdZ_3d%NR, dbdZ_3d%NPHI, &
                  dbdZ_3d%NZ,&
                  dbdZ_3d%BCSR, dbdZ_3d%BCSPHI, dbdZ_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             dbdZ_3d%PHI%x1 = F%X%R
             ! dbdZ_3d%PHI%x2 = F%X%PHI
             dbdZ_3d%PHI%x3 = F%X%Z

             call EZspline_setup(dbdZ_3d%PHI, F%dBdZ_3D%PHI, ezerr)
             call EZspline_error(ezerr)

             !write(output_unit_write,*) dbdZ_3d%PHI%x2


             ! Initializing Z component of interpolant
             call EZspline_init(dbdZ_3d%Z, dbdZ_3d%NR, dbdZ_3d%NPHI, &
                  dbdZ_3d%NZ,&
                  dbdZ_3d%BCSR, dbdZ_3d%BCSPHI, dbdZ_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             dbdZ_3d%Z%x1 = F%X%R
             ! dbdZ_3d%Z%x2 = F%X%PHI
             dbdZ_3d%Z%x3 = F%X%Z

             call EZspline_setup(dbdZ_3d%Z, F%dBdZ_3D%Z, ezerr)
             call EZspline_error(ezerr)

          end if
       end if

       fields_domain%Ro = F%X%R(1)
       fields_domain%Zo = F%X%Z(1)

       ! * * * * * * * * ELECTRIC FIELD * * * * * * * * !
       if (F%Efield.AND.(params%field_eval.eq.'interp')) then
          if (F%axisymmetric_fields) then

             if (params%mpi_params%rank .EQ. 0) then
                write(output_unit_write,*) '2D electric fields'
             endif

             if(F%Dim2x1t) then

                if (.not.(EZspline_allocated(efield_2d%PHI))) then

                   efield_2d%NR = F%dims(1)
                   efield_2d%NZ = F%dims(3)

                   ! Initializing R component
                   call EZspline_init(efield_2d%R,efield_2d%NR,efield_2d%NZ, &
                        efield_2d%BCSR,efield_2d%BCSZ,ezerr)
                   call EZspline_error(ezerr)

                   efield_2d%R%x1 = F%X%R
                   efield_2d%R%x2 = F%X%Z

                   ! Initializing PHI component
                   call EZspline_init(efield_2d%PHI,efield_2d%NR,efield_2d%NZ, &
                        efield_2d%BCSR,efield_2d%BCSZ,ezerr)
                   call EZspline_error(ezerr)

                   efield_2d%PHI%x1 = F%X%R
                   efield_2d%PHI%x2 = F%X%Z

                   ! Initializing Z component
                   call EZspline_init(efield_2d%Z,efield_2d%NR,efield_2d%NZ, &
                        efield_2d%BCSR,efield_2d%BCSZ,ezerr)
                   call EZspline_error(ezerr)

                   efield_2d%Z%x1 = F%X%R
                   efield_2d%Z%x2 = F%X%Z

                end if

                call EZspline_setup(efield_2d%R, F%E_3D%R(:,F%ind_2x1t,:), &
                     ezerr, .TRUE.)
                call EZspline_error(ezerr)

                call EZspline_setup(efield_2d%PHI, F%E_3D%PHI(:,F%ind_2x1t,:), &
                     ezerr, .TRUE.)
                call EZspline_error(ezerr)

                call EZspline_setup(efield_2d%Z, F%E_3D%Z(:,F%ind_2x1t,:), &
                     ezerr, .TRUE.)
                call EZspline_error(ezerr)

!                write(output_unit_write,'("efield_2d%PHI: ",E17.10)') efield_2d%PHI%fspl(1,:,:)

             else

                efield_2d%NR = F%dims(1)
                efield_2d%NZ = F%dims(3)

                ! Initializing R component
                call EZspline_init(efield_2d%R,efield_2d%NR,efield_2d%NZ, &
                     efield_2d%BCSR,efield_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                efield_2d%R%x1 = F%X%R
                efield_2d%R%x2 = F%X%Z

                call EZspline_setup(efield_2d%R, F%E_2D%R, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                ! Initializing PHI component
                call EZspline_init(efield_2d%PHI,efield_2d%NR,efield_2d%NZ, &
                     efield_2d%BCSR,efield_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                efield_2d%PHI%x1 = F%X%R
                efield_2d%PHI%x2 = F%X%Z

                call EZspline_setup(efield_2d%PHI, F%E_2D%PHI, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                ! Initializing Z component
                call EZspline_init(efield_2d%Z,efield_2d%NR,efield_2d%NZ, &
                     efield_2d%BCSR,efield_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                efield_2d%Z%x1 = F%X%R
                efield_2d%Z%x2 = F%X%Z

                call EZspline_setup(efield_2d%Z, F%E_2D%Z, ezerr, .TRUE.)
                call EZspline_error(ezerr)

             end if
          else

             write(output_unit_write,*) '3D electric fields'

             efield_3d%NR = F%dims(1)
             efield_3d%NPHI = F%dims(2)
             efield_3d%NZ = F%dims(3)

             ! Initializing R component of interpolant
             call EZspline_init(efield_3d%R, efield_3d%NR, efield_3d%NPHI, efield_3d%NZ,&
                  efield_3d%BCSR, efield_3d%BCSPHI, efield_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             efield_3d%R%x1 = F%X%R
             ! efield_3d%R%x2 = F%X%PHI
             efield_3d%R%x3 = F%X%Z

             call EZspline_setup(efield_3d%R, F%E_3D%R, ezerr)
             call EZspline_error(ezerr)

             ! Initializing PHI component of interpolant
             call EZspline_init(efield_3d%PHI, efield_3d%NR, efield_3d%NPHI, &
                  efield_3d%NZ,efield_3d%BCSR, efield_3d%BCSPHI, efield_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             efield_3d%PHI%x1 = F%X%R
             ! efield_3d%PHI%x2 = F%X%PHI
             efield_3d%PHI%x3 = F%X%Z

             call EZspline_setup(efield_3d%PHI, F%E_3D%PHI, ezerr)
             call EZspline_error(ezerr)

             ! Initializing Z component of interpolant
             call EZspline_init(efield_3d%Z, efield_3d%NR, efield_3d%NPHI, efield_3d%NZ,&
                  efield_3d%BCSR, efield_3d%BCSPHI, efield_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             efield_3d%Z%x1 = F%X%R
             ! efield_3d%Z%x2 = F%X%PHI
             efield_3d%Z%x3 = F%X%Z

             call EZspline_setup(efield_3d%Z, F%E_3D%Z, ezerr)
             call EZspline_error(ezerr)
          end if
       end if

       if(F%E1field) then

          write(output_unit_write,*) '2D AORSA electric fields'
          flush(output_unit_write)

          e1Refield_2dx%NR = F%dims(1)
          e1Refield_2dx%NZ = F%dims(3)

          ! Initializing BX1Re interpolant
          call EZspline_init(e1Refield_2dx%X,e1Refield_2dx%NR, &
               e1Refield_2dx%NZ,e1Refield_2dx%BCSR,e1Refield_2dx%BCSZ,ezerr)

          call EZspline_error(ezerr)

          e1Refield_2dx%X%x1 = F%X%R
          e1Refield_2dx%X%x2 = F%X%Z

          !write(output_unit_write,'("R",E17.10)') F%X%R
          !write(output_unit_write,'("Z",E17.10)') F%X%Z

          call EZspline_setup(e1Refield_2dx%X, F%E1Re_2DX%X, ezerr, .TRUE.)
          call EZspline_error(ezerr)

          !write(6,*) 'size',size(F%E1Re_2D%R(:,200))
          !write(6,*) 'E1Re_2D%R',F%E1Re_2D%R(:,200)*params%cpp%Bo


          ! Initializing BY1Re interpolant
          call EZspline_init(e1Refield_2dx%Y,e1Refield_2dx%NR, &
               e1Refield_2dx%NZ,e1Refield_2dx%BCSR,e1Refield_2dx%BCSZ,ezerr)

          call EZspline_error(ezerr)

          e1Refield_2dx%Y%x1 = F%X%R
          e1Refield_2dx%Y%x2 = F%X%Z

          !write(output_unit_write,'("R",E17.10)') F%X%R
          !write(output_unit_write,'("Z",E17.10)') F%X%Z

          call EZspline_setup(e1Refield_2dx%Y, F%E1Re_2DX%Y, ezerr, .TRUE.)
          call EZspline_error(ezerr)

          !write(6,*) 'size',size(F%E1Re_2D%R(:,200))
          !write(6,*) 'E1Re_2D%R',F%E1Re_2D%R(:,200)*params%cpp%Bo

          ! Initializing BZ1Re interpolant
          call EZspline_init(e1Refield_2dx%Z,e1Refield_2dx%NR, &
               e1Refield_2dx%NZ,e1Refield_2dx%BCSR,e1Refield_2dx%BCSZ,ezerr)

          call EZspline_error(ezerr)

          e1Refield_2dx%Z%x1 = F%X%R
          e1Refield_2dx%Z%x2 = F%X%Z

          !write(output_unit_write,'("R",E17.10)') F%X%R
          !write(output_unit_write,'("Z",E17.10)') F%X%Z

          call EZspline_setup(e1Refield_2dx%Z, F%E1Re_2DX%Z, ezerr, .TRUE.)
          call EZspline_error(ezerr)

          !write(6,*) 'size',size(F%E1Re_2D%R(:,200))
          !write(6,*) 'E1Re_2D%R',F%E1Re_2D%R(:,200)*params%cpp%Bo

          e1Imfield_2dx%NR = F%dims(1)
          e1Imfield_2dx%NZ = F%dims(3)

          ! Initializing BX1Im interpolant
          call EZspline_init(e1Imfield_2dx%X,e1Imfield_2dx%NR, &
               e1Imfield_2dx%NZ,e1Imfield_2dx%BCSR,e1Imfield_2dx%BCSZ,ezerr)

          call EZspline_error(ezerr)

          e1Imfield_2dx%X%x1 = F%X%R
          e1Imfield_2dx%X%x2 = F%X%Z

          !write(output_unit_write,'("R",E17.10)') F%X%R
          !write(output_unit_write,'("Z",E17.10)') F%X%Z

          call EZspline_setup(e1Imfield_2dx%X, F%E1Im_2DX%X, ezerr, .TRUE.)
          call EZspline_error(ezerr)

          !write(6,*) 'size',size(F%E1Im_2D%R(:,200))
          !write(6,*) 'E1Im_2D%R',F%E1Im_2D%R(:,200)*params%cpp%Bo


          ! Initializing BY1Im interpolant
          call EZspline_init(e1Imfield_2dx%Y,e1Imfield_2dx%NR, &
               e1Imfield_2dx%NZ,e1Imfield_2dx%BCSR,e1Imfield_2dx%BCSZ,ezerr)

          call EZspline_error(ezerr)

          e1Imfield_2dx%Y%x1 = F%X%R
          e1Imfield_2dx%Y%x2 = F%X%Z

          !write(output_unit_write,'("R",E17.10)') F%X%R
          !write(output_unit_write,'("Z",E17.10)') F%X%Z

          call EZspline_setup(e1Imfield_2dx%Y, F%E1Im_2DX%Y, ezerr, .TRUE.)
          call EZspline_error(ezerr)

          !write(6,*) 'size',size(F%E1Im_2D%R(:,200))
          !write(6,*) 'E1Im_2D%R',F%E1Im_2D%R(:,200)*params%cpp%Bo

          ! Initializing BZ1Im interpolant
          call EZspline_init(e1Imfield_2dx%Z,e1Imfield_2dx%NR, &
               e1Imfield_2dx%NZ,e1Imfield_2dx%BCSR,e1Imfield_2dx%BCSZ,ezerr)

          call EZspline_error(ezerr)

          e1Imfield_2dx%Z%x1 = F%X%R
          e1Imfield_2dx%Z%x2 = F%X%Z

          !write(output_unit_write,'("R",E17.10)') F%X%R
          !write(output_unit_write,'("Z",E17.10)') F%X%Z

          call EZspline_setup(e1Imfield_2dx%Z, F%E1Im_2DX%Z, ezerr, .TRUE.)
          call EZspline_error(ezerr)

          !write(6,*) 'size',size(F%E1Im_2D%R(:,200))
          !write(6,*) 'E1Im_2D%R',F%E1Im_2D%R(:,200)*params%cpp%Bo


       end if

       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'("* * * * * * INTERPOLANT INITIALIZED * * * * * *",/)')
       end if
    else if (params%field_model(1:10) .EQ. 'ANALYTICAL') then
       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'("* * * * USING ANALYTICAL MAGNETIC FIELD * * * *",/)')
       end if
    else if (params%field_model .EQ. 'UNIFORM') then
       if (params%mpi_params%rank .EQ. 0) then
          write(output_unit_write,'("* * * * USING UNIFORM MAGNETIC FIELD * * * *",/)')
       end if
    end if
  end subroutine initialize_fields_interpolant

  subroutine initialize_Hollmann_interpolant(params,Nrho,NE,rho_axis,g_axis,fRE_E,fRE_pitch)
    !! @note Subroutine that initializes fields interpolants. @endnote
    !! This subroutine initializes either 2-D or 3-D PSPLINE interpolants
    !! using the data of fields in the KORC-dervied-type variable F.
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    !! Core KORC simulation parameters.
    INTEGER, INTENT(IN)       :: Nrho,NE
    REAL(rp),DIMENSION(Nrho), INTENT(IN)  :: rho_axis
    REAL(rp),DIMENSION(NE), INTENT(IN)  :: g_axis
    REAL(rp),DIMENSION(Nrho,NE), INTENT(IN)  :: fRE_E,fRE_pitch


    hollmann_2d%NRHO = Nrho
    hollmann_2d%NE = NE

    ! initialize fRE_E
    call EZspline_init(hollmann_2d%fRE_E,hollmann_2d%NRHO, &
         hollmann_2d%NE,hollmann_2d%BCRHO,hollmann_2d%BCE,ezerr)
    call EZspline_error(ezerr)

    hollmann_2d%fRE_E%x1 = rho_axis
    hollmann_2d%fRE_E%x2 = g_axis

    call EZspline_setup(hollmann_2d%fRE_E, fRE_E, ezerr, .TRUE.)
    call EZspline_error(ezerr)


    ! initialize fRE_pitch
    call EZspline_init(hollmann_2d%fRE_pitch,hollmann_2d%NRHO, &
         hollmann_2d%NE,hollmann_2d%BCRHO,hollmann_2d%BCE,ezerr)
    call EZspline_error(ezerr)

    hollmann_2d%fRE_pitch%x1 = rho_axis
    hollmann_2d%fRE_pitch%x2 = g_axis


    call EZspline_setup(hollmann_2d%fRE_pitch, fRE_pitch, &
         ezerr, .TRUE.)
    call EZspline_error(ezerr)


  end subroutine initialize_Hollmann_interpolant

  subroutine initialize_SC1D_field_interpolant(params,F)
    !! @note Subroutine that initializes fields interpolants. @endnote
    !! This subroutine initializes either 2-D or 3-D PSPLINE interpolants
    !! using the data of fields in the KORC-dervied-type variable F.
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)       :: F
    !! An instance of KORC's derived type FIELDS containing all the information
    !! about the fields used in the simulation.
    !! See [[korc_types]] and [[korc_fields]].
    integer :: ii,jj


!    if (EZspline_allocated(efield_SC1d%PHI)) &
!         call Ezspline_free(efield_SC1d%PHI, ezerr)

    if (.not.(EZspline_allocated(efield_SC1d%PHI))) then
       efield_SC1d%Nrm = F%dim_1D

       call EZspline_init(efield_SC1d%PHI,efield_SC1d%Nrm, &
            efield_SC1d%BCSrm,ezerr)
       call EZspline_error(ezerr)

       efield_SC1d%PHI%x1 = F%r_1D/params%cpp%length
    end if

    call EZspline_setup(efield_SC1d%PHI, F%E_SC_1D%PHI, ezerr, .TRUE.)
    call EZspline_error(ezerr)

    if (.not.ALLOCATED(fields_domain%FLAG1D)) &
         ALLOCATE(fields_domain%FLAG1D(efield_SC1d%Nrm))

    fields_domain%Drm = ABS(F%r_1D(2) - F%r_1D(1))

  end subroutine initialize_SC1D_field_interpolant

  subroutine initialize_SC1D_field_interpolant_FS(params,F)
    !! @note Subroutine that initializes fields interpolants. @endnote
    !! This subroutine initializes either 2-D or 3-D PSPLINE interpolants
    !! using the data of fields in the KORC-dervied-type variable F.
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    !! Core KORC simulation parameters.
    TYPE(FIELDS), INTENT(IN)       :: F
    !! An instance of KORC's derived type FIELDS containing all the information
    !! about the fields used in the simulation.
    !! See [[korc_types]] and [[korc_fields]].
    integer :: ii,jj


!    if (EZspline_allocated(efield_SC1d%PHI)) &
!         call Ezspline_free(efield_SC1d%PHI, ezerr)

    if (.not.(EZspline_allocated(efield_SC1d%PHI))) then
       efield_SC1d%NPSIP = F%dim_1D

       call EZspline_init(efield_SC1d%PHI,efield_SC1d%NPSIP, &
            efield_SC1d%BCSPSIP,ezerr)
       call EZspline_error(ezerr)

       efield_SC1d%PHI%x1 = F%PSIP_1D/(params%cpp%Bo*params%cpp%length**2)
    end if

    call EZspline_setup(efield_SC1d%PHI, F%E_SC_1D%PHI, ezerr, .TRUE.)
    call EZspline_error(ezerr)

    if (.not.ALLOCATED(fields_domain%FLAG1D)) &
         ALLOCATE(fields_domain%FLAG1D(efield_SC1d%Nrm))

    fields_domain%DPSIP = ABS(F%PSIP_1D(2) - F%PSIP_1D(1))

  end subroutine initialize_SC1D_field_interpolant_FS

  subroutine check_if_in_fields_domain(F,Y,flag)
    !! @note Subrotuine that checks if particles in the simulation are within
    !! the spatial domain where interpolants and fields are known. @endnote
    !! External fields and interpolants can have different spatial domains where
    !! they are defined. Therefore, it is necessary to
    !! check if a given particle has left these spatial domains to stop
    !! following it, otherwise this will cause an error in the simulation.
    TYPE(FIELDS), INTENT(IN)                                   :: F
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
    !! Particles' position in cylindrical coordinates,
    !! Y(1,:) = \(R\), Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
    INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
    !! Flag that determines whether particles are followed in the
    !! simulation (flag=1), or not (flag=0).
    INTEGER                                                :: IR
    !! Variable used to localize the grid cell in the \((R,\phi,Z)\)
    !! or \((R,Z)\) grid containing the fields data that corresponds
    !! to the radial position of the particles.
    INTEGER                                                :: IPHI
    !! Variable used to localize the grid cell in the \((R,\phi,Z)\)
    !! or \((R,Z)\) grid containing the fields data that corresponds
    !! to the azimuthal position of the particles.
    INTEGER                                                :: IZ
    !! Variable used to localize the grid cell in the \((R,\phi,Z)\)
    !! or \((R,Z)\) grid containing the fields data that corresponds
    !! to the vertical position of the particles.
    INTEGER(ip)                                            :: pp
    !! Particle iterator.
    INTEGER(ip)                                            :: ss
    !! Species iterator.
    REAL(rp) :: Rwall,xtmp,ytmp,ztmp,DiMESsurf,DiMESrad
    REAL(rp),DIMENSION(3) :: DiMESloc_cyl,DiMESloc_cart
    REAL(rp),DIMENSION(2) :: DiMESdims

    if (size(Y,1).eq.1) then
       ss = size(Y,1)
    else
       if (Y(2,1).eq.0) then
          ss=1_idef
       else
          ss = size(Y,1)
       end if
    endif

!    write(output_unit_write,'("R: ",E15.10)') Y(1,1)
!    write(output_unit_write,'("PHI: ",E15.10)') Y(2,1)
!    write(output_unit_write,'("Z: ",E15.10)') Y(1,3)

!    write(output_unit_write,*) 'Flag',flag(1)

    if (ALLOCATED(fields_domain%FLAG3D)) then
       if (F%Dim2x1t) then

          !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,IR,IZ) &
          !$OMP& SHARED(Y,flag,fields_domain,bfield_2X1T)
          do pp=1_idef,ss
             IR = INT(FLOOR((Y(pp,1)  - fields_domain%Ro + 0.5_rp* &
                  fields_domain%DR)/fields_domain%DR) + 1.0_rp,idef)

             IZ = INT(FLOOR((Y(pp,3)  + ABS(fields_domain%Zo) + 0.5_rp* &
                  fields_domain%DZ)/fields_domain%DZ) + 1.0_rp,idef)

             if ((fields_domain%FLAG3D(IR,1,IZ).NE.1_is).OR. &
                  ((IR.GT.bfield_2X1T%NR).OR.(IZ.GT.bfield_2X1T%NZ))) then
                if (F%Analytic_IWL.eq.'NONE') then
                   flag(pp) = 0_is
                else if (F%Analytic_IWL.eq.'D3D') then
                   if ((IR.lt.floor(bfield_2d%NR/6._rp)).and. &
                        (IZ.gt.floor(bfield_2d%NZ/5._rp)).and. &
                        (IZ.lt.floor(4._rp*bfield_2d%NZ/5._rp))) then

                      Rwall=F%circumradius*cos(C_PI/F%ntiles)/ &
                           (cos((modulo(Y(pp,2),2*C_PI/F%ntiles))-C_PI/F%ntiles))
                      if (Y(pp,1).lt.Rwall) flag(pp) = 0_is

                   else
                      flag(pp) = 0_is
                   endif
                endif
                !write(output_unit_write,'("YR:",E17.10)') Y(1,1)
                !write(output_unit_write,'("YZ:",E17.10)') Y(1,3)
                !write(output_unit_write,'("IR: ",I16)') IR
                !write(output_unit_write,'("IZ: ",I16)') IZ
             end if
          end do
          !$OMP END PARALLEL DO
       else
          !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,IR,IPHI,IZ) &
          !$OMP& SHARED(Y,flag,fields_domain,bfield_3d)
          do pp=1_idef,ss
             IR = INT(FLOOR((Y(pp,1)  - fields_domain%Ro + 0.5_rp* &
                  fields_domain%DR)/fields_domain%DR) + 1.0_rp,idef)
             IPHI = INT(FLOOR((Y(pp,2)  + 0.5_rp*fields_domain%DPHI)/ &
                  fields_domain%DPHI) + 1.0_rp,idef)
             IZ = INT(FLOOR((Y(pp,3)  + ABS(fields_domain%Zo) + 0.5_rp* &
                  fields_domain%DZ)/fields_domain%DZ) + 1.0_rp,idef)

             if ((fields_domain%FLAG3D(IR,IPHI,IZ).NE.1_is).OR. &
                  ((IR.GT.bfield_3d%NR).OR.(IZ.GT.bfield_3d%NZ))) then
                flag(pp) = 0_is
                !write(output_unit_write,'("YR:",E17.10)') Y(1,1)
                !write(output_unit_write,'("YZ:",E17.10)') Y(1,3)
                !write(output_unit_write,'("IR: ",I16)') IR
                !write(output_unit_write,'("IZ: ",I16)') IZ
             end if
          end do
          !$OMP END PARALLEL DO
       end if
    else if (ALLOCATED(fields_domain%FLAG2D)) then

       !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,IR,IZ) &
       !$OMP& SHARED(Y,flag,fields_domain,bfield_2d)
       do pp=1_idef,ss

          IR = INT(FLOOR((Y(pp,1)  - fields_domain%Ro + 0.5_rp* &
               fields_domain%DR)/fields_domain%DR) + 1.0_rp,idef)
          IZ = INT(FLOOR((Y(pp,3)  + ABS(fields_domain%Zo) + 0.5_rp* &
               fields_domain%DZ)/fields_domain%DZ) + 1.0_rp,idef)

#if DBG_CHECK
          if ((IR.lt.0).or.(IZ.lt.0).or.(IR.GT. &
               bfield_2d%NR).OR.(IZ.GT.bfield_2d%NZ)) then
             !write(6,'("YR:",E17.10)') Y(1,1)
             !write(6,'("YZ:",E17.10)') Y(1,3)
             !write(6,'("IR: ",I16)') IR
             !write(6,'("IZ: ",I16)') IZ
             !call KORC_ABORT(23)
          end if
#endif

          !write(output_unit_write,'("IR: ",I16)') IR
          !write(output_unit_write,'("IZ: ",I16)') IZ

          if (((IR.lt.0).or.(IZ.lt.0).or. &
               (IR.GT.bfield_2d%NR).OR.(IZ.GT.bfield_2d%NZ)).or. &
               (fields_domain%FLAG2D(IR,IZ).NE.1_is)) then
             !write(output_unit_write,*) 'here'
             if (F%Analytic_IWL.eq.'NONE') then
                flag(pp) = 0_is
             else if (F%Analytic_IWL.eq.'D3D') then
                if ((IR.lt.floor(bfield_2d%NR/6._rp)).and. &
                     (IZ.gt.floor(bfield_2d%NZ/5._rp)).and. &
                     (IZ.lt.floor(4._rp*bfield_2d%NZ/5._rp))) then

                   Rwall=F%circumradius*cos(C_PI/F%ntiles)/ &
                        (cos((modulo(Y(pp,2),2*C_PI/F%ntiles))-C_PI/F%ntiles))
                   if (Y(pp,1).lt.Rwall) flag(pp) = 0_is

                else
                   flag(pp) = 0_is
                endif
             endif
          end if

          if (F%useDiMES) THEN

             DiMESloc_cyl=F%DiMESloc
             DiMESloc_cyl(2)=C_PI*DiMESloc_cyl(2)/180._rp
             DiMESdims=F%DiMESdims

             if ((abs(Y(pp,1)-DiMESloc_cyl(1)).le.DiMESdims(1)).and. &
                  (abs(Y(pp,3)-DiMESloc_cyl(3)).le.DiMESdims(2)).and.&
                  (abs(Y(pp,2)-DiMESloc_cyl(2)).le.asin(DiMESdims(1)/DiMESloc_cyl(1)))) THEN

                xtmp=Y(pp,1)*cos(Y(pp,2))
                ytmp=Y(pp,1)*sin(Y(pp,2))
                ztmp=Y(pp,3)

                DiMESloc_cart(1)=DiMESloc_cyl(1)*cos(DiMESloc_cyl(2))
                DiMESloc_cart(2)=DiMESloc_cyl(1)*sin(DiMESloc_cyl(2))
                DiMESloc_cart(3)=DiMESloc_cyl(3)

                DiMESrad=DiMESdims(1)**2-(xtmp-DiMESloc_cart(1))**2-(ytmp-DiMESloc_cart(2))**2

                if (DiMESrad.le.0._rp) THEN
                   cycle
                end if

                DiMESsurf=DiMESloc_cart(3)+(DiMESdims(2)/DiMESdims(1))*sqrt(DiMESrad)

                if (ztmp.le.DiMESsurf) flag(pp)=0_is

             end if

          end if

       end do
       !$OMP END PARALLEL DO

    end if
  end subroutine check_if_in_fields_domain

  subroutine check_if_in_fields_domain_p(pchunk,F,Y_R,Y_PHI,Y_Z,flag)
    !! @note Subrotuine that checks if particles in the simulation are within
    !! the spatial domain where interpolants and fields are known. @endnote
    !! External fields and interpolants can have different spatial domains where
    !! they are defined. Therefore, it is necessary to
    !! check if a given particle has left these spatial domains to
    !! stop following it, otherwise this will cause an error in the simulation.
    INTEGER, INTENT(IN)  :: pchunk
    TYPE(FIELDS), INTENT(IN)                                   :: F
    REAL(rp), DIMENSION(pchunk),  INTENT(IN)      :: Y_R,Y_PHI,Y_Z
    INTEGER(is), DIMENSION(pchunk), INTENT(INOUT)  :: flag
    !! Flag that determines whether particles are followed in the
    !! simulation (flag=1), or not (flag=0).
    INTEGER                                                :: IR
    !! Variable used to localize the grid cell in the \((R,\phi,Z)\)
    !! or \((R,Z)\) grid containing the fields data that corresponds
    !! to the radial position of the particles.
    INTEGER                                                :: IPHI
    !! Variable used to localize the grid cell in the \((R,\phi,Z)\)
    !! or \((R,Z)\) grid containing the fields data that corresponds
    !! to the azimuthal position of the particles.
    INTEGER                                                :: IZ
    !! Variable used to localize the grid cell in the \((R,\phi,Z)\)
    !! or \((R,Z)\) grid containing the fields data that corresponds
    !! to the vertical position of the particles.
    INTEGER(ip)                                            :: pp
    !! Particle iterator.
    INTEGER(ip)                                            :: ss
    !! Species iterator.
    REAL(rp) :: Rwall,lscale,xtmp,ytmp,ztmp,DiMESsurf,DiMESrad
    REAL(rp),DIMENSION(3) :: DiMESloc_cyl,DiMESloc_cart
    REAL(rp),DIMENSION(2) :: DiMESdims

    lscale=7.983650528137727e-04

    !write(output_unit_write,'("YR:",E17.10)') Y_R
    !write(output_unit_write,'("YPHI:",E17.10)') Y_PHI
    !write(output_unit_write,'("YZ:",E17.10)') Y_Z

    !write(output_unit_write,'("Ro:",E17.10)') fields_domain%Ro
    !write(output_unit_write,'("Zo:",E17.10)') fields_domain%Zo
    !write(output_unit_write,'("DR:",E17.10)') fields_domain%DR
    !write(output_unit_write,'("DZ:",E17.10)') fields_domain%DZ
    !write(output_unit_write,'("DT:",E17.10)') fields_domain%DT


    if (ALLOCATED(fields_domain%FLAG3D)) then
       if (F%Dim2x1t) then
          !$OMP SIMD
          !       !$OMP&  aligned(IR,IPHI,IZ)
          do pp=1_idef,pchunk

             IR = INT(FLOOR((Y_R(pp)  - fields_domain%Ro + &
                  0.5_rp*fields_domain%DR)/fields_domain%DR) + 1.0_rp,idef)

             IZ = INT(FLOOR((Y_Z(pp)  + ABS(fields_domain%Zo) + &
                  0.5_rp*fields_domain%DZ)/fields_domain%DZ) + 1.0_rp,idef)

!             write(output_unit_write,'("IR: ",I16)') IR
!             write(output_unit_write,'("IPHI: ",I16)') IPHI
!             write(output_unit_write,'("IZ: ",I16)') IZ

             if ((fields_domain%FLAG3D(IR,1,IZ).NE.1_is).OR. &
                  ((IR.GT.bfield_2X1T%NR).OR.(IZ.GT.bfield_2X1T%NZ))) then
                if (F%Analytic_IWL.eq.'NONE') then
                   flag(pp) = 0_is
                else if (F%Analytic_IWL.eq.'D3D') then
                   if ((IR.lt.floor(bfield_2d%NR/6._rp)).and. &
                        (IZ.gt.floor(bfield_2d%NZ/5._rp)).and. &
                        (IZ.lt.floor(4._rp*bfield_2d%NZ/5._rp))) then

                      Rwall=F%circumradius*cos(C_PI/F%ntiles)/ &
                           (cos((modulo(Y_PHI(pp),2*C_PI/F%ntiles))-C_PI/F%ntiles))
                      if (Y_R(pp).lt.Rwall) flag(pp) = 0_is

                   else
                      flag(pp) = 0_is
                   endif
                endif

                !write(output_unit_write,'("YR:",E17.10)') Y_R(pp)
                !write(output_unit_write,'("YPHI:",E17.10)') Y_PHI(pp)
                !write(output_unit_write,'("YZ:",E17.10)') Y_Z(pp)

                !write(output_unit_write,'("IR: ",I16)') IR
                !write(output_unit_write,'("IPHI: ",I16)') IPHI
                !write(output_unit_write,'("IZ: ",I16)') IZ

                !call KORC_ABORT()

             end if

             !write(output_unit_write,'("IPHI: ",I16)') IPHI
             !write(output_unit_write,'("flag: ",I16)') flag(pp)


          end do
          !$OMP END SIMD
       else
          !$OMP SIMD
          !       !$OMP&  aligned(IR,IPHI,IZ)
          do pp=1_idef,pchunk

             IR = INT(FLOOR((Y_R(pp)  - fields_domain%Ro + &
                  0.5_rp*fields_domain%DR)/fields_domain%DR) + 1.0_rp,idef)
             IPHI = INT(FLOOR((Y_PHI(pp)  + 0.5_rp*fields_domain%DPHI)/ &
                  fields_domain%DPHI) + 1.0_rp,idef)
             IZ = INT(FLOOR((Y_Z(pp)  + ABS(fields_domain%Zo) + &
                  0.5_rp*fields_domain%DZ)/fields_domain%DZ) + 1.0_rp,idef)

             if ((fields_domain%FLAG3D(IR,IPHI,IZ).NE.1_is).OR. &
                  ((IR.GT.bfield_3d%NR).OR.(IZ.GT.bfield_3d%NZ))) then
                flag(pp) = 0_is

                !write(output_unit_write,'("YR:",E17.10)') Y_R
                !write(output_unit_write,'("YPHI:",E17.10)') Y_PHI
                !write(output_unit_write,'("YZ:",E17.10)') Y_Z

                !write(output_unit_write,'("IR: ",I16)') IR
                !write(output_unit_write,'("IPHI: ",I16)') IPHI
                !write(output_unit_write,'("IZ: ",I16)') IZ

                !call KORC_ABORT()

             end if

             !write(output_unit_write,'("IPHI: ",I16)') IPHI
             !write(output_unit_write,'("flag: ",I16)') flag(pp)


          end do
          !$OMP END SIMD
       end if
    else if (ALLOCATED(fields_domain%FLAG2D)) then
       !$OMP SIMD
!       !$OMP& aligned(IR,IZ)
       do pp=1_idef,pchunk

          !write(output_unit_write,*) Y_R(pp),Y_Z(pp),fields_domain%Ro,fields_domain%DR,fields_domain%Zo,fields_domain%DZ

          IR = INT(FLOOR((Y_R(pp)  - fields_domain%Ro + &
               0.5_rp*fields_domain%DR)/fields_domain%DR) + 1.0_rp,idef)
          IZ = INT(FLOOR((Y_Z(pp)  + ABS(fields_domain%Zo) + &
               0.5_rp*fields_domain%DZ)/fields_domain%DZ) + 1.0_rp,idef)

!          write(output_unit_write,*) pp

!          write(output_unit_write,'("Size of fields_domain R: ",I16)') &
!               size(fields_domain%FLAG2D,1)
!          write(output_unit_write,'("Size of fields_domain Z: ",I16)') &
!               size(fields_domain%FLAG2D,2)

!          if ((IR.lt.0).or.(IZ.lt.0)) then
!             write(output_unit_write,'("YR:",E17.10)') Y_R(pp)
!             write(output_unit_write,'("YZ:",E17.10)') Y_Z(pp)
!             write(output_unit_write,'("IR: ",I16)') IR
!             write(output_unit_write,'("IZ: ",I16)') IZ
!          end if

#if DBG_CHECK
          if ((IR.lt.1).or.(IZ.lt.1).or.(IR.GT. &
               bfield_2d%NR).OR.(IZ.GT.bfield_2d%NZ)) then

             !write(6,'("YR:",E17.10)') Y_R(pp)
             !write(6,'("YZ:",E17.10)') Y_Z(pp)
             !write(6,'("IR: ",I16)') IR
             !write(6,'("IZ: ",I16)') IZ
             !write(6,*) 'NR',bfield_2d%NR
             !write(6,*) 'NZ',bfield_2d%NZ
             !write(6,*) 'FLAG',fields_domain%FLAG2D(IR,IZ)
             !call KORC_ABORT(23)
          end if
#endif

          if ((IR.lt.1).or.(IZ.lt.1).or. &
               (IR.GT.bfield_2d%NR).OR.(IZ.GT.bfield_2d%NZ).or. &
               (fields_domain%FLAG2D(IR,IZ).NE.1_is)) then


             !write(6,*) 'Y',Y_R(pp)*lscale,Y_PHI(pp),Y_Z(pp)*lscale
             !write(6,*) 'flagCon',flag(pp)
             !write(6,*) 'IFLAG',IR,IZ
             !write(6,*) 'NFLAG',bfield_2d%NR,bfield_2d%NZ
             !write(6,*) 'BFLAG',floor(bfield_2d%NR/6._rp),floor(bfield_2d%NR/5._rp),floor(4._rp*bfield_2d%NR/5._rp)
             !write(6,*) 'FLAGR',fields_domain%Ro*lscale,fields_domain%DR*lscale
             !write(6,*) 'FLAGZ',fields_domain%Zo*lscale,fields_domain%DZ*lscale

             if (F%Analytic_IWL.eq.'NONE') then
                flag(pp) = 0_is
             else if (F%Analytic_IWL.eq.'D3D') then
                if ((IR.lt.floor(bfield_2d%NR/6._rp)).and. &
                     (IZ.gt.floor(bfield_2d%NZ/5._rp)).and. &
                     (IZ.lt.floor(4._rp*bfield_2d%NZ/5._rp))) then

                   Rwall=F%circumradius*cos(C_PI/F%ntiles)/ &
                        (cos(modulo(Y_PHI(pp),2*C_PI/F%ntiles)-C_PI/F%ntiles))

                   !write(6,*) 'Rc,nt',F%circumradius*lscale,F%ntiles
                   !write(6,*) 'Rwall',Rwall*lscale
                   !write(6,*) 'mod',modulo(Y_PHI(pp),2*C_PI/F%ntiles)

                   if (Y_R(pp).lt.Rwall) flag(pp) = 0_is

                else
                   flag(pp) = 0_is
                endif
             endif
!             write(output_unit_write,'("Shit''s fucked.")')
          end if

          if (F%useDiMES) THEN

            DiMESloc_cyl=F%DiMESloc
            DiMESloc_cyl(2)=C_PI*DiMESloc_cyl(2)/180._rp
            DiMESdims=F%DiMESdims

            if ((abs(Y_R(pp)-DiMESloc_cyl(1)).le.DiMESdims(1)).and. &
                 (abs(Y_Z(pp)-DiMESloc_cyl(3)).le.DiMESdims(2)).and.&
                 (abs(Y_PHI(pp)-DiMESloc_cyl(2)).le.asin(DiMESdims(1)/DiMESloc_cyl(1)))) THEN

               xtmp=Y_R(pp)*cos(Y_PHI(pp))
               ytmp=Y_R(pp)*sin(Y_PHI(pp))
               ztmp=Y_Z(pp)

               DiMESloc_cart(1)=DiMESloc_cyl(1)*cos(DiMESloc_cyl(2))
               DiMESloc_cart(2)=DiMESloc_cyl(1)*sin(DiMESloc_cyl(2))
               DiMESloc_cart(3)=DiMESloc_cyl(3)

               DiMESrad=DiMESdims(1)**2-(xtmp-DiMESloc_cart(1))**2-(ytmp-DiMESloc_cart(2))**2

               if (DiMESrad.le.0._rp) THEN
                  cycle
               end if

               DiMESsurf=DiMESloc_cart(3)+(DiMESdims(2)/DiMESdims(1))*sqrt(DiMESrad)

               if (ztmp.le.DiMESsurf) flag(pp)=0_is

            end if

         end if

       end do
       !$OMP END SIMD
!       write(output_unit_write,'("Shit''s not fucked.")')
    end if
  end subroutine check_if_in_fields_domain_p

  subroutine check_if_in_LCFS(F,Y,inLCFS)
    !! @note Subrotuine that checks if particles in the simulation are within
    !! the spatial domain where interpolants and fields are known. @endnote
    !! External fields and interpolants can have different spatial domains where
    !! they are defined. Therefore, it is necessary to
    !! check if a given particle has left these spatial domains to stop
    !! following it, otherwise this will cause an error in the simulation.
    TYPE(FIELDS), INTENT(IN)                                   :: F
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
    !! Particles' position in cylindrical coordinates,
    !! Y(1,:) = \(R\), Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
    INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: inLCFS
    !! Flag that determines whether particles are followed in the
    !! simulation (flag=1), or not (flag=0).
    INTEGER                                                :: IR
    !! Variable used to localize the grid cell in the \((R,\phi,Z)\)
    !! or \((R,Z)\) grid containing the fields data that corresponds
    !! to the radial position of the particles.
    INTEGER                                                :: IPHI
    !! Variable used to localize the grid cell in the \((R,\phi,Z)\)
    !! or \((R,Z)\) grid containing the fields data that corresponds
    !! to the azimuthal position of the particles.
    INTEGER                                                :: IZ
    !! Variable used to localize the grid cell in the \((R,\phi,Z)\)
    !! or \((R,Z)\) grid containing the fields data that corresponds
    !! to the vertical position of the particles.
    INTEGER(ip)                                            :: pp
    !! Particle iterator.
    INTEGER(ip)                                            :: ss
    !! Species iterator.
    REAL(rp) :: Rwall

    if (size(Y,1).eq.1) then
       ss = size(Y,1)
    else
       if (Y(2,1).eq.0) then
          ss=1_idef
       else
          ss = size(Y,1)
       end if
    endif

!    write(output_unit_write,'("R: ",E15.10)') Y(1,1)
!    write(output_unit_write,'("PHI: ",E15.10)') Y(2,1)
!    write(output_unit_write,'("Z: ",E15.10)') Y(1,3)

    !    write(output_unit_write,*) 'Flag',flag(1)

    !write(6,*) 'ss',ss
    !write(6,*) 'inLCFS',inLCFS



    !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,IR,IZ) &
    !$OMP& SHARED(Y,inLCFS,fields_domain,bfield_2d)
    do pp=1_idef,ss

       IR = INT(FLOOR((Y(pp,1)  - fields_domain%Ro + 0.5_rp* &
            fields_domain%DR)/fields_domain%DR) + 1.0_rp,idef)
       IZ = INT(FLOOR((Y(pp,3)  + ABS(fields_domain%Zo) + 0.5_rp* &
            fields_domain%DZ)/fields_domain%DZ) + 1.0_rp,idef)

#if DBG_CHECK
       if ((IR.lt.0).or.(IZ.lt.0).or.(IR.GT. &
            bfield_2d%NR).OR.(IZ.GT.bfield_2d%NZ)) then
          !write(6,'("YR:",E17.10)') Y(1,1)
          !write(6,'("YZ:",E17.10)') Y(1,3)
          !write(6,'("IR: ",I16)') IR
          !write(6,'("IZ: ",I16)') IZ
          !call KORC_ABORT(23)
       end if
#endif

       !write(output_unit_write,'("IR: ",I16)') IR
       !write(output_unit_write,'("IZ: ",I16)') IZ

       if (((IR.lt.0).or.(IZ.lt.0).or. &
            (IR.GT.bfield_2d%NR).OR.(IZ.GT.bfield_2d%NZ)).or. &
            (fields_domain%LCFS2D(IR,IZ).NE.1_is)) then
          !write(output_unit_write,*) 'here'

          inLCFS(pp) = 0_is

       end if
    end do
    !$OMP END PARALLEL DO


  end subroutine check_if_in_LCFS

  subroutine initialize_profiles_interpolant(params,P)
    !! @note Subroutine that initializes plasma profiles interpolants. @endnote
    !! This subroutine initializes either 2-D or 3-D PSPLINE interpolants
    !! using the data of plasma profiles in the KORC-dervied-type variable P.
    TYPE(KORC_PARAMS), INTENT(IN) :: params
    !! Core KORC simulation parameters.
    TYPE(PROFILES), INTENT(INOUT) :: P
    !! An instance of KORC's derived type PROFILES containing
    !! all the information about the plasma profiles used in the simulation.
    !! See [[korc_types]] and [[korc_profiles]].

!#ifdef FIO
!    P%FIO_ne = -1
!    P%FIO_te = -1
!    P%FIO_zeff = -1
!#endif

    if (params%collisions.or.params%radiation) then
       if (params%profile_model(1:8) .EQ. 'EXTERNAL') then

          if (params%mpi_params%rank .EQ. 0) then
             write(output_unit_write,'("* * * * INITIALIZING PROFILES INTERPOLANT * * * *")')
          end if

          if (P%axisymmetric) then

             if (params%mpi_params%rank .EQ. 0) then
                write(output_unit_write,*) '2D ne, Te, Zeff'
                flush(output_unit_write)
             end if

             profiles_2d%NR = P%dims(1)
             profiles_2d%NZ = P%dims(3)

!             write(output_unit_write,'("NR",I15)') profiles_2d%NR
!             write(output_unit_write,'("NZ",I15)') profiles_2d%NR

             ! Initializing ne
             call EZspline_init(profiles_2d%ne,profiles_2d%NR,profiles_2d%NZ, &
                  profiles_2d%BCSR,profiles_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             profiles_2d%ne%x1 = P%X%R
             profiles_2d%ne%x2 = P%X%Z

             call EZspline_setup(profiles_2d%ne, P%ne_2D, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             ! Initializing Te
             call EZspline_init(profiles_2d%Te,profiles_2d%NR,profiles_2d%NZ, &
                  profiles_2d%BCSR,profiles_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             profiles_2d%Te%x1 = P%X%R
             profiles_2d%Te%x2 = P%X%Z

!             write(output_unit_write,'("Te_interp_R",E17.10)') profiles_2d%Te%x1
!             write(output_unit_write,'("Te_interp_Z",E17.10)') profiles_2d%Te%x2

!             write(output_unit_write,'("Te",E17.10)') P%Te_2D(10,:)

             call EZspline_setup(profiles_2d%Te, P%Te_2D, ezerr, .TRUE.)
             call EZspline_error(ezerr)


             ! Initializing Zeff
             call EZspline_init(profiles_2d%Zeff,profiles_2d%NR, &
                  profiles_2d%NZ,profiles_2d%BCSR,profiles_2d%BCSZ,ezerr)
             call EZspline_error(ezerr)

             profiles_2d%Zeff%x1 = P%X%R
             profiles_2d%Zeff%x2 = P%X%Z

             call EZspline_setup(profiles_2d%Zeff, P%Zeff_2D, ezerr, .TRUE.)
             call EZspline_error(ezerr)

             if (params%profile_model(10:10) .EQ. 'H') then

                if (params%mpi_params%rank .EQ. 0) then
                   write(output_unit_write,*) '2D Hollmann impurities'
                   flush(output_unit_write)
                end if

                ! Initializing RHON
                call EZspline_init(profiles_2d%RHON,profiles_2d%NR, &
                     profiles_2d%NZ,profiles_2d%BCSR,profiles_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                profiles_2d%RHON%x1 = P%X%R
                profiles_2d%RHON%x2 = P%X%Z

                call EZspline_setup(profiles_2d%RHON, P%RHON, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                !write(output_unit_write,'("profiles_2d%RHON: ",E17.10)') profiles_2d%RHON%fspl(1,:,:)

                ! Initializing nRE
                call EZspline_init(profiles_2d%nRE,profiles_2d%NR, &
                     profiles_2d%NZ,profiles_2d%BCSR,profiles_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                profiles_2d%nRE%x1 = P%X%R
                profiles_2d%nRE%x2 = P%X%Z

                call EZspline_setup(profiles_2d%nRE, P%nRE_2D, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                ! Initializing nAr0
                call EZspline_init(profiles_2d%nAr0,profiles_2d%NR, &
                     profiles_2d%NZ,profiles_2d%BCSR,profiles_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                profiles_2d%nAr0%x1 = P%X%R
                profiles_2d%nAr0%x2 = P%X%Z

                call EZspline_setup(profiles_2d%nAr0, P%nAr0_2D, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                ! Initializing nAr1
                call EZspline_init(profiles_2d%nAr1,profiles_2d%NR, &
                     profiles_2d%NZ,profiles_2d%BCSR,profiles_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                profiles_2d%nAr1%x1 = P%X%R
                profiles_2d%nAr1%x2 = P%X%Z

                call EZspline_setup(profiles_2d%nAr1, P%nAr1_2D, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                ! Initializing nAr2
                call EZspline_init(profiles_2d%nAr2,profiles_2d%NR, &
                     profiles_2d%NZ,profiles_2d%BCSR,profiles_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                profiles_2d%nAr2%x1 = P%X%R
                profiles_2d%nAr2%x2 = P%X%Z

                call EZspline_setup(profiles_2d%nAr2, P%nAr2_2D, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                ! Initializing nAr3
                call EZspline_init(profiles_2d%nAr3,profiles_2d%NR, &
                     profiles_2d%NZ,profiles_2d%BCSR,profiles_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                profiles_2d%nAr3%x1 = P%X%R
                profiles_2d%nAr3%x2 = P%X%Z

                call EZspline_setup(profiles_2d%nAr3, P%nAr3_2D, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                ! Initializing nD
                call EZspline_init(profiles_2d%nD,profiles_2d%NR, &
                     profiles_2d%NZ,profiles_2d%BCSR,profiles_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                profiles_2d%nD%x1 = P%X%R
                profiles_2d%nD%x2 = P%X%Z

                call EZspline_setup(profiles_2d%nD, P%nD_2D, ezerr, .TRUE.)
                call EZspline_error(ezerr)

                ! Initializing nD1
                call EZspline_init(profiles_2d%nD1,profiles_2d%NR, &
                     profiles_2d%NZ,profiles_2d%BCSR,profiles_2d%BCSZ,ezerr)
                call EZspline_error(ezerr)

                profiles_2d%nD1%x1 = P%X%R
                profiles_2d%nD1%x2 = P%X%Z

                call EZspline_setup(profiles_2d%nD1, P%nD1_2D, ezerr, .TRUE.)
                call EZspline_error(ezerr)

             end if

             ALLOCATE(profiles_domain%FLAG2D(profiles_2d%NR,profiles_2d%NZ))
             profiles_domain%FLAG2D = P%FLAG2D

             profiles_domain%DR = ABS(P%X%R(2) - P%X%R(1))
             profiles_domain%DZ = ABS(P%X%Z(2) - P%X%Z(1))
          else

             if (params%mpi_params%rank .EQ. 0) then
                write(output_unit_write,*) '3D ne, Te, Zeff'
             end if

             profiles_3d%NR = P%dims(1)
             profiles_3d%NPHI = P%dims(2)
             profiles_3d%NZ = P%dims(3)

             ! Initializing ne
             call EZspline_init(profiles_3d%ne, profiles_3d%NR, &
                  profiles_3d%NPHI, profiles_3d%NZ,&
                  profiles_3d%BCSR, profiles_3d%BCSPHI, profiles_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             profiles_3d%ne%x1 = P%X%R
             ! profiles_3d%ne%x2 = P%X%PHI
             profiles_3d%ne%x3 = P%X%Z

             call EZspline_setup(profiles_3d%ne, P%ne_3D, ezerr)
             call EZspline_error(ezerr)

             ! Initializing Te
             call EZspline_init(profiles_3d%Te, profiles_3d%NR, &
                  profiles_3d%NPHI, profiles_3d%NZ,&
                  profiles_3d%BCSR, profiles_3d%BCSPHI, profiles_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             profiles_3d%Te%x1 = P%X%R
             ! profiles_3d%Te%x2 = P%X%PHI
             profiles_3d%Te%x3 = P%X%Z

             call EZspline_setup(profiles_3d%Te, P%Te_3D, ezerr)
             call EZspline_error(ezerr)

             ! Initializing Zeff
             call EZspline_init(profiles_3d%Zeff, profiles_3d%NR, &
                  profiles_3d%NPHI, profiles_3d%NZ,&
                  profiles_3d%BCSR, profiles_3d%BCSPHI, profiles_3d%BCSZ, ezerr)
             call EZspline_error(ezerr)

             profiles_3d%Zeff%x1 = P%X%R
             ! profiles_3d%Zeff%x2 = P%X%PHI
             profiles_3d%Zeff%x3 = P%X%Z

             call EZspline_setup(profiles_3d%Zeff, P%Zeff_3D, ezerr)
             call EZspline_error(ezerr)

             ALLOCATE(profiles_domain%FLAG3D(profiles_3d%NR,profiles_3d%NPHI, &
                  profiles_3d%NZ))
             profiles_domain%FLAG3D = P%FLAG3D

             profiles_domain%DR = ABS(P%X%R(2) - P%X%R(1))
             profiles_domain%DPHI = 2.0_rp*C_PI/profiles_3d%NPHI
             profiles_domain%DZ = ABS(P%X%Z(2) - P%X%Z(1))
          end if

          profiles_domain%Ro = P%X%R(1)
          profiles_domain%Zo = P%X%Z(1)

          if (params%mpi_params%rank .EQ. 0) then
             write(output_unit_write,'("* * * * * * INTERPOLANT   INITIALIZED * * * * * *",/)')
          end if
       else if (params%profile_model(1:10) .EQ. 'ANALYTICAL') then
          if (params%mpi_params%rank .EQ. 0) then
             write(output_unit_write,'("* * * * USING ANALYTICAL PROFILES * * * *",/)')
          end if
       else if (params%profile_model .EQ. 'UNIFORM') then
          if (params%mpi_params%rank .EQ. 0) then
             write(output_unit_write,'("* * * * UNIFORM PLASMA: NO PROFILES USED * * * *",/)')
          end if
       end if
    end if
  end subroutine initialize_profiles_interpolant


  subroutine check_if_in_profiles_domain(Y,flag)
    !! @note Subrotuine that checks if particles in the simulation are
    !! within the spatial domain where interpolants and plasma profiles
    !! are known. @endnote
    !!External plasma profiles and interpolants can have different spatial
    !! domains where they are defined. Therefore, it is necessary to check
    !! if a given particle has left these spatial domains to stop following
    !! it, otherwise this will cause an error in the simulation.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
    !! Particles' position in cylindrical coordinates,
    !! Y(1,:) = \(R\), Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
    INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
    !! Flag that determines whether particles are followed
    !! in the simulation (flag=1), or not (flag=0).
    INTEGER                                                :: IR
    !! @param IR Variable used to localize the grid cell in
    !! the \((R,\phi,Z)\) or \((R,Z)\) grid containing the fields data that
    !! corresponds to the radial position of the particles.
    INTEGER                                                :: IPHI
    !! @param IPHI Variable used to localize the grid cell in
    !! the \((R,\phi,Z)\) or \((R,Z)\) grid containing the fields data that
    !! corresponds to the azimuthal position of the particles.
    INTEGER                                                :: IZ
    !! @param IZ Variable used to localize the grid cell in the
    !! \((R,\phi,Z)\) or \((R,Z)\) grid containing the fields data that
    !! corresponds to the vertical position of the particles.
    INTEGER(ip)                                            :: pp
    !! @param pp Particle iterator.
    INTEGER(ip)                                            :: ss
    !! @param ss Species iterator.

    if (size(Y,1).eq.1) then
       ss = size(Y,1)
    else
       if (Y(2,1).eq.0) then
          ss=1_idef
       else
          ss = size(Y,1)
       end if
    endif

    if (ALLOCATED(profiles_domain%FLAG3D)) then
       !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,IR,IPHI,IZ) &
       !$OMP& SHARED(Y,flag,profiles_domain,profiles_3d)
       do pp=1_idef,ss
          IR = INT(FLOOR((Y(pp,1)  - profiles_domain%Ro + &
               0.5_rp*profiles_domain%DR)/profiles_domain%DR) + 1.0_rp,idef)
          IPHI = INT(FLOOR((Y(pp,2)  + 0.5_rp*profiles_domain%DPHI)/ &
               profiles_domain%DPHI) + 1.0_rp,idef)
          IZ = INT(FLOOR((Y(pp,3)  + ABS(profiles_domain%Zo) + &
               0.5_rp*profiles_domain%DZ)/profiles_domain%DZ) + 1.0_rp,idef)

          if ((profiles_domain%FLAG3D(IR,IPHI,IZ).NE.1_is).OR. &
               ((IR.GT.profiles_3d%NR).OR.(IZ.GT.profiles_3d%NZ))) then
             flag(pp) = 0_is
          end if
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,IR,IZ) &
       !$OMP& SHARED(Y,flag,profiles_domain,profiles_2d)
       do pp=1_idef,ss
          IR = INT(FLOOR((Y(pp,1)  - profiles_domain%Ro + &
               0.5_rp*profiles_domain%DR)/profiles_domain%DR) + 1.0_rp,idef)
          IZ = INT(FLOOR((Y(pp,3)  + ABS(profiles_domain%Zo) + &
               0.5_rp*profiles_domain%DZ)/profiles_domain%DZ) + 1.0_rp,idef)

          if ((profiles_domain%FLAG2D(IR,IZ).NE.1_is).OR. &
               ((IR.GT.profiles_2d%NR).OR.(IZ.GT.profiles_2d%NZ))) then
             flag(pp) = 0_is
          end if
       end do
       !$OMP END PARALLEL DO
    end if
  end subroutine check_if_in_profiles_domain

  subroutine check_if_in_profiles_domain_p(pchunk,Y_R,Y_PHI,Y_Z,flag)
    INTEGER, INTENT(IN)  :: pchunk
    REAL(rp), DIMENSION(pchunk),  INTENT(IN)      :: Y_R,Y_PHI,Y_Z
    INTEGER(is), DIMENSION(pchunk), INTENT(INOUT)  :: flag
    !! Flag that determines whether particles are followed
    !! in the simulation (flag=1), or not (flag=0).
    INTEGER                                                :: IR
    !! @param IR Variable used to localize the grid cell in
    !! the \((R,\phi,Z)\) or \((R,Z)\) grid containing the fields data that
    !! corresponds to the radial position of the particles.
    INTEGER                                                :: IPHI
    !! @param IPHI Variable used to localize the grid cell in
    !! the \((R,\phi,Z)\) or \((R,Z)\) grid containing the fields data that
    !! corresponds to the azimuthal position of the particles.
    INTEGER                                                :: IZ
    !! @param IZ Variable used to localize the grid cell in the
    !! \((R,\phi,Z)\) or \((R,Z)\) grid containing the fields data that
    !! corresponds to the vertical position of the particles.
    INTEGER(ip)                                            :: pp
    !! @param pp Particle iterator.
    INTEGER(ip)                                            :: ss
    !! @param ss Species iterator.


    if (ALLOCATED(profiles_domain%FLAG3D)) then
       !$OMP SIMD
!       !$OMP& aligned(IR,IPHI,IZ)
       do pp=1_idef,pchunk
          IR = INT(FLOOR((Y_R(pp)  - profiles_domain%Ro + &
               0.5_rp*profiles_domain%DR)/profiles_domain%DR) + 1.0_rp,idef)
          IPHI = INT(FLOOR((Y_PHI(pp)  + 0.5_rp*profiles_domain%DPHI)/ &
               profiles_domain%DPHI) + 1.0_rp,idef)
          IZ = INT(FLOOR((Y_Z(pp)  + ABS(profiles_domain%Zo) + &
               0.5_rp*profiles_domain%DZ)/profiles_domain%DZ) + 1.0_rp,idef)

          if ((profiles_domain%FLAG3D(IR,IPHI,IZ).NE.1_is).OR. &
               ((IR.GT.profiles_3d%NR).OR.(IZ.GT.profiles_3d%NZ))) then
             flag(pp) = 0_is
          end if
       end do
       !$OMP END SIMD
    else
       !$OMP SIMD
!       !$OMP& aligned(IR,IZ)
       do pp=1_idef,pchunk
          IR = INT(FLOOR((Y_R(pp)  - profiles_domain%Ro + &
               0.5_rp*profiles_domain%DR)/profiles_domain%DR) + 1.0_rp,idef)
          IZ = INT(FLOOR((Y_Z(pp)  + ABS(profiles_domain%Zo) + &
               0.5_rp*profiles_domain%DZ)/profiles_domain%DZ) + 1.0_rp,idef)

          if ((profiles_domain%FLAG2D(IR,IZ).NE.1_is).OR. &
               ((IR.GT.profiles_2d%NR).OR.(IZ.GT.profiles_2d%NZ))) then
             flag(pp) = 0_is

!             write(output_unit_write,'("Shit''s fucked.")')
          end if
       end do
       !$OMP END SIMD
!       write(output_unit_write,'("Shit''s not fucked.")')
    end if
  end subroutine check_if_in_profiles_domain_p

subroutine interp_2D_bfields(params,Y,B,flag)
  !! @note Subroutine for interpolating the pre-computed, axisymmetric magnetic
  !! field to the particles' position. @endnote
  TYPE(KORC_PARAMS), INTENT(IN)      :: params
  !! Core KORC simulation parameters.
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
  !! Particles' position in cylindrical coordinates, Y(1,:) = \(R\),
  !! Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: B
  !! Cartesian components of interpolated magnetic field components.
  !! B(1,:)=\(B_x\), B(2,:)=\(B_y\), and B(3,:)=\(B_z\).
  REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: F
  !! Cylindrical components of interpolated magnetic field components.
  !! F(1,:)=\(B_R\), F(2,:)=\(B_\phi\), and F(3,:)=\(B_Z\).
  INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
  !! Flag that indicates whether particles are followed in the simulation
  !! (flag=1), or not (flag=0).
  INTEGER                                                :: pp
  !! Particle iterator.
  INTEGER                                                :: ss
  !! Species iterator.

  ss = size(Y,1)

  ALLOCATE(F(ss,3))
  !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
  !$OMP& SHARED(F,Y,B,flag,bfield_2d)
  do pp=1_idef,ss
     if ( flag(pp) .EQ. 1_is ) then
        call EZspline_interp(bfield_2d%R, Y(pp,1), Y(pp,3), F(pp,1), ezerr)
        call EZspline_error(ezerr)

        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        end if

        call EZspline_interp(bfield_2d%PHI, Y(pp,1), Y(pp,3), F(pp,2), ezerr)
        call EZspline_error(ezerr)

        call EZspline_interp(bfield_2d%Z, Y(pp,1), Y(pp,3), F(pp,3), ezerr)
        call EZspline_error(ezerr)

        if (.not.params%GC_coords) then
           B(pp,1) = F(pp,1)*COS(Y(pp,2)) - F(pp,2)*SIN(Y(pp,2))
           B(pp,2) = F(pp,1)*SIN(Y(pp,2)) + F(pp,2)*COS(Y(pp,2))
           B(pp,3) = F(pp,3)
        else
           B(pp,1) = F(pp,1)
           B(pp,2) = F(pp,2)
           B(pp,3) = F(pp,3)
        end if
     end if
  end do
  !$OMP END PARALLEL DO
  DEALLOCATE(F)
end subroutine interp_2D_bfields

subroutine gradient_2D_Bfields(Y,BR,BPHI,BZ,flag)
  !! @note Subroutine for interpolating the pre-computed, axisymmetric
  !! gradient of the magnitude of themagnetic field to the particles'
  !! position. Stored as cylindrical components of field. @endnote
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
  !! Particles' position in cylindrical coordinates, Y(1,:) = \(R\),
  !! Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: BR
  !! Cylindrical components of gradient of R-component of magnetic field.
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: BPHI
  !! Cylindrical components of gradient of R-component of magnetic field.
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: BZ
  !! Cylindrical components of gradient of R-component of magnetic field.
  REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: F
  !! Cylindrical components of interpolated magnetic field components.
  !! F(1,:)=\(B_R\), F(2,:)=\(B_\phi\), and F(3,:)=\(B_Z\).
  INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
  !! Flag that indicates whether particles are followed in the simulation
  !! (flag=1), or not (flag=0).
  INTEGER                                                :: pp
  !! Particle iterator.
  INTEGER                                                :: ss
  !! Species iterator.

  ss = size(Y,1)

  ALLOCATE(F(2,ss))
  !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
  !$OMP& SHARED(Y,BR,BPHI,BZ,flag,bfield_2d)
  do pp=1_idef,ss
     if ( flag(pp) .EQ. 1_is ) then
        call EZspline_gradient(bfield_2d%R, Y(pp,1), Y(pp,3), F(:,pp), ezerr)
        call EZspline_error(ezerr)

        BR(pp,1) = F(pp,1)
        BR(pp,2) = 0._rp
        BR(pp,3) = F(pp,2)

        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        end if

        call EZspline_gradient(bfield_2d%PHI, Y(pp,1), Y(pp,3), F(:,pp), &
             ezerr)
        call EZspline_error(ezerr)

        BPHI(pp,1) = F(pp,1)
        BPHI(pp,2) = 0._rp
        BPHI(pp,3) = F(pp,2)

        call EZspline_gradient(bfield_2d%Z, Y(pp,1), Y(pp,3), F(:,pp), ezerr)
        call EZspline_error(ezerr)

        BZ(pp,1) = F(pp,1)
        BZ(pp,2) = 0._rp
        BZ(pp,3) = F(pp,2)

     end if
  end do
  !$OMP END PARALLEL DO
  DEALLOCATE(F)
end subroutine gradient_2D_Bfields

subroutine interp_2D_gradBfields(Y,gradB,flag)
  !! @note Subroutine for interpolating the pre-computed, axisymmetric
  !! gradient of the magnitude of themagnetic field to the particles'
  !! position. Stored as cylindrical components of field. @endnote
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
  !! Particles' position in cylindrical coordinates, Y(1,:) = \(R\),
  !! Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: gradB
  !! Cylindirical components of interpolated gradient of magnitude of
  !! magnetic field.
  REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: F
  !! Cylindrical components of interpolated magnetic field components.
  !! F(1,:)=\(B_R\), F(2,:)=\(B_\phi\), and F(3,:)=\(B_Z\).
  INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
  !! Flag that indicates whether particles are followed in the simulation
  !! (flag=1), or not (flag=0).
  INTEGER                                                :: pp
  !! Particle iterator.
  INTEGER                                                :: ss
  !! Species iterator.

  ss = size(Y,1)

  ALLOCATE(F(ss,3))
  !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
  !$OMP& SHARED(F,Y,gradB,flag,gradB_2d)
  do pp=1_idef,ss
     if ( flag(pp) .EQ. 1_is ) then
        call EZspline_interp(gradB_2d%R, Y(pp,1), Y(pp,3), F(pp,1), ezerr)
        call EZspline_error(ezerr)

        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        end if

        call EZspline_interp(gradB_2d%PHI, Y(pp,1), Y(pp,3), F(pp,2), ezerr)
        call EZspline_error(ezerr)

        call EZspline_interp(gradB_2d%Z, Y(pp,1), Y(pp,3), F(pp,3), ezerr)
        call EZspline_error(ezerr)

!        write(output_unit_write,'("PS R-gradB ",E17.10)') F(pp,1)
!        write(output_unit_write,'("PS PHI-gradB ",E17.10)') F(pp,2)
!        write(output_unit_write,'("PS Z-gradB ",E17.10)') F(pp,3)

        gradB(pp,1) = F(pp,1)
        gradB(pp,2) = F(pp,2)
        gradB(pp,3) = F(pp,3)

!        write(output_unit_write,'("PHI-gradB ",E17.10)') gradB(2,1)

     end if
  end do
  !$OMP END PARALLEL DO
  DEALLOCATE(F)
end subroutine interp_2D_gradBfields

subroutine interp_2D_curlbfields(Y,curlb,flag)
  !! @note Subroutine for interpolating the pre-computed, axisymmetric
  !! curl of the magnetic field unit vector to the particles'
  !! position. Stored as cylindrical components of field. @endnote
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
  !! Particles' position in cylindrical coordinates, Y(1,:) = \(R\),
  !! Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: curlb
  !! Cylindirical components of interpolated curl of direction of
  !! magnetic field.
  REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: F
  !! Cylindrical components of interpolated magnetic field components.
  !! F(1,:)=\(B_R\), F(2,:)=\(B_\phi\), and F(3,:)=\(B_Z\).
  INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
  !! Flag that indicates whether particles are followed in the simulation
  !! (flag=1), or not (flag=0).
  INTEGER                                                :: pp
  !! Particle iterator.
  INTEGER                                                :: ss
  !! Species iterator.

  ss = size(Y,1)

  ALLOCATE(F(ss,3))
  !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
  !$OMP& SHARED(F,Y,curlb,flag,curlb_2d)
  do pp=1_idef,ss
     if ( flag(pp) .EQ. 1_is ) then
        call EZspline_interp(curlb_2d%R, Y(pp,1), Y(pp,3), F(pp,1), ezerr)
        call EZspline_error(ezerr)

        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        end if

        call EZspline_interp(curlb_2d%PHI, Y(pp,1), Y(pp,3), F(pp,2), ezerr)
        call EZspline_error(ezerr)

        call EZspline_interp(curlb_2d%Z, Y(pp,1), Y(pp,3), F(pp,3), ezerr)
        call EZspline_error(ezerr)

        curlb(pp,1) = F(pp,1)
        curlb(pp,2) = F(pp,2)
        curlb(pp,3) = F(pp,3)
     end if
  end do
  !$OMP END PARALLEL DO
  DEALLOCATE(F)
end subroutine interp_2D_curlbfields

subroutine interp_FOfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_X,B_Y,B_Z, &
     E_X,E_Y,E_Z,PSIp,flag_cache)
  INTEGER, INTENT(IN)  :: pchunk
  TYPE(FIELDS), INTENT(IN)                               :: F
  REAL(rp),DIMENSION(pchunk),INTENT(IN)   :: Y_R,Y_PHI,Y_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: B_X,B_Y,B_Z
  REAL(rp),DIMENSION(pchunk)   :: B_R,B_PHI
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: E_X,E_Y,E_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: PSIp
  REAL(rp),DIMENSION(pchunk)   :: E_R,E_PHI
  REAL(rp),DIMENSION(pchunk)   :: cP,sP
  !  INTEGER(ip) :: ezerr
  INTEGER                                      :: cc
  !! Particle chunk iterator.
  INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flag_cache

  call check_if_in_fields_domain_p(pchunk,F,Y_R,Y_PHI,Y_Z,flag_cache)

  call EZspline_interp(bfield_2d%A,pchunk,Y_R, Y_Z,PSIp, ezerr)
  call EZspline_error(ezerr)

  call EZspline_interp(bfield_2d%R,bfield_2d%PHI,bfield_2d%Z,efield_2d%R, &
       efield_2d%PHI,efield_2d%Z,pchunk,Y_R,Y_Z,B_R,B_PHI,B_Z, &
       E_R,E_PHI,E_Z,ezerr)
  call EZspline_error(ezerr)

  !$OMP SIMD
!  !$OMP& aligned (cP,sP,B_X,B_Y,E_X,E_Y,Y_PHI,B_R,B_PHI,E_R,E_PHI)
  do cc=1_idef,pchunk
     cP(cc)=cos(Y_PHI(cc))
     sP(cc)=sin(Y_PHI(cc))

     B_X(cc) = B_R(cc)*cP(cc) - B_PHI(cc)*sP(cc)
     B_Y(cc) = B_R(cc)*sP(cc) + B_PHI(cc)*cP(cc)

     E_X(cc) = E_R(cc)*cP(cc) - E_PHI(cc)*sP(cc)
     E_Y(cc) = E_R(cc)*sP(cc) + E_PHI(cc)*cP(cc)

  end do
  !$OMP END SIMD

end subroutine interp_FOfields_p

subroutine interp_FO3Dfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_X,B_Y,B_Z, &
     E_X,E_Y,E_Z,flag_cache)
  INTEGER, INTENT(IN)  :: pchunk
  TYPE(FIELDS), INTENT(IN)                               :: F
  REAL(rp),DIMENSION(pchunk),INTENT(IN)   :: Y_R,Y_PHI,Y_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: B_X,B_Y,B_Z
  REAL(rp),DIMENSION(pchunk)   :: B_R,B_PHI
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: E_X,E_Y,E_Z
  REAL(rp),DIMENSION(pchunk)   :: E_R,E_PHI
  REAL(rp),DIMENSION(pchunk)   :: cP,sP
  !  INTEGER(ip) :: ezerr
  INTEGER                                      :: cc
  !! Particle chunk iterator.
  INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flag_cache

  call check_if_in_fields_domain_p(pchunk,F,Y_R,Y_PHI,Y_Z,flag_cache)

  call EZspline_interp(bfield_3d%R,bfield_3d%PHI,bfield_3d%Z,efield_3d%R, &
       efield_3d%PHI,efield_3d%Z,pchunk,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
       E_R,E_PHI,E_Z,ezerr)
  call EZspline_error(ezerr)

  !write(6,*) 'B_R',B_R
  !write(6,*) 'B_PHI',B_PHI
  !write(6,*) 'B_Z',B_Z

  !$OMP SIMD
!  !$OMP& aligned (cP,sP,B_X,B_Y,E_X,E_Y,Y_PHI,B_R,B_PHI,E_R,E_PHI)
  do cc=1_idef,pchunk
     cP(cc)=cos(Y_PHI(cc))
     sP(cc)=sin(Y_PHI(cc))



     B_X(cc) = B_R(cc)*cP(cc) - B_PHI(cc)*sP(cc)
     B_Y(cc) = B_R(cc)*sP(cc) + B_PHI(cc)*cP(cc)

     E_X(cc) = E_R(cc)*cP(cc) - E_PHI(cc)*sP(cc)
     E_Y(cc) = E_R(cc)*sP(cc) + E_PHI(cc)*cP(cc)

  end do
  !$OMP END SIMD

end subroutine interp_FO3Dfields_p

subroutine interp_FOfields1_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_X,B_Y,B_Z, &
     E_X,E_Y,E_Z,PSIp,flag_cache)
  INTEGER, INTENT(IN)  :: pchunk
  TYPE(FIELDS), INTENT(IN)                               :: F
  REAL(rp),DIMENSION(pchunk),INTENT(IN)   :: Y_R,Y_PHI,Y_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: B_X,B_Y,B_Z
  REAL(rp),DIMENSION(pchunk)   :: B_R,B_PHI
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: E_X,E_Y,E_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: PSIp
  REAL(rp),DIMENSION(pchunk)   :: E_R,E_PHI
  REAL(rp),DIMENSION(pchunk)   :: cP,sP
  !  INTEGER(ip) :: ezerr
  INTEGER                                      :: cc
  !! Particle chunk iterator.
  INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flag_cache

  call check_if_in_fields_domain_p(pchunk,F,Y_R,Y_PHI,Y_Z,flag_cache)

  call EZspline_interp(bfield_2d%A,pchunk,Y_R, Y_Z,PSIp, ezerr)
  call EZspline_error(ezerr)

  call calculate_magnetic_field_p(pchunk,F,Y_R,Y_Z,B_R,B_PHI,B_Z)

  call EZspline_interp(efield_2d%R,efield_2d%PHI,efield_2d%Z,pchunk,Y_R,Y_Z, &
       E_R,E_PHI,E_Z,ezerr)
  call EZspline_error(ezerr)

  !$OMP SIMD
!  !$OMP& aligned (cP,sP,B_X,B_Y,E_X,E_Y,Y_PHI,B_R,B_PHI,E_R,E_PHI)
  do cc=1_idef,pchunk
     cP(cc)=cos(Y_PHI(cc))
     sP(cc)=sin(Y_PHI(cc))

     B_X(cc) = B_R(cc)*cP(cc) - B_PHI(cc)*sP(cc)
     B_Y(cc) = B_R(cc)*sP(cc) + B_PHI(cc)*cP(cc)

     E_X(cc) = E_R(cc)*cP(cc) - E_PHI(cc)*sP(cc)
     E_Y(cc) = E_R(cc)*sP(cc) + E_PHI(cc)*cP(cc)

  end do
  !$OMP END SIMD

end subroutine interp_FOfields1_p

subroutine interp_FOfields_mars(prtcls, F, params)
  TYPE(KORC_PARAMS), INTENT(IN)                              :: params
  TYPE(PARTICLES), INTENT(INOUT) :: prtcls
  TYPE(FIELDS), INTENT(IN)       :: F
  REAL(rp),DIMENSION(1)   :: B_R,B_PHI,Y_R,Y_Z,Y_PHI
  REAL(rp),DIMENSION(1)   :: B0_R,B0_PHI,B0_Z
  REAL(rp),DIMENSION(1)   :: B1_R,B1_PHI,B1_Z
  REAL(rp),DIMENSION(1)   :: B1Re_R,B1Re_PHI,B1Re_Z
  REAL(rp),DIMENSION(1)   :: B1Im_R,B1Im_PHI,B1Im_Z
  REAL(rp),DIMENSION(1)   :: cP,sP,cPshift,sPshift
  REAL(rp), DIMENSION(1,3)  :: A
  !  INTEGER(ip) :: ezerr
  INTEGER                                      :: pp,ss
  !! Particle chunk iterator.
  REAL(rp) :: psip_conv
  REAL(rp) :: amp,phase



  if (size(prtcls%Y,1).eq.1) then
     ss = size(prtcls%Y,1)
  else
     if (prtcls%Y(2,1).eq.0) then
        ss=1_idef
     else
        ss = size(prtcls%Y,1)
     end if
  endif

  psip_conv=F%psip_conv
  amp=F%AMP
  phase=F%MARS_phase


  !$OMP PARALLEL DO &
  !$OMP& PRIVATE(pp,ezerr,A,B1Re_R,B1Re_PHI,B1Re_Z,B1Im_R,B1Im_PHI,B1Im_Z, &
  !$OMP& B0_R,B0_PHI,B0_Z,B1_R,B1_PHI,B1_Z,cP,sP,B_R,B_PHI,Y_R,Y_Z,Y_PHI) &
  !$OMP& FIRSTPRIVATE(psip_conv,amp) &
  !$OMP& SHARED(prtcls,F,params)
  do pp = 1,ss

     Y_R(1)=prtcls%Y(pp,1)
     Y_PHI(1)=prtcls%Y(pp,2)
     Y_Z(1)=prtcls%Y(pp,3)

     call EZspline_interp(bfield_2d%A,b1Refield_2d%R,b1Refield_2d%PHI, &
          b1Refield_2d%Z,b1Imfield_2d%R,b1Imfield_2d%PHI,b1Imfield_2d%Z, &
          1,Y_R,Y_Z,A,B1Re_R,B1Re_PHI,B1Re_Z, &
          B1Im_R,B1Im_PHI,B1Im_Z,ezerr)
     call EZspline_error(ezerr)

     prtcls%PSI_P=A(1,1)

     B0_R = psip_conv*A(1,3)/Y_R(1)
     B0_PHI = -F%Bo*F%Ro/Y_R(1)
     B0_Z = -psip_conv*A(1,2)/Y_R(1)

     cP=cos(Y_PHI(1))
     sP=sin(Y_PHI(1))
     cPshift=cos(Y_PHI(1)-phase)
     sPshift=sin(Y_PHI(1)-phase)

     B1_R = amp*(B1Re_R*cPshift-B1Im_R*sPshift)
     B1_PHI = amp*(B1Re_PHI*cPshift-B1Im_PHI*sPshift)
     B1_Z = amp*(B1Re_Z*cPshift-B1Im_Z*sPshift)

     B_R = B0_R+B1_R
     B_PHI = B0_PHI+B1_PHI
     prtcls%B(pp,3) = B0_Z(1)+B1_Z(1)

     prtcls%B(pp,1) = B_R(1)*cP(1) - B_PHI(1)*sP(1)
     prtcls%B(pp,2) = B_R(1)*sP(1) + B_PHI(1)*cP(1)

     !write(6,*) '(R,PHI,Z)',Y_R*params%cpp%length,Y_PHI,Y_Z*params%cpp%length
     !write(6,*) 'amp',amp,'cP,sP',cP,sP,'cnP,snP',cnP,snP
     !write(6,*) 'psi',PSIp*params%cpp%Bo*params%cpp%length**2
     !write(6,*) 'dpsidR',A(:,2)*params%cpp%Bo*params%cpp%length
     !write(6,*) 'dpsidZ',A(:,3)*params%cpp%Bo*params%cpp%length
     !write(6,*) 'B0',B0_R*params%cpp%Bo,B0_PHI*params%cpp%Bo,B0_Z*params%cpp%Bo
     !write(6,*) 'AMP',amp
     !write(6,*) 'B1Re',B1Re_R*params%cpp%Bo,B1Re_PHI*params%cpp%Bo,B1Re_Z*params%cpp%Bo
     !write(6,*) 'B1Im',B1Im_R*params%cpp%Bo,B1Im_PHI*params%cpp%Bo,B1Im_Z*params%cpp%Bo
     !write(6,*) 'B1',B1_R*params%cpp%Bo,B1_PHI*params%cpp%Bo,B1_Z*params%cpp%Bo
     !write(6,*) 'B',B_R*params%cpp%Bo,B_PHI*params%cpp%Bo,prtcls%B(pp,3)*params%cpp%Bo

  end do
  !$OMP END PARALLEL DO

end subroutine interp_FOfields_mars

subroutine interp_FOfields_mars_p(params,pchunk,F,Y_R,Y_PHI,Y_Z,B_X,B_Y,B_Z, &
     PSIp,flag_cache)
  TYPE(KORC_PARAMS), INTENT(IN)                              :: params
  INTEGER, INTENT(IN)  :: pchunk
  TYPE(FIELDS), INTENT(IN)                               :: F
  REAL(rp),DIMENSION(pchunk),INTENT(IN)   :: Y_R,Y_PHI,Y_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: B_X,B_Y,B_Z
  REAL(rp),DIMENSION(pchunk)   :: B_R,B_PHI
  REAL(rp),DIMENSION(pchunk)   :: B0_R,B0_PHI,B0_Z
  REAL(rp),DIMENSION(pchunk)   :: B1_R,B1_PHI,B1_Z
  REAL(rp),DIMENSION(pchunk)   :: B1Re_R,B1Re_PHI,B1Re_Z
  REAL(rp),DIMENSION(pchunk)   :: B1Im_R,B1Im_PHI,B1Im_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: PSIp
  REAL(rp),DIMENSION(pchunk)   :: cP,sP,cPshift,sPshift
  REAL(rp), DIMENSION(pchunk,3)  :: A
  !  INTEGER(ip) :: ezerr
  INTEGER                                      :: cc
  !! Particle chunk iterator.
  INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flag_cache
  REAL(rp) :: psip_conv
  REAL(rp) :: amp,phase

  psip_conv=F%psip_conv
  amp=F%AMP
  phase=F%MARS_phase

  call check_if_in_fields_domain_p(pchunk,F,Y_R,Y_PHI,Y_Z,flag_cache)

  call EZspline_interp(bfield_2d%A,b1Refield_2d%R,b1Refield_2d%PHI, &
       b1Refield_2d%Z,b1Imfield_2d%R,b1Imfield_2d%PHI,b1Imfield_2d%Z, &
       pchunk,Y_R,Y_Z,A,B1Re_R,B1Re_PHI,B1Re_Z,B1Im_R,B1Im_PHI,B1Im_Z,ezerr)
  call EZspline_error(ezerr)

  !$OMP SIMD
  do cc=1_idef,pchunk
     PSIp(cc)=A(cc,1)

     B0_R(cc) = psip_conv*A(cc,3)/Y_R(cc)
     B0_PHI(cc) = -F%Bo*F%Ro/Y_R(cc)
     B0_Z(cc) = -psip_conv*A(cc,2)/Y_R(cc)

     cP(cc)=cos(Y_PHI(cc))
     sP(cc)=sin(Y_PHI(cc))
     cPshift(cc)=cos(Y_PHI(cc)-phase)
     sPshift(cc)=sin(Y_PHI(cc)-phase)

     B1_R(cc) = amp*(B1Re_R(cc)*cPshift(cc)-B1Im_R(cc)*sPshift(cc))
     B1_PHI(cc) = amp*(B1Re_PHI(cc)*cPshift(cc)-B1Im_PHI(cc)*sPshift(cc))
     B1_Z(cc) = amp*(B1Re_Z(cc)*cPshift(cc)-B1Im_Z(cc)*sPshift(cc))

     B_R(cc) = B0_R(cc)+B1_R(cc)
     B_PHI(cc) = B0_PHI(cc)+B1_PHI(cc)
     B_Z(cc) = B0_Z(cc)+B1_Z(cc)

     B_X(cc) = B_R(cc)*cP(cc) - B_PHI(cc)*sP(cc)
     B_Y(cc) = B_R(cc)*sP(cc) + B_PHI(cc)*cP(cc)

  end do
  !$OMP END SIMD

  !write(6,*) '(R,PHI,Z)',Y_R*params%cpp%length,Y_PHI*params%cpp%length, &
  !     Y_Z*params%cpp%length
  !write(6,*) 'psi',PSIp*params%cpp%Bo*params%cpp%length**2
  !write(6,*) 'dpsidR',A(:,2)*params%cpp%Bo*params%cpp%length
  !write(6,*) 'dpsidZ',A(:,3)*params%cpp%Bo*params%cpp%length
  !write(6,*) 'B0',B0_R*params%cpp%Bo,B0_PHI*params%cpp%Bo,B0_Z*params%cpp%Bo
  !write(6,*) 'AMP',amp
  !write(6,*) 'B1Re',B1Re_R*params%cpp%Bo,B1Re_PHI*params%cpp%Bo,B1Re_Z*params%cpp%Bo
  !write(6,*) 'B1Im',B1Im_R*params%cpp%Bo,B1Im_PHI*params%cpp%Bo,B1Im_Z*params%cpp%Bo
  !write(6,*) 'B1',B1_R*params%cpp%Bo,B1_PHI*params%cpp%Bo,B1_Z*params%cpp%Bo
  !write(6,*) 'B',B_X*params%cpp%Bo,B_Y*params%cpp%Bo,B_Z*params%cpp%Bo


end subroutine interp_FOfields_mars_p

subroutine interp_FOfields_aorsa(prtcls, F, params)
  TYPE(KORC_PARAMS), INTENT(IN)                              :: params
  TYPE(PARTICLES), INTENT(INOUT) :: prtcls
  TYPE(FIELDS), INTENT(IN)       :: F
  REAL(rp),DIMENSION(1)   :: Y_R,Y_PHI,Y_Z
  REAL(rp),DIMENSION(1)   :: B0_R,B0_PHI,B0_Z
  REAL(rp),DIMENSION(1)   :: B0_X,B0_Y
  REAL(rp),DIMENSION(1)   :: B1_X,B1_Y,B1_Z
  REAL(rp),DIMENSION(1)   :: B1Re_X,B1Re_Y,B1Re_Z
  REAL(rp),DIMENSION(1)   :: B1Im_X,B1Im_Y,B1Im_Z
  REAL(rp),DIMENSION(1)   :: E1Re_X,E1Re_Y,E1Re_Z
  REAL(rp),DIMENSION(1)   :: E1Im_X,E1Im_Y,E1Im_Z
  REAL(rp),DIMENSION(1)   :: cP,sP,cnP,snP
  REAL(rp), DIMENSION(1,3)  :: A
  !  INTEGER(ip) :: ezerr
  INTEGER                                      :: pp,ss
  !! Particle chunk iterator.
  REAL(rp) :: psip_conv
  REAL(rp) :: amp,nmode

  if (size(prtcls%Y,1).eq.1) then
     ss = size(prtcls%Y,1)
  else
     if (prtcls%Y(2,1).eq.0) then
        ss=1_idef
     else
        ss = size(prtcls%Y,1)
     end if
  endif

  psip_conv=F%psip_conv
  amp=F%AMP
  nmode=F%AORSA_nmode


  !$OMP PARALLEL DO &
  !$OMP& PRIVATE(pp,ezerr,A,B1Re_X,B1Re_Y,B1Re_Z,B1Im_X,B1Im_Y,B1Im_Z, &
  !$OMP& E1Re_X,E1Re_Y,E1Re_Z,E1Im_X,E1Im_Y,E1Im_Z, &
  !$OMP& B0_R,B0_PHI,B0_X,B0_Y,B0_Z,B1_X,B1_Y,B1_Z,cP,sP,cnP,snP,Y_R,Y_Z,Y_PHI)&
  !$OMP& FIRSTPRIVATE(psip_conv,amp,nmode) &
  !$OMP& SHARED(prtcls,F,params)
  do pp = 1,ss

     Y_R(1)=prtcls%Y(pp,1)
     Y_PHI(1)=prtcls%Y(pp,2)
     Y_Z(1)=prtcls%Y(pp,3)

     call EZspline_interp(bfield_2d%A,&
          b1Refield_2dx%X,b1Refield_2dx%Y,b1Refield_2dx%Z, &
          b1Imfield_2dx%X,b1Imfield_2dx%Y,b1Imfield_2dx%Z, &
          e1Refield_2dx%X,e1Refield_2dx%Y,e1Refield_2dx%Z, &
          e1Imfield_2dx%X,e1Imfield_2dx%Y,e1Imfield_2dx%Z, &
          1,Y_R,Y_Z,A,B1Re_X,B1Re_Y,B1Re_Z,B1Im_X,B1Im_Y,B1Im_Z, &
          E1Re_X,E1Re_Y,E1Re_Z,E1Im_X,E1Im_Y,E1Im_Z,ezerr)
     call EZspline_error(ezerr)

     prtcls%PSI_P=A(1,1)

     B0_R = psip_conv*A(1,3)/prtcls%Y(pp,1)
     B0_PHI = -F%Bo*F%Ro/prtcls%Y(pp,1)
     B0_Z = -psip_conv*A(1,2)/prtcls%Y(pp,1)

     cnP=cos(nmode*Y_PHI)
     snP=sin(nmode*Y_PHI)

     cP=cos(Y_PHI)
     sP=sin(Y_PHI)

     B1_X = amp*(B1Re_X*cnP-B1Im_X*snP)
     B1_Y = amp*(B1Re_Y*cnP-B1Im_Y*snP)
     B1_Z = amp*(B1Re_Z*cnP-B1Im_Z*snP)

     prtcls%E(pp,1) = amp*(E1Re_X(1)*cnP(1)-E1Im_X(1)*snP(1))
     prtcls%E(pp,2) = amp*(E1Re_Y(1)*cnP(1)-E1Im_Y(1)*snP(1))
     prtcls%E(pp,3) = amp*(E1Re_Z(1)*cnP(1)-E1Im_Z(1)*snP(1))

     B0_X = B0_R*cP - B0_PHI*sP
     B0_Y = B0_R*sP + B0_PHI*cP
     prtcls%B(pp,3) = B0_Z(1)+B1_Z(1)

     prtcls%B(pp,1) = B0_X(1)+B1_X(1)
     prtcls%B(pp,2) = B0_Y(1)+B1_Y(1)

     !write(6,*) '(R,PHI,Z)',Y_R*params%cpp%length,Y_PHI, &
     !     Y_Z*params%cpp%length
     !write(6,*) 'amp',amp,'cP,sP',cP,sP,'cnP,snP',cnP,snP
     !write(6,*) 'psi',PSIp*params%cpp%Bo*params%cpp%length**2
     !write(6,*) 'dpsidR',A(:,2)*params%cpp%Bo*params%cpp%length
     !write(6,*) 'dpsidZ',A(:,3)*params%cpp%Bo*params%cpp%length
     !write(6,*) 'B0',B0_R*params%cpp%Bo,B0_PHI*params%cpp%Bo,B0_Z*params%cpp%Bo
     !write(6,*) 'AMP',amp
     !write(6,*) 'B1Re',B1Re_X*params%cpp%Bo,B1Re_Y*params%cpp%Bo,B1Re_Z*params%cpp%Bo
     !write(6,*) 'B1Im',B1Im_X*params%cpp%Bo,B1Im_Y*params%cpp%Bo,B1Im_Z*params%cpp%Bo
     !write(6,*) 'B1',B1_X*params%cpp%Bo,B1_Y*params%cpp%Bo,B1_Z*params%cpp%Bo
     !write(6,*) 'B',B_X*params%cpp%Bo,B_Y*params%cpp%Bo,B_Z*params%cpp%Bo

  end do
  !$OMP END PARALLEL DO



end subroutine interp_FOfields_aorsa

subroutine interp_FOfields_aorsa_p(time,params,pchunk,F,Y_R,Y_PHI,Y_Z, &
     B_X,B_Y,B_Z,E_X,E_Y,E_Z,PSIp,flag_cache)
  TYPE(KORC_PARAMS), INTENT(IN)                              :: params
  REAL(rp), INTENT(IN)  :: time
  INTEGER, INTENT(IN)  :: pchunk
  TYPE(FIELDS), INTENT(IN)                               :: F
  REAL(rp),DIMENSION(pchunk),INTENT(IN)   :: Y_R,Y_PHI,Y_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: B_X,B_Y,B_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: E_X,E_Y,E_Z
  REAL(rp),DIMENSION(pchunk)   :: B0_R,B0_PHI,B0_Z,B0_X,B0_Y
  REAL(rp),DIMENSION(pchunk)   :: B1_X,B1_Y,B1_Z
  REAL(rp),DIMENSION(pchunk)   :: B1Re_X,B1Re_Y,B1Re_Z
  REAL(rp),DIMENSION(pchunk)   :: B1Im_X,B1Im_Y,B1Im_Z
  REAL(rp),DIMENSION(pchunk)   :: E1Re_X,E1Re_Y,E1Re_Z
  REAL(rp),DIMENSION(pchunk)   :: E1Im_X,E1Im_Y,E1Im_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: PSIp
  REAL(rp),DIMENSION(pchunk)   :: cP,sP,cnP,snP
  REAL(rp), DIMENSION(pchunk,3)  :: A
  !  INTEGER(ip) :: ezerr
  INTEGER                                      :: cc
  !! Particle chunk iterator.
  INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flag_cache
  REAL(rp) :: psip_conv
  REAL(rp) :: amp,nmode,omega

  psip_conv=F%psip_conv
  amp=F%AMP
  nmode=F%AORSA_nmode
  omega=2*C_PI*F%AORSA_freq

  call check_if_in_fields_domain_p(pchunk,F,Y_R,Y_PHI,Y_Z,flag_cache)

  call EZspline_interp(bfield_2d%A,b1Refield_2dx%X,b1Refield_2dx%Y, &
       b1Refield_2dx%Z,b1Imfield_2dx%X,b1Imfield_2dx%Y,b1Imfield_2dx%Z, &
       e1Refield_2dx%X,e1Refield_2dx%Y,e1Refield_2dx%Z, &
       e1Imfield_2dx%X,e1Imfield_2dx%Y,e1Imfield_2dx%Z, &
       pchunk,Y_R,Y_Z,A,B1Re_X,B1Re_Y,B1Re_Z,B1Im_X,B1Im_Y,B1Im_Z, &
       E1Re_X,E1Re_Y,E1Re_Z,E1Im_X,E1Im_Y,E1Im_Z,ezerr)
  call EZspline_error(ezerr)


  !$OMP SIMD
  do cc=1_idef,pchunk
     PSIp(cc)=A(cc,1)

     B0_R(cc) = psip_conv*A(cc,3)/Y_R(cc)
     B0_PHI(cc) = -F%Bo*F%Ro/Y_R(cc)
     B0_Z(cc) = -psip_conv*A(cc,2)/Y_R(cc)

     cP(cc)=cos(Y_PHI(cc))
     sP(cc)=sin(Y_PHI(cc))

     B0_X(cc) = B0_R(cc)*cP(cc) - B0_PHI(cc)*sP(cc)
     B0_Y(cc) = B0_R(cc)*sP(cc) + B0_PHI(cc)*cP(cc)


     cnP(cc)=cos(omega*time+nmode*Y_PHI(cc))
     snP(cc)=sin(omega*time+nmode*Y_PHI(cc))

     B1_X(cc) = amp*(B1Re_X(cc)*cnP(cc)-B1Im_X(cc)*snP(cc))
     B1_Y(cc) = amp*(B1Re_Y(cc)*cnP(cc)-B1Im_Y(cc)*snP(cc))
     B1_Z(cc) = amp*(B1Re_Z(cc)*cnP(cc)-B1Im_Z(cc)*snP(cc))

     B_X(cc) = B0_X(cc)+B1_X(cc)
     B_Y(cc) = B0_Y(cc)+B1_Y(cc)
     B_Z(cc) = B0_Z(cc)+B1_Z(cc)

     E_X(cc) = amp*(E1Re_X(cc)*cnP(cc)-E1Im_X(cc)*snP(cc))
     E_Y(cc) = amp*(E1Re_Y(cc)*cnP(cc)-E1Im_Y(cc)*snP(cc))
     E_Z(cc) = amp*(E1Re_Z(cc)*cnP(cc)-E1Im_Z(cc)*snP(cc))


  end do
  !$OMP END SIMD

#if DBG_CHECK
  !write(6,*) '(R,PHI,Z,time)',Y_R*params%cpp%length,Y_PHI, &
  !     Y_Z*params%cpp%length,time
  !write(6,*) 'psi',PSIp*params%cpp%Bo*params%cpp%length**2
  !write(6,*) 'dpsidR',A(:,2)*params%cpp%Bo*params%cpp%length
  !write(6,*) 'dpsidZ',A(:,3)*params%cpp%Bo*params%cpp%length
  !write(6,*) 'B0',B0_R*params%cpp%Bo,B0_PHI*params%cpp%Bo,B0_Z*params%cpp%Bo
  !write(6,*) 'AMP',amp
  !write(6,*) 'B1Re',B1Re_X*params%cpp%Bo,B1Re_Y*params%cpp%Bo,B1Re_Z*params%cpp%Bo
  !write(6,*) 'B1Im',B1Im_X*params%cpp%Bo,B1Im_Y*params%cpp%Bo,B1Im_Z*params%cpp%Bo
  !write(6,*) 'B1',B1_X*params%cpp%Bo,B1_Y*params%cpp%Bo,B1_Z*params%cpp%Bo
  !write(6,*) 'B',B_X*params%cpp%Bo,B_Y*params%cpp%Bo,B_Z*params%cpp%Bo
#endif

end subroutine interp_FOfields_aorsa_p

subroutine interp_FOcollision_p(pchunk,Y_R,Y_PHI,Y_Z,ne,Te,Zeff,flag_cache)
  INTEGER, INTENT(IN)  :: pchunk
  REAL(rp),DIMENSION(pchunk),INTENT(IN)   :: Y_R,Y_PHI,Y_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: ne,Te,Zeff
  INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flag_cache

  call check_if_in_profiles_domain_p(pchunk,Y_R,Y_PHI,Y_Z,flag_cache)
!  write(output_unit_write,'("YR: ",E17.10)') Y_R(1)
!  write(output_unit_write,'("YPHI: ",E17.10)') Y_PHI(1)
!  write(output_unit_write,'("YZ: ",E17.10)') Y_Z(1)

!  write(output_unit_write,'("Te_interp_R",E17.10)') profiles_2d%Te%x1
!  write(output_unit_write,'("Te_interp_Z",E17.10)') profiles_2d%Te%x2

  call EZspline_interp(profiles_2d%ne,profiles_2d%Te, &
       profiles_2d%Zeff,pchunk,Y_R,Y_Z,ne,Te,Zeff,ezerr)
  ! this will call PSPLINE routine EZspline_interp2_bmag_cloud as there
  ! is the same number of entries
  call EZspline_error(ezerr)


end subroutine interp_FOcollision_p


subroutine interp_fields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z, &
     curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,flag_cache)
  INTEGER, INTENT(IN)  :: pchunk
  TYPE(FIELDS), INTENT(IN)                               :: F
  REAL(rp),DIMENSION(pchunk),INTENT(IN)   :: Y_R,Y_PHI,Y_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: B_R,B_PHI,B_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: gradB_R,gradB_PHI,gradB_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: curlB_R,curlB_PHI,curlB_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: E_R,E_PHI,E_Z
  INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flag_cache

  !write(output_unit_write,*) Y_R,Y_Z,flag_cache

  call check_if_in_fields_domain_p(pchunk,F,Y_R,Y_PHI,Y_Z,flag_cache)

  call EZspline_interp(bfield_2d%R,bfield_2d%PHI,bfield_2d%Z,efield_2d%R, &
       efield_2d%PHI,efield_2d%Z,gradB_2d%R,gradB_2d%PHI,gradB_2d%Z, &
       curlb_2d%R,curlb_2d%PHI,curlb_2d%Z,pchunk,Y_R,Y_Z,B_R,B_PHI,B_Z, &
       E_R,E_PHI,E_Z,gradB_R,gradB_PHI,gradB_Z,curlb_R,curlb_PHI,curlb_Z, &
       ezerr)
  call EZspline_error(ezerr)

end subroutine interp_fields_p

subroutine interp_fields_3D_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z, &
     curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,flag_cache)
  INTEGER, INTENT(IN)  :: pchunk
  TYPE(FIELDS), INTENT(IN)                               :: F
  REAL(rp),DIMENSION(pchunk),INTENT(IN)   :: Y_R,Y_PHI,Y_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: B_R,B_PHI,B_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: gradB_R,gradB_PHI,gradB_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: curlB_R,curlB_PHI,curlB_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: E_R,E_PHI,E_Z
  INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flag_cache
  REAL(rp), DIMENSION(pchunk)     :: Y_PHI_mod

  Y_PHI_mod=modulo(Y_PHI,2._rp*C_PI)
!  write(output_unit_write,*) Y_PHI(1)
!  write(output_unit_write,*) Y_PHI_mod(1)

  call check_if_in_fields_domain_p(pchunk,F,Y_R,Y_PHI_mod,Y_Z,flag_cache)

  call EZspline_interp(bfield_3d%R,bfield_3d%PHI,bfield_3d%Z,efield_3d%R, &
       efield_3d%PHI,efield_3d%Z,gradB_3d%R,gradB_3d%PHI,gradB_3d%Z, &
       curlb_3d%R,curlb_3d%PHI,curlb_3d%Z,pchunk,Y_R,Y_PHI_mod,Y_Z,B_R,B_PHI,B_Z, &
       E_R,E_PHI,E_Z,gradB_R,gradB_PHI,gradB_Z,curlb_R,curlb_PHI,curlb_Z, &
       ezerr)
  call EZspline_error(ezerr)

end subroutine interp_fields_3D_p

subroutine interp_collision_p(pchunk,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z, &
     ne,Te,Zeff,flag_cache)
  INTEGER, INTENT(IN)  :: pchunk
  REAL(rp),DIMENSION(pchunk),INTENT(IN)   :: Y_R,Y_PHI,Y_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: B_R,B_PHI,B_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: E_R,E_PHI,E_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: ne,Te,Zeff
  INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flag_cache
!  INTEGER(ip) :: ezerr

  call check_if_in_profiles_domain_p(pchunk,Y_R,Y_PHI,Y_Z,flag_cache)

  call EZspline_interp(bfield_2d%R,bfield_2d%PHI,bfield_2d%Z,efield_2d%R, &
       efield_2d%PHI,efield_2d%Z,profiles_2d%ne,profiles_2d%Te, &
       profiles_2d%Zeff,pchunk,Y_R,Y_Z,B_R,B_PHI,B_Z, &
       E_R,E_PHI,E_Z,ne,Te,Zeff,ezerr)
  call EZspline_error(ezerr)


end subroutine interp_collision_p

subroutine interp_bmag_p(pchunk,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z)
  INTEGER, INTENT(IN)  :: pchunk
  REAL(rp),DIMENSION(pchunk),INTENT(IN)   :: Y_R,Y_PHI,Y_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: B_R,B_PHI,B_Z
!  INTEGER(ip) :: ezerr

  call EZspline_interp(bfield_2d%R,bfield_2d%PHI,bfield_2d%Z, &
       pchunk,Y_R,Y_Z,B_R,B_PHI,B_Z,ezerr)
  call EZspline_error(ezerr)


end subroutine interp_bmag_p

!> @brief Subroutine for interpolating the pre-computed, 3-D magnetic field to the particles' position.
!!
!! @param[in] Y Particles' position in cylindrical coordinates, Y(1,:) = \(R\), Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
!! @param[in,out] B Cartesian components of interpolated magnetic field components. B(1,:)=\(B_x\), B(2,:)=\(B_y\), and B(3,:)=\(B_z\).
!! @param F Cylindrical components of interpolated magnetic field components. F(1,:)=\(B_R\), F(2,:)=\(B_\phi\), and F(3,:)=\(B_Z\).
!! @param flag Flag that indicates whether particles are followed in the simulation (flag=1), or not (flag=0).
!! @param pp Particle iterator.
!! @param ss Species iterator.
subroutine interp_3D_bfields(params,Y,B,flag)
  TYPE(KORC_PARAMS), INTENT(IN)      :: params
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: B
  REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: F
  INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
  INTEGER                                                :: pp
  INTEGER                                                :: ss

  ss = size(Y,1)

  ALLOCATE(F(ss,3))
  !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
  !$OMP& SHARED(F,Y,B,flag,bfield_3d)
  do pp=1_idef,ss
     if ( flag(pp) .EQ. 1_is ) then
        call EZspline_interp(bfield_3d%R, Y(pp,1), Y(pp,2), Y(pp,3), &
             F(pp,1), ezerr)
        call EZspline_error(ezerr)

        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        end if

        call EZspline_interp(bfield_3d%PHI, Y(pp,1), Y(pp,2), Y(pp,3), &
             F(pp,2), ezerr)
        call EZspline_error(ezerr)

        call EZspline_interp(bfield_3d%Z, Y(pp,1), Y(pp,2), Y(pp,3), &
             F(pp,3), ezerr)
        call EZspline_error(ezerr)

        if (.not.params%GC_coords) then
           B(pp,1) = F(pp,1)*COS(Y(pp,2)) - F(pp,2)*SIN(Y(pp,2))
           B(pp,2) = F(pp,1)*SIN(Y(pp,2)) + F(pp,2)*COS(Y(pp,2))
           B(pp,3) = F(pp,3)
        else
           B(pp,1) = F(pp,1)
           B(pp,2) = F(pp,2)
           B(pp,3) = F(pp,3)
        end if


     end if
  end do
  !$OMP END PARALLEL DO
  DEALLOCATE(F)
end subroutine interp_3D_bfields


!> @brief Subroutine that calculates the axisymmetric magnetic field to the particles' position using the poloidal magnetic flux.
!! @details When the poloidal magnetic flux \(\Psi(R,Z)\) is used in a KORC simulation, the magnetic field components are calculated as it follows:
!!
!!
!! $$B_R = \frac{1}{R}\frac{\partial \Psi}{\partial Z},$$
!! $$B_\phi = \frac{RoBo}{R}\),$$
!! $$B_Z = -\frac{1}{R}\frac{\partial \Psi}{\partial R},$$
!!
!!
!! where \(Ro\) and \(Bo\) are the radial position of the magnetic axis and the magnetic field as measured at the magnetic axis, respectively.
!! First, the derivatives of the poloidal magnetic flux are calculated at the particles' position using the PSPLINE interpolant of
!! the poloidal magnetic flux. Then, we calculate the cylindrical components of the magnetic field, and finally we calculate its Cartesian
!! components that will be used in the particle pusher.
!!
!! @param[in] Y Particles' position in cylindrical coordinates, Y(1,:) = \(R\), Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
!! @param[in] F An instance of KORC's derived type FIELDS containing all the information about the fields used in the simulation.
!! See korc_types.f90 and korc_fields.f90.
!! @param[in,out] B Cartesian components of interpolated magnetic field components. B(1,:)=\(B_x\), B(2,:)=\(B_y\), and B(3,:)=\(B_x\).
!! @param flag Flag that indicates whether particles are followed in the simulation (flag=1), or not (flag=0).
!! @param A Variable containing the partial derivatives of the poloidal magnetic flux \(\Psi(R,Z)\) and the cylindrical components
!! of the magnetic field (its value changes through the subroutine).
!! @param pp Particle iterator.
!! @param ss Species iterator.
subroutine calculate_magnetic_field(params,Y,F,B,E,PSI_P,flag)
  TYPE(KORC_PARAMS), INTENT(IN)      :: params
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
  TYPE(FIELDS), INTENT(IN)                               :: F
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: B
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: E
  REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)   :: PSI_P
  INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
  REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: A
  INTEGER                                                :: pp
  INTEGER                                                :: ss
  REAL(rp)  :: psip_conv

  psip_conv=F%psip_conv

  if (size(Y,1).eq.1) then
     ss = size(Y,1)
  else
     if (Y(2,1).eq.0) then
        ss=1_idef
     else
        ss = size(Y,1)
     end if
  endif

  ALLOCATE(A(ss,3))
  A=0._rp

  if(F%Dim2x1t.and.(.not.F%ReInterp_2x1t)) then
     !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
     !$OMP& SHARED(F,Y,A,B,flag,bfield_2X1T,PSI_P)
     do pp=1_idef,ss


        !        write(output_unit_write,'("pp: ",I16)') pp

        !        write(output_unit_write,'("Y_R: ",E17.10)') Y(:,1)
        !      write(output_unit_write,'("Y_PHI: ",E17.10)') Y(:,2)
        !      write(output_unit_write,'("Y_Z: ",E17.10)') Y(:,3)


        call EZspline_interp(bfield_2X1T%A, Y(pp,1), F%t0_2x1t, Y(pp,3), &
             PSI_P(pp), ezerr)
        call EZspline_error(ezerr)

        !write(output_unit_write,'("PSI_P: ",E17.10)') PSI_P(1)

        ! FR = (dA/dZ)/R
        call EZspline_derivative(bfield_2X1T%A, 0, 0, 1, Y(pp,1), F%t0_2x1t, &
             Y(pp,3), A(pp,1), ezerr)
        !			call EZspline_error(ezerr)

        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        else

           !           write(output_unit_write,'("R*B_R: ",E17.10)') A(pp,1)

           !if(params%SC_E) A(pp,1)=A(pp,1)/(2*C_PI)

           A(pp,1) = psip_conv*A(pp,1)/Y(pp,1)

           ! FPHI = Fo*Ro/R
           A(pp,2) = -F%Bo*F%Ro/Y(pp,1)

           ! FR = -(dA/dR)/R
           call EZspline_derivative(bfield_2X1T%A, 1, 0, 0, Y(pp,1), &
                F%t0_2x1t, Y(pp,3), A(pp,3), ezerr)
           call EZspline_error(ezerr)

           !           write(output_unit_write,'("R*B_Z: ",E17.10)') A(pp,3)

           !if(params%SC_E) A(pp,3)=A(pp,3)/(2*C_PI)

           call EZspline_derivative(bfield_2X1T%A, 0, 1, 0, Y(pp,1), &
                F%t0_2x1t, Y(pp,3), E(pp,2), ezerr)

           E(pp,2) = E(pp,2)/(2*C_PI*Y(pp,1))

           A(pp,3) = -psip_conv*A(pp,3)/Y(pp,1)

           if (.not.params%GC_coords) then
              B(pp,1) = A(pp,1)*COS(Y(pp,2)) - A(pp,2)*SIN(Y(pp,2))
              B(pp,2) = A(pp,1)*SIN(Y(pp,2)) + A(pp,2)*COS(Y(pp,2))
              B(pp,3) = A(pp,3)

              E(pp,1) = -E(pp,2)*sin(Y(pp,2))
              E(pp,2) = E(pp,2)*cos(Y(pp,2))
           else
              B(pp,1) = A(pp,1)
              B(pp,2) = A(pp,2)
              B(pp,3) = A(pp,3)
           end if


        end if
     end do
     !$OMP END PARALLEL DO


  else
     !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
     !$OMP& SHARED(F,Y,A,B,flag,bfield_2d,PSI_P)
     do pp=1_idef,ss


        !        write(output_unit_write,'("pp: ",I16)') pp

        !        write(output_unit_write,'("Y_R: ",E17.10)') Y(:,1)
        !      write(output_unit_write,'("Y_PHI: ",E17.10)') Y(:,2)
        !      write(output_unit_write,'("Y_Z: ",E17.10)') Y(:,3)


        call EZspline_interp(bfield_2d%A, Y(pp,1), Y(pp,3), &
             PSI_P(pp), ezerr)
        call EZspline_error(ezerr)

        !        write(output_unit_write,'("PSI_P: ",E17.10)') PSI_P(1)

        ! FR = (dA/dZ)/R
        call EZspline_derivative(bfield_2d%A, 0, 1, Y(pp,1), Y(pp,3), &
             A(pp,1), ezerr)
        !			call EZspline_error(ezerr)

        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        else

           !           write(output_unit_write,'("R*B_R: ",E17.10)') A(pp,1)

           !if(params%SC_E) A(pp,1)=A(pp,1)/(2*C_PI)

           A(pp,1) = psip_conv*A(pp,1)/Y(pp,1)

           ! FPHI = Fo*Ro/R
           A(pp,2) = -F%Bo*F%Ro/Y(pp,1)

           ! FR = -(dA/dR)/R
           call EZspline_derivative(bfield_2d%A, 1, 0, Y(pp,1), Y(pp,3), &
                A(pp,3), ezerr)
           call EZspline_error(ezerr)

           !           write(output_unit_write,'("R*B_Z: ",E17.10)') A(pp,3)

           !if(params%SC_E) A(pp,3)=A(pp,3)/(2*C_PI)

           A(pp,3) = -psip_conv*A(pp,3)/Y(pp,1)

           if (.not.params%GC_coords) then
              B(pp,1) = A(pp,1)*COS(Y(pp,2)) - A(pp,2)*SIN(Y(pp,2))
              B(pp,2) = A(pp,1)*SIN(Y(pp,2)) + A(pp,2)*COS(Y(pp,2))
              B(pp,3) = A(pp,3)
           else
              B(pp,1) = A(pp,1)
              B(pp,2) = A(pp,2)
              B(pp,3) = A(pp,3)
           end if


        end if
     end do
     !$OMP END PARALLEL DO
  end if

!  write(output_unit_write,'("calculate_fields")')

!  write(output_unit_write,'("B_R: ",E17.10)') A(:,1)
!  write(output_unit_write,'("B_PHI: ",E17.10)') A(:,2)
!  write(output_unit_write,'("B_Z: ",E17.10)') A(:,3)

!  write(output_unit_write,'("B_X: ",E17.10)') B(:,1)
!  write(output_unit_write,'("B_Y: ",E17.10)') B(:,2)
!  write(output_unit_write,'("B_Z: ",E17.10)') B(:,3)

  DEALLOCATE(A)
end subroutine calculate_magnetic_field


subroutine calculate_magnetic_field_p(pchunk,F,Y_R,Y_Z,B_R,B_PHI,B_Z)
  INTEGER, INTENT(IN)  :: pchunk
  REAL(rp), DIMENSION(pchunk), INTENT(IN)      :: Y_R,Y_Z
  TYPE(FIELDS), INTENT(IN)                               :: F
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: B_R,B_PHI,B_Z
  INTEGER                                                :: pp
  REAL(rp), DIMENSION(pchunk)  :: PSIp
  REAL(rp), DIMENSION(pchunk,2)  :: A
  REAL(rp)  :: psip_conv

  psip_conv=F%psip_conv

  call EZspline_interp(bfield_2d%A, pchunk, Y_R, Y_Z, &
       PSIp, ezerr)
  call EZspline_error(ezerr)

  ! FR = (dA/dZ)/R
  call EZspline_gradient(bfield_2d%A, pchunk, Y_R, Y_Z, &
       A, ezerr)
  call EZspline_error(ezerr)

  !write(output_unit_write,'("dPSIp/dR: ",E17.10)') A(:,1)
  !write(output_unit_write,'("dPSIp/dZ: ",E17.10)') A(:,2)
  !write(output_unit_write,'("Y_R: ",E17.10)') Y_R

  B_R = psip_conv*A(:,2)/Y_R

  ! FPHI = Fo*Ro/R
  B_PHI = -F%Bo*F%Ro/Y_R

  ! FR = -(dA/dR)/R

  !     write(output_unit_write,'("R*B_Z: ",E17.10)') B_Z(1)

  B_Z= -psip_conv*A(:,1)/Y_R


!  write(output_unit_write,'("PSIp: ",E17.10)') PSIp

!  write(output_unit_write,'("Y_R: ",E17.10)') Y_R(1)
!  write(output_unit_write,'("Y_Z: ",E17.10)') Y_Z(1)

!  write(output_unit_write,'("B_R: ",E17.10)') B_R
!  write(output_unit_write,'("B_PHI: ",E17.10)') B_PHI
!  write(output_unit_write,'("B_Z: ",E17.10)') B_Z

end subroutine calculate_magnetic_field_p

subroutine calculate_2DBdBfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
     E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
     flag_cache,PSIp)
  INTEGER, INTENT(IN)  :: pchunk
  REAL(rp), DIMENSION(pchunk), INTENT(IN)      :: Y_R,Y_PHI,Y_Z
  TYPE(FIELDS), INTENT(IN)                               :: F
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: B_R,B_PHI,B_Z
  REAL(rp), DIMENSION(pchunk,3)   :: BR,BPHI,BZ
  REAL(rp), DIMENSION(pchunk)   :: dBRdR,dBPHIdR,dBZdR
  REAL(rp), DIMENSION(pchunk)   :: dBRdPHI,dBPHIdPHI,dBZdPHI
  REAL(rp), DIMENSION(pchunk)   :: dBRdZ,dBPHIdZ,dBZdZ
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: gradB_R,gradB_PHI,gradB_Z
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: curlb_R,curlb_PHI,curlb_Z
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: E_R,E_PHI,E_Z
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: PSIp
  REAL(rp), DIMENSION(pchunk)   :: Bmag
  INTEGER                                                :: cc
  REAL(rp), DIMENSION(pchunk,6)  :: A
  INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flag_cache

  call check_if_in_fields_domain_p(pchunk,F,Y_R,Y_PHI,Y_Z,flag_cache)

  call EZspline_interp(bfield_2d%R,bfield_2d%PHI,bfield_2d%Z, &
       efield_2d%R,efield_2d%PHI,efield_2d%Z,bfield_2d%A, &
       pchunk,Y_R,Y_Z,BR,BPHI,BZ,E_R,E_PHI,E_Z,PSIp,ezerr)
  call EZspline_error(ezerr)


!  call EZspline_interp(bfield_2d%R,bfield_2d%PHI,bfield_2d%Z, &
!       dbdR_2d%R,dbdR_2d%PHI,dBdR_2d%Z, &
!       dbdPHI_2d%R,dbdPHI_2d%PHI,dbdPHI_2d%Z, &
!       dbdZ_2d%R,dbdZ_2d%PHI,dbdZ_2d%Z, &
!       efield_2d%R,efield_2d%PHI,efield_2d%Z,8,Y_R,Y_Z,B_R,B_PHI,B_Z, &
!       dBRdR,dBPHIdR,dBZdR,dBRdPHI,dBPHIdPHI,dBZdPHI,dBRdZ,dBPHIdZ,dBZdZ, &
!       E_R,E_PHI,E_Z,ezerr)
!  call EZspline_error(ezerr)

  !$OMP SIMD
!    !$OMP& aligned(PSIp,A,B_R,Y_R,B_PHI,B_Z,Bmag,gradB_R,gradB_PHI,gradB_Z, &
!    !$OMP& curlb_R,curlb_PHI,curlb_Z,E_R,E_PHI,E_Z)
  do cc=1_idef,pchunk

     B_R(cc)=BR(cc,1)
     B_PHI(cc)=BPHI(cc,1)
     B_Z(cc)=BZ(cc,1)

     dBRdR(cc)=BR(cc,2)
     dBRdPHI(cc)=0._rp
     dBRdZ(cc)=BR(cc,3)

     dBPHIdR(cc)=BPHI(cc,2)
     dBPHIdPHI(cc)=0._rp
     dBPHIdZ(cc)=BPHI(cc,3)

     dBZdR(cc)=BZ(cc,2)
     dBZdPHI(cc)=0._rp
     dBZdZ(cc)=BZ(cc,3)

     Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

     gradB_R(cc)=(B_R(cc)*dBRdR(cc)+B_PHI(cc)*dBPHIdR(cc)+ &
          B_Z(cc)*dBZdR(cc))/Bmag(cc)
     gradB_PHI(cc)=(B_R(cc)*dBRdPHI(cc)+B_PHI(cc)*dBPHIdPHI(cc)+ &
          B_Z(cc)*dBZdPHI(cc))/(Y_R(cc)*Bmag(cc))
     gradB_Z(cc)=(B_R(cc)*dBRdZ(cc)+B_PHI(cc)*dBPHIdZ(cc)+ &
          B_Z(cc)*dBZdZ(cc))/Bmag(cc)

     curlb_R(cc)=(Bmag(cc)*dBZdPHI(cc)/Y_R(cc)-B_Z(cc)*gradB_PHI(cc)- &
          Bmag(cc)*dBPHIdZ(cc)+B_PHI(cc)*gradB_Z(cc))/(Bmag(cc)*bmag(cc))
     curlb_PHI(cc)=(Bmag(cc)*dBRdZ(cc)-B_R(cc)*gradB_Z(cc)- &
          Bmag(cc)*dBZdR(cc)+B_Z(cc)*gradB_R(cc))/(Bmag(cc)*bmag(cc))
     curlb_Z(cc)=(Bmag(cc)*B_PHI(cc)/Y_R(cc)+Bmag(cc)*dBPHIdR(cc)- &
          B_PHI(cc)*gradB_R(cc)-Bmag(cc)*dBRdPHI(cc)/Y_R(cc)+ &
          B_R(cc)*gradB_PHI(cc))/(Bmag(cc)*bmag(cc))

  end do
  !$OMP END SIMD


!  write(output_unit_write,'("PSIp: ",E17.10)') PSIp

!  write(output_unit_write,'("Y_R: ",E17.10)') Y_R
!  write(output_unit_write,'("Y_Z: ",E17.10)') Y_Z

!  write(output_unit_write,'("B_R: ",E17.10)') B_R
!  write(output_unit_write,'("B_PHIinterp: ",E17.10)') B_PHI
!  write(output_unit_write,'("B_Z: ",E17.10)') B_Z

end subroutine calculate_2DBdBfields_p

subroutine calculate_3DBdBfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
     E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
     flag_cache)
  INTEGER, INTENT(IN)  :: pchunk
  REAL(rp), DIMENSION(pchunk), INTENT(IN)      :: Y_R,Y_PHI,Y_Z
  real(rp), DIMENSION(pchunk) :: Y_PHI_mod
  TYPE(FIELDS), INTENT(IN)                               :: F
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: B_R,B_PHI,B_Z
  REAL(rp), DIMENSION(pchunk)   :: dBRdR,dBPHIdR,dBZdR
  REAL(rp), DIMENSION(pchunk)   :: dBRdPHI,dBPHIdPHI,dBZdPHI
  REAL(rp), DIMENSION(pchunk)   :: dBRdZ,dBPHIdZ,dBZdZ
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: gradB_R,gradB_PHI,gradB_Z
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: curlb_R,curlb_PHI,curlb_Z
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: E_R,E_PHI,E_Z
  REAL(rp), DIMENSION(pchunk)   :: Bmag
  INTEGER                                                :: cc
  INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flag_cache

  Y_PHI_mod=modulo(Y_PHI,2._rp*C_PI)

  call check_if_in_fields_domain_p(pchunk,F,Y_R,Y_PHI_mod,Y_Z,flag_cache)


  call EZspline_interp(bfield_2d%R,bfield_2d%PHI,bfield_2d%Z, &
       dbdR_2d%R,dbdR_2d%PHI,dBdR_2d%Z, &
       dbdPHI_2d%R,dbdPHI_2d%PHI,dbdPHI_2d%Z, &
       dbdZ_2d%R,dbdZ_2d%PHI,dbdZ_2d%Z, &
       efield_2d%R,efield_2d%PHI,efield_2d%Z,pchunk,Y_R,Y_Z,B_R,B_PHI,B_Z, &
       dBRdR,dBPHIdR,dBZdR,dBRdPHI,dBPHIdPHI,dBZdPHI,dBRdZ,dBPHIdZ,dBZdZ, &
       E_R,E_PHI,E_Z,ezerr)
  call EZspline_error(ezerr)

  !$OMP SIMD
!    !$OMP& aligned(PSIp,A,B_R,Y_R,B_PHI,B_Z,Bmag,gradB_R,gradB_PHI,gradB_Z, &
!    !$OMP& curlb_R,curlb_PHI,curlb_Z,E_R,E_PHI,E_Z)
  do cc=1_idef,pchunk

     Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

     gradB_R(cc)=(B_R(cc)*dBRdR(cc)+B_PHI(cc)*dBPHIdR(cc)+ &
          B_Z(cc)*dBZdR(cc))/Bmag(cc)
     gradB_PHI(cc)=(B_R(cc)*dBRdPHI(cc)+B_PHI(cc)*dBPHIdPHI(cc)+ &
          B_Z(cc)*dBZdPHI(cc))/(Y_R(cc)*Bmag(cc))
     gradB_Z(cc)=(B_R(cc)*dBRdZ(cc)+B_PHI(cc)*dBPHIdZ(cc)+ &
          B_Z(cc)*dBZdZ(cc))/Bmag(cc)

     curlb_R(cc)=(Bmag(cc)*dBZdPHI(cc)/Y_R(cc)-B_Z(cc)*gradB_PHI(cc)- &
          Bmag(cc)*dBPHIdZ(cc)+B_PHI(cc)*gradB_Z(cc))/(Bmag(cc)*Bmag(cc))
     curlb_PHI(cc)=(Bmag(cc)*dBRdZ(cc)-B_R(cc)*gradB_Z(cc)- &
          Bmag(cc)*dBZdR(cc)+B_Z(cc)*gradB_R(cc))/(Bmag(cc)*Bmag(cc))
     curlb_Z(cc)=(Bmag(cc)*B_PHI(cc)/Y_R(cc)+Bmag(cc)*dBPHIdR(cc)- &
          B_PHI(cc)*gradB_R(cc)-Bmag(cc)*dBRdPHI(cc)/Y_R(cc)+ &
          B_R(cc)*gradB_PHI(cc))/(Bmag(cc)*Bmag(cc))

  end do
  !$OMP END SIMD


!  write(output_unit_write,'("PSIp: ",E17.10)') PSIp

!  write(output_unit_write,'("Y_R: ",E17.10)') Y_R
!  write(output_unit_write,'("Y_Z: ",E17.10)') Y_Z

!  write(output_unit_write,'("B_R: ",E17.10)') B_R
!  write(output_unit_write,'("B_PHIinterp: ",E17.10)') B_PHI
!  write(output_unit_write,'("B_Z: ",E17.10)') B_Z

end subroutine calculate_3DBdBfields_p

subroutine calculate_3DBdBfields1_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
     E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
     flag_cache,PSIp)
  INTEGER, INTENT(IN)  :: pchunk
  REAL(rp), DIMENSION(pchunk), INTENT(IN)      :: Y_R,Y_PHI,Y_Z
  real(rp), DIMENSION(pchunk) :: Y_PHI_mod
  TYPE(FIELDS), INTENT(IN)                               :: F
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: B_R,B_PHI,B_Z
  REAL(rp), DIMENSION(pchunk,4)   :: BR,BPHI,BZ
  REAL(rp), DIMENSION(pchunk)   :: dBRdR,dBPHIdR,dBZdR
  REAL(rp), DIMENSION(pchunk)   :: dBRdPHI,dBPHIdPHI,dBZdPHI
  REAL(rp), DIMENSION(pchunk)   :: dBRdZ,dBPHIdZ,dBZdZ
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: gradB_R,gradB_PHI,gradB_Z
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: curlb_R,curlb_PHI,curlb_Z
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: E_R,E_PHI,E_Z
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: PSIp
  REAL(rp), DIMENSION(pchunk)   :: Bmag
  INTEGER                                                :: cc
  REAL(rp), DIMENSION(pchunk,6)  :: A
  INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flag_cache

  Y_PHI_mod=modulo(Y_PHI,2._rp*C_PI)

  call check_if_in_fields_domain_p(pchunk,F,Y_R,Y_PHI_mod,Y_Z,flag_cache)

  call EZspline_interp(bfield_3d%R,bfield_3d%PHI,bfield_3d%Z, &
       efield_3d%R,efield_3d%PHI,efield_3d%Z,bfield_3d%A, &
       pchunk,Y_R,Y_PHI_mod,Y_Z,BR,BPHI,BZ,E_R,E_PHI,E_Z,PSIp,ezerr)
  call EZspline_error(ezerr)


!  call EZspline_interp(bfield_2d%R,bfield_2d%PHI,bfield_2d%Z, &
!       dbdR_2d%R,dbdR_2d%PHI,dBdR_2d%Z, &
!       dbdPHI_2d%R,dbdPHI_2d%PHI,dbdPHI_2d%Z, &
!       dbdZ_2d%R,dbdZ_2d%PHI,dbdZ_2d%Z, &
!       efield_2d%R,efield_2d%PHI,efield_2d%Z,8,Y_R,Y_Z,B_R,B_PHI,B_Z, &
!       dBRdR,dBPHIdR,dBZdR,dBRdPHI,dBPHIdPHI,dBZdPHI,dBRdZ,dBPHIdZ,dBZdZ, &
!       E_R,E_PHI,E_Z,ezerr)
!  call EZspline_error(ezerr)

  !$OMP SIMD
!    !$OMP& aligned(PSIp,A,B_R,Y_R,B_PHI,B_Z,Bmag,gradB_R,gradB_PHI,gradB_Z, &
!    !$OMP& curlb_R,curlb_PHI,curlb_Z,E_R,E_PHI,E_Z)
  do cc=1_idef,pchunk

     B_R(cc)=BR(cc,1)
     B_PHI(cc)=BPHI(cc,1)
     B_Z(cc)=BZ(cc,1)

     dBRdR(cc)=BR(cc,2)
     dBRdPHI(cc)=BR(cc,3)
     dBRdZ(cc)=BR(cc,4)

     dBPHIdR(cc)=BPHI(cc,2)
     dBPHIdPHI(cc)=BPHI(cc,3)
     dBPHIdZ(cc)=BPHI(cc,4)

     dBZdR(cc)=BZ(cc,2)
     dBZdPHI(cc)=BZ(cc,3)
     dBZdZ(cc)=BZ(cc,4)

     Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

     gradB_R(cc)=(B_R(cc)*dBRdR(cc)+B_PHI(cc)*dBPHIdR(cc)+ &
          B_Z(cc)*dBZdR(cc))/Bmag(cc)
     gradB_PHI(cc)=(B_R(cc)*dBRdPHI(cc)+B_PHI(cc)*dBPHIdPHI(cc)+ &
          B_Z(cc)*dBZdPHI(cc))/(Y_R(cc)*Bmag(cc))
     gradB_Z(cc)=(B_R(cc)*dBRdZ(cc)+B_PHI(cc)*dBPHIdZ(cc)+ &
          B_Z(cc)*dBZdZ(cc))/Bmag(cc)

     curlb_R(cc)=(Bmag(cc)*dBZdPHI(cc)/Y_R(cc)-B_Z(cc)*gradB_PHI(cc)- &
          Bmag(cc)*dBPHIdZ(cc)+B_PHI(cc)*gradB_Z(cc))/(Bmag(cc)*Bmag(cc))
     curlb_PHI(cc)=(Bmag(cc)*dBRdZ(cc)-B_R(cc)*gradB_Z(cc)- &
          Bmag(cc)*dBZdR(cc)+B_Z(cc)*gradB_R(cc))/(Bmag(cc)*Bmag(cc))
     curlb_Z(cc)=(Bmag(cc)*B_PHI(cc)/Y_R(cc)+Bmag(cc)*dBPHIdR(cc)- &
          B_PHI(cc)*gradB_R(cc)-Bmag(cc)*dBRdPHI(cc)/Y_R(cc)+ &
          B_R(cc)*gradB_PHI(cc))/(Bmag(cc)*Bmag(cc))

  end do
  !$OMP END SIMD


!  write(output_unit_write,'("PSIp: ",E17.10)') PSIp

!  write(output_unit_write,'("Y_R: ",E17.10)') Y_R
!  write(output_unit_write,'("Y_Z: ",E17.10)') Y_Z

!  write(output_unit_write,'("B_R: ",E17.10)') B_R
!  write(output_unit_write,'("B_PHIinterp: ",E17.10)') B_PHI
!  write(output_unit_write,'("B_Z: ",E17.10)') B_Z

end subroutine calculate_3DBdBfields1_p

subroutine calculate_GCfieldswE_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z,E_R,E_PHI,E_Z, &
     curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,flag_cache,PSIp)
  INTEGER, INTENT(IN)  :: pchunk
  REAL(rp), DIMENSION(pchunk), INTENT(IN)      :: Y_R,Y_PHI,Y_Z
  TYPE(FIELDS), INTENT(IN)                               :: F
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: B_R,B_PHI,B_Z
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: gradB_R,gradB_PHI,gradB_Z
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: curlb_R,curlb_PHI,curlb_Z
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: E_R,E_PHI,E_Z
  REAL(rp), DIMENSION(pchunk)   :: Bmag,EPHI
  INTEGER                                                :: cc
  REAL(rp), DIMENSION(pchunk),INTENT(OUT)  :: PSIp
  REAL(rp), DIMENSION(pchunk,6)  :: A
  INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flag_cache
  REAL(rp) :: psip_conv

  psip_conv=F%psip_conv

  call check_if_in_fields_domain_p(pchunk,F,Y_R,Y_PHI,Y_Z,flag_cache)

  call EZspline_derivative(bfield_2d%A, efield_2d%PHI, pchunk, Y_R, Y_Z, A, &
       EPHI, ezerr)
  call EZspline_error(ezerr)

  !A(:,1) = PSIp
  !A(:,2) = dPSIp/dR
  !A(:,3) = dPSIp/dZ
  !A(:,4) = d^2PSIp/dR^2
  !A(:,5) = d^2PSIp/dZ^2
  !A(:,6) = d^2PSIp/dRdZ

  !$OMP SIMD
!    !$OMP& aligned(PSIp,A,B_R,Y_R,B_PHI,B_Z,Bmag,gradB_R,gradB_PHI,gradB_Z, &
!    !$OMP& curlb_R,curlb_PHI,curlb_Z,E_R,E_PHI,E_Z)
  do cc=1_idef,pchunk
     PSIp(cc)=A(cc,1)

     B_R(cc) = psip_conv*A(cc,3)/Y_R(cc)
     ! BR = (dA/dZ)/R
     B_PHI(cc) = -F%Bo*F%Ro/Y_R(cc)
     ! BPHI = Fo*Ro/R
     B_Z(cc) = -psip_conv*A(cc,2)/Y_R(cc)
     ! BR = -(dA/dR)/R



     Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

     gradB_R(cc)=(B_R(cc)*psip_conv*A(cc,6)-B_Z(cc)*psip_conv*A(cc,4)- &
          Bmag(cc)*Bmag(cc))/(Y_R(cc)*Bmag(cc))
     gradB_PHI(cc)=0._rp
     gradB_Z(cc)=(B_R(cc)*psip_conv*A(cc,5)-B_Z(cc)*psip_conv*A(cc,6))/ &
          (Y_R(cc)*Bmag(cc))

     curlb_R(cc)=B_PHI(cc)*gradB_Z(cc)/(Bmag(cc)*Bmag(cc))
     curlb_PHI(cc)=(Bmag(cc)/Y_R(cc)*(B_Z(cc)+psip_conv*A(cc,4)+ &
          psip_conv*A(cc,5))-B_R(cc)*gradB_Z(cc)+B_Z(cc)*gradB_R(cc))/ &
          (Bmag(cc)*Bmag(cc))
     curlb_Z(cc)=-B_PHI(cc)*gradB_R(cc)/(Bmag(cc)*Bmag(cc))

     if (F%E_2x1t) then
        E_R(cc) = 0._rp
        E_PHI(cc) = EPHI(cc)
        E_Z(cc) = 0._rp
     else
        E_R(cc) = 0._rp
        E_PHI(cc) = 0._rp
        E_Z(cc) = 0._rp
     end if


  end do

end subroutine calculate_GCfieldswE_p

subroutine calculate_GCfields_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
     E_R,E_PHI,E_Z, &
     curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z,flag_cache,PSIp)
  INTEGER, INTENT(IN)  :: pchunk
  REAL(rp), DIMENSION(pchunk), INTENT(IN)      :: Y_R,Y_PHI,Y_Z
  TYPE(FIELDS), INTENT(IN)                               :: F
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: B_R,B_PHI,B_Z
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: gradB_R,gradB_PHI,gradB_Z
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: curlb_R,curlb_PHI,curlb_Z
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: E_R,E_PHI,E_Z
  REAL(rp), DIMENSION(pchunk)   :: Bmag
  INTEGER                                                :: cc
  REAL(rp), DIMENSION(pchunk),INTENT(OUT)  :: PSIp
  REAL(rp), DIMENSION(pchunk,6)  :: A
  INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flag_cache
  REAL(rp) :: psip_conv

  psip_conv=F%psip_conv

  call check_if_in_fields_domain_p(pchunk,F,Y_R,Y_PHI,Y_Z,flag_cache)

  call EZspline_derivative(bfield_2d%A, pchunk, Y_R, Y_Z, A, ezerr)
  call EZspline_error(ezerr)

  !A(:,1) = PSIp
  !A(:,2) = dPSIp/dR
  !A(:,3) = dPSIp/dZ
  !A(:,4) = d^2PSIp/dR^2
  !A(:,5) = d^2PSIp/dZ^2
  !A(:,6) = d^2PSIp/dRdZ

  !$OMP SIMD
!    !$OMP& aligned(PSIp,A,B_R,Y_R,B_PHI,B_Z,Bmag,gradB_R,gradB_PHI,gradB_Z, &
!    !$OMP& curlb_R,curlb_PHI,curlb_Z,E_R,E_PHI,E_Z)
  do cc=1_idef,pchunk
     PSIp(cc)=A(cc,1)

     B_R(cc) = psip_conv*A(cc,3)/Y_R(cc)
     ! BR = (dA/dZ)/R
     B_PHI(cc) = -F%Bo*F%Ro/Y_R(cc)
     ! BPHI = Fo*Ro/R
     B_Z(cc) = -psip_conv*A(cc,2)/Y_R(cc)
     ! BR = -(dA/dR)/R



     Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

     gradB_R(cc)=(B_R(cc)*psip_conv*A(cc,6)-B_Z(cc)*psip_conv*A(cc,4)- &
          Bmag(cc)*Bmag(cc))/(Y_R(cc)*Bmag(cc))
     gradB_PHI(cc)=0._rp
     gradB_Z(cc)=(B_R(cc)*psip_conv*A(cc,5)-B_Z(cc)*psip_conv*A(cc,6))/ &
          (Y_R(cc)*Bmag(cc))

     curlb_R(cc)=B_PHI(cc)*gradB_Z(cc)/(Bmag(cc)*Bmag(cc))
     curlb_PHI(cc)=(Bmag(cc)/Y_R(cc)*(B_Z(cc)+psip_conv*A(cc,4)+ &
          psip_conv*A(cc,5))-B_R(cc)*gradB_Z(cc)+B_Z(cc)*gradB_R(cc))/ &
          (Bmag(cc)*Bmag(cc))
     curlb_Z(cc)=-B_PHI(cc)*gradB_R(cc)/(Bmag(cc)*Bmag(cc))

     E_R(cc) = 0._rp
     E_PHI(cc) = F%Eo*F%Ro/Y_R(cc)
     E_Z(cc) = 0._rp



  end do

end subroutine calculate_GCfields_p

subroutine calculate_GCfields_2x1t_p(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
     E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
     flag_cache,PSIp,time)
  INTEGER, INTENT(IN)  :: pchunk
  REAL(rp), DIMENSION(pchunk), INTENT(IN)      :: Y_R,Y_PHI,Y_Z
  TYPE(FIELDS), INTENT(IN)                               :: F
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: B_R,B_PHI,B_Z
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: gradB_R,gradB_PHI,gradB_Z
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: curlb_R,curlb_PHI,curlb_Z
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: E_R,E_PHI,E_Z
  REAL(rp), DIMENSION(pchunk)   :: Bmag
  INTEGER                                                :: cc
  REAL(rp), DIMENSION(pchunk),INTENT(OUT)  :: PSIp
  REAL(rp), DIMENSION(pchunk,7)  :: A
  INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flag_cache
  REAL(rp), INTENT(IN) :: time
  REAL(rp), DIMENSION(pchunk) :: Y_T
  REAL(rp)  :: psip_conv

  psip_conv=F%psip_conv

  !$OMP SIMD
  do cc=1_idef,pchunk
     Y_T(cc)=F%t0_2x1t+time
  end do
  !$OMP END SIMD

  !write(output_unit_write,*) 't0',F%t0_2x1t,'time',time,'Y_T',Y_T(1)

  call check_if_in_fields_domain_p(pchunk,F,Y_R,Y_T,Y_Z,flag_cache)

  call EZspline_derivative(bfield_2X1T%A, pchunk, Y_R, Y_T, Y_Z, A, ezerr)
  call EZspline_error(ezerr)

  !A(:,1) = PSIp
  !A(:,2) = dPSIp/dR
  !A(:,3) = dPSIp/dT
  !A(:,4) = dPSIp/dZ
  !A(:,5) = d^2PSIp/dR^2
  !A(:,6) = d^2PSIp/dZ^2
  !A(:,7) = d^2PSIp/dRdZ

  !$OMP SIMD
!    !$OMP& aligned(PSIp,A,B_R,Y_R,B_PHI,B_Z,Bmag,gradB_R,gradB_PHI,gradB_Z, &
!    !$OMP& curlb_R,curlb_PHI,curlb_Z,E_R,E_PHI,E_Z)
  do cc=1_idef,pchunk
     PSIp(cc)=A(cc,1)

     B_R(cc) = psip_conv*A(cc,4)/Y_R(cc)
     ! BR = (dA/dZ)/R
     B_PHI(cc) = -F%Bo*F%Ro/Y_R(cc)
     ! BPHI = Fo*Ro/R
     B_Z(cc) = -psip_conv*A(cc,2)/Y_R(cc)
     ! BR = -(dA/dR)/R


     Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

     gradB_R(cc)=(B_R(cc)*psip_conv*A(cc,7)-B_Z(cc)*psip_conv*A(cc,5)- &
          Bmag(cc)*Bmag(cc))/(Y_R(cc)*Bmag(cc))
     gradB_PHI(cc)=0._rp
     gradB_Z(cc)=(B_R(cc)*psip_conv*A(cc,6)-B_Z(cc)*psip_conv*A(cc,7))/ &
          (Y_R(cc)*Bmag(cc))

     curlb_R(cc)=B_PHI(cc)*gradB_Z(cc)/(Bmag(cc)*Bmag(cc))
     curlb_PHI(cc)=(Bmag(cc)/Y_R(cc)*(B_Z(cc)+psip_conv*A(cc,5)+ &
          psip_conv*A(cc,6))-B_R(cc)*gradB_Z(cc)+B_Z(cc)*gradB_R(cc))/ &
          (Bmag(cc)*Bmag(cc))
     curlb_Z(cc)=-B_PHI(cc)*gradB_R(cc)/(Bmag(cc)*Bmag(cc))

     if (F%E_2x1t) then
        E_R(cc) = 0._rp
        E_PHI(cc) = A(cc,3)/(2._rp*C_PI*Y_R(cc))
        E_Z(cc) = 0._rp
     else
        E_R(cc) = 0._rp
        E_PHI(cc) = 0._rp
        E_Z(cc) = 0._rp
     end if

  end do
  !$OMP END SIMD

end subroutine calculate_GCfields_2x1t_p

subroutine calculate_GCfields_p_FS(pchunk,F,Y_R,Y_PHI,Y_Z,B_R,B_PHI,B_Z, &
     E_R,E_PHI,E_Z,curlb_R,curlb_PHI,curlb_Z,gradB_R,gradB_PHI,gradB_Z, &
     flag_cache,PSIp)
  INTEGER, INTENT(IN)  :: pchunk
  REAL(rp), DIMENSION(pchunk), INTENT(IN)      :: Y_R,Y_PHI,Y_Z
  TYPE(FIELDS), INTENT(IN)                               :: F
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: B_R,B_PHI,B_Z
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: gradB_R,gradB_PHI,gradB_Z
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: curlb_R,curlb_PHI,curlb_Z
  REAL(rp), DIMENSION(pchunk),  INTENT(OUT)   :: E_R,E_PHI,E_Z
  REAL(rp), DIMENSION(pchunk)   :: Bmag
  INTEGER                                                :: cc
  REAL(rp), DIMENSION(pchunk),INTENT(OUT)  :: PSIp
  REAL(rp), DIMENSION(pchunk,6)  :: A
  INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flag_cache
  REAL(rp)  :: psip_conv

  psip_conv=F%psip_conv

  call check_if_in_fields_domain_p(pchunk,F,Y_R,Y_PHI,Y_Z,flag_cache)

  call EZspline_derivative(bfield_2d%A, pchunk, Y_R, Y_Z, A, ezerr)
  call EZspline_error(ezerr)

  !write (output_unit_write,*) A(1,1),A(1,2)

  !A(:,1) = PSIp
  !A(:,2) = dPSIp/dR
  !A(:,3) = dPSIp/dZ
  !A(:,4) = d^2PSIp/dR^2
  !A(:,5) = d^2PSIp/dZ^2
  !A(:,6) = d^2PSIp/dRdZ

  !$OMP SIMD
!    !$OMP& aligned(PSIp,A,B_R,Y_R,B_PHI,B_Z,Bmag,gradB_R,gradB_PHI,gradB_Z, &
!    !$OMP& curlb_R,curlb_PHI,curlb_Z,E_R,E_PHI,E_Z)
  do cc=1_idef,pchunk
     PSIp(cc)=A(cc,1)

     A(cc,2)=psip_conv*A(cc,2)/(2*C_PI)
     A(cc,3)=psip_conv*A(cc,3)/(2*C_PI)
     A(cc,4)=psip_conv*A(cc,4)/(2*C_PI)
     A(cc,5)=psip_conv*A(cc,5)/(2*C_PI)
     A(cc,6)=psip_conv*A(cc,6)/(2*C_PI)


     B_R(cc) = A(cc,3)/Y_R(cc)
     ! BR = (dA/dZ)/R
     B_PHI(cc) = -F%Bo*F%Ro/Y_R(cc)
     ! BPHI = Fo*Ro/R
     B_Z(cc) = -A(cc,2)/Y_R(cc)
     ! BR = -(dA/dR)/R



     Bmag(cc)=sqrt(B_R(cc)*B_R(cc)+B_PHI(cc)*B_PHI(cc)+B_Z(cc)*B_Z(cc))

     gradB_R(cc)=(B_R(cc)*A(cc,6)-B_Z(cc)*A(cc,4)-Bmag(cc)*Bmag(cc))/ &
          (Y_R(cc)*Bmag(cc))
     gradB_PHI(cc)=0._rp
     gradB_Z(cc)=(B_R(cc)*A(cc,5)-B_Z(cc)*A(cc,6))/ &
          (Y_R(cc)*Bmag(cc))

     curlb_R(cc)=B_PHI(cc)*gradB_Z(cc)/(Bmag(cc)*Bmag(cc))
     curlb_PHI(cc)=(Bmag(cc)/Y_R(cc)*(B_Z(cc)+A(cc,4)+A(cc,5))- &
          B_R(cc)*gradB_Z(cc)+B_Z(cc)*gradB_R(cc))/ &
          (Bmag(cc)*Bmag(cc))
     curlb_Z(cc)=-B_PHI(cc)*gradB_R(cc)/(Bmag(cc)*Bmag(cc))

     E_R(cc) = 0._rp
     E_PHI(cc) = F%Eo*F%Ro/Y_R(cc)
     E_Z(cc) = 0._rp

  end do

end subroutine calculate_GCfields_p_FS

subroutine add_interp_SCE_p(params,F,Y_R,Y_PHI,Y_Z,E_PHI)
  TYPE(KORC_PARAMS), INTENT(IN)                              :: params
  TYPE(FIELDS), INTENT(IN)                               :: F
  REAL(rp), DIMENSION(params%pchunk), INTENT(IN)      :: Y_R,Y_PHI,Y_Z
  REAL(rp), DIMENSION(params%pchunk), INTENT(INOUT)      :: E_PHI

  REAL(rp),DIMENSION(params%pchunk) :: rm,E_SC_PHI
  REAL(rp) :: R0,Z0
  INTEGER :: cc,pchunk

  pchunk=params%pchunk
  R0=F%Ro
  Z0=F%Zo

  !$OMP SIMD
  do cc=1_idef,pchunk
     rm(cc)=sqrt((Y_R(cc)-R0)*(Y_R(cc)-R0)+(Y_Z(cc)-Z0)*(Y_Z(cc)-Z0))
  end do
  !$OMP END SIMD

  call EZspline_interp(efield_SC1d%PHI,pchunk, rm, E_SC_PHI, ezerr)
  call EZspline_error(ezerr)

  !$OMP SIMD
  do cc=1_idef,pchunk
     E_PHI(cc)=E_PHI(cc)+E_SC_PHI(cc)
  end do
  !$OMP END SIMD

end subroutine add_interp_SCE_p

subroutine add_interp_SCE_p_FS(params,F,PSIp,E_PHI)
  TYPE(KORC_PARAMS), INTENT(IN)                              :: params
  TYPE(FIELDS), INTENT(IN)                               :: F
  REAL(rp), DIMENSION(params%pchunk), INTENT(IN)      :: PSIp
  REAL(rp), DIMENSION(params%pchunk), INTENT(INOUT)      :: E_PHI

  REAL(rp),DIMENSION(params%pchunk) :: E_SC_PHI
  INTEGER :: cc,pchunk

  pchunk=params%pchunk

  call EZspline_interp(efield_SC1d%PHI,pchunk, PSIp, E_SC_PHI, ezerr)
  call EZspline_error(ezerr)

  !$OMP SIMD
  do cc=1_idef,pchunk
     E_PHI(cc)=E_PHI(cc)+E_SC_PHI(cc)
  end do
  !$OMP END SIMD

end subroutine add_interp_SCE_p_FS

subroutine calculate_initial_magnetic_field(F)

  TYPE(FIELDS), INTENT(INOUT)                               :: F
  REAL(rp),dimension(F%dims(1),F%dims(3),2)                  :: gradA
  INTEGER                                                :: ii
  INTEGER                                                :: jj
  REAL(rp) :: psip_conv

  psip_conv=F%psip_conv

  call EZspline_interp(bfield_2d%A,F%dims(1),F%dims(3),F%X%R, F%X%Z, &
       F%PSIp, ezerr)
  call EZspline_error(ezerr)

  ! FR = (dA/dZ)/R
  call EZspline_gradient(bfield_2d%A,F%dims(1),F%dims(3),F%X%R, F%X%Z, &
       gradA, ezerr)
  call EZspline_error(ezerr)

  do ii=1,F%dims(1)
     F%B_2D%R(ii,:) = psip_conv*gradA(ii,:,2)/F%X%R(ii)
     F%B_2D%PHI(ii,:) = -F%Bo*F%Ro/F%X%R(ii)
     F%B_2D%Z(ii,:) = -psip_conv*gradA(ii,:,1)/F%X%R(ii)
  end do

  !        write(output_unit_write,'("AR",E17.10)') gradA(1)
  !        write(output_unit_write,'("AZ",E17.10)') gradA(2)

end subroutine calculate_initial_magnetic_field

subroutine sample_poloidal_flux(F)

  TYPE(FIELDS), INTENT(INOUT)                               :: F

  ! FR = (dA/dZ)/R
  call EZspline_interp(bfield_2d%A,F%dims(1),F%dims(3),F%X%R, F%X%Z, &
       F%PSIp, ezerr)
  call EZspline_error(ezerr)


  !        write(output_unit_write,'("AR",E17.10)') gradA(1)
  !        write(output_unit_write,'("AZ",E17.10)') gradA(2)

end subroutine sample_poloidal_flux

!> @brief Subroutine for interpolating the pre-computed, axisymmetric electric field to the particles' position.
!!
!! @param[in] Y Particles' position in cylindrical coordinates, Y(1,:) = \(R\), Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
!! @param[in,out] E Cartesian components of interpolated electric field components. E(1,:)=\(E_x\), E(2,:)=\(E_y\), and E(3,:)=\(E_z\).
!! @param F Cylindrical components of interpolated magnetic field components. F(1,:)=\(E_R\), F(2,:)=\(E_\phi\), and F(3,:)=\(E_Z\).
!! @param flag Flag that indicates whether particles are followed in the simulation (flag=1), or not (flag=0).
!! @param pp Particle iterator.
!! @param ss Species iterator.
subroutine interp_2D_efields(params,Y,E,flag)
  TYPE(KORC_PARAMS), INTENT(IN)      :: params
  !! Core KORC simulation parameters.
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: E
  REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: F
  INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
  INTEGER                                                :: pp
  INTEGER                                                :: ss

!  write(output_unit_write,*) 'interp E fields'

  if (size(Y,1).eq.1) then
     ss = size(Y,1)
  else
     if (Y(2,1).eq.0) then
        ss=1_idef
     else
        ss = size(Y,1)
     end if
  endif

  ALLOCATE(F(ss,3))
  !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
  !$OMP& SHARED(F,Y,E,flag,efield_2d)
  do pp=1_idef,ss
     if ( flag(pp) .EQ. 1_is ) then
        call EZspline_interp(efield_2d%R, Y(pp,1), Y(pp,3), F(pp,1), ezerr)
        call EZspline_error(ezerr)

        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        end if

        call EZspline_interp(efield_2d%PHI, Y(pp,1), Y(pp,3), F(pp,2), ezerr)
        call EZspline_error(ezerr)

        call EZspline_interp(efield_2d%Z, Y(pp,1), Y(pp,3), F(pp,3), ezerr)
        call EZspline_error(ezerr)

        if (.not.params%GC_coords) then
           E(pp,1) = F(pp,1)*COS(Y(pp,2)) - F(pp,2)*SIN(Y(pp,2))
           E(pp,2) = F(pp,1)*SIN(Y(pp,2)) + F(pp,2)*COS(Y(pp,2))
           E(pp,3) = F(pp,3)
        else
           E(pp,1) = F(pp,1)
           E(pp,2) = F(pp,2)
           E(pp,3) = F(pp,3)
        end if

        !write(output_unit_write,*) 'EPHI',E(pp,2)
     end if

  end do
  !$OMP END PARALLEL DO
  DEALLOCATE(F)
end subroutine interp_2D_efields


!> @brief Subroutine for interpolating the pre-computed 3-D electric field to the particles' position.
!!
!! @param[in] Y Particles' position in cylindrical coordinates, Y(1,:) = \(R\), Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
!! @param[in,out] E Cartesian components of interpolated electric field components. E(1,:)=\(E_x\), E(2,:)=\(E_y\), and E(3,:)=\(E_z\).
!! @param F Cylindrical components of interpolated magnetic field components. F(1,:)=\(E_R\), F(2,:)=\(E_\phi\), and F(3,:)=\(E_Z\).
!! @param flag Flag that indicates whether particles are followed in the simulation (flag=1), or not (flag=0).
!! @param pp Particle iterator.
!! @param ss Species iterator.
subroutine interp_3D_efields(params,Y,E,flag)
  TYPE(KORC_PARAMS), INTENT(IN)      :: params
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)   :: E
  REAL(rp), DIMENSION(:,:), ALLOCATABLE                  :: F
  INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
  INTEGER                                                :: pp
  INTEGER                                                :: ss

  ss = size(Y,1)

  ALLOCATE(F(ss,3))
  !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
  !$OMP& SHARED(F,Y,E,flag,efield_3d)
  do pp=1_idef,ss
     if ( flag(pp) .EQ. 1_is ) then
        call EZspline_interp(efield_3d%R, Y(pp,1), Y(pp,2), Y(pp,3), &
             F(pp,1), ezerr)
        call EZspline_error(ezerr)

        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        end if

        call EZspline_interp(efield_3d%PHI, Y(pp,1), Y(pp,2), Y(pp,3), &
             F(pp,2), ezerr)
        call EZspline_error(ezerr)

        call EZspline_interp(efield_3d%Z, Y(pp,1), Y(pp,2), Y(pp,3), &
             F(pp,3), ezerr)
        call EZspline_error(ezerr)

        if (.not.params%GC_coords) then
           E(pp,1) = F(pp,1)*COS(Y(pp,2)) - F(pp,2)*SIN(Y(pp,2))
           E(pp,2) = F(pp,1)*SIN(Y(pp,2)) + F(pp,2)*COS(Y(pp,2))
           E(pp,3) = F(pp,3)
        else
           E(pp,1) = F(pp,1)
           E(pp,2) = F(pp,2)
           E(pp,3) = F(pp,3)
        end if

     end if
  end do
  !$OMP END PARALLEL DO
  DEALLOCATE(F)
end subroutine interp_3D_efields
#endif

subroutine interp_fields(params,prtcls,F)
  !! @note Subroutine that works as an interface for calling the
  !! appropriate subroutines for interpolating or calculating the
  !! electric and magnetic fields. @endnote
  TYPE(KORC_PARAMS), INTENT(IN)      :: params
    !! Core KORC simulation parameters.
  TYPE(PARTICLES), INTENT(INOUT) :: prtcls
    !! An instance of PARTICLES containing the variables of a given species.
  TYPE(FIELDS), INTENT(IN)       :: F
    !! An instance of KORC's derived type FIELDS containing all the
    !! information about the fields used in the simulation.
  !! See [[korc_types]] and [[korc_fields]].
  integer :: pp

  if (.not.params%GC_coords) call cart_to_cyl(prtcls%X,prtcls%Y)

!  write(output_unit_write,'("BR: ",E17.10)') prtcls%BR(:,1)
!  write(output_unit_write,'("Y: ",E17.10)') prtcls%X(2,1)
!  write(output_unit_write,'("Z: ",E17.10)') prtcls%X(3,1)

#ifdef PSPLINE
  call check_if_in_fields_domain(F,prtcls%Y, prtcls%flagCon)

  if (F%useLCFS) call check_if_in_LCFS(F,prtcls%Y, prtcls%initLCFS)
#endif

  !write(output_unit_write,*) 'checked domain'

#ifdef FIO
  if (TRIM(params%field_model) .eq. 'M3D_C1'.or. &
       TRIM(params%field_model) .eq. 'NIMROD') then

     if (F%FIO_B .ge. 0) then

        !write(6,*) 'interp_fields'

        call get_fio_magnetic_fields(prtcls, F, params)

     end if

     if (F%FIO_E .ge. 0) then
        call get_fio_electric_fields(prtcls, F, params)
     end if

     if (F%FIO_A .ge. 0) then
        call get_fio_vector_potential(prtcls, F, params)
     end if

     do pp=1,sizeof(prtcls%flagCon)
        if (prtcls%flagCon(pp)==0.and. &
             (.not.(params%restart.OR.params%proceed))) then
           write(6,*) 'RE initialized outside of computational domain!!!'
           call KORC_ABORT(15)
        end if
     end do

  end if
#endif

#ifdef PSPLINE
  if (params%field_model(10:13).eq.'MARS') then
     call interp_FOfields_mars(prtcls, F, params)
  end if

  if (params%field_model(10:14).eq.'AORSA') then
     call interp_FOfields_aorsa(prtcls, F, params)
  end if

  if ((ALLOCATED(F%PSIp).and.F%Bflux).or. &
       (F%ReInterp_2x1t.and.(.not.((TRIM(params%field_model).eq.'M3D_C1').or. &
       (params%field_model(10:13).eq.'MARS').or. &
       (TRIM(params%field_model).eq.'NIMROD'))))) then

     !     write(output_unit_write,'("3 size of PSI_P: ",I16)') size(prtcls%PSI_P)

     !     write(output_unit_write,'("B_X: ",E17.10)') prtcls%B(:,1)
     !     write(output_unit_write,'("B_Z: ",E17.10)') prtcls%B(:,3)
     !     write(output_unit_write,'("B_Y: ",E17.10)') prtcls%B(:,2)
     !     write(output_unit_write,'("PSI_P: ",E17.10)') prtcls%PSI_P

     call calculate_magnetic_field(params,prtcls%Y,F,prtcls%B,prtcls%E, &
          prtcls%PSI_P,prtcls%flagCon)

     !write(output_unit_write,*) 'interp PSIp'

     !     write(output_unit_write,'("interp_fields")')
     !     write(output_unit_write,'("B_X: ",E17.10)') prtcls%B(:,1)
     !     write(output_unit_write,'("B_Z: ",E17.10)') prtcls%B(:,3)
     !     write(output_unit_write,'("B_Y: ",E17.10)') prtcls%B(:,2)

  end if

  if (ALLOCATED(F%PSIp3D).and.F%Bflux3D) then

!     write(output_unit_write,'("3 size of PSI_P: ",I16)') size(prtcls%PSI_P)

!     write(output_unit_write,'("B_X: ",E17.10)') prtcls%B(:,1)
!     write(output_unit_write,'("B_Z: ",E17.10)') prtcls%B(:,3)
!     write(output_unit_write,'("B_Y: ",E17.10)') prtcls%B(:,2)
!     write(output_unit_write,'("PSI_P: ",E17.10)') prtcls%PSI_P

     call calculate_magnetic_field(params,prtcls%Y,F,prtcls%B,prtcls%E, &
          prtcls%PSI_P,prtcls%flagCon)

!     write(output_unit_write,'("interp_fields")')
!     write(output_unit_write,'("B_X: ",E17.10)') prtcls%B(:,1)
!     write(output_unit_write,'("B_Z: ",E17.10)') prtcls%B(:,3)
!     write(output_unit_write,'("B_Y: ",E17.10)') prtcls%B(:,2)

  end if

  if (ALLOCATED(F%B_2D%R).and.F%Bfield) then
     call interp_2D_bfields(params,prtcls%Y,prtcls%B,prtcls%flagCon)
  end if

  if (ALLOCATED(F%B_3D%R).and.F%Bfield) then
     call interp_3D_bfields(params,prtcls%Y,prtcls%B,prtcls%flagCon)
  end if

!  if (ALLOCATED(F%E_2D%R).and.F%Efield) then
!     call interp_2D_efields(params,prtcls%Y,prtcls%E,prtcls%flagCon)
!  end if

  if (ALLOCATED(F%E_3D%R).and.F%Efield.and.(.not.F%Dim2x1t)) then
     call interp_3D_efields(params,prtcls%Y,prtcls%E,prtcls%flagCon)

  end if

  if (ALLOCATED(F%E_3D%R).and.F%Efield.and.F%Dim2x1t) then
     call interp_2D_efields(params,prtcls%Y,prtcls%E,prtcls%flagCon)

!     write(output_unit_write,*) 'interpolated efield'

  end if

  if (params%GC_coords.and.ALLOCATED(F%gradB_2D%R).and.F%Bfield) then
     call interp_2D_gradBfields(prtcls%Y,prtcls%gradB,prtcls%flagCon)

  end if

  if (params%GC_coords.and.ALLOCATED(F%gradB_2D%R).and.F%Bfield) then
     call interp_2D_curlbfields(prtcls%Y,prtcls%curlb,prtcls%flagCon)
  end if

  if(params%GC_coords.and.params%orbit_model(3:6)=='grad') then
     call gradient_2D_bfields(prtcls%Y,prtcls%BR,prtcls%BPHI, &
          prtcls%BZ,prtcls%flagCon)
  end if
#endif

end subroutine interp_fields

#ifdef PSPLINE
subroutine interp_Hcollision_p(pchunk,Y_R,Y_PHI,Y_Z,ne,Te,Zeff, &
     nAr0,nAr1,nAr2,nAr3,nD,nD1,flag_cache)
  INTEGER, INTENT(IN)  :: pchunk
  REAL(rp),DIMENSION(pchunk),INTENT(IN)   :: Y_R,Y_PHI,Y_Z
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: ne,Te,Zeff
  REAL(rp),DIMENSION(pchunk),INTENT(OUT)   :: nAr0,nAr1,nAr2,nAr3,nD,nD1
  INTEGER(is),DIMENSION(pchunk),INTENT(INOUT)   :: flag_cache

  call check_if_in_profiles_domain_p(pchunk,Y_R,Y_PHI,Y_Z,flag_cache)
!  write(output_unit_write,'("YR: ",E17.10)') Y_R(1)
!  write(output_unit_write,'("YPHI: ",E17.10)') Y_PHI(1)
!  write(output_unit_write,'("YZ: ",E17.10)') Y_Z(1)

!  write(output_unit_write,'("Te_interp_R",E17.10)') profiles_2d%Te%x1
!  write(output_unit_write,'("Te_interp_Z",E17.10)') profiles_2d%Te%x2

  call EZspline_interp(profiles_2d%ne,profiles_2d%Te,profiles_2d%Zeff, &
       profiles_2d%nAr0,profiles_2d%nAr1,profiles_2d%nAr2,profiles_2d%nAr3, &
       profiles_2d%nD,profiles_2d%nD1,pchunk,Y_R,Y_Z,ne,Te,Zeff, &
       nAr0,nAr1,nAr2,nAr3,nD,nD1,ezerr)
  ! this will call PSPLINE routine EZspline_interp2_bmag_cloud as there
  ! is the same number of entries
  call EZspline_error(ezerr)

end subroutine interp_Hcollision_p

subroutine interp_nRE(params,Y_R,Y_PHI,Y_Z,PSIp,EPHI,ne,Te,nRE, &
     nAr0,nAr1,nAr2,nAr3,nD,nD1,g_test,fRE_out,rho1D)
  TYPE(KORC_PARAMS), INTENT(IN) 	:: params
  REAL(rp),INTENT(IN)   :: Y_R,Y_PHI,Y_Z,g_test
  REAL(rp),INTENT(OUT)   :: PSIp,EPHI,ne,Te,nRE,fRE_out
  REAL(rp),INTENT(OUT)   :: nAr0,nAr1,nAr2,nAr3,nD,nD1
  REAL(rp)   :: RHON
  INTEGER(is),DIMENSION(1)   :: flag_cache
  REAL(rp), INTENT(IN),optional  :: rho1D

  flag_cache=1_is
  call check_if_in_profiles_domain_p(1,(/Y_R/),(/Y_PHI/),(/Y_Z/), &
       flag_cache)

  if (flag_cache(1).ne.0_is) then
     call EZspline_interp(bfield_2d%A,efield_2d%PHI, &
          profiles_2d%ne,profiles_2d%Te, &
          profiles_2d%nRE,profiles_2d%nAr0,profiles_2d%nAr1, &
          profiles_2d%nAr2,profiles_2d%nAr3,profiles_2d%nD,profiles_2d%nD1, &
          profiles_2d%RHON,Y_R,Y_Z, &
          PSIp,EPHI,ne,Te,nRE,nAr0,nAr1,nAr2,nAr3,nD,nD1,RHON,ezerr)
     call EZspline_error(ezerr)

     !write(6,*) 'RHON',RHON

     if ((RHON.le.1).and.(RHON.gt.0)) then
        if(present(rho1D)) RHON=rho1D

        call EZspline_interp(hollmann_2d%fRE_E,RHON,g_test,fRE_out,ezerr)
        call EZspline_error(ezerr)
        if (ezerr .NE. 0) then ! We flag the particle as lost
           write(6,*) 'R,Z',Y_R*params%cpp%length,Y_Z*params%cpp%length
           write(6,*) 'RHON',RHON
        end if

     else
        fRE_out=0._rp
     endif

     if (fRE_OUT.lt.0) fRE_OUT=0._rp
  else
     fRE_OUT=0._rp
     nRE=0._rp
  endif

end subroutine interp_nRE

subroutine interp_2D_profiles(Y,ne,Te,Zeff,flag)
  !! @note Subroutine for interpolating the pre-computed, axisymmetric
  !! plasma profiles to the particles' position. @endnote
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
  !! Particles' position in cylindrical coordinates,
  !! Y(1,:) = \(R\), Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
  REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: ne
  !! Interpolated background electron density !!\(n_e(R,Z)\).
  REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: Te
  !! Interpolated background electron temperature \(T_e(R,Z)\).
  REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: Zeff
  !! Interpolated effective charge number \(Z_{eff}(R,Z)\).
  INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
  !! Flag that indicates whether particles are followed in the
  !! simulation (flag=1), or not (flag=0).
  INTEGER                                                :: pp
  !! Particle iterator.
  INTEGER                                                :: ss
  !! Species iterator.

  if (size(Y,1).eq.1) then
     ss = size(Y,1)
  else
     if (Y(2,1).eq.0) then
        ss=1_idef
     else
        ss = size(Y,1)
     end if
  endif

!  write(output_unit_write,'("Also R_buffer: ",E17.10)') Y(1,ss)

  !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
  !$OMP& SHARED(Y,ne,Te,Zeff,flag,profiles_2d)
  do pp=1_idef,ss
     if ( flag(pp) .EQ. 1_is ) then
        call EZspline_interp(profiles_2d%ne, Y(pp,1), Y(pp,3), ne(pp), ezerr)
        call EZspline_error(ezerr)

!        write(output_unit_write,'("Also R_buffer: ",E17.10)') Y(pp,1)

        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        end if

        call EZspline_interp(profiles_2d%Te, Y(pp,1), Y(pp,3), Te(pp), ezerr)
        call EZspline_error(ezerr)

        call EZspline_interp(profiles_2d%Zeff, Y(pp,1), Y(pp,3), Zeff(pp), ezerr)
        call EZspline_error(ezerr)
     end if
  end do
  !$OMP END PARALLEL DO
end subroutine interp_2D_profiles


subroutine interp_3D_profiles(Y,ne,Te,Zeff,flag)
  !! @note Subroutine for interpolating the pre-computed,
  !! 3-D plasma profiles to the particles' position. @endnote
  REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: Y
  !! Particles' position in cylindrical coordinates,
  !! Y(1,:) = \(R\), Y(2,:) = \(\phi\), and Y(3,:) = \(Z\).
  REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: ne
  !! Interpolated background electron density \(n_e(R,\phi,Z)\).
  REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: Te
  !! Interpolated background electron temperature \(T_e(R,\phi,Z)\).
  REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: Zeff
  !! Interpolated effective charge number \(Z_{eff}(R,\phi,Z)\).
  INTEGER(is), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: flag
  !! Flag that indicates whether particles are followed in
  !! the simulation (flag=1), or not (flag=0).
  INTEGER                                                :: pp
  !! Particle iterator.
  INTEGER                                                :: ss
  !! Species iterator.

  ss = size(Y,1)

  !$OMP PARALLEL DO FIRSTPRIVATE(ss) PRIVATE(pp,ezerr) &
  !$OMP& SHARED(Y,ne,Te,Zeff,flag,profiles_2d)
  do pp=1_idef,ss
     if ( flag(pp) .EQ. 1_is ) then
        call EZspline_interp(profiles_3d%ne, Y(pp,1), Y(pp,2), Y(pp,3), &
             ne(pp), ezerr)
        call EZspline_error(ezerr)

        if (ezerr .NE. 0) then ! We flag the particle as lost
           flag(pp) = 0_is
        end if

        call EZspline_interp(profiles_3d%Te, Y(pp,1), Y(pp,2), Y(pp,3), &
             Te(pp), ezerr)
        call EZspline_error(ezerr)

        call EZspline_interp(profiles_3d%Zeff, Y(pp,1), Y(pp,2), Y(pp,3), &
             Zeff(pp), ezerr)
        call EZspline_error(ezerr)
     end if
  end do
  !$OMP END PARALLEL DO
end subroutine interp_3D_profiles
#endif


subroutine interp_profiles(params,prtcls,P)
  !! @note Subroutine that calls the appropriate subroutines for
  !! interpolating the 2-D or 3-D plasma profiles to the particles'
  !! position. @endnote
  TYPE(KORC_PARAMS), INTENT(IN)      :: params
  !! Core KORC simulation parameters.
  TYPE(PARTICLES), INTENT(INOUT) :: prtcls
  !! An instance of PARTICLES containing the variables of a
  !! given species. Call to this subroutine generally passes spp%vars.
  TYPE(PROFILES), INTENT(IN)     :: P
  !! An instance of KORC's derived type PROFILES containing all the
  !! information about the plasma profiles used in the simulation.
  !! See[[ korc_types]] and [[korc_profiles]].

  if (.not.params%GC_coords) call cart_to_cyl(prtcls%X,prtcls%X)

  !write(output_unit_write,'("Also R_buffer: ",E17.10)') prtcls%Y(1,1)
#ifdef PSPLINE
  call check_if_in_profiles_domain(prtcls%Y, prtcls%flagCon)
#endif

  if (ALLOCATED(P%ne_2D)) then
#ifdef PSPLINE
!     write(output_unit_write,'("Also R_buffer: ",E17.10)') prtcls%X(1,1)
     call interp_2D_profiles(prtcls%Y,prtcls%ne,prtcls%Te,prtcls%Zeff, &
          prtcls%flagCon)
  else if (ALLOCATED(P%ne_3D)) then
     call interp_3D_profiles(prtcls%Y,prtcls%ne,prtcls%Te,prtcls%Zeff, &
          prtcls%flagCon)
#endif
#ifdef FIO
  else if (P%FIO_ne   .ge. 0 .or.     &
       P%FIO_te   .ge. 0 .or.         &
       P%FIO_zeff .ge. 0) then
     call get_fio_profile(prtcls, P, params)

#endif
  else
     write(output_unit_write,'("Error: NO PROFILES ALLOCATED")')
     call KORC_ABORT(16)
  end if
end subroutine interp_profiles


#ifdef PSPLINE
!> @brief Subroutine that frees memory allocated for PSPLINE interpolants.
!!
!! @param[in] params Core KORC simulation parameters.
subroutine finalize_interpolants(params)
  TYPE(KORC_PARAMS), INTENT(IN) :: params

  if ((params%field_model(1:8) .EQ. 'EXTERNAL').or. &
       (params%field_eval.eq.'interp')) then
     if (params%mpi_params%rank .EQ. 0) then
        write(output_unit_write,'("* * * * FINALIZING FIELD INTERPOLANT * * * *")')
     end if

     if (EZspline_allocated(bfield_3d%R)) call Ezspline_free(bfield_3d%R, ezerr)
     if (EZspline_allocated(bfield_3d%PHI)) &
          call Ezspline_free(bfield_3d%PHI,ezerr)

     if (EZspline_allocated(bfield_3d%Z)) call Ezspline_free(bfield_3d%Z, ezerr)
     if (EZspline_allocated(bfield_2d%A)) call Ezspline_free(bfield_2d%A, ezerr)
     if (EZspline_allocated(bfield_2d%R)) call Ezspline_free(bfield_2d%R, ezerr)
     if (EZspline_allocated(bfield_2d%PHI)) &
          call Ezspline_free(bfield_2d%PHI,ezerr)

     if (EZspline_allocated(bfield_2d%Z)) call Ezspline_free(bfield_2d%Z, ezerr)

     if (EZspline_allocated(gradB_2d%R)) call Ezspline_free(gradB_2d%R, ezerr)
     if (EZspline_allocated(gradB_2d%PHI)) &
          call Ezspline_free(gradB_2d%PHI, ezerr)

     if (EZspline_allocated(gradB_2d%Z)) call Ezspline_free(gradB_2d%Z, ezerr)

     if (EZspline_allocated(curlb_2d%R)) call Ezspline_free(curlb_2d%R, ezerr)
     if (EZspline_allocated(curlb_2d%PHI)) &
          call Ezspline_free(curlb_2d%PHI, ezerr)

     if (EZspline_allocated(gradB_3d%R)) call Ezspline_free(gradB_3d%R, ezerr)
     if (EZspline_allocated(gradB_3d%PHI)) &
          call Ezspline_free(gradB_3d%PHI, ezerr)

     if (EZspline_allocated(gradB_3d%Z)) call Ezspline_free(gradB_3d%Z, ezerr)

     if (EZspline_allocated(curlb_3d%R)) call Ezspline_free(curlb_3d%R, ezerr)
     if (EZspline_allocated(curlb_3d%PHI)) &
          call Ezspline_free(curlb_3d%PHI, ezerr)

     if (EZspline_allocated(curlb_3d%Z)) call Ezspline_free(curlb_3d%Z, ezerr)

     if (ALLOCATED(profiles_domain%FLAG2D)) DEALLOCATE(profiles_domain%FLAG2D)
     if (ALLOCATED(profiles_domain%FLAG3D)) DEALLOCATE(profiles_domain%FLAG3D)

     if (params%mpi_params%rank .EQ. 0) then
        write(output_unit_write,'("* * * * FIELD INTERPOLANT FINALIZED * * * *")')
     end if
  end if

  if (params%profile_model(1:8) .EQ. 'EXTERNAL') then
     if (params%mpi_params%rank .EQ. 0) then
        write(output_unit_write,'("* * * * FINALIZING PROFILE INTERPOLANT * * * *")')
     end if

     if (EZspline_allocated(profiles_3d%ne)) &
          call Ezspline_free(profiles_3d%ne,ezerr)

     if (EZspline_allocated(profiles_3d%Te)) &
          call Ezspline_free(profiles_3d%Te,ezerr)

     if (EZspline_allocated(profiles_3d%Zeff)) call Ezspline_free( &
          profiles_3d%Zeff, ezerr)
     if (EZspline_allocated(profiles_2d%ne)) &
          call Ezspline_free(profiles_2d%ne,ezerr)

     if (EZspline_allocated(profiles_2d%Te)) &
          call Ezspline_free(profiles_2d%Te,ezerr)

     if (EZspline_allocated(profiles_2d%Zeff)) call Ezspline_free( &
          profiles_2d%Zeff, ezerr)

     if (ALLOCATED(profiles_domain%FLAG2D)) DEALLOCATE(profiles_domain%FLAG2D)
     if (ALLOCATED(profiles_domain%FLAG3D)) DEALLOCATE(profiles_domain%FLAG3D)

     if (params%mpi_params%rank .EQ. 0) then
        write(output_unit_write,'("* * * * PROFILE INTERPOLANT FINALIZED * * * *")')
     end if
  end if
end subroutine finalize_interpolants
#endif

#ifdef FIO
  !!  @note FIXME Add documentation
subroutine get_fio_magnetic_fields(prtcls, F, params)
  USE omp_lib
  IMPLICIT NONE

    TYPE(PARTICLES), INTENT(INOUT) :: prtcls
    TYPE(FIELDS), INTENT(IN)       :: F
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    INTEGER (C_INT)                :: status
    INTEGER                        :: pp
    REAL(rp), DIMENSION(3)         :: x
    REAL(rp), DIMENSION(3)         :: Btmp
    TYPE(C_PTR), DIMENSION(size(prtcls%hint)) :: hint
    INTEGER             :: thread_num

!    write(output_unit_write,*) 'in m3dc1 B'

    if (prtcls%cart) then
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(pp,status,x)
       do pp = 1, SIZE(prtcls%hint)
          if (prtcls%flagCon(pp) .EQ. 1_is) then
             x = prtcls%X(pp,:)*params%cpp%length
             status = fio_eval_field(F%FIO_B, x(1),   &
                  prtcls%B(pp,1),                        &
                  prtcls%hint(pp))

             if (status .eq. FIO_NO_DATA) then
                prtcls%B(pp,:) = 0
                prtcls%flagCon(pp) = 0_is
             else if (status .ne. FIO_SUCCESS) then
                prtcls%flagCon(pp) = 0_is
             end if
          end if
       end do
       !$OMP END PARALLEL DO
    else

!       write(output_unit_write,*) 'in cart false'
       !hint=prtcls%hint
       !write(output_unit_write,*) 'hint: ',hint

       Btmp=0._rp

       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP& SHARED(prtcls,params,F) &
       !$OMP& PRIVATE(pp,status,x,thread_num) &
       !$OMP& FIRSTPRIVATE(Btmp)
       do pp = 1, SIZE(prtcls%hint)

          thread_num = OMP_GET_THREAD_NUM()

          if (prtcls%flagCon(pp) .EQ. 1_is) then
             x(1) = prtcls%Y(pp,1)*params%cpp%length
             x(2) = prtcls%Y(pp,2)
             x(3) = prtcls%Y(pp,3)*params%cpp%length

             !prtcls%hint(pp)=c_null_ptr

             !write(6,*) 'thread',thread_num,'X',x

             !             prtcls%hint(pp)=c_null_ptr

             !write(output_unit_write,*) 'thread',thread_num,'before interpolating B'

             status = fio_eval_field(F%FIO_B, x(1),                      &
                  Btmp(1),prtcls%hint(pp))


             if (status .eq. FIO_NO_DATA) then
                prtcls%B(pp,:) = 0
                prtcls%flagCon(pp) = 0_is
             else if (status .ne. FIO_SUCCESS) then
                prtcls%flagCon(pp) = 0_is
             end if


             if (.not.params%GC_coords) then

                prtcls%B(pp,1)=(Btmp(1)*cos(x(2))-Btmp(2)*sin(x(2)))/ &
                     params%cpp%Bo
                prtcls%B(pp,2)=(Btmp(1)*sin(x(2))+Btmp(2)*cos(x(2)))/ &
                     params%cpp%Bo
                prtcls%B(pp,3)=Btmp(3)/params%cpp%Bo

             else

                prtcls%B(pp,1)=Btmp(1)/params%cpp%Bo
                prtcls%B(pp,2)=Btmp(2)/params%cpp%Bo
                prtcls%B(pp,3)=Btmp(3)/params%cpp%Bo
             end if


          end if
       end do
       !$OMP END PARALLEL DO

    end if
  end subroutine get_fio_magnetic_fields

  subroutine get_fio_FOmagnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       B_X,B_Y,B_Z,flag,hint)
    TYPE(FIELDS), INTENT(IN)       :: F
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    REAL(rp), DIMENSION(params%pchunk), INTENT(IN)  :: Y_R,Y_PHI,Y_Z
    REAL(rp), DIMENSION(params%pchunk), INTENT(INOUT)  :: B_X,B_Y,B_Z
    INTEGER(is), DIMENSION(params%pchunk), INTENT(INOUT)  :: flag
    TYPE(C_PTR), DIMENSION(params%pchunk), INTENT(INOUT)  :: hint
    INTEGER (C_INT)                :: status
    INTEGER                        :: pp,pchunk
    REAL(rp), DIMENSION(3)         :: x
    REAL(rp), DIMENSION(3)         :: Btmp

    pchunk=params%pchunk

    do pp = 1,pchunk
       if (flag(pp) .EQ. 1_is) then
          x(1) = Y_R(pp)*params%cpp%length
          x(2) = Y_PHI(pp)
          x(3) = Y_Z(pp)*params%cpp%length

          !             prtcls%hint(pp)=c_null_ptr

          status = fio_eval_field(F%FIO_B, x(1),                      &
               Btmp(1),hint(pp))

          if (status .eq. FIO_SUCCESS) then
             B_X(pp)=(Btmp(1)*cos(x(2))-Btmp(2)*sin(x(2)))/ &
                  params%cpp%Bo
             B_Y(pp)=(Btmp(1)*sin(x(2))+Btmp(2)*cos(x(2)))/ &
                  params%cpp%Bo
             B_Z(pp)=Btmp(3)/params%cpp%Bo
          else if (status .eq. FIO_NO_DATA) then
             flag(pp) = 0_is
          else if (status .ne. FIO_SUCCESS) then
             flag(pp) = 0_is
          end if

       end if
    end do

  end subroutine get_fio_FOmagnetic_fields_p

  subroutine get_fio_GCmagnetic_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       B_R,B_PHI,B_Z,gradB_R,gradB_PHI,gradB_Z, &
       curlb_R,curlb_PHI,curlb_Z,flag,hint)
    TYPE(FIELDS), INTENT(IN)       :: F
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    REAL(rp), DIMENSION(params%pchunk), INTENT(IN)  :: Y_R,Y_PHI,Y_Z
    REAL(rp), DIMENSION(params%pchunk), INTENT(INOUT)  :: B_R,B_PHI,B_Z
    REAL(rp), DIMENSION(params%pchunk), INTENT(INOUT)  :: gradB_R,gradB_PHI,gradB_Z
    REAL(rp), DIMENSION(params%pchunk), INTENT(INOUT)  :: curlb_R,curlb_PHI,curlb_Z
    INTEGER(is), DIMENSION(params%pchunk), INTENT(INOUT)  :: flag
    TYPE(C_PTR), DIMENSION(params%pchunk), INTENT(INOUT)  :: hint
    INTEGER (C_INT)                :: status
    INTEGER                        :: pp,pchunk
    REAL(rp), DIMENSION(3)         :: x
    REAL(rp), DIMENSION(3)         :: Btmp
    REAL(rp), DIMENSION(9)         :: dBtmp
    REAL(rp)  :: Bmag,dBRdR,dBPHIdR,dBZdR,dBRdPHI,dBPHIdPHI,dBZdPHI,dBRdZ,dBPHIdZ,dBZdZ

    pchunk=params%pchunk

    do pp = 1,pchunk
       if (flag(pp) .EQ. 1_is) then
          x(1) = Y_R(pp)*params%cpp%length
          x(2) = Y_PHI(pp)
          x(3) = Y_Z(pp)*params%cpp%length

          !             prtcls%hint(pp)=c_null_ptr

          !if (pp.eq.1) write(output_unit_write,*) 'Yinterp',x

          status = fio_eval_field(F%FIO_B, x(1), &
               Btmp(1),hint(pp))

          if (status .eq. FIO_SUCCESS) then

             !if (pp.eq.1) write(output_unit_write,*) 'interp success!'

             B_R(pp)=Btmp(1)/params%cpp%Bo
             B_PHI(pp)=Btmp(2)/params%cpp%Bo
             B_Z(pp)=Btmp(3)/params%cpp%Bo

             Bmag=sqrt(B_R(pp)*B_R(pp)+B_PHI(pp)*B_PHI(pp)+B_Z(pp)*B_Z(pp))

             status = fio_eval_field_deriv(F%FIO_B, x(1),dBtmp(1),hint(pp))

             !dBRdR=dBtmp(FIO_DR_R)*(params%cpp%length/params%cpp%Bo)
             !dBPHIdR=dBtmp(FIO_DR_PHI)*(params%cpp%length/params%cpp%Bo)
             !dBZdR=dBtmp(FIO_DR_Z)*(params%cpp%length/params%cpp%Bo)
             !dBRdPHI=dBtmp(FIO_DPHI_R)*(params%cpp%length/params%cpp%Bo)
             !dBPHIdPHI=dBtmp(FIO_DPHI_PHI)*(params%cpp%length/params%cpp%Bo)
             !dBZdPHI=dBtmp(FIO_DPHI_Z)*(params%cpp%length/params%cpp%Bo)
             !dBRdZ=dBtmp(FIO_DZ_R)*(params%cpp%length/params%cpp%Bo)
             !dBPHIdZ=dBtmp(FIO_DZ_PHI)*(params%cpp%length/params%cpp%Bo)
             !dBZdZ=dBtmp(FIO_DZ_Z)*(params%cpp%length/params%cpp%Bo)

             dBRdR=dBtmp(1)*(params%cpp%length/params%cpp%Bo)
             dBPHIdR=dBtmp(2)*(params%cpp%length/params%cpp%Bo)
             dBZdR=dBtmp(3)*(params%cpp%length/params%cpp%Bo)
             dBRdPHI=dBtmp(4)/params%cpp%Bo
             dBPHIdPHI=dBtmp(5)/params%cpp%Bo
             dBZdPHI=dBtmp(6)/params%cpp%Bo
             dBRdZ=dBtmp(7)*(params%cpp%length/params%cpp%Bo)
             dBPHIdZ=dBtmp(8)*(params%cpp%length/params%cpp%Bo)
             dBZdZ=dBtmp(9)*(params%cpp%length/params%cpp%Bo)

             !write(6,*) 'dBRdR',dBRdR
             !write(6,*) 'dBPHIdR',dBPHIdR
             !write(6,*) 'dBZdR',dBZdR
             !write(6,*) 'dBRdPHI',dBRdPHI
             !write(6,*) 'dBPHIdPHI',dBPHIdPHI
             !write(6,*) 'dBZdPHI',dBZdPHI
             !write(6,*) 'dBRdZ',dBRdZ
             !write(6,*) 'dBPHIdZ',dBPHIdZ
             !write(6,*) 'dBZdZ',dBZdZ

             gradB_R(pp)=(B_R(pp)*dBRdR+B_PHI(pp)*dBPHIdR+B_Z(pp)*dBZdR)/ &
                  Bmag
             gradB_PHI(pp)=(B_R(pp)*dBRdPHI+B_PHI(pp)*dBPHIdPHI+ &
                  B_Z(pp)*dBZdPHI)/(Y_R(pp)*Bmag)
             gradB_Z(pp)=(B_R(pp)*dBRdZ+B_PHI(pp)*dBPHIdZ+B_Z(pp)*dBZdZ)/ &
                  Bmag

             curlb_R(pp)=(Bmag/Y_R(pp)*dBZdPHI-B_Z(pp)*gradB_PHI(pp)- &
                  Bmag*dBPHIdZ+B_PHI(pp)*gradB_Z(pp))/(Bmag*Bmag)
             curlb_PHI(pp)=(Bmag*dBRdZ-B_R(pp)*gradB_Z(pp)- &
                  Bmag*dBZdR+B_Z(pp)*gradB_R(pp))/(Bmag*Bmag)
             curlb_Z(pp)=(Bmag/Y_R(pp)*B_PHI(pp)+Bmag*dBPHIdR- &
                  B_PHI(pp)*gradB_R(pp)- &
                  Bmag/Y_R(pp)*dBRdPHI+B_R(pp)*gradB_PHI(pp))/(Bmag*Bmag)

          else if (status .eq. FIO_NO_DATA) then
             flag(pp) = 0_is
          else if (status .ne. FIO_SUCCESS) then
             flag(pp) = 0_is
          end if

       end if
    end do

  end subroutine get_fio_GCmagnetic_fields_p

  subroutine get_fio_vector_potential(prtcls, F, params)
    TYPE(PARTICLES), INTENT(INOUT) :: prtcls
    TYPE(FIELDS), INTENT(IN)       :: F
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    INTEGER (C_INT)                :: status
    INTEGER                        :: pp
    REAL(rp), DIMENSION(3)         :: x
    REAL(rp), DIMENSION(3)         :: Atmp
    integer(ip)  ::  ss

    if (size(prtcls%Y,1).eq.1) then
       ss = size(prtcls%Y,1)
    else
       if (prtcls%Y(2,1).eq.0) then
          ss=1_idef
       else
          ss = size(prtcls%Y,1)
       end if
    endif

    Atmp=0._rp

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP& PRIVATE(pp,status,x) &
    !$OMP& FIRSTPRIVATE(Atmp)
    do pp = 1,ss
       if (prtcls%flagCon(pp) .EQ. 1_is) then
          x(1) = prtcls%Y(pp,1)*params%cpp%length
          x(2) = prtcls%Y(pp,2)
          x(3) = prtcls%Y(pp,3)*params%cpp%length

          !prtcls%hint(pp)=c_null_ptr

          !write(output_unit_write,*) F%FIO_A,x,Atmp

          status = fio_eval_field(F%FIO_A, x(1),                      &
               Atmp(1),prtcls%hint(pp))

          if (status .eq. FIO_SUCCESS) then
             prtcls%PSI_P(pp)=-Atmp(2)*x(1)
          else if (status .eq. FIO_NO_DATA) then
             prtcls%PSI_P(pp) = 100._rp
             prtcls%flagCon(pp) = 0_is
          else if (status .ne. FIO_SUCCESS) then
             prtcls%flagCon(pp) = 0_is
             prtcls%PSI_P(pp) = 100._rp
          end if


       end if
    end do
    !$OMP END PARALLEL DO

  end subroutine get_fio_vector_potential

  subroutine get_fio_vector_potential_p(params,F,Y_R,Y_PHI,Y_Z, &
       PSIp,flag,hint)
    TYPE(FIELDS), INTENT(IN)       :: F
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    REAL(rp), DIMENSION(params%pchunk), INTENT(IN)  :: Y_R,Y_PHI,Y_Z
    REAL(rp), DIMENSION(params%pchunk), INTENT(INOUT)  :: PSIp
    INTEGER(is), DIMENSION(params%pchunk), INTENT(INOUT)  :: flag
    TYPE(C_PTR), DIMENSION(params%pchunk), INTENT(INOUT)  :: hint
    INTEGER (C_INT)                :: status
    INTEGER                        :: pp,pchunk
    REAL(rp), DIMENSION(3)         :: x
    REAL(rp), DIMENSION(3)         :: Atmp
    integer(ip)  ::  ss

    pchunk=params%pchunk

    do pp = 1,pchunk
       if (flag(pp) .EQ. 1_is) then
          x(1) = Y_R(pp)*params%cpp%length
          x(2) = Y_PHI(pp)
          x(3) = Y_Z(pp)*params%cpp%length

          !             prtcls%hint(pp)=c_null_ptr


          status = fio_eval_field(F%FIO_A, x(1),                      &
               Atmp(1),hint(pp))

          if (status .eq. FIO_SUCCESS) then
             PSIp(pp)=-Atmp(2)*x(1)
          else if (status .eq. FIO_NO_DATA) then
             PSIp(pp) = 100._rp
             flag(pp) = 0_is
          else if (status .ne. FIO_SUCCESS) then
             PSIp(pp) = 100._rp
             flag(pp) = 0_is
          end if

       end if
    end do

  end subroutine get_fio_vector_potential_p

  !!  @note FIXME Add documentation
  subroutine get_fio_electric_fields(prtcls, F, params)
    TYPE(PARTICLES), INTENT(INOUT) :: prtcls
    TYPE(FIELDS), INTENT(IN)       :: F
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    INTEGER (C_INT)                :: status
    INTEGER                        :: pp
    REAL(rp), DIMENSION(3)         :: x
    REAL(rp), DIMENSION(3)         :: Etmp

    if (prtcls%cart) then
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(pp,status,x)
       do pp = 1, SIZE(prtcls%hint)
          if (prtcls%flagCon(pp) .EQ. 1_is) then
             x = prtcls%X(pp,:)*params%cpp%length
             status = fio_eval_field(F%FIO_E, x(1),                      &
                  prtcls%E(pp,1),                        &
                  prtcls%hint(pp))

             if (status .eq. FIO_NO_DATA) then
                prtcls%E(pp,:) = 0
             else if (status .ne. FIO_SUCCESS) then
                prtcls%flagCon(pp) = 0_is
             end if
          end if
       end do
       !$OMP END PARALLEL DO
    else

       Etmp=0._rp

       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP& SHARED(prtcls,params,F) &
       !$OMP& PRIVATE(pp,status,x) &
       !$OMP& FIRSTPRIVATE(Etmp)
       do pp = 1, SIZE(prtcls%hint)
          if (prtcls%flagCon(pp) .EQ. 1_is) then
             x(1) = prtcls%Y(pp,1)*params%cpp%length
             x(2) = prtcls%Y(pp,2)
             x(3) = prtcls%Y(pp,3)*params%cpp%length

             status = fio_eval_field(F%FIO_E, x(1),                      &
                  Etmp(1),prtcls%hint(pp))

             if (status .eq. FIO_NO_DATA) then
                prtcls%E(pp,:) = 0
                prtcls%flagCon(pp) = 0_is
             else if (status .ne. FIO_SUCCESS) then
                prtcls%flagCon(pp) = 0_is
             end if

             if (.not.params%GC_coords) then
                prtcls%E(pp,1)=(Etmp(1)*cos(x(2))-Etmp(2)*sin(x(2)))/ &
                     params%cpp%Eo
                prtcls%E(pp,2)=(Etmp(1)*sin(x(2))+Etmp(2)*cos(x(2)))/ &
                     params%cpp%Eo
                prtcls%E(pp,3)=Etmp(3)/params%cpp%Eo
             else
                prtcls%E(pp,1)=Etmp(1)/params%cpp%Eo
                prtcls%E(pp,2)=Etmp(2)/params%cpp%Eo
                prtcls%E(pp,3)=Etmp(3)/params%cpp%Eo
             end if

          end if
       end do
       !$OMP END PARALLEL DO
    end if
  end subroutine get_fio_electric_fields

  subroutine get_fio_FOelectric_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       E_X,E_Y,E_Z,flag,hint)
    TYPE(FIELDS), INTENT(IN)       :: F
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    REAL(rp), DIMENSION(params%pchunk), INTENT(IN)  :: Y_R,Y_PHI,Y_Z
    REAL(rp), DIMENSION(params%pchunk), INTENT(OUT)  :: E_X,E_Y,E_Z
    INTEGER(is), DIMENSION(params%pchunk), INTENT(INOUT)  :: flag
    TYPE(C_PTR), DIMENSION(params%pchunk), INTENT(INOUT)  :: hint
    INTEGER (C_INT)                :: status
    INTEGER                        :: pp,pchunk
    REAL(rp), DIMENSION(3)         :: x
    REAL(rp), DIMENSION(3)         :: Etmp

    pchunk=params%pchunk

    do pp = 1,pchunk
       if (flag(pp) .EQ. 1_is) then
          x(1) = Y_R(pp)*params%cpp%length
          x(2) = Y_PHI(pp)
          x(3) = Y_Z(pp)*params%cpp%length

          !             prtcls%hint(pp)=c_null_ptr


          status = fio_eval_field(F%FIO_E, x(1),                      &
               Etmp(1),hint(pp))

          if (status .eq. FIO_SUCCESS) then
             E_X(pp)=(Etmp(1)*cos(x(2))-Etmp(2)*sin(x(2)))/ &
                  params%cpp%Eo
             E_Y(pp)=(Etmp(1)*sin(x(2))+Etmp(2)*cos(x(2)))/ &
                  params%cpp%Eo
             E_Z(pp)=Etmp(3)/params%cpp%Eo
          else if (status .eq. FIO_NO_DATA) then
             E_X(pp) = 0
             E_Y(pp) = 0
             E_Z(pp) = 0
             flag(pp) = 0_is
          else if (status .ne. FIO_SUCCESS) then
             flag(pp) = 0_is
          end if

       end if
    end do

  end subroutine get_fio_FOelectric_fields_p

  subroutine get_fio_GCelectric_fields_p(params,F,Y_R,Y_PHI,Y_Z, &
       E_R,E_PHI,E_Z,flag,hint)
    TYPE(FIELDS), INTENT(IN)       :: F
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    REAL(rp), DIMENSION(params%pchunk), INTENT(IN)  :: Y_R,Y_PHI,Y_Z
    REAL(rp), DIMENSION(params%pchunk), INTENT(OUT)  :: E_R,E_PHI,E_Z
    INTEGER(is), DIMENSION(params%pchunk), INTENT(INOUT)  :: flag
    TYPE(C_PTR), DIMENSION(params%pchunk), INTENT(INOUT)  :: hint
    INTEGER (C_INT)                :: status
    INTEGER                        :: pp,pchunk
    REAL(rp), DIMENSION(3)         :: x
    REAL(rp), DIMENSION(3)         :: Etmp

    pchunk=params%pchunk

    do pp = 1,pchunk
       if (flag(pp) .EQ. 1_is) then
          x(1) = Y_R(pp)*params%cpp%length
          x(2) = Y_PHI(pp)
          x(3) = Y_Z(pp)*params%cpp%length

          status = fio_eval_field(F%FIO_E, x(1),                      &
               Etmp(1),hint(pp))

          if (status .eq. FIO_SUCCESS) then
             E_R(pp)=Etmp(1)/params%cpp%Eo
             E_PHI(pp)=Etmp(2)/params%cpp%Eo
             E_Z(pp)=Etmp(3)/params%cpp%Eo
          else if (status .eq. FIO_NO_DATA) then
             E_R(pp) = 0
             E_PHI(pp) = 0
             E_Z(pp) = 0
             flag(pp) = 0_is
          else if (status .ne. FIO_SUCCESS) then
             flag(pp) = 0_is
          end if

!          write(6,*) E_R,E_PHI,E_Z

       end if
    end do

  end subroutine get_fio_GCelectric_fields_p

  subroutine get_fio_profile(prtcls, P, params)
    TYPE(PARTICLES), INTENT(INOUT) :: prtcls
    TYPE(PROFILES), INTENT(IN)     :: P
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    INTEGER (C_INT)                :: status
    INTEGER                        :: pp
    REAL(rp), DIMENSION(3)         :: x
    integer(ip)  ::  ss
    REAL(rp)        :: netmp=-1._rp,Tetmp=-1._rp

    if (size(prtcls%Y,1).eq.1) then
       ss = size(prtcls%Y,1)
    else
       if (prtcls%Y(2,1).eq.0) then
          ss=1_idef
       else
          ss = size(prtcls%Y,1)
       end if
    endif

    if (prtcls%cart) then
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(pp,status,x)
       do pp = 1,ss
          if (prtcls%flagCon(pp) .EQ. 1_is) then
             x = prtcls%X(pp,:)*params%cpp%length
             status = fio_eval_field(P%FIO_ne, x(1),                     &
                  netmp,                         &
                  prtcls%hint(pp))

             if (status .eq. FIO_SUCCESS) then

                if(netmp.le.0) netmp=1._rp

                prtcls%ne(pp) = netmp/params%cpp%density

             else if (status .eq. FIO_NO_DATA) then
                prtcls%ne(pp) = 0
             else if (status .ne. FIO_SUCCESS) then
                prtcls%flagCon(pp) = 0_is
                CYCLE
             end if

             status = fio_eval_field(P%FIO_te, x(1),                     &
                  Tetmp,                         &
                  prtcls%hint(pp))

             if (status .eq. FIO_SUCCESS) then
                if(Tetmp.le.0) Tetmp=0.1_rp

                prtcls%te(pp) = Tetmp/(params%cpp%temperature/C_E)
             else if (status .eq. FIO_NO_DATA) then
                prtcls%te(pp) = 0
             end if

             status = fio_eval_field(P%FIO_zeff, x(1), &
                  prtcls%Zeff(pp),prtcls%hint(pp))

             if (status .eq. FIO_NO_DATA) then
                prtcls%Zeff(pp) = 1
             end if
          end if
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(pp,status,x)
       do pp = 1,ss
          if (prtcls%flagCon(pp) .EQ. 1_is) then
             x(1) = prtcls%Y(1,pp)*params%cpp%length
             x(2) = prtcls%Y(2,pp)
             x(3) = prtcls%Y(3,pp)*params%cpp%length
             status = fio_eval_field(P%FIO_ne, x(1), &
                  netmp,prtcls%hint(pp))

             if (status .eq. FIO_SUCCESS) then

                if(netmp.le.0) netmp=1._rp

                prtcls%ne(pp) = netmp/params%cpp%density

             else if (status .eq. FIO_NO_DATA) then
                prtcls%ne(pp) = 0
             else if (status .ne. FIO_SUCCESS) then
                prtcls%flagCon(pp) = 0_is
                CYCLE
             end if

             status = fio_eval_field(P%FIO_te, x(1),&
                  tetmp,prtcls%hint(pp))

             if (status .eq. FIO_SUCCESS) then
                if(Tetmp.le.0) Tetmp=0.1_rp

                prtcls%te(pp) = Tetmp/(params%cpp%temperature/C_E)
             else if (status .eq. FIO_NO_DATA) then
                prtcls%te(pp) = 0
             end if

             status = fio_eval_field(P%FIO_zeff, x(1), &
                  prtcls%Zeff(pp),prtcls%hint(pp))

             if (status .eq. FIO_NO_DATA) then
                prtcls%Zeff(pp) = 1
             end if
          end if
       end do
       !$OMP END PARALLEL DO
    end if
  end subroutine get_fio_profile

  subroutine get_fio_profile_p(params,P,Y_R,Y_PHI,Y_Z, &
       n_e,T_e,flag,hint)
    TYPE(PROFILES), INTENT(IN)       :: P
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    REAL(rp), DIMENSION(params%pchunk), INTENT(IN)  :: Y_R,Y_PHI,Y_Z
    REAL(rp), DIMENSION(params%pchunk), INTENT(INOUT)  :: n_e,T_e
    INTEGER(is), DIMENSION(params%pchunk), INTENT(INOUT)  :: flag
    TYPE(C_PTR), DIMENSION(params%pchunk), INTENT(INOUT)  :: hint
    INTEGER (C_INT)                :: status
    INTEGER                        :: pp,pchunk
    REAL(rp), DIMENSION(3)         :: x
    REAL(rp)        :: netmp=-1._rp,Tetmp=-1._rp

    pchunk=params%pchunk

    do pp = 1,pchunk
       if (flag(pp) .EQ. 1_is) then
          x(1) = Y_R(pp)*params%cpp%length
          x(2) = Y_PHI(pp)
          x(3) = Y_Z(pp)*params%cpp%length

          !write(6,*) P%FIO_ne,x

          status = fio_eval_field(P%FIO_ne, x(1), &
               netmp,hint(pp))

          if (status .eq. FIO_SUCCESS) then

             if(netmp.le.0) netmp=1._rp

             n_e(pp) = netmp/params%cpp%density
          else if (status .eq. FIO_NO_DATA) then
             flag(pp) = 0_is
          else if (status .ne. FIO_SUCCESS) then
             flag(pp) = 0_is
             CYCLE
          end if

          status = fio_eval_field(P%FIO_te, x(1), &
               Tetmp,hint(pp))

          if (status .eq. FIO_SUCCESS) then

             if(Tetmp.le.0) Tetmp=0.1_rp

             T_e(pp) = Tetmp/(params%cpp%temperature/C_E)
          end if

!          write(6,*) E_R,E_PHI,E_Z

       end if
    end do

  end subroutine get_fio_profile_p

  subroutine get_fio_ion_p(params,P,Y_R,Y_PHI,Y_Z, &
       n_e,n_i,nimp,Zeff,flag,hint)
    TYPE(PROFILES), INTENT(IN)       :: P
    TYPE(KORC_PARAMS), INTENT(IN)  :: params
    REAL(rp), DIMENSION(params%pchunk), INTENT(IN)  :: Y_R,Y_PHI,Y_Z
    REAL(rp), DIMENSION(params%pchunk,params%num_impurity_species),&
         & INTENT(INOUT)  :: nimp
    REAL(rp), DIMENSION(params%pchunk), INTENT(INOUT)  :: n_i,Zeff
    REAL(rp), DIMENSION(params%pchunk), INTENT(IN)  :: n_e
    INTEGER(is), DIMENSION(params%pchunk), INTENT(INOUT)  :: flag
    TYPE(C_PTR), DIMENSION(params%pchunk), INTENT(INOUT)  :: hint
    INTEGER (C_INT)                :: status
    INTEGER                        :: ii,pp,pchunk,num_imp
    REAL(rp), DIMENSION(3)         :: x
    REAL(rp)        :: nimptmp=-1._rp,nitmp=-1._rp

    pchunk=params%pchunk
    num_imp=params%num_impurity_species

    do pp = 1,pchunk
       if (flag(pp) .EQ. 1_is) then
          x(1) = Y_R(pp)*params%cpp%length
          x(2) = modulo(Y_PHI(pp),2*C_PI)
          x(3) = Y_Z(pp)*params%cpp%length

          !write(6,*) 'X',x

          do ii = 1,num_imp
             status = fio_eval_field(P%FIO_nimp(ii), x(1), &
               nimptmp,hint(pp))

!             write(6,*) P%FIO_nimp(ii)
!             write(6,*) 'nimp_',ii,nimptmp

             if (status .eq. FIO_SUCCESS) then

                if(nimptmp.le.0) nimptmp=1._rp

                nimp(pp,ii) = nimptmp/params%cpp%density
             else if (status .eq. FIO_NO_DATA) then
                flag(pp) = 0_is
             else if (status .ne. FIO_SUCCESS) then
                flag(pp) = 0_is
                CYCLE
             end if
          end do

          status = fio_eval_field(P%FIO_ni, x(1), &
               nitmp,hint(pp))

          if (status .eq. FIO_SUCCESS) then
             n_i(pp) = nitmp/params%cpp%density
          end if

          Zeff(pp)=n_i(pp)
          do ii=1,params%num_impurity_species
             Zeff(pp)=Zeff(pp)+nimp(pp,ii)*params%Zj(ii)**2
          end do
          Zeff(pp)=Zeff(pp)/n_e(pp)

          if ((Zeff(pp).gt.100).or.(Zeff(pp).lt.0.01).or. &
               (isnan(Zeff(pp)))) then
             Zeff(pp)=1._rp
             !write(6,*) 'Zeff is ',Zeff(pp),'!!'
             !write(6,*) 'R:',x(1)
             !write(6,*) 'PHI:',x(2)
             !write(6,*) 'Z:',x(3)
             !do ii = 1,num_imp
             !   write(6,*) 'n_imp_',ii,':',nimp(pp,ii)
             !end do
             !write(6,*) 'n_e:',n_e(pp)
             !write(6,*) 'n_i:',n_i(pp)
             !call korc_abort(21)
          endif

       end if
    end do

  end subroutine get_fio_ion_p

#endif


end module korc_interp
