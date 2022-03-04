!> @brief Module containing physical and mathematical constants to be used in KORC.
!! @details In this module we define the constant parameters to be used in KORC. Notice that the numerical precision of these quantities is '_rp', see korc_types.f90.
!! Any new constant needs to be compliant with the numerical precision used in KORC.
module korc_constants
  USE korc_types

  IMPLICIT NONE

  REAL(rp), PARAMETER :: C_PI = 4.0_rp*ATAN(1.0_rp) !< Definition of @f$\pi@f$
  REAL(rp), PARAMETER :: C_E = 1.602176E-19_rp !< Absolute value of electron charge in Coulombs (C).
  REAL(rp), PARAMETER :: C_ME = 9.109382E-31_rp !< Electron mass in kg
  REAL(rp), PARAMETER :: C_MP = 1.672621E-27_rp !< Proton mass in kg
  REAL(rp), PARAMETER :: C_U = 1.660538E-27_rp !< Atomic mass unit in kg
  REAL(rp), PARAMETER :: C_KB = 1.380650E-23_rp !< Boltzmann constant in Joules/Kelvin
  REAL(rp), PARAMETER :: C_C = 299792458.0_rp !< Light speed in m/s
  REAL(rp), PARAMETER :: C_MU = 4.0_rp*C_PI*1E-7_rp !< Vacuum permeability in N/A^2
  REAL(rp), PARAMETER :: C_E0 = 1.0_rp/(C_MU*C_C**2) !< Vacuum permittivity in C^2/(N*m^2)
  REAL(rp), PARAMETER :: C_Ke = 1.0_rp/(4.0_rp*C_PI*C_E0) !< Coulomb constant in N*m^2/C^2
  REAL(rp), PARAMETER :: C_RE = C_E**2/( 4.0_rp*C_PI*C_E0*C_ME*C_C**2 ) !< Classical electron radius
  REAL(rp), PARAMETER :: C_h = 6.6261E-34_rp !< Planck constant in Joules*s
  REAL(rp), PARAMETER :: C_a = 1._rp/137._rp !< Fine-structure constant
  
  
end module korc_constants
