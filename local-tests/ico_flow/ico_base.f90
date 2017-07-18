module ico_base_module
  use amrex_base_module
  implicit none

  public

  ! component indices into variable and flux multifabs
  integer, parameter :: U_i = 1
  integer, parameter :: V_i = 2
  integer, parameter :: P_i = 3

end module ico_base_module
