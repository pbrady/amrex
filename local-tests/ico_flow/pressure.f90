module pressure_module

  use amrex_base_module
  use hypre_module
  implicit none

  private

  public :: init_pressure

contains

  subroutine init_pressure(P)
    implicit none
    type(amrex_multifab), intent(inout) :: P
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: dp(:,:,:,:)


    print *, "init_pressure"
    call amrex_mfiter_build(mfi, P)
    do while(mfi%next())
       bx = mfi%tilebox()
       call amrex_print(bx)

       dp => P%dataPtr(mfi)
       dp = 0.0d0
       print *, lbound(dp)
       print *, ubound(dp)
    end do
    call amrex_mfiter_destroy(mfi)

    call init_hypre_pres(P)

  end subroutine init_pressure

end module pressure_module
