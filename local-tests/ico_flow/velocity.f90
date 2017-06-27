module velocity_module

  use amrex_base_module
  use hypre_module
  implicit none

  private

  public :: init_velocity
  protected :: rey

  real(amrex_real) :: rey

contains

  ! hard-coded initial condition of zero velocity
  subroutine init_velocity(U, V)
    implicit none
    type(amrex_multifab), intent(inout) :: U, V
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: dp(:,:,:,:)
    type(amrex_parmparse) :: pp

    call amrex_parmparse_build(pp, "velocity")
    call pp%get("reynolds", rey)
    call amrex_parmparse_destroy(pp)


    print *, "init_u"
    call amrex_mfiter_build(mfi, U)
    do while(mfi%next())
       bx = mfi%tilebox()
       call amrex_print(bx)
       dp => U%dataPtr(mfi)
       dp = 0.0d0
       print *, lbound(dp)
       print *, ubound(dp)
!!$
!!$       dp => U%dataPtr(mfi)
!!$       dp = 0.0d0
!!$       dp => V%dataPtr(mfi)
!!$       dp = 0.0d0
    end do
    call amrex_mfiter_destroy(mfi)


    print *, "init_v"
    call amrex_mfiter_build(mfi, V)
    do while(mfi%next())
       bx = mfi%tilebox()
       call amrex_print(bx)
       dp => V%dataPtr(mfi)
       dp = 0.0d0
       print *, lbound(dp)
       print *, ubound(dp)
!!$
!!$       dp => U%dataPtr(mfi)
!!$       dp = 0.0d0
!!$       dp => V%dataPtr(mfi)
!!$       dp = 0.0d0
    end do
    call amrex_mfiter_destroy(mfi)

    !print *, 'init_hypre_u'
    call init_hypre_u(U)
    !print *, 'init_hypre_v'
    call init_hypre_v(V)

  end subroutine init_velocity

end module velocity_module
