module timestep_module
  use ico_base_module
  use velocity_module, only: rey
  implicit none

  private

  public :: init_timestep, predict_timestep, update_timestep, done_timestep

  real(amrex_real) :: dt, time, maxtime
  integer :: step, maxsteps

  real(amrex_real) :: cfl_v, cfl_c

contains

  subroutine init_timestep
    implicit none
    type(amrex_parmparse) :: pp

    cfl_v = 2.0d0
    cfl_c = 0.8d0
    ! linker does not properly zero out module variables
    maxtime = 0.0d0
    time = 0.0d0
    dt = 0.0d0
    step = 0
    maxsteps = 0

    call amrex_parmparse_build(pp, "timestep")
    call pp%query("maxtime", maxtime)
    call pp%query("maxsteps", maxsteps)
    call pp%query("viscous_cfl", cfl_v)
    call pp%query("convection_cfl", cfl_c)
    call amrex_parmparse_destroy(pp)

    if (maxtime.le.0.0d0.and.maxsteps.le.0) &
         stop "We shant run forever"

    if (cfl_v.le.0.0d0.or.cfl_c.le.0.0d0) &
         stop "positive cfl's only"

  end subroutine init_timestep

  subroutine predict_timestep(field, dx, time_n, timestep, currentstep)
    implicit none
    type(amrex_multifab), intent(in) :: field
    real(amrex_real), intent(in) :: dx
    real(amrex_real), intent(out) :: time_n, timestep
    integer, intent(out) :: currentstep
    !-
    real(amrex_real) :: dtc, dtv, umax, vmax

    ! viscous timestep restriction
    dtv = cfl_v * 0.25d0*rey*dx**2
    ! convective timestep restriction
    umax = field%norm0(U_i)
    vmax = field%norm0(V_i)
    dtc = cfl_c * dx/sqrt(umax**2+vmax**2)

    dt = min(dtc, dtv)
    time_n = time
    timestep = dt
    currentstep = step

  end subroutine predict_timestep

  subroutine update_timestep()
    implicit none
    step = step+1
    time = time+dt
  end subroutine update_timestep

  function done_timestep() result (p)
    implicit none
    logical p

    p = .false.
    if (maxsteps.gt.0.and.step.ge.maxsteps) p=.true.
    if (maxtime.gt.0d0.and.time.ge.maxtime) p=.true.

  end function done_timestep

end module timestep_module
