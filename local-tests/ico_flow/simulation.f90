module simulation_module
  use ico_base_module
  use hypre_module
  use pressure_module
  use velocity_module
  use timestep_module
  use io_module, only: check_io
  implicit none

  private

  public :: step_simulation

contains
  ! --
  ! Advance the solution from the n timestep to n+1
  !--
  ! flux_n is H^n-1 on input
  ! flux_o will become H^n-1 and H^n will be computed and stored in flux_n
  ! field_n is the solution at the n timestep with physical boundary conditions
  !   properly set and internal communication boundaries properly set
  ! field_n will become the solution at the n+1 level
  ! field_o will become the solution at the n level
  ! mtm_rhs eliminates the need for memory allocations
  subroutine step_simulation (field_n, field_o, flux_n, flux_o, &
       mtm_rhs, geom)
    implicit none
    type(amrex_multifab), intent(inout) :: field_n, field_o, flux_n, flux_o, &
         mtm_rhs
    type(amrex_geometry), intent(in) :: geom
    !-
    real(amrex_real) :: dx, dt, time, umax, umin, vmax, vmin, pmax, pmin
    integer :: step

    ! For simplicity, we assume dx(1) == dx(2) == dx(3)
    dx = geom%dx(1)

    call predict_timestep(field_n, dx, time, dt, step)

    if ( amrex_parallel_IOProcessor() ) then
       print*,'Advancing time step',step,'with dt=',dt, 'to time=', time
    end if

    if (step.eq.0) then
       call check_io(0, 0.0d0, 0.0d0, field_n, geom)
       ! store flux in flux_n so the swap stores the correct flux in flux_o
       call velocity_conv_flux_u(field_n, flux_n, geom)
       call velocity_conv_flux_v(field_n, flux_n, geom)
    end if

    ! convective fluxes at n and n-1 (o)
    call amrex_multifab_swap(flux_n, flux_o)
    call velocity_conv_flux_u(field_n, flux_n, geom)
    call velocity_conv_flux_v(field_n, flux_n, geom)

    ! store current field in _o and update velocity components of field_n to *
    call amrex_multifab_swap(field_n, field_o)

    call velocity_rhs_u(field_o, flux_n, flux_o, mtm_rhs, geom, dt)
    call velocity_solve_u(field_n, field_o, mtm_rhs, geom, dt)

    call velocity_rhs_v(field_o, flux_n, flux_o, mtm_rhs, geom, dt)
    call velocity_solve_v(field_n, field_o, mtm_rhs, geom, dt)

    ! ensure solution is consistent across procs
    call velocity_bc(field_n, geom, time+dt)
    call field_n%fill_boundary(geom)

    call pressure_rhs(field_n, geom, dt)
    call pressure_solve(field_n, field_o, geom)

    ! ensure pressure solution is consistent across procs
    call field_n%fill_boundary(geom)

    call pressure_project(field_n, geom, dt)
    ! ensure velocity solution is consistent across procs
    call velocity_bc(field_n, geom, time+dt)
    call field_n%fill_boundary(geom)

    umin = field_n%min(U_i)
    vmin = field_n%min(V_i)
    pmin = field_n%min(P_i)
    umax = field_n%max(U_i)
    vmax = field_n%max(V_i)
    pmax = field_n%max(P_i)

    if ( amrex_parallel_IOProcessor() ) then
       write(*,'(6(es13.5))') umin, umax, vmin, vmax, pmin, pmax
    end if


    call check_io(step+1, dt, time+dt, field_n, geom)
    call update_timestep()

  end subroutine step_simulation

end module simulation_module
