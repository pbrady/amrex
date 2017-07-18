! A single-phase incompressible flow solver based on:
! "Application of a Fractional-Step Method to the Incompressible Navier-Stokes Equations"
! - Kim and Moin 1985
!
! The solver utilizes a staggered mesh (although built on amrex cell centered)
! Velocity is evaluated at the cell faces with pressure at the cell centers.
! The indexing convention is:
!
!               ^ v(i,j+1)
!        -------|-------
!        |             |
!        |             |
!        |             |
! u(i,j) ->   P(i,j)   -> u(i+1,j)
!        |             |
!        |      ^      |
!        |------|------|
!               v(i,j)
!
! A 2nd order Adams-Bashforth scheme is used for the convective terms
! with a 2nd order Crank-Nicholson for the viscous terms
!
! The discretization for the x-momentum can be written as
!
! (1-A_1)(1-A_2)(u*-u^n) = dt/2 (3H^n - 3H^n-1) + dP/dx + 2(A_1+A_2)u^n
!
! where H^k = -d(u^k*u^k)dx - d(u^k*v^k)
! and A_i = dt/(2*Rey) * d^2/dx_i^2
!
! Second order interpolation and derivative stencils are used
!
! The pressure poisson system is
!
! div(grad (P^n+1 - P^n)) = (1/dt) div u*
!
! P^n+1 is constrained to have a zero mean
!
! The projection step is then
!
! u^n+1 = u* - dt * dP^n+1/dx
! v^n+1 = v* - dt * dP^n+1/dy
subroutine amrex_fmain () bind(c)
  use ico_base_module
  use pressure_module, only : init_pressure
  use velocity_module, only : init_velocity
  use io_module, only : init_io, io_write_centerlines
  use timestep_module, only : init_timestep, done_timestep
  use simulation_module, only : step_simulation

  implicit none

  integer :: n_cell, max_grid_size
  integer, parameter :: nghost = 1 ! one ghost cell for second order stencils
  type(amrex_parmparse) :: pp
  type(amrex_box) :: domain
  type(amrex_boxarray)  :: ba
  type(amrex_distromap) :: dm
  type(amrex_geometry)  :: geom
  type(amrex_multifab)  :: vars_new, vars_old ! primitive variables
  type(amrex_multifab)  :: convf_new, convf_old ! convective fluxes
  type(amrex_multifab)  :: mtm_rhs

  ! amrex_parmparse is way of reading inputs from the inputs file
  ! "get" means it must be set in the inputs file, whereas
  ! "query" means it may not may not be in the inputs file
  call amrex_parmparse_build(pp)

  call pp%get("n_cell", n_cell)  ! # of cells in each dimension

  max_grid_size = 32   ! default max grid size
  call pp%query("max_grid_size", max_grid_size)

  call amrex_parmparse_destroy(pp)

  ! Define a single box covering the domain
  domain = amrex_box([0,0], [n_cell-1, n_cell-1])

  ! Initialize the boxarrays
  call amrex_boxarray_build(ba, domain)

  ! Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
  call ba%maxSize(max_grid_size)

  ! Build a DistributionMapping for the boxarray
  call amrex_distromap_build(dm, ba)

  ! This defines a amrex_geometry object.
  call amrex_geometry_build(geom, domain)

  ! multifabs
  call amrex_multifab_build(vars_new, ba, dm, 3, nghost)
  call amrex_multifab_build(vars_old, ba, dm, 3, nghost)
  call amrex_multifab_build(convf_new, ba, dm, 2, nghost)
  call amrex_multifab_build(convf_old, ba, dm, 2, nghost)
  call amrex_multifab_build(mtm_rhs, ba, dm, 2, nghost)

  ! initialize all multifabs to zero
  call vars_new%setval(0.0d0)
  call vars_old%setval(0.0d0)
  call convf_old%setval(0.0d0)
  call convf_new%setval(0.0d0)
  call mtm_rhs%setval(0.0d0)

  ! initialize io multifab and module
  call init_io(ba, dm, 3)

  ! Intialize data
  call init_pressure(vars_new, geom)
  call init_velocity(vars_new, convf_new, geom)

  call init_timestep()

  ! some cleanup before main running loop
  call amrex_distromap_destroy(dm)
  call amrex_boxarray_destroy(ba)

  do while(.not.done_timestep())
     call step_simulation(vars_new, vars_old, convf_new, convf_old, &
          mtm_rhs, geom)
  end do
  call io_write_centerlines(vars_new, geom)

end subroutine amrex_fmain
