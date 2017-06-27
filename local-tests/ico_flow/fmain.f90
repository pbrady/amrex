subroutine amrex_fmain () bind(c)
  use amrex_base_module
  use pressure_module, only : init_pressure
  use velocity_module, only : init_velocity
  use io_module, only : init_io
  use timestep_module, only : init_timestep, done_timestep
  use simulation_module, only : step_simulation

  implicit none

  integer :: n_cell, max_grid_size
  integer, parameter :: ncomp = 2, nghost = 1  ! two components per multifab, one ghost
  type(amrex_parmparse) :: pp
  type(amrex_box) :: domain
  type(amrex_boxarray)  :: ba
  type(amrex_distromap) :: dm
  type(amrex_geometry)  :: geom
  type(amrex_multifab)  :: U, V, P

  ! amrex_parmparse is way of reading inputs from the inputs file
  ! "get" means it must be set in the inputs file, whereas
  ! "query" means it may not may not be in the inputs file
  call amrex_parmparse_build(pp)

  call pp%get("n_cell", n_cell)  ! # of cells in each dimension

  max_grid_size = 32   ! default max grid size
  call pp%query("max_grid_size", max_grid_size)

  call amrex_parmparse_destroy(pp)

  ! Define a single box covering the domain
  domain = amrex_box([0,0,0], [n_cell-1, n_cell-1, n_cell-1])

  ! Initialize the boxarrays
  call amrex_boxarray_build(ba, domain)

  ! Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
  call ba%maxSize(max_grid_size)

  ! Build a DistributionMapping for the boxarray
  call amrex_distromap_build(dm, ba)

  ! This defines a amrex_geometry object.
  call amrex_geometry_build(geom, domain)

  ! multifabs
  call amrex_multifab_build(P, ba, dm, ncomp, nghost, [.false., .false.])
  call amrex_multifab_build(U, ba, dm, ncomp, nghost, [.true., .false.])
  call amrex_multifab_build(V, ba, dm, ncomp, nghost, [.false., .true.])

  ! initialize io multifab and module
  call init_io(ba, dm, 3)

  ! Intialize data
  call init_pressure(P)
  call init_velocity(U, V)

  call init_timestep()

  ! some cleanup before main running loop
  call amrex_distromap_destroy(dm)
  call amrex_boxarray_destroy(ba)

  do while(.not.done_timestep())
     call step_simulation(P, U, V, geom)
     exit
  end do

end subroutine amrex_fmain
