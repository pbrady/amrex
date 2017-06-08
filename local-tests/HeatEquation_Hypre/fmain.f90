
subroutine amrex_fmain () bind(c)

  use amrex_base_module
  use init_phi_module, only : init_phi
  use advance_module, only : advance, writeplotfile, lo_bc, hi_bc !, init_lu

  implicit none

  integer :: n_cell, max_grid_size, nsteps, plot_int
  integer, parameter :: ncomp = 1, nghost = 1  ! one component, one ghost
  integer :: istep, plot_step
  real(amrex_real) :: dt, time
  type(amrex_parmparse) :: pp
  type(amrex_box) :: domain
  type(amrex_boxarray)  :: ba
  type(amrex_distromap) :: dm
  type(amrex_geometry)  :: geom
  type(amrex_multifab)  :: new_phi, old_phi
  character(len=127) :: plot_file  = "plt"

  ! amrex_parmparse is way of reading inputs from the inputs file
  ! "get" means it must be set in the inputs file, whereas
  ! "query" means it may not may not be in the inputs file
  call amrex_parmparse_build(pp)

  call pp%get("n_cell", n_cell);  ! # of cells in each dimension
  call pp%get("nsteps", nsteps) ! # of steps

  max_grid_size = 32   ! default max grid size
  call pp%query("max_grid_size", max_grid_size);

  plot_int = -1 ! default to no plotfiles
  call pp%query("plot_int", plot_int);
  call pp%query("plot_file", plot_file);

  call amrex_parmparse_destroy(pp)

  ! Define a single box covering the domain
  domain = amrex_box((/0,0,0/), (/n_cell-1, n_cell-1, n_cell-1/))

  ! Initialize the boxarray "ba" from the single box "bx"
  call amrex_boxarray_build(ba, domain)

  ! Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
  call ba%maxSize(max_grid_size)

  ! Build a DistributionMapping for the boxarray
  call amrex_distromap_build(dm, ba)

  ! This defines a amrex_geometry object.
  call amrex_geometry_build(geom, domain)

  ! Build data multifabs
  call amrex_multifab_build(new_phi, ba, dm, ncomp, nghost)
  call amrex_multifab_build(old_phi, ba, dm, ncomp, nghost)

  call amrex_distromap_destroy(dm)
  call amrex_boxarray_destroy(ba)

  ! Intialize data
  call init_phi(new_phi, geom)
  istep = 0
  plot_step = 0
  time = 0.d0
  print *, 'dx', geom%dx
  ! initialize boundary conditions
  lo_bc = amrex_bc_reflect_odd
  hi_bc = amrex_bc_reflect_odd

  ! write initialdata
  call writeplotfile(plot_file, new_phi, geom, time, plot_step)

  ! choose a time step with a diffusive CFL of 2.0
  dt = 2.0d0*geom%dx(1)**2/(2.d0*amrex_spacedim)
  !dt = 1.0d0

  do istep = 1, nsteps

     if ( amrex_parallel_IOProcessor() ) then
        print*,'Advancing time step',istep,'with dt=',dt
     end if

     ! Swap the guts of multifabs so we don't have to allocate and de-allocate data
     call amrex_multifab_swap(new_phi, old_phi)

     ! advance phi
     call advance(old_phi, new_phi, geom, time, dt)

     time = time + dt

     if (plot_int.gt.0 .and. mod(istep, plot_int).eq.0) then
        plot_step = plot_step+1
        call writeplotfile(plot_file, new_phi, geom, time, plot_step)
     endif
  end do

  call amrex_multifab_destroy(new_phi)
  call amrex_multifab_destroy(old_phi)

  call amrex_geometry_destroy(geom)

end subroutine amrex_fmain
