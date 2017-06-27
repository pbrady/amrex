module simulation_module

  use amrex_base_module
  use hypre_module
  use pressure_module
  use velocity_module
  use timestep_module
  use io_module, only: check_io
  implicit none

  private

  public :: step_simulation

contains


  subroutine step_simulation (P, U, V, geom)
    implicit none
    type(amrex_multifab), intent(inout) :: P, U, V
    type(amrex_geometry), intent(in) :: geom
    !-
    real(amrex_real) :: dx, dt, time
    integer :: step
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: dp
    integer :: i,j

    ! For simplicity, we assume dx(1) == dx(2) == dx(3)
    dx = geom%dx(1)

    call predict_timestep(U, V, time, dt, step)

    if ( amrex_parallel_IOProcessor() ) then
       print*,'Advancing time step',step,'with dt=',dt
    end if

    if (step.eq.0) &
         call check_io(step, dt, time, P, U, V, geom)
!!$
!!$
!!$    ! This fills periodic ghost cells and ghost cells from neighboring grids
!!$    call old_phi%fill_boundary(geom)
!!$    ! update physcial boundary conditions on new_phi for implicit solve
!!$    call fill_physbc(new_phi, geom, time+dt)
!!$
!!$    call hypre_update(old_phi, new_phi, geom, 0.5d0*dt)
!!$    call hypre_solve(new_phi)
!!$
!!$    ! transform wall values of new_phi into ghost values based on solution
!!$    call amrex_mfiter_build(mfi, new_phi, tiling=.true.)
!!$
!!$    do while (mfi%next())
!!$
!!$       bx = mfi%tilebox()
!!$
!!$       p => new_phi%dataptr(mfi)
!!$
!!$       if (bx%lo(1).eq.geom%domain%lo(1)) then ! xmin wall
!!$          i = bx%lo(1)-1
!!$          do j=bx%lo(2),bx%hi(2)
!!$             !write(*,*) i,j,p(i,j,1,1),p(i+1,j,1,1)
!!$             p(i,j,1,1) = 2.0d0*p(i,j,1,1)-p(i+1,j,1,1)
!!$          end do
!!$       end if
!!$       if (bx%hi(1).eq.geom%domain%hi(1)) then ! xmax wall
!!$          i = bx%hi(1)+1
!!$          do j=bx%lo(2),bx%hi(2)
!!$             p(i,j,1,1) = 2.0d0*p(i,j,1,1)-p(i-1,j,1,1)
!!$          end do
!!$       end if
!!$       if (bx%lo(2).eq.geom%domain%lo(2)) then ! ymin wall
!!$          j = bx%lo(2)-1
!!$          do i=bx%lo(1),bx%hi(1)
!!$             p(i,j,1,1) = 2.0d0*p(i,j,1,1)-p(i,j+1,1,1)
!!$          end do
!!$       end if
!!$       if (bx%hi(2).eq.geom%domain%hi(2)) then ! ymax wall
!!$          j = bx%hi(2)+1
!!$          do i=bx%lo(1),bx%hi(1)
!!$             p(i,j,1,1) = 2.0d0*p(i,j,1,1)-p(i,j-1,1,1)
!!$          end do
!!$       end if
!!$
!!$    end do
!!$
!!$    call amrex_mfiter_destroy(mfi)

  end subroutine step_simulation



!!$  subroutine fill_physbc (mf, geom, time)
!!$    use amrex_filcc_module, only : amrex_filcc
!!$    type(amrex_geometry) :: geom
!!$    type(amrex_multifab) :: mf
!!$    real(amrex_real), intent(in) :: time
!!$    !-
!!$
!!$    type(amrex_mfiter) :: mfi
!!$    type(amrex_box) :: bx
!!$    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: p
!!$    integer :: i,j
!!$    real(amrex_real) :: x,y,lx,ly,xm,xp,ym,yp,t
!!$
!!$    xm = amrex_problo(1)
!!$    xp = amrex_probhi(1)
!!$    ym = amrex_problo(2)
!!$    yp = amrex_probhi(2)
!!$    lx = xp-xm
!!$    ly = yp-ym
!!$    t = time*tpi
!!$
!!$    !$omp parallel private(mfi,p,plo,phi)
!!$    call amrex_mfiter_build(mfi, mf, tiling=.true.)
!!$    do while(mfi%next())
!!$
!!$       bx = mfi%tilebox()
!!$       p => mf%dataptr(mfi)
!!$
!!$       if (bx%lo(1).eq.geom%domain%lo(1)) then ! xmin wall
!!$          ! phi = sin(t) * sin(2pi*y/ly)
!!$          i = bx%lo(1)-1
!!$          do j=bx%lo(2),bx%hi(2)
!!$             y = amrex_problo(2)+(dble(j)+0.5d0)*geom%dx(2)
!!$             p(i,j,1,1) = sin(t) * sin(tpi*(y-ym)/ly)
!!$          end do
!!$       end if
!!$       if (bx%hi(1).eq.geom%domain%hi(1)) then ! xmax wall
!!$          ! phi = sin(t) * sin(4pi*(y)/ly)
!!$          i = bx%hi(1)+1
!!$          do j=bx%lo(2),bx%hi(2)
!!$             y = amrex_problo(2)+(dble(j)+0.5d0)*geom%dx(2)
!!$             p(i,j,1,1) = sin(t) * sin(2.0d0*tpi*(y-ym)/ly)
!!$          end do
!!$       end if
!!$       if (bx%lo(2).eq.geom%domain%lo(2)) then ! ymin wall
!!$          ! phi = sin(t) * sin(3pi/2 * x/lx)
!!$          j = bx%lo(2)-1
!!$          do i=bx%lo(1),bx%hi(1)
!!$             x = amrex_problo(1)+(dble(i)+0.5d0)*geom%dx(1)
!!$             p(i,j,1,1) = sin(t) * sin(0.75d0*tpi*(x-xm)/lx)
!!$          end do
!!$       end if
!!$       if (bx%hi(2).eq.geom%domain%hi(2)) then ! ymax wall
!!$          ! phi = - sin(t) * sin(3pi/2 * x/lx)
!!$          j = bx%hi(2)+1
!!$          do i=bx%lo(1),bx%hi(1)
!!$             x = amrex_problo(1)+(dble(i)+0.5d0)*geom%dx(1)
!!$             p(i,j,1,1) = sin(t) * sin(0.75d0*tpi*(x-xm)/lx)
!!$          end do
!!$       end if
!!$
!!$    end do
!!$    !$omp end parallel
!!$
!!$  end subroutine fill_physbc

end module simulation_module
