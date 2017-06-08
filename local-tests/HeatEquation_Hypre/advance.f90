
module advance_module

  use amrex_base_module
  use hypre_module
  implicit none

  private

  public :: advance, writeplotfile

  integer, public :: lo_bc(amrex_spacedim,1), hi_bc(amrex_spacedim,1)
  real(amrex_real), parameter :: tpi = 6.283185307179586d0

contains

  subroutine writeplotfile (plot_file, phi, geom, time, step)
    use amrex_base_module

    character(*), intent(in) :: plot_file
    type(amrex_multifab), intent(in) :: phi
    type(amrex_geometry), intent(in) :: geom
    real(amrex_real), intent(in) :: time
    integer, intent(in) :: step
    !
    integer :: nlevs, rr(1), stepno(1)
    character(len=127) :: name
    character(len=16)  :: current_step
    type(amrex_string) :: varname(1)
    type(amrex_multifab) :: phi_(1)
    type(amrex_geometry) :: geom_(1)

    phi_(1) = phi
    geom_(1) = geom
    stepno(1) = step
    nlevs = 1
    rr(1) = 1

    if      (step .lt. 1000000) then
       write(current_step,fmt='(i5.5)') step
    else if (step .lt. 10000000) then
       write(current_step,fmt='(i6.6)') step
    else if (step .lt. 100000000) then
       write(current_step,fmt='(i7.7)') step
    else if (step .lt. 1000000000) then
       write(current_step,fmt='(i8.8)') step
    else
       write(current_step,fmt='(i15.15)') step
    end if
    name = trim(plot_file) // current_step

    call amrex_string_build(varname(1), "phi")

    call amrex_write_plotfile(name, nlevs, phi_, varname, geom_, &
         time, stepno, rr)

  end subroutine writeplotfile


  subroutine advance (old_phi, new_phi, geom, time, dt)
    type(amrex_multifab), intent(inout) :: old_phi, new_phi
    type(amrex_geometry), intent(in) :: geom
    real(amrex_real), intent(in) :: time, dt
    !-
    real(amrex_real) :: dx
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: p
    integer :: i,j

    ! For simplicity, we assume dx(1) == dx(2) == dx(3)
    dx = geom%dx(1)


    ! This fills periodic ghost cells and ghost cells from neighboring grids
    call old_phi%fill_boundary(geom)
    ! update physcial boundary conditions on new_phi for implicit solve
    call fill_physbc(new_phi, geom, time+dt)

    call hypre_update(old_phi, new_phi, geom, 0.5d0*dt)
    call hypre_solve(new_phi)

    ! transform wall values of new_phi into ghost values based on solution
    call amrex_mfiter_build(mfi, new_phi, tiling=.true.)

    do while (mfi%next())

       bx = mfi%tilebox()

       p => new_phi%dataptr(mfi)

       if (bx%lo(1).eq.geom%domain%lo(1)) then ! xmin wall
          i = bx%lo(1)-1
          do j=bx%lo(2),bx%hi(2)
             !write(*,*) i,j,p(i,j,1,1),p(i+1,j,1,1)
             p(i,j,1,1) = 2.0d0*p(i,j,1,1)-p(i+1,j,1,1)
          end do
       end if
       if (bx%hi(1).eq.geom%domain%hi(1)) then ! xmax wall
          i = bx%hi(1)+1
          do j=bx%lo(2),bx%hi(2)
             p(i,j,1,1) = 2.0d0*p(i,j,1,1)-p(i-1,j,1,1)
          end do
       end if
       if (bx%lo(2).eq.geom%domain%lo(2)) then ! ymin wall
          j = bx%lo(2)-1
          do i=bx%lo(1),bx%hi(1)
             p(i,j,1,1) = 2.0d0*p(i,j,1,1)-p(i,j+1,1,1)
          end do
       end if
       if (bx%hi(2).eq.geom%domain%hi(2)) then ! ymax wall
          j = bx%hi(2)+1
          do i=bx%lo(1),bx%hi(1)
             p(i,j,1,1) = 2.0d0*p(i,j,1,1)-p(i,j-1,1,1)
          end do
       end if

    end do

    call amrex_mfiter_destroy(mfi)

  end subroutine advance



  subroutine fill_physbc (mf, geom, time)
    use amrex_filcc_module, only : amrex_filcc
    type(amrex_geometry) :: geom
    type(amrex_multifab) :: mf
    real(amrex_real), intent(in) :: time
    !-

    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: p
    integer :: i,j
    real(amrex_real) :: x,y,lx,ly,xm,xp,ym,yp,t

    xm = amrex_problo(1)
    xp = amrex_probhi(1)
    ym = amrex_problo(2)
    yp = amrex_probhi(2)
    lx = xp-xm
    ly = yp-ym
    t = time*tpi

    !$omp parallel private(mfi,p,plo,phi)
    call amrex_mfiter_build(mfi, mf, tiling=.true.)
    do while(mfi%next())

       bx = mfi%tilebox()
       p => mf%dataptr(mfi)

       if (bx%lo(1).eq.geom%domain%lo(1)) then ! xmin wall
          ! phi = sin(t) * sin(2pi*y/ly)
          i = bx%lo(1)-1
          do j=bx%lo(2),bx%hi(2)
             y = amrex_problo(2)+(dble(j)+0.5d0)*geom%dx(2)
             p(i,j,1,1) = sin(t) * sin(tpi*(y-ym)/ly)
          end do
       end if
       if (bx%hi(1).eq.geom%domain%hi(1)) then ! xmax wall
          ! phi = sin(t) * sin(4pi*(y)/ly)
          i = bx%hi(1)+1
          do j=bx%lo(2),bx%hi(2)
             y = amrex_problo(2)+(dble(j)+0.5d0)*geom%dx(2)
             p(i,j,1,1) = sin(t) * sin(2.0d0*tpi*(y-ym)/ly)
          end do
       end if
       if (bx%lo(2).eq.geom%domain%lo(2)) then ! ymin wall
          ! phi = sin(t) * sin(3pi/2 * x/lx)
          j = bx%lo(2)-1
          do i=bx%lo(1),bx%hi(1)
             x = amrex_problo(1)+(dble(i)+0.5d0)*geom%dx(1)
             p(i,j,1,1) = sin(t) * sin(0.75d0*tpi*(x-xm)/lx)
          end do
       end if
       if (bx%hi(2).eq.geom%domain%hi(2)) then ! ymax wall
          ! phi = - sin(t) * sin(3pi/2 * x/lx)
          j = bx%hi(2)+1
          do i=bx%lo(1),bx%hi(1)
             x = amrex_problo(1)+(dble(i)+0.5d0)*geom%dx(1)
             p(i,j,1,1) = sin(t) * sin(0.75d0*tpi*(x-xm)/lx)
          end do
       end if

    end do
    !$omp end parallel

  end subroutine fill_physbc

end module advance_module
