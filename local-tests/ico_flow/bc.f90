module boundary_condition_module
  use ico_base_module
  implicit none

  private

  public :: velocity_bc

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

end module boundary_condition_module
