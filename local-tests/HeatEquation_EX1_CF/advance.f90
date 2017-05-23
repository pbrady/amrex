
module advance_module

  use amrex_base_module

  implicit none

  private

  public :: advance, writeplotfile !, init_lu

  integer, public :: lo_bc(amrex_spacedim,1), hi_bc(amrex_spacedim,1)

  real(amrex_real), allocatable :: ddxlu(:,:), ddylu(:,:)
  integer, allocatable :: ipivx(:), ipivy(:)

contains

!!$  subroutine init_lu(geom, dt)
!!$    type(amrex_geometry), intent(in) :: geom
!!$    integer :: nx, ny
!!$
!!$    associate( d => geom%domain )
!!$      nx = d%hi(1)-d%lo(1)+1
!!$      ny = d%hi(2)-d%lo(2)+1
!!$    end associate
!!$
!!$    call tridiag_lufactor(dt, nx, ddxlu, ipivx)
!!$    call tridiag_lufactor(dt, ny, ddylu, ipivy)
!!$
!!$  contains
!!$
!!$    subroutine tridiag_lufactor(dt, nx, lu, ipiv)
!!$      real(amrex_real), intent(in) :: dt
!!$      integer, intent(in) :: nx
!!$      real(amrex_real), allocatable, intent(out) :: lu(:,:)
!!$      integer, allocatable, intent(out) :: ipiv(:)
!!$      !
!!$      integer :: ldab, ku, kl, ldab, info
!!$
!!$      ! tridiagonal
!!$      ku = 1
!!$      kl = 1
!!$      ldab = 2*kl+ku+1
!!$      allocate(lu(ldab,nx))
!!$      allocate(ipiv(nx))
!!$
!!$      lu(2,2:nx) = -dt/2.0d0
!!$      lu(3,:) = 1.0d0+dt
!!$      lu(4,1:nx-1) = -dt/2.0d0
!!$
!!$      call dgbtrf(nx, nx, kl, ku, lu, ldab, ipiv, info)
!!$      if (info.ne.0) then
!!$         write(2, *) "dgbtrf failed for lu with info: ", info
!!$         stop
!!$      end if
!!$
!!$    end subroutine tridiag_lufactor
!!$
!!$  end subroutine init_lu

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


  subroutine advance (old_phi, new_phi, geom, dt)
    type(amrex_multifab), intent(inout) :: old_phi, new_phi
    type(amrex_geometry), intent(in) :: geom
    real(amrex_real), intent(in) :: dt

    integer :: plo(4), phi(4)
    real(amrex_real) :: dx
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: po, pn

    ! For simplicity, we assume dx(1) == dx(2) == dx(3)
    dx = geom%dx(1)

    ! This fills periodic ghost cells and ghost cells from neighboring grids
    call old_phi%fill_boundary(geom)
    ! update physcial boundary conditions
    call fill_physbc(old_phi, geom)

    !$omp parallel private(mfi,bx,po,pn,plo,phi)
    call amrex_mfiter_build(mfi, old_phi, tiling=.true.)

    do while (mfi%next())

       bx = mfi%tilebox()

       po => old_phi%dataptr(mfi)
       pn => new_phi%dataptr(mfi)

       plo = lbound(po)
       phi = ubound(pn)

       select case (amrex_spacedim)
       case (2)
          call update_phi_2d(bx%lo, bx%hi, po, pn, plo, phi, dx, dt)
       case (3)
          call update_phi_3d(bx%lo, bx%hi, po, pn, plo, phi, dx, dt)
       end select
    end do

    call amrex_mfiter_destroy(mfi)
    !$omp end parallel

  end subroutine advance

  subroutine update_phi_2d (lo, hi, pold, pnew, plo, phi, dx, dt)
    integer, intent(in) :: lo(2), hi(2), plo(2), phi(2)
    real(amrex_real), intent(in   ) :: pold(plo(1):phi(1), plo(2):phi(2))
    real(amrex_real), intent(inout) :: pnew(plo(1):phi(1), plo(2):phi(2))
    real(amrex_real), intent(in) :: dx, dt

    integer :: i,j
    real(amrex_real) :: dxinv, dtdx
    real(amrex_real) :: fx(lo(1):hi(1)+1,lo(2):hi(2)  )
    real(amrex_real) :: fy(lo(1):hi(1)  ,lo(2):hi(2)+1)

    dxinv = 1.d0/dx
    dtdx = dt*dxinv

    ! x-fluxes
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          fx(i,j) = ( pold(i-1,j) - pold(i,j) ) * dxinv
       end do
    end do

    ! y-fluxes
    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          fy(i,j) = ( pold(i,j-1) - pold(i,j) ) * dxinv
       end do
    end do

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          pnew(i,j) = pold(i,j) - dtdx * &
               ( fx(i+1,j)-fx(i,j) &
               + fy(i,j+1)-fy(i,j) )

       end do
    end do

  end subroutine update_phi_2d

  subroutine update_phi_3d (lo, hi, pold, pnew, plo, phi, dx, dt)
    integer, intent(in) :: lo(3), hi(3), plo(3), phi(3)
    real(amrex_real), intent(in   ) :: pold(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
    real(amrex_real), intent(inout) :: pnew(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
    real(amrex_real), intent(in) :: dx, dt

    integer :: i,j,k
    real(amrex_real) :: dxinv, dtdx
    real(amrex_real) :: fx(lo(1):hi(1)+1,lo(2):hi(2)  ,lo(3):hi(3))
    real(amrex_real) :: fy(lo(1):hi(1)  ,lo(2):hi(2)+1,lo(3):hi(3))
    real(amrex_real) :: fz(lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)+1)

    dxinv = 1.d0/dx
    dtdx = dt*dxinv

    ! x-fluxes
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             fx(i,j,k) = ( pold(i-1,j,k) - pold(i,j,k) ) * dxinv
          end do
       end do
    end do

    ! y-fluxes
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             fy(i,j,k) = ( pold(i,j-1,k) - pold(i,j,k) ) * dxinv
          end do
       end do
    end do

    ! z-fluxes
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             fz(i,j,k) = ( pold(i,j,k-1) - pold(i,j,k) ) * dxinv
          end do
       end do
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             pnew(i,j,k) = pold(i,j,k) - dtdx * &
                  ( fx(i+1,j,k)-fx(i,j,k) &
                  + fy(i,j+1,k)-fy(i,j,k) &
                  + fz(i,j,k+1)-fz(i,j,k) )

          end do
       end do
    end do

  end subroutine update_phi_3d


  subroutine fill_physbc (mf, geom)
    use amrex_filcc_module, only : amrex_filcc
    type(amrex_geometry) :: geom
    type(amrex_multifab) :: mf

    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: p
    integer :: plo(4), phi(4)

    !$omp parallel private(mfi,p,plo,phi)
    call amrex_mfiter_build(mfi, mf, tiling=.false.)
    do while(mfi%next())
       p => mf%dataptr(mfi)
       if (.not. geom%domain%contains(p)) then ! part of this box is outside the domain
          plo = lbound(p)
          phi = ubound(p)
          call amrex_filcc(p, plo, phi,         & ! fortran array and bounds
               geom%domain%lo, geom%domain%hi,  & ! index extent of whole problem domain
               geom%dx,                         & ! cell size in real
               geom%get_physical_location(plo), & ! physical location of lower left corner
               lo_bc, hi_bc)                      ! bc types for each component

          ! amrex_filcc doesn't fill EXT_DIR (see amrex_bc_types_module for a list of bc types
          ! In that case, the user needs to fill it.
       end if
    end do
    !$omp end parallel

  end subroutine fill_physbc

end module advance_module
