
module init_phi_module

  use amrex_base_module
  use hypre_module
  implicit none

  private

  public :: init_phi

contains

  subroutine init_phi(phi,geom)

    type(amrex_multifab), intent(inout) :: phi
    type(amrex_geometry), intent(in   ) :: geom

    ! local variables
    integer :: lo(4), hi(4)
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, pointer :: p(:,:,:,:)
    real(amrex_real) :: lx

    lx = amrex_probhi(1)-amrex_problo(1)
    !$omp parallel private(bx,p,lo,hi)
    call amrex_mfiter_build(mfi, phi, tiling=.true.)
    do while(mfi%next())
       bx = mfi%tilebox()
       p => phi%dataPtr(mfi)
       lo = lbound(p)
       hi = ubound(p)
       select case (amrex_spacedim)
       case (2)
          call init_phi_2d(bx%lo(1:2), bx%hi(1:2), p(:,:,1,1), lo(1:2), hi(1:2), &
               amrex_problo(1:2), geom%dx(1:2), lx)
       case (3)
          call init_phi_3d(bx%lo, bx%hi, p(:,:,:,1), lo(1:3), hi(1:3), &
               amrex_problo, geom%dx)
       end select
    end do
    call amrex_mfiter_destroy(mfi)
    !$omp end parallel

    call init_hypre(phi)

  end subroutine init_phi

  subroutine init_phi_2d(lo, hi, phi, dlo, dhi, prob_lo, dx, lx)
    integer          :: lo(2), hi(2), dlo(2), dhi(2)
    real(amrex_real) :: phi(dlo(1):dhi(1),dlo(2):dhi(2))
    real(amrex_real) :: prob_lo(2)
    real(amrex_real) :: dx(2), lx

    ! local varables
    integer          :: i,j
    real(amrex_real) :: x,y,r2

    ! handle boundary conditions
    phi = 0.0d0

    do j=lo(2),hi(2)
       y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

          r2 = ((x)**2 + (y)**2) / 0.1d0
          !r2 = ((x-0.125d0*lx)**2 + (y-0.125d0*lx)**2) / 0.1d0
          phi(i,j) = exp(-r2)

       end do
    end do

  end subroutine init_phi_2d

  subroutine init_phi_3d(lo, hi, phi, dlo, dhi, prob_lo, dx)
    integer          :: lo(3), hi(3), dlo(3), dhi(3)
    real(amrex_real) :: phi(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
    real(amrex_real) :: prob_lo(3)
    real(amrex_real) :: dx(3)

    ! local varables
    integer          :: i,j,k
    real(amrex_real) :: x,y,z,r2

    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

             r2 = ((x-0.25d0)**2 + (y-0.25d0)**2 + (z-0.25d0)**2) / 0.01d0
             phi(i,j,k) = 1.d0 + exp(-r2)

          end do
       end do
    end do

  end subroutine init_phi_3d

end module init_phi_module