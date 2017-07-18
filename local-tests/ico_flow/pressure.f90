module pressure_module
  use ico_base_module
  use amrex_fi_mpi
  use hypre_module
  implicit none

  private

  public :: init_pressure, pressure_rhs, pressure_solve, pressure_project

contains

  subroutine init_pressure(vars, geom)
    implicit none
    type(amrex_multifab), intent(inout) :: vars
    type(amrex_geometry), intent(in) :: geom
    !-

    ! default initialization of pressure field
    ! handled in fmain

    call init_hypre_pres(vars, geom)
  end subroutine init_pressure

  ! rhs = div(vel)/dt
  ! scaled by -dx**2 so system is SPD
  subroutine pressure_rhs(field, geom, dt)
    implicit none
    type(amrex_multifab), intent(inout) :: field
    type(amrex_geometry), intent(in) :: geom
    real(amrex_real), intent(in) :: dt
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: dp(:,:,:,:), &
         u(:,:), v(:,:), rhs(:,:)
    real(amrex_real) :: dti, dxi
    integer :: i,j,lb(4)

    dti = 1.0d0/dt
    dxi = 1.0d0/geom%dx(1)

    call amrex_mfiter_build(mfi, field)
    do while(mfi%next())
       bx = mfi%tilebox()

       dp => field%dataPtr(mfi)
       lb = lbound(dp)

       u(lb(1):,lb(2):) => dp(:,:,1,U_i)
       v(lb(1):,lb(2):) => dp(:,:,1,V_i)
       ! store rhs of pressure equation in pressure slot of field
       rhs(lb(1):,lb(2):) => dp(:,:,1,P_i)

       do j=bx%lo(2),bx%hi(2)
          do i=bx%lo(1),bx%hi(1)
             rhs(i,j) = dti*dxi*(u(i+1,j)-u(i,j)+v(i,j+1)-v(i,j))
          end do
       end do
    end do
    call amrex_mfiter_destroy(mfi)

  end subroutine pressure_rhs


  ! -dx**2 * div(grad (P^n+1 - P^n)) = rhs_p
  subroutine pressure_solve(field_n, field_o, geom)
    implicit none
    type(amrex_multifab), intent(inout) :: field_n
    type(amrex_multifab), intent(in) :: field_o
    type(amrex_geometry), intent(in) :: geom
    !-

    call hypre_prepare_p(field_n, geom)
    call hypre_solve_p()

    call hypre_update_p(field_n, field_o, geom)
  end subroutine pressure_solve


  ! u^n+1 = u* - dt * dP^n+1/dx
  ! v^n+1 = v* - dt * dP^n+1/dy
  subroutine pressure_project(field_n, geom, dt)
    implicit none
    type(amrex_multifab), intent(inout) :: field_n
    type(amrex_geometry), intent(in) :: geom
    real(amrex_real), intent(in) :: dt
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    integer :: i,j,imin,jmin,lb(4)
    real(amrex_real), contiguous, pointer :: dp(:,:,:,:), &
         pn(:,:), un(:,:), vn(:,:)
    real(amrex_real) :: f

    f = dt/geom%dx(1)
    imin = geom%domain%lo(1)
    jmin = geom%domain%lo(2)

    call amrex_mfiter_build(mfi, field_n)
    do while(mfi%next())

       bx = mfi%tilebox()
       ! need un-uo on boundary
       dp => field_n%dataPtr(mfi)
       lb = lbound(dp)

       pn(lb(1):,lb(2):) => dp(:,:,1,P_i)
       un(lb(1):,lb(2):) => dp(:,:,1,U_i)
       vn(lb(1):,lb(2):) => dp(:,:,1,V_i)

       do j=bx%lo(2),bx%hi(2)
          do i=bx%lo(1),bx%hi(1)

             if (i.gt.imin) un(i,j) = un(i,j)-f*(pn(i,j)-pn(i-1,j))
             if (j.gt.jmin) vn(i,j) = vn(i,j)-f*(pn(i,j)-pn(i,j-1))
          end do
       end do
    end do
    call amrex_mfiter_destroy(mfi)

  end subroutine pressure_project
end module pressure_module
