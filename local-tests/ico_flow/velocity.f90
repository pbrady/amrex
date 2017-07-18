module velocity_module
  use ico_base_module
  use amrex_fi_mpi
  use hypre_module
  implicit none

  real(amrex_real), protected :: rey

contains

  ! hard-coded initial condition of zero velocity
  subroutine init_velocity(vars, fluxes, geom)
    implicit none
    type(amrex_multifab), intent(inout) :: vars, fluxes
    type(amrex_geometry), intent(in) :: geom
    !-
    type(amrex_parmparse) :: pp

    call amrex_parmparse_build(pp, "velocity")
    call pp%get("reynolds", rey)
    call amrex_parmparse_destroy(pp)

    ! default initialization of zero velocity handled in fmain

    call init_hypre_u(vars, geom)
    call init_hypre_v(vars, geom)

    ! update boundaries
    call velocity_bc(vars, geom, 0.0d0)
    call vars%fill_boundary(geom)

  end subroutine init_velocity

  ! rhs_u = dt/2 (3H^n - 3H^n-1) - dP/dx + 2(A_1+A_2)u^n
  subroutine velocity_rhs_u(field, flux_n, flux_o, mtm_rhs, geom, dt)
    implicit none
    type(amrex_multifab), intent(in) :: field, flux_n, flux_o
    type(amrex_geometry), intent(in) :: geom
    type(amrex_multifab), intent(inout) :: mtm_rhs
    real(amrex_real), intent(in) :: dt
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: dp(:,:,:,:), &
         u(:,:), v(:,:), p(:,:), rhs(:,:), fn(:,:), fo(:,:)
    integer :: i, j, lb(4)
    real(amrex_real) :: dxi, reyi

    dxi = 1.0d0/geom%dx(1)
    reyi = 1.0d0/rey

    call amrex_mfiter_build(mfi, field)
    do while(mfi%next())
       bx = mfi%tilebox()
       ! grab component data
       dp => field%dataPtr(mfi)
       lb = lbound(dp)
       u(lb(1):,lb(2):) => dp(:,:,1,U_i)
       v(lb(1):,lb(2):) => dp(:,:,1,V_i)
       p(lb(1):,lb(2):) => dp(:,:,1,P_i)

       dp => flux_n%dataPtr(mfi)
       fn(lb(1):,lb(2):) => dp(:,:,1,U_i)
       dp => flux_o%dataPtr(mfi)
       fo(lb(1):,lb(2):) => dp(:,:,1,U_i)

       dp => mtm_rhs%dataPtr(mfi)
       rhs(lb(1):,lb(2):) => dp(:,:,1,U_i)

       do j=bx%lo(2),bx%hi(2)
          do i=bx%lo(1),bx%hi(1)
             rhs(i,j) = dt*(&
                  -dxi*(p(i,j)-p(i-1,j)) &
                  +0.5d0*(3*fn(i,j)-fo(i,j)) &
                  +reyi*(dxi**2*(u(i-1,j)-2.0d0*u(i,j)+u(i+1,j))) &
                  +reyi*(dxi**2*(u(i,j-1)-2.0d0*u(i,j)+u(i,j+1))))
          end do
       end do
    end do
    call amrex_mfiter_destroy(mfi)
  end subroutine velocity_rhs_u

  ! (1-A_1)(1-A_2)(u*-u^n) = rhs_u
  subroutine velocity_solve_u(field_n, field_o, mtm_rhs, geom, dt)
    implicit none
    type(amrex_multifab), intent(in) :: field_o, mtm_rhs
    type(amrex_geometry), intent(in) :: geom
    type(amrex_multifab), intent(inout) :: field_n
    real(amrex_real), intent(in) :: dt
    !-
    real(amrex_real) :: f

    f = dt/(2.0d0*rey)

    ! solve:  (1-A_1) X = rhs_u
    call hypre_prepare_u_x(mtm_rhs, geom, f)
    call hypre_solve_u_x()

    ! solve:  (1-A_2)(u*-u^n) = X
    call hypre_prepare_u_y(field_n, geom, f)
    call hypre_solve_u_y()

    ! update:  u* = u^n + (u*-u^n)
    call hypre_update_u(field_n, field_o, geom)
  end subroutine velocity_solve_u

  ! rhs_v = dt/2 (3H^n - 3H^n-1) + dP/dy + 2(A_1+A_2)v^n
  subroutine velocity_rhs_v(field, flux_n, flux_o, mtm_rhs, geom, dt)
    implicit none
    type(amrex_multifab), intent(in) :: field, flux_n, flux_o
    type(amrex_geometry), intent(in) :: geom
    type(amrex_multifab), intent(inout) :: mtm_rhs
    real(amrex_real), intent(in) :: dt
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: dp(:,:,:,:), &
         u(:,:), v(:,:), p(:,:), rhs(:,:), fn(:,:), fo(:,:)
    integer :: i, j, lb(4)
    real(amrex_real) :: dxi, reyi

    dxi = 1.0d0/geom%dx(1)
    reyi = 1.0d0/rey

    call amrex_mfiter_build(mfi, field)
    do while(mfi%next())
       bx = mfi%tilebox()
       ! grab component data
       dp => field%dataPtr(mfi)
       lb = lbound(dp)
       u(lb(1):,lb(2):) => dp(:,:,1,U_i)
       v(lb(1):,lb(2):) => dp(:,:,1,V_i)
       p(lb(1):,lb(2):) => dp(:,:,1,P_i)

       dp => flux_n%dataPtr(mfi)
       fn(lb(1):,lb(2):) => dp(:,:,1,V_i)
       dp => flux_o%dataPtr(mfi)
       fo(lb(1):,lb(2):) => dp(:,:,1,V_i)

       dp => mtm_rhs%dataPtr(mfi)
       rhs(lb(1):,lb(2):) => dp(:,:,1,V_i)

       do j=bx%lo(2),bx%hi(2)
          do i=bx%lo(1),bx%hi(1)
             rhs(i,j) = dt*(&
                  -dxi*(p(i,j)-p(i,j-1)) &
                  +0.5d0*(3*fn(i,j)-fo(i,j)) &
                  +reyi*(dxi**2*(v(i-1,j)-2.0d0*v(i,j)+v(i+1,j))) &
                  +reyi*(dxi**2*(v(i,j-1)-2.0d0*v(i,j)+v(i,j+1))))
          end do
       end do
    end do
    call amrex_mfiter_destroy(mfi)
  end subroutine velocity_rhs_v

  ! (1-A_1)(1-A_2)(v*-v^n) = rhs_v
  subroutine velocity_solve_v(field_n, field_o, mtm_rhs, geom, dt)
    implicit none
    type(amrex_multifab), intent(in) :: field_o, mtm_rhs
    type(amrex_geometry), intent(in) :: geom
    type(amrex_multifab), intent(inout) :: field_n
    real(amrex_real), intent(in) :: dt
    !-
    real(amrex_real) :: f

    f = dt/(2.0d0*rey)

    ! solve: (1-A_1) X = rhs_v
    call hypre_prepare_v_x(mtm_rhs, geom, f)
    call hypre_solve_v_x()

    ! solve: (1-A_2)(v*-v^n) = X
    call hypre_prepare_v_y(field_n, geom, f)
    call hypre_solve_v_y()

    ! update: v* = v^n + (v*-v^n)
    call hypre_update_v(field_n, field_o, geom)
  end subroutine velocity_solve_v


  ! hardcoded dirichlet boundary conditions for lid driven cavity
  subroutine velocity_bc(field, geom, time)
    implicit none
    type(amrex_geometry) :: geom
    type(amrex_multifab) :: field
    real(amrex_real), intent(in) :: time
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: dp
    real(amrex_real), contiguous, pointer, dimension(:,:) :: u, v
    integer :: i,j,ng,c, lb(4)

    ng = field%nghost()

    call amrex_mfiter_build(mfi, field)
    do while (mfi%next())

       bx = mfi%tilebox()
       dp => field%dataptr(mfi)
       lb = lbound(dp)

       u(lb(1):,lb(2):) => dp(:,:,1,U_i)
       v(lb(1):,lb(2):) => dp(:,:,1,V_i)

       if (bx%lo(1).eq.geom%domain%lo(1)) then ! xmin wall
          c = bx%lo(1)
          ! u(c,:,:) located at x-face
          ! v(c,:,:) located at x-center
          u(c,:) = 0.0d0
          do i = 1,ng
             u(c-i,:) = -u(c+i,:)
             v(c-i,:) = -v(c+i-1,:)
          end do
       end if
       if (bx%hi(1).eq.geom%domain%hi(1)) then ! xmax wall
          c = bx%hi(1)
          ! u(c+1,:) located at x-face
          ! v(c,:) located at x-center
          u(c+1,:) = 0.0d0
          do i=2,ng
             u(c+i,:) = -u(c+2-i,:)
          end do
          do i = 1,ng
             v(c+i,:) = -v(c-i+1,:)
          end do
       end if
       if (bx%lo(2).eq.geom%domain%lo(2)) then ! ymin wall
          c = bx%lo(2)
          ! u(:,c) located at y-center
          ! v(:,c) located at y-face
          v(:,c) = 0.0d0
          do j = 1,ng
             v(:,c-j) = -v(:,c+j)
             u(:,c-j) = -u(:,c+j-1)
          end do
       end if
       if (bx%hi(2).eq.geom%domain%hi(2)) then ! ymax wall
          c = bx%hi(2)
          ! u(:,c) located at y-center
          ! v(:,c+1) located at y-face
          v(:,c+1) = 0.0d0
          do j = 2,ng
             v(:,c+j) = -v(:,c+2-j)
          end do
          ! hardcoded top wall velocity of (1,0)
          do j = 1,ng
             u(:,c+j) = 2.0d0-u(:,c-j+1)
          end do
       end if
    end do
    call amrex_mfiter_destroy(mfi)
  end subroutine velocity_bc

  ! store -d(u*u)/dx-d(u*v)/dy in flux
  ! flux values live at same location as u-velocity
  !
  !     v(i-1,j+1)         v(i,j+1)
  !          |------yp-----|
  !          |             |
  !          |             |
  !         xm    u(i,j)   xp
  !          |             |
  !          |             |
  !          |------ym-----|
  !      v(i-1,j)          v(i,j)
  !
  ! d(u*u)/dx = (u_xp^2-u_xm^2)/dx
  ! d(u*v)/dy = (u_yp*v_yp-u_ym*v_ym)/dx
  !
  ! where:
  !  u_xp = (u(i+1,j)+u(i,j))/2
  !  u_xm = (u(i-1,j)+u(i,j))/2
  !  u_yp = (u(i,j+1)+u(i,j))/2
  !  u_ym = (u(i,j-1)+u(i,j))/2
  !  v_yp = (v(i-1,j+1)+v(i,j+1))/2
  !  v_ym = (v(i-1,j)+v(i,j))/2
  subroutine velocity_conv_flux_u(field, flux, geom)
    implicit none
    type(amrex_multifab), intent(in) :: field
    type(amrex_multifab), intent(inout) :: flux
    type(amrex_geometry), intent(in) :: geom
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: dp
    real(amrex_real), contiguous, pointer, dimension(:,:) :: u, v, f
    integer :: i,j, lb(4)
    real(amrex_real) :: dxi, u_xm, u_xp, u_ym, u_yp, v_ym, v_yp

    ! For simplicity, we assume dx(1) == dx(2) == dx(3)
    dxi = 1.0d0/geom%dx(1)

    call amrex_mfiter_build(mfi, field)
    do while (mfi%next())

       bx = mfi%tilebox()
       dp => field%dataptr(mfi)
       lb = lbound(dp)

       u(lb(1):,lb(2):) => dp(:,:,1,U_i)
       v(lb(1):,lb(2):) => dp(:,:,1,V_i)
       dp => flux%dataptr(mfi)
       f(lb(1):,lb(2):) => dp(:,:,1,U_i)

       do j = bx%lo(2),bx%hi(2)
          do i = bx%lo(1),bx%hi(1)
             u_xp = 0.5d0*(u(i+1,j)+u(i,j))
             u_xm = 0.5d0*(u(i-1,j)+u(i,j))
             u_yp = 0.5d0*(u(i,j+1)+u(i,j))
             u_ym = 0.5d0*(u(i,j-1)+u(i,j))
             v_yp = 0.5d0*(v(i-1,j+1)+v(i,j+1))
             v_ym = 0.5d0*(v(i-1,j)+v(i,j))

             f(i,j) = -dxi*(u_xp*u_xp-u_xm*u_xm) &
                  -dxi*(u_yp*v_yp-u_ym*v_ym)
          end do
       end do
    end do
    call amrex_mfiter_destroy(mfi)
  end subroutine velocity_conv_flux_u

  ! store -d(u*v)/dx-d(v*v)/dy in flux
  ! flux values live at same location as v-velocity
  !
  !
  !        u(i,j)        u(i+1,j)
  !          |------yp-----|
  !          |             |
  !          |             |
  !          |             |
  !         xm    v(i,j)   xp
  !          |             |
  !          |             |
  !          |------ym-----|
  !        u(i,j-1)      u(i+1,j-1)
  !
  ! d(u*v)/dx = (u_xp*v_xp-u_xm*v_xm)/dx
  ! d(v*v)/dy = (v_yp*v_yp-v_ym*v_ym)/dx
  !
  ! where:
  !  u_xp = (u(i+1,j)+u(i+1,j-1))/2
  !  u_xm = (u(i,j)+u(i,j-1))/2
  !  v_xp = (v(i+1,j)+v(i,j))/2
  !  v_xm = (v(i-1,j)+v(i,j))/2
  !  v_yp = (v(i,j+1)+v(i,j))/2
  !  v_ym = (v(i,j-1)+v(i,j))/2
  subroutine velocity_conv_flux_v(field, flux, geom)
    implicit none
    type(amrex_multifab), intent(in) :: field
    type(amrex_multifab), intent(inout) :: flux
    type(amrex_geometry), intent(in) :: geom
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: dp
    real(amrex_real), contiguous, pointer, dimension(:,:) :: u, v, f
    integer :: i,j,lb(4)
    real(amrex_real) :: dxi, u_xp, u_xm, v_xp, v_xm, v_yp, v_ym

    ! For simplicity, we assume dx(1) == dx(2) == dx(3)
    dxi = 1.0d0/geom%dx(1)

    call amrex_mfiter_build(mfi, field)
    do while (mfi%next())

       bx = mfi%tilebox()
       dp => field%dataptr(mfi)
       lb = lbound(dp)

       u(lb(1):,lb(2):) => dp(:,:,1,U_i)
       v(lb(1):,lb(2):) => dp(:,:,1,V_i)
       dp => flux%dataptr(mfi)
       f(lb(1):,lb(2):) => dp(:,:,1,V_i)

       do j = bx%lo(2),bx%hi(2)
          do i = bx%lo(1),bx%hi(1)
             u_xp = 0.5d0*(u(i+1,j)+u(i+1,j-1))
             u_xm = 0.5d0*(u(i,j)+u(i,j-1))
             v_xp = 0.5d0*(v(i+1,j)+v(i,j))
             v_xm = 0.5d0*(v(i-1,j)+v(i,j))
             v_yp = 0.5d0*(v(i,j+1)+v(i,j))
             v_ym = 0.5d0*(v(i,j-1)+v(i,j))

             f(i,j) = -dxi*(u_xp*v_xp-u_xm*v_xm) &
                  -dxi*(v_yp*v_yp-v_ym*v_ym)
          end do
       end do
    end do
    call amrex_mfiter_destroy(mfi)
  end subroutine velocity_conv_flux_v

end module velocity_module
