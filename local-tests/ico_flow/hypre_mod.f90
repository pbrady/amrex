module hypre_module
  use ico_base_module
  use amrex_fi_mpi
  use iso_fortran_env, only: int64
  implicit none

  private

  public :: init_hypre_pres, init_hypre_u, init_hypre_v, &
       hypre_prepare_u_x, hypre_solve_u_x, &
       hypre_prepare_u_y, hypre_solve_u_y, &
       hypre_prepare_v_x, hypre_solve_v_x, &
       hypre_prepare_v_y, hypre_solve_v_y, &
       hypre_prepare_p, hypre_solve_p, &
       hypre_update_u, hypre_update_v, hypre_update_p
       !, hypre_update_v, hypre_solve_v

  include 'HYPREf.h'

  !character(len=64) :: hypre_solver
  real(amrex_real) :: pressure_tol
  integer :: pressure_iter

  ! need separate hypre structures for cell-centered and face based data
  integer(int64) :: grid_p, rhs_p, sol_p, stencil_p, matrix_p, solver_p
  integer(int64) :: stencil_x, stencil_y
  integer(int64) :: grid_u, rhs_u, sol_u, matrix_u_x, matrix_u_y, solver_u
  integer(int64) :: grid_v, rhs_v, sol_v, matrix_v_x, matrix_v_y, solver_v

  integer :: comm, ierr

contains

  subroutine init_hypre_pres(vars, g)
    implicit none
    type(amrex_multifab), intent(in) :: vars
    type(amrex_geometry), intent(in) :: g
    !-
    type(amrex_parmparse) :: pp

    pressure_tol=1.0e-6
    pressure_iter=20

    call amrex_parmparse_build(pp, "hypre")
    call pp%query("tol", pressure_tol)
    call pp%query("maxiter", pressure_iter)
    call amrex_parmparse_destroy(pp)

    comm = MPI_COMM_WORLD
    call grid_init_p(vars, grid_p, amrex_spacedim)
    call stencil_init_p(stencil_p, amrex_spacedim)
    call stencil_init_x(stencil_x, amrex_spacedim)
    call stencil_init_y(stencil_y, amrex_spacedim)

    call matrix_init_p(vars, g, grid_p, stencil_p, matrix_p)

    ! simple initialization of rhs, and solution
    call HYPRE_StructVectorCreate(comm, grid_p, rhs_p, ierr)
    call HYPRE_StructVectorInitialize(rhs_p, ierr)

    call HYPRE_StructVectorCreate(comm, grid_p, sol_p, ierr)
    call HYPRE_StructVectorInitialize(sol_p, ierr)

  end subroutine init_hypre_pres


  subroutine init_hypre_u(vars, g)
    implicit none
    type(amrex_multifab), intent(in) :: vars
    type(amrex_geometry), intent(in) :: g

    comm = MPI_COMM_WORLD
    call grid_init_u(vars, grid_u, g, amrex_spacedim)

    ! simple initialization of matrix, rhs, and solution
    call HYPRE_StructMatrixCreate(comm, grid_u, stencil_x, matrix_u_x, ierr)
    call HYPRE_StructMatrixInitialize(matrix_u_x, ierr)
    call HYPRE_StructMatrixCreate(comm, grid_u, stencil_y, matrix_u_y, ierr)
    call HYPRE_StructMatrixInitialize(matrix_u_y, ierr)

    call HYPRE_StructVectorCreate(comm, grid_u, rhs_u, ierr)
    call HYPRE_StructVectorInitialize(rhs_u, ierr)

    call HYPRE_StructVectorCreate(comm, grid_u, sol_u, ierr)
    call HYPRE_StructVectorInitialize(sol_u, ierr)

  end subroutine init_hypre_u


  subroutine init_hypre_v(vars, g)
    implicit none
    type(amrex_multifab), intent(in) :: vars
    type(amrex_geometry), intent(in) :: g

    comm = MPI_COMM_WORLD
    call grid_init_v(vars, grid_v, g, amrex_spacedim)

    ! simple initialization of matrix, rhs, and solution
    call HYPRE_StructMatrixCreate(comm, grid_v, stencil_x, matrix_v_x, ierr)
    call HYPRE_StructMatrixInitialize(matrix_v_x, ierr)
    call HYPRE_StructMatrixCreate(comm, grid_v, stencil_y, matrix_v_y, ierr)
    call HYPRE_StructMatrixInitialize(matrix_v_y, ierr)

    call HYPRE_StructVectorCreate(comm, grid_v, rhs_v, ierr)
    call HYPRE_StructVectorInitialize(rhs_v, ierr)

    call HYPRE_StructVectorCreate(comm, grid_v, sol_v, ierr)
    call HYPRE_StructVectorInitialize(sol_v, ierr)

  end subroutine init_hypre_v



  ! same grid initialization routine can be used for cell-centered
  ! and face based data since mfi%tilebox() will do the right thing
  subroutine grid_init_p(phi, grid, n)
    type(amrex_multifab), intent(in) :: phi
    integer(int64), intent(out) :: grid
    integer, intent(in) :: n
    !-
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi

    call HYPRE_StructGridCreate(comm, n, grid, ierr)

    call amrex_mfiter_build(mfi, phi)
    do while (mfi%next())
       bx = mfi%tilebox()
       call HYPRE_StructGridSetExtents(grid, bx%lo(1:n), bx%hi(1:n), ierr)
    end do
    call amrex_mfiter_destroy(mfi)

    call HYPRE_StructGridAssemble(grid, ierr)

  end subroutine grid_init_p


  subroutine grid_init_u(phi, grid, geom, n)
    type(amrex_multifab), intent(in) :: phi
    integer(int64), intent(out) :: grid
    type(amrex_geometry), intent(in) :: geom
    integer, intent(in) :: n
    !-
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    integer :: lo(n), hi(n)

    call HYPRE_StructGridCreate(comm, n, grid, ierr)

    call amrex_mfiter_build(mfi, phi)
    do while (mfi%next())
       bx = mfi%tilebox()
       lo = bx%lo(1:n)
       hi = bx%hi(1:n)
       if (lo(1).eq.geom%domain%lo(1)) lo(1) = lo(1)+1
       call HYPRE_StructGridSetExtents(grid, lo, hi, ierr)
    end do
    call amrex_mfiter_destroy(mfi)

    call HYPRE_StructGridAssemble(grid, ierr)

  end subroutine grid_init_u

  subroutine grid_init_v(phi, grid, geom, n)
    type(amrex_multifab), intent(in) :: phi
    integer(int64), intent(out) :: grid
    type(amrex_geometry), intent(in) :: geom
    integer, intent(in) :: n
    !-
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    integer :: lo(n), hi(n)

    call HYPRE_StructGridCreate(comm, n, grid, ierr)

    call amrex_mfiter_build(mfi, phi)
    do while (mfi%next())
       bx = mfi%tilebox()
       lo = bx%lo(1:n)
       hi = bx%hi(1:n)
       if (lo(2).eq.geom%domain%lo(2)) lo(2) = lo(2)+1
       call HYPRE_StructGridSetExtents(grid, lo, hi, ierr)
    end do
    call amrex_mfiter_destroy(mfi)

    call HYPRE_StructGridAssemble(grid, ierr)

  end subroutine grid_init_v


  ! multidimensional stencil for pressure
  subroutine stencil_init_p(stencil, n)
    integer, intent(in) :: n
    integer(int64), intent(out) :: stencil
    !
    integer :: offsets(n, 2*n+1), i, j, o

    call HYPRE_StructStencilCreate(n, 2*n+1, stencil, ierr)

    offsets = 0
    o = 2
    do i=1,n
       do j=-1,1,2
          offsets(i, o) = j
          o = o+1
       end do
    end do
    do i=1,2*n+1
       call HYPRE_StructStencilSetElement(stencil, i-1, offsets(1,i), ierr)
    end do

  end subroutine stencil_init_p

  ! tridiagonal in x stencil
  subroutine stencil_init_x(stencil, n)
    integer, intent(in) :: n
    integer(int64), intent(out) :: stencil
    !
    integer :: offsets(n, 3), i, n2

    n2 = size(offsets, dim=2)

    call HYPRE_StructStencilCreate(n, n2, stencil, ierr)

    offsets(:, 1) = [0, 0]
    offsets(:, 2) = [-1,0]
    offsets(:, 3) = [1, 0]
    do i=1,n2
       call HYPRE_StructStencilSetElement(stencil, i-1, offsets(1,i), ierr)
    end do
  end subroutine stencil_init_x

  ! tridiagonal in y stencil
  subroutine stencil_init_y(stencil, n)
    integer, intent(in) :: n
    integer(int64), intent(out) :: stencil
    !
    integer :: offsets(n, 3), i, n2

    n2 = size(offsets, dim=2)

    call HYPRE_StructStencilCreate(n, n2, stencil, ierr)

    offsets(:, 1) = [0, 0]
    offsets(:, 2) = [0,-1]
    offsets(:, 3) = [0, 1]
    do i=1,n2
       call HYPRE_StructStencilSetElement(stencil, i-1, offsets(1,i), ierr)
    end do

  end subroutine stencil_init_y

  subroutine matrix_init_p(field, geom, grid, stencil, matrix)
    implicit none
    type(amrex_multifab), intent(in) :: field
    type(amrex_geometry), intent(in) :: geom
    integer(int64), intent(in) :: grid, stencil
    integer(int64), intent(out) :: matrix
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real) :: v(5)
    integer :: i,j,ij(2),st_ind(5)

    call HYPRE_StructMatrixCreate(comm, grid, stencil, matrix, ierr)
    call HYPRE_StructMatrixInitialize(matrix, ierr)

    st_ind = [0,1,2,3,4]

    call amrex_mfiter_build(mfi, field)
    do while(mfi%next())
       bx = mfi%tilebox()

       do j=bx%lo(2),bx%hi(2)
          do i=bx%lo(1),bx%hi(1)

             ij = [i,j]

             v = 0.0d0

             if (i.gt.geom%domain%lo(1)) v(2) = -1.0d0
             if (i.lt.geom%domain%hi(1)) v(3) = -1.0d0
             if (j.gt.geom%domain%lo(2)) v(4) = -1.0d0
             if (j.lt.geom%domain%hi(2)) v(5) = -1.0d0

             v(1) = -sum(v)
             call HYPRE_StructMatrixSetValues(matrix, ij, 5, st_ind, v, ierr)
          end do
       end do
    end do
    call amrex_mfiter_destroy(mfi)
    call HYPRE_StructMatrixAssemble(matrix, ierr)
  end subroutine matrix_init_p


  subroutine hypre_prepare_u_x(mtm_rhs, geom, c)
    implicit none
    type(amrex_multifab), intent(in) :: mtm_rhs
    type(amrex_geometry), intent(in) :: geom
    real(amrex_real), intent(in) :: c
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real) :: lhs(3), v, coeff(3), dxi2
    integer :: i,j,ij(2),lo(2),hi(2),st_ind(3),lb(4)
    real(amrex_real), contiguous, pointer :: dp(:,:,:,:), rhs(:,:)

    st_ind = [0,1,2]
    dxi2 = 1.0d0/geom%dx(1)**2
    ! lhs coefficients in st_ind order
    coeff = [1.0d0+2.0d0*c*dxi2, -c*dxi2, -c*dxi2]

    call amrex_mfiter_build(mfi, mtm_rhs)
    do while(mfi%next())

       bx = mfi%tilebox()

       dp => mtm_rhs%dataPtr(mfi)
       lb = lbound(dp)

       rhs(lb(1):,lb(2):) => dp(:,:,1,U_i)

       lo = bx%lo(1:amrex_spacedim)
       hi = bx%hi(1:amrex_spacedim)
       if (bx%lo(1).eq.geom%domain%lo(1)) lo(1) = lo(1)+1

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ij = [i,j]

             lhs = coeff
             v = rhs(i,j)

             ! x-derivatives reaching into boundary.. adjust for
             ! time independent dirichlet conditions
             if (i.eq.geom%domain%lo(1)+1) lhs(2) = 0.0d0
             if (i.eq.geom%domain%hi(1))   lhs(3) = 0.0d0

             call HYPRE_StructMatrixSetValues(matrix_u_x, ij, 3, st_ind, &
                  lhs, ierr)
             call HYPRE_StructVectorSetValues(rhs_u, ij, v, ierr)
          end do
       end do
    end do
    call amrex_mfiter_destroy(mfi)
    call HYPRE_StructMatrixAssemble(matrix_u_x, ierr)
    call HYPRE_StructVectorAssemble(rhs_u, ierr)
    call HYPRE_StructVectorAssemble(sol_u, ierr) ! not sure if required
  end subroutine hypre_prepare_u_x

  subroutine hypre_solve_u_x()
    implicit none


    call HYPRE_StructCycRedCreate(comm, solver_u, ierr)
    call HYPRE_StructCycRedSetup(solver_u, matrix_u_x, rhs_u, sol_u, ierr)
    call HYPRE_StructCycRedSolve(solver_u, matrix_u_x, rhs_u, sol_u, ierr)
    call HYPRE_StructCycRedDestroy(solver_u, ierr)

!!$    call HYPRE_StructMatrixPrint(matrix_u_x, 1, ierr)
!!$    call execute_command_line("mv HYPRE_StructMatrix.out.00000 HYPRE_matrix_u_x")
!!$    call HYPRE_StructVectorPrint(rhs_u, 1, ierr)
!!$    call execute_command_line("mv HYPRE_StructVector.out.00000 HYPRE_rhs_u_x")
!!$    call HYPRE_StructVectorPrint(sol_u, 1, ierr)
!!$    call execute_command_line("mv HYPRE_StructVector.out.00000 HYPRE_sol_u_x")
!!$    stop
  end subroutine hypre_solve_u_x


  subroutine hypre_prepare_u_y(field, geom, c)
    implicit none
    type(amrex_multifab), intent(in) :: field
    type(amrex_geometry), intent(in) :: geom
    real(amrex_real), intent(in) :: c
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real) :: lhs(3), coeff(3), dxi2
    integer :: i,j,ij(2),lo(2),hi(2),st_ind(3)

    st_ind = [0,1,2]
    dxi2 = 1.0d0/geom%dx(1)**2
    ! lhs coefficients in st_ind order
    coeff = [1.0d0+2.0d0*c*dxi2, -c*dxi2, -c*dxi2]

    call amrex_mfiter_build(mfi, field)
    do while(mfi%next())

       bx = mfi%tilebox()

       lo = bx%lo(1:amrex_spacedim)
       hi = bx%hi(1:amrex_spacedim)
       if (bx%lo(1).eq.geom%domain%lo(1)) lo(1) = lo(1)+1

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ij = [i,j]

             lhs = coeff

             ! y-derivatives reaching into boundary.. adjust for
             ! time independent dirichlet conditions
             if (j.eq.geom%domain%lo(2)) then
                lhs(1) = lhs(1)-lhs(2)
                lhs(2) = 0.0d0
             end if
             if (j.eq.geom%domain%hi(2)) then
                lhs(1) = lhs(1)-lhs(3)
                lhs(3) = 0.0d0
             end if

             call HYPRE_StructMatrixSetValues(matrix_u_y, ij, 3, st_ind, &
                  lhs, ierr)
          end do
       end do
    end do
    call amrex_mfiter_destroy(mfi)
    call HYPRE_StructMatrixAssemble(matrix_u_y, ierr)
    call HYPRE_StructVectorAssemble(rhs_u, ierr)
    call HYPRE_StructVectorAssemble(sol_u, ierr)

  end subroutine hypre_prepare_u_y

  subroutine hypre_solve_u_y()
    implicit none

    call HYPRE_StructCycRedCreate(comm, solver_u, ierr)
    call HYPRE_StructCycRedSetTDim(solver_u, 1, ierr)
    call HYPRE_StructCycRedSetup(solver_u, matrix_u_y, sol_u, rhs_u, ierr)
    call HYPRE_StructCycRedSolve(solver_u, matrix_u_y, sol_u, rhs_u, ierr)
    call HYPRE_StructCycRedDestroy(solver_u, ierr)

!!$    call HYPRE_StructMatrixPrint(matrix_u_y, 1, ierr)
!!$    call execute_command_line("mv HYPRE_StructMatrix.out.00000 HYPRE_matrix_u_y")
!!$    call HYPRE_StructVectorPrint(sol_u, 1, ierr)
!!$    call execute_command_line("mv HYPRE_StructVector.out.00000 HYPRE_rhs_u_y")
!!$    call HYPRE_StructVectorPrint(rhs_u, 1, ierr)
!!$    call execute_command_line("mv HYPRE_StructVector.out.00000 HYPRE_sol_u_y")

  end subroutine hypre_solve_u_y

  subroutine hypre_update_u(field_n, field_o, geom)
    implicit none
    type(amrex_multifab), intent(in) :: field_o
    type(amrex_multifab), intent(inout) :: field_n
    type(amrex_geometry), intent(in) :: geom
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    integer :: i,j,ij(2),lo(2),hi(2),lb(4)
    real(amrex_real), contiguous, pointer :: dp(:,:,:,:), &
         un(:,:), uo(:,:)
    real(amrex_real) :: v


    call amrex_mfiter_build(mfi, field_n)
    do while(mfi%next())

       bx = mfi%tilebox()
       ! need un-uo on boundary
       dp => field_n%dataPtr(mfi)
       lb = lbound(dp)
       un(lb(1):,lb(2):) => dp(:,:,1,U_i)
       dp => field_o%dataPtr(mfi)
       uo(lb(1):,lb(2):) => dp(:,:,1,U_i)

       lo = bx%lo(1:amrex_spacedim)
       hi = bx%hi(1:amrex_spacedim)
       if (bx%lo(1).eq.geom%domain%lo(1)) lo(1) = lo(1)+1

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ij = [i,j]
             call HYPRE_StructVectorGetValues(rhs_u, ij, v, ierr)
             un(i,j) = uo(i,j)+v
          end do
       end do
    end do
    call amrex_mfiter_destroy(mfi)
  end subroutine hypre_update_u


  subroutine hypre_prepare_v_x(mtm_rhs, geom, c)
    implicit none
    type(amrex_multifab), intent(in) :: mtm_rhs
    type(amrex_geometry), intent(in) :: geom
    real(amrex_real), intent(in) :: c
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real) :: lhs(3), v, coeff(3), dxi2
    integer :: i,j,ij(2),lo(2),hi(2),st_ind(3),lb(4)
    real(amrex_real), contiguous, pointer :: dp(:,:,:,:), rhs(:,:)

    st_ind = [0,1,2]
    dxi2 = 1.0d0/geom%dx(1)**2
    ! lhs coefficients in st_ind order
    coeff = [1.0d0+2.0d0*c*dxi2, -c*dxi2, -c*dxi2]

    call amrex_mfiter_build(mfi, mtm_rhs)
    do while(mfi%next())

       bx = mfi%tilebox()

       dp => mtm_rhs%dataPtr(mfi)
       lb = lbound(dp)
       rhs(lb(1):,lb(2):) => dp(:,:,1,V_i)

       lo = bx%lo(1:amrex_spacedim)
       hi = bx%hi(1:amrex_spacedim)
       if (bx%lo(2).eq.geom%domain%lo(2)) lo(2) = lo(2)+1

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ij = [i,j]

             lhs = coeff
             v = rhs(i,j)

             ! x-derivatives reaching into boundary ... adjust for
             ! time independent dirichlet conditions
             if (i.eq.geom%domain%lo(1)) then
                lhs(1) = lhs(1)-lhs(2)
                lhs(2) = 0.0d0
             end if
             if (i.eq.geom%domain%hi(1))   then
                lhs(1) = lhs(1)-lhs(3)
                lhs(3) = 0.0d0
             end if

             call HYPRE_StructMatrixSetValues(matrix_v_x, ij, 3, st_ind, &
                  lhs, ierr)
             call HYPRE_StructVectorSetValues(rhs_v, ij, v, ierr)
          end do
       end do
    end do
    call amrex_mfiter_destroy(mfi)
    call HYPRE_StructMatrixAssemble(matrix_v_x, ierr)
    call HYPRE_StructVectorAssemble(rhs_v, ierr)
    call HYPRE_StructVectorAssemble(sol_v, ierr) ! not sure if required

  end subroutine hypre_prepare_v_x

  subroutine hypre_solve_v_x()
    implicit none

    call HYPRE_StructCycRedCreate(comm, solver_v, ierr)
    call HYPRE_StructCycRedSetup(solver_v, matrix_v_x, rhs_v, sol_v, ierr)
    call HYPRE_StructCycRedSolve(solver_v, matrix_v_x, rhs_v, sol_v, ierr)
    call HYPRE_StructCycRedDestroy(solver_v, ierr)

  end subroutine hypre_solve_v_x


  subroutine hypre_prepare_v_y(field, geom, c)
    implicit none
    type(amrex_multifab), intent(in) :: field
    type(amrex_geometry), intent(in) :: geom
    real(amrex_real), intent(in) :: c
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real) :: lhs(3), coeff(3), dxi2
    integer :: i,j,ij(2),lo(2),hi(2),st_ind(3)

    st_ind = [0,1,2]
    dxi2 = 1.0d0/geom%dx(1)**2
    ! lhs coefficients in st_ind order
    coeff = [1.0d0+2.0d0*c*dxi2, -c*dxi2, -c*dxi2]

    call amrex_mfiter_build(mfi, field)
    do while(mfi%next())

       bx = mfi%tilebox()

       lo = bx%lo(1:amrex_spacedim)
       hi = bx%hi(1:amrex_spacedim)
       if (bx%lo(2).eq.geom%domain%lo(2)) lo(2) = lo(2)+1

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ij = [i,j]

             lhs = coeff

             if (j.eq.geom%domain%lo(2)+1) lhs(2) = 0.0d0
             if (j.eq.geom%domain%hi(2))   lhs(3) = 0.0d0

             call HYPRE_StructMatrixSetValues(matrix_v_y, ij, 3, st_ind, &
                  lhs, ierr)
          end do
       end do
    end do
    call amrex_mfiter_destroy(mfi)
    call HYPRE_StructMatrixAssemble(matrix_v_y, ierr)
    call HYPRE_StructVectorAssemble(rhs_v, ierr)
    call HYPRE_StructVectorAssemble(sol_v, ierr)

  end subroutine hypre_prepare_v_y

  subroutine hypre_solve_v_y()
    implicit none

    call HYPRE_StructCycRedCreate(comm, solver_v, ierr)
    call HYPRE_StructCycRedSetTDim(solver_v, 1, ierr)
    call HYPRE_StructCycRedSetup(solver_v, matrix_v_y, sol_v, rhs_v, ierr)
    call HYPRE_StructCycRedSolve(solver_v, matrix_v_y, sol_v, rhs_v, ierr)
    call HYPRE_StructCycRedDestroy(solver_v, ierr)

  end subroutine hypre_solve_v_y

  subroutine hypre_update_v(field_n, field_o, geom)
    implicit none
    type(amrex_multifab), intent(in) :: field_o
    type(amrex_multifab), intent(inout) :: field_n
    type(amrex_geometry), intent(in) :: geom
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    integer :: i,j,ij(2),lo(2),hi(2),lb(4)
    real(amrex_real), contiguous, pointer :: dp(:,:,:,:), &
         un(:,:), uo(:,:)
    real(amrex_real) :: v


    call amrex_mfiter_build(mfi, field_n)
    do while(mfi%next())

       bx = mfi%tilebox()
       ! need un-uo on boundary
       dp => field_n%dataPtr(mfi)
       lb = lbound(dp)

       un(lb(1):,lb(2):) => dp(:,:,1,V_i)
       dp => field_o%dataPtr(mfi)
       uo(lb(1):,lb(2):) => dp(:,:,1,V_i)

       lo = bx%lo(1:amrex_spacedim)
       hi = bx%hi(1:amrex_spacedim)
       if (bx%lo(2).eq.geom%domain%lo(2)) lo(2) = lo(2)+1

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ij = [i,j]
             call HYPRE_StructVectorGetValues(rhs_v, ij, v, ierr)
             un(i,j) = uo(i,j)+v
          end do
       end do
    end do
    call amrex_mfiter_destroy(mfi)

  end subroutine hypre_update_v


  subroutine hypre_prepare_p(field, geom)
    implicit none
    type(amrex_multifab), intent(in) :: field
    type(amrex_geometry), intent(in) :: geom
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    integer :: i,j,ij(2),lb(4)
    real(amrex_real), contiguous, pointer :: dp(:,:,:,:), rhs(:,:)
    real(amrex_real) :: d

    d = geom%dx(1)**2

    call amrex_mfiter_build(mfi, field)
    do while(mfi%next())

       bx = mfi%tilebox()

       dp => field%dataPtr(mfi)
       lb = lbound(dp)

       rhs(lb(1):,lb(2):) => dp(:,:,1,P_i)

       do j=bx%lo(2),bx%hi(2)
          do i=bx%lo(1),bx%hi(1)
             ij = [i,j]
             call HYPRE_StructVectorSetValues(rhs_p, ij, -d*rhs(i,j), ierr)
          end do
       end do
    end do
    call amrex_mfiter_destroy(mfi)
    call HYPRE_StructVectorAssemble(rhs_p, ierr)
    call HYPRE_StructVectorAssemble(sol_p, ierr)

  end subroutine hypre_prepare_p


  subroutine hypre_solve_p
    implicit none
    integer :: it
    real(amrex_real) :: tol

    call HYPRE_StructSMGCreate(comm, solver_p, ierr)
    call HYPRE_StructSMGSetTol(solver_p, pressure_tol, ierr)
    call HYPRE_StructSMGSetMaxIter(solver_p, pressure_iter, ierr)
    call HYPRE_StructSMGSetLogging(solver_p, 1, ierr)
    call HYPRE_StructSMGSetup(solver_p, matrix_p, rhs_p, sol_p, ierr)

    call HYPRE_StructSMGSolve(solver_p, matrix_p, rhs_p, sol_p, ierr)
    call HYPRE_StructSMGGetNumIterations(solver_p, it, ierr)
    call HYPRE_StructSMGGetFinalRelative(solver_p, tol, ierr)
    call HYPRE_StructSMGDestroy(solver_p, ierr)

    if ( amrex_parallel_IOProcessor() ) then
       print *, 'pressure residual=', tol,'with iterations=', it
    end if

  end subroutine hypre_solve_p


  subroutine hypre_update_p(field_n, field_o, geom)
    implicit none
    type(amrex_multifab), intent(in) :: field_o
    type(amrex_multifab), intent(inout) :: field_n
    type(amrex_geometry), intent(in) :: geom
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    integer :: i,j,ij(2),lb(4)
    real(amrex_real), contiguous, pointer :: dp(:,:,:,:), &
         pn(:,:), po(:,:)
    real(amrex_real) :: v, sum, avg, d

    sum = 0.0d0
    d = geom%dx(1)**2

    call amrex_mfiter_build(mfi, field_n)
    do while(mfi%next())

       bx = mfi%tilebox()
       ! need un-uo on boundary
       dp => field_n%dataPtr(mfi)
       lb = lbound(dp)

       pn(lb(1):,lb(2):) => dp(:,:,1,P_i)
       dp => field_o%dataPtr(mfi)
       po(lb(1):,lb(2):) => dp(:,:,1,P_i)

       do j=bx%lo(2),bx%hi(2)
          do i=bx%lo(1),bx%hi(1)

             ij = [i,j]
             call HYPRE_StructVectorGetValues(sol_p, ij, v, ierr)
             pn(i,j) = po(i,j)+v
             sum = sum + d*pn(i,j)
          end do
       end do
    end do
    call amrex_mfiter_destroy(mfi)

    ! demean pressure field
    call MPI_AllReduce(sum, avg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    call amrex_mfiter_build(mfi, field_n)
    do while(mfi%next())

       bx = mfi%tilebox()
       ! need un-uo on boundary
       dp => field_n%dataPtr(mfi)
       lb = lbound(dp)
       pn(lb(1):,lb(2):) => dp(:,:,1,P_i)

       do j=bx%lo(2),bx%hi(2)
          do i=bx%lo(1),bx%hi(1)
             pn(i,j) = pn(i,j)-avg
          end do
       end do
    end do
    call amrex_mfiter_destroy(mfi)
  end subroutine hypre_update_p

end module hypre_module
