module hypre_module
  use amrex_base_module
  use amrex_fi_mpi
  use iso_fortran_env, only: int64
  implicit none

  private

  public :: init_hypre, hypre_update, hypre_solve

  include 'HYPREf.h'

  integer(int64) :: grid, stencil, matrix, rhs, sol, solver
  integer :: comm
  contains

    subroutine init_hypre(phi)
      type(amrex_multifab), intent(in) :: phi
      integer :: ierr

      comm = MPI_COMM_WORLD !parallel_communicator()
      call grid_init(phi, amrex_spacedim)
      call stencil_init(amrex_spacedim)


      ! simple initialization of matrix, rhs, and solution
      call HYPRE_StructMatrixCreate(comm, grid, stencil, matrix, ierr)
      call HYPRE_StructMatrixInitialize(matrix, ierr)

      call HYPRE_StructVectorCreate(comm, grid, rhs, ierr)
      call HYPRE_StructVectorInitialize(rhs, ierr)

      call HYPRE_StructVectorCreate(comm, grid, sol, ierr)
      call HYPRE_StructVectorInitialize(sol, ierr)

    end subroutine init_hypre

    ! prepare the hypre data structures for the solve step
    ! requires that phi_old and phi_new have the appropriate
    ! boundary conditions set
    subroutine hypre_update(phi_old, phi_new, geom, c)
      implicit none
      type(amrex_multifab), intent(in) :: phi_old, phi_new
      type(amrex_geometry), intent(in) :: geom
      real(amrex_real), intent(in) :: c
      !-
      integer :: ierr

      call matrix_update(phi_new, geom, c)
      call rhs_update(phi_old, phi_new, geom, c)
      call sol_update(phi_old)

      ! debugging
!!$      call HYPRE_StructMatrixPrint(matrix, 1, ierr)
!!$      call HYPRE_StructVectorPrint(rhs, 1, ierr)
!!$      call execute_command_line("mv HYPRE_StructVector.out.00000 HYPRE_RHS.out.00000")
!!$      call HYPRE_StructVectorPrint(sol, 1, ierr)

!!$      call HYPRE_StructSMGCreate(comm, solver, ierr)
!!$      call HYPRE_StructSMGSetMaxIter(solver, 100, ierr)
!!$      !call HYPRE_StructSMGSetTol(solver, 1e-10, ierr)
!!$      call HYPRE_StructSMGSetLogging(solver, 2, ierr)
!!$      call HYPRE_StructSMGSetPrintLevel(solver, 2, ierr)
!!$      call HYPRE_StructSMGSetup(solver, matrix, rhs, sol, ierr)
!!$      print *, 'hypre smgsetup return code: ', ierr

      call HYPRE_StructPCGCreate(comm, solver, ierr)
      call HYPRE_StructPCGSetTol(solver, 1e-6, ierr)
      call HYPRE_StructPCGSetRelChange(solver, 1e-10, ierr)
!!$      call HYPRE_StructPCGSetPrintLevel(solver, 2, ierr)
      call HYPRE_StructPCGSetup(solver, matrix, rhs, sol, ierr)
    end subroutine hypre_update

    ! solve for the solution following a successful call to hypre_update
    subroutine hypre_solve(phi_new)
      implicit none
      type(amrex_multifab), intent(inout) :: phi_new
      !-
      integer :: ierr, it, i, j, ij(2)
      type(amrex_box) :: bx
      type(amrex_mfiter) :: mfi
      real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: p
      real(amrex_real) :: v, tol

      call HYPRE_StructPCGSolve(solver, matrix, rhs, sol, ierr)
!!$      print *, 'hypre solver return code: ', ierr
      call HYPRE_StructPCGGetNumIterations(solver, it, ierr)
      call HYPRE_StructPCGGetFinalRelative(solver, tol, ierr)

      print *, 'iterations: ', it
      print *, 'residual  : ', tol


      ! transfer solution
      call amrex_mfiter_build(mfi, phi_new, tiling=.true.)

      do while (mfi%next())

         bx = mfi%tilebox()
         p => phi_new%dataptr(mfi)

         do j=bx%lo(2),bx%hi(2)
            do i=bx%lo(1),bx%hi(1)

               ij = [i,j]
               call HYPRE_StructVectorGetValues(sol, ij, v, ierr)
               p(i,j,1,1) = v

            end do
         end do
      end do

      call amrex_mfiter_destroy(mfi)

      call HYPRE_StructPCGDestroy(solver, ierr)
      !stop
    end subroutine hypre_solve


    subroutine grid_init(phi, n)
      type(amrex_multifab), intent(in) :: phi
      integer, intent(in) :: n
      !-
      type(amrex_box) :: bx
      type(amrex_mfiter) :: mfi
      integer :: ierr

      call HYPRE_StructGridCreate(comm, n, grid, ierr)
      if (ierr.ne.0) then
         write(*,*) "HYPRE_StructGridCreate failed with ", ierr
         stop
      end if

      call amrex_mfiter_build(mfi, phi, tiling=.true.)

      do while (mfi%next())

         bx = mfi%tilebox()

         call HYPRE_StructGridSetExtents(grid, bx%lo(1:n), bx%hi(1:n), &
              ierr)
         if (ierr.ne.0) then
            write(*,*) "HYPRE_StructGridSetExtents failed with ", ierr
            stop
         end if

      end do

      call amrex_mfiter_destroy(mfi)

      call HYPRE_StructGridAssemble(grid, ierr)
      if (ierr.ne.0) then
         write(*,*) "HYPRE_StructGridAssemble failed with ", ierr
         stop
      end if

    end subroutine grid_init


    subroutine stencil_init(n)
      integer, intent(in) :: n
      !
      integer :: offsets(n, 2*n+1), ierr, i, j, o

      call HYPRE_StructStencilCreate(n, 2*n+1, stencil, ierr)
      if (ierr.ne.0) then
         write(*,*) "HYPRE_StructStencilCreate failed with ", ierr
         stop
      end if

      offsets = 0
      o = 2
      do i=1,n
         do j=-1,1,2
            offsets(i, o) = j
            o = o+1
         end do
      end do
      do i=1,2*n+1
         print *, offsets(:,i)
      end do

      do i=1,2*n+1
         call HYPRE_StructStencilSetElement(stencil, i-1, offsets(1,i), &
              ierr)
      end do

    end subroutine stencil_init


    subroutine matrix_update(phi, geom, c)
      implicit none

      type(amrex_multifab), intent(in) :: phi
      type(amrex_geometry), intent(in) :: geom
      real(amrex_real), intent(in) :: c
      !-
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(amrex_real) :: values(5),r(3)
      integer :: i,j,ij(2),st_ind(5),ierr

      r=1.0d0/geom%dx**2
      st_ind = [0,1,2,3,4]

      call amrex_mfiter_build(mfi, phi, tiling=.true.)
      do while(mfi%next())

         bx = mfi%tilebox()

         do j=bx%lo(2),bx%hi(2)
            do i=bx%lo(1),bx%hi(1)

               ij = [i,j]

               values(1) = 1.0d0+c*(2.0d0*sum(r(1:amrex_spacedim)))

               if (i.eq.geom%domain%lo(1)) then
                  values(1) = values(1) + c*r(1)
                  values(2) = 0.0d0
                  values(3) = -c*r(1)
               elseif (i.eq.geom%domain%hi(1)) then
                  values(1) = values(1) + c*r(1)
                  values(2) = -c*r(1)
                  values(3) = 0.0d0
               else
                  values(2) = -c*r(1)
                  values(3) = -c*r(1)
               endif

               if (j.eq.geom%domain%lo(2)) then
                  values(1) = values(1) + c*r(2)
                  values(4) = 0.0d0
                  values(5) = -c*r(2)
               elseif (j.eq.geom%domain%hi(2)) then
                  values(1) = values(1) + c*r(2)
                  values(4) = -c*r(2)
                  values(5) = 0.0d0
               else
                  values(4) = -c*r(2)
                  values(5) = -c*r(2)
               endif

               call HYPRE_StructMatrixSetValues(matrix, ij, 5, st_ind, values, &
                    ierr)
            end do
         end do
      end do
      call HYPRE_StructMatrixAssemble(matrix, ierr)

      if (ierr.ne.0) then
         write(*,*) "HYPRE_StructMatrixAssemble failed with ", ierr
         stop
      end if

    end subroutine matrix_update

    ! make sure boundary cells are set in phi
    subroutine rhs_update(phi_old, phi_new, geom, c)
      implicit none

      type(amrex_multifab), intent(in) :: phi_old, phi_new
      type(amrex_geometry), intent(in) :: geom
      real(amrex_real), intent(in) :: c
      !-
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(amrex_real) :: r(3),v
      real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: po, pn
      integer :: i,j,ij(2),ierr

      r=1.0d0/geom%dx**2

      call amrex_mfiter_build(mfi, phi_old, tiling=.true.)
      do while(mfi%next())

         bx = mfi%tilebox()
         po => phi_old%dataptr(mfi)
         pn => phi_new%dataptr(mfi)

         do j=bx%lo(2),bx%hi(2)
            do i=bx%lo(1),bx%hi(1)

               ij = [i,j]

               ! initialize rhs based on 'n' timestep data
               v = po(i,j,1,1)*(1.0d0-c*(2.0d0*sum(r(1:amrex_spacedim)))) &
                    +(po(i-1,j,1,1)+po(i+1,j,1,1))*(c*r(1)) &
                    +(po(i,j-1,1,1)+po(i,j+1,1,1))*(c*r(2))

               ! correct rhs for specified boundaries at 'n+1'
               ! assume the ghost cell value is the dirichlet wall value
               if (i.eq.geom%domain%lo(1)) then
                  v = v + c*2.0d0*pn(i-1,j,1,1)*r(1)
               elseif (i.eq.geom%domain%hi(1)) then
                  v = v + c*2.0d0*pn(i+1,j,1,1)*r(1)
               endif

               if (j.eq.geom%domain%lo(2)) then
                  v = v + c*2.0d0*pn(i,j-1,1,1)*r(2)
               elseif (j.eq.geom%domain%hi(2)) then
                  v = v + c*2.0d0*pn(i,j+1,1,1)*r(2)
               endif

               call HYPRE_StructVectorSetValues(rhs, ij, v, ierr)
            end do
         end do
      end do
      call HYPRE_StructVectorAssemble(rhs, ierr)

      if (ierr.ne.0) then
         write(*,*) "HYPRE_StructVectorAssemble failed with ", ierr
         stop
      end if

    end subroutine rhs_update

    ! initial guess -> old solution
    subroutine sol_update(phi)
      implicit none

      type(amrex_multifab), intent(in) :: phi
      !-
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(amrex_real) :: v
      real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: p
      integer :: i,j,ij(2),ierr

      call amrex_mfiter_build(mfi, phi, tiling=.true.)
      do while(mfi%next())

         bx = mfi%tilebox()
         p => phi%dataptr(mfi)

         do j=bx%lo(2),bx%hi(2)
            do i=bx%lo(1),bx%hi(1)

               ij = [i,j]
               v = p(i,j,1,1)

               call HYPRE_StructVectorSetValues(sol, ij, v, ierr)
            end do
         end do
      end do
      call HYPRE_StructVectorAssemble(sol, ierr)

      if (ierr.ne.0) then
         write(*,*) "HYPRE_StructVectorAssemble (sol) failed with ", ierr
         stop
      end if

    end subroutine sol_update



end module hypre_module
