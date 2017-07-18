module io_module
  use ico_base_module
  use amrex_fi_mpi
  implicit none

  private

  public :: init_io, check_io, io_write_centerlines

  type(amrex_multifab) :: io
  integer :: plot_step_interval, plot_step
  real(amrex_real) :: plot_dt_interval, plot_time
  character(len=127) :: plot_file

contains

  ! initialize module number of plotting components
  subroutine init_io(ba, dm, ncomp)
    implicit none
    type(amrex_boxarray), intent(in) :: ba
    type(amrex_distromap), intent(in) :: dm
    integer, intent(in) :: ncomp
    !-
    type(amrex_parmparse) :: pp

    plot_file = "plt"
    plot_step_interval = 0
    plot_dt_interval = 0.0d0

    call amrex_parmparse_build(pp, "io")
    call pp%query("plot_dt_interval", plot_dt_interval)
    call pp%query("plot_step_interval", plot_step_interval)
    call pp%query("plot_file", plot_file)
    call amrex_parmparse_destroy(pp)

    if (plot_dt_interval.gt.0.0d0.and.plot_step_interval.gt.0) &
       stop "You can't do that"

    plot_step = 0
    plot_time = 0.0d0

    call amrex_multifab_build(io, ba, dm, ncomp, 0)

  end subroutine init_io


  subroutine check_io(step, dt, time, vars, geom)
    implicit none
    integer, intent(in) :: step
    real(amrex_real), intent(in) :: time, dt
    type(amrex_multifab), intent(in) :: vars
    type(amrex_geometry), intent(in) :: geom

    plot_time = plot_time+dt

    if (dumpp(step, time, dt)) then
       call data_to_io(vars)
       call write_io(geom, time)
    end if

  contains

    ! predicate function returning true if it's time to dump data
    function dumpp(step, time, dt) result (p)
      implicit none
      real(amrex_real), intent(in) :: time, dt
      integer, intent(in) :: step
      logical :: p


      p = .false.

      if (plot_step_interval.le.0.and.plot_dt_interval.le.0.0d0) &
           return

      ! dump initial data
      if (step.eq.0) then
         p = .true.
         return
      end if

      ! logic if using step_interval
      if (plot_step_interval.gt.0.and.mod(step, plot_step_interval).eq.0) then
         plot_step = plot_step+1
         p = .true.
         return
      end if

      ! logic if using dt_interval
      if (plot_dt_interval.gt.0.0d0.and.&
           plot_dt_interval-0.5d0*dt.le.plot_time) then
         p = .true.
         plot_step = plot_step+1
         plot_time = plot_time-plot_dt_interval
         return
      end if
    end function dumpp
  end subroutine check_io


  ! write 1st component of multifabs to io
  subroutine data_to_io(vars)
    implicit none
    type(amrex_multifab), intent(in) :: vars
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: io_dp, dp
    integer :: i, j

    call amrex_mfiter_build(mfi, io)
    do while (mfi%next())
       bx = mfi%tilebox()
       io_dp => io%dataptr(mfi)
       dp => vars%dataptr(mfi)

       do j=bx%lo(2),bx%hi(2)
          do i=bx%lo(1),bx%hi(1)
             io_dp(i,j,1,P_i) = dp(i,j,1,P_i)
             io_dp(i,j,1,U_i) = dp(i,j,1,U_i) !0.5d0*(dp(i,j,1,U_i)+dp(i+1,j,1,U_i))
             io_dp(i,j,1,V_i) = dp(i,j,1,V_i) !0.5d0*(dp(i,j,1,V_i)+dp(i,j+1,1,V_i))
          end do
       end do
    end do
    call amrex_mfiter_destroy(mfi)

  end subroutine data_to_io


  subroutine write_io(geom, time)
    implicit none
    type(amrex_geometry), intent(in) :: geom
    real(amrex_real), intent(in) :: time
    !-
    integer :: nlevs, rr(1), stepno(1)
    character(len=145) :: name
    character(len=16)  :: current_step
    type(amrex_string) :: varname(3)
    type(amrex_multifab) :: phi(1)
    type(amrex_geometry) :: geom_(1)

    phi(1) = io
    geom_(1) = geom
    stepno(1) = plot_step
    nlevs = 1
    rr(1) = 1

    if      (plot_step .lt. 1000000) then
       write(current_step,fmt='(i5.5)') plot_step
    else if (plot_step .lt. 10000000) then
       write(current_step,fmt='(i6.6)') plot_step
    else if (plot_step .lt. 100000000) then
       write(current_step,fmt='(i7.7)') plot_step
    else if (plot_step .lt. 1000000000) then
       write(current_step,fmt='(i8.8)') plot_step
    else
       write(current_step,fmt='(i15.15)') plot_step
    end if
    name = trim(plot_file) // current_step

    call amrex_string_build(varname(U_i), "U")
    call amrex_string_build(varname(V_i), "V")
    call amrex_string_build(varname(P_i), "P")

    call amrex_write_plotfile(name, nlevs, phi, varname, geom_, &
         time, stepno, rr)
  end subroutine write_io


  ! write y positions and u-velocity to u_centerline_$nx
  ! write x positions and v-velocity to v_centerline_$nx
  subroutine io_write_centerlines(field, geom)
    implicit none
    type(amrex_multifab), intent(in) :: field
    type(amrex_geometry), intent(in) :: geom
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: dp
    real(amrex_real), allocatable, dimension(:) :: u, v, pos
    integer :: i, j, i_cl, j_cl, nx, n_, lo, hi, status(MPI_STATUS_SIZE)
    integer(kind=MPI_OFFSET_KIND) :: pos_off, vel_off
    real(amrex_real) :: dx
    character(len=32) :: u_name, v_name, cnx
    integer :: u_unit, v_unit, ierr

    nx = geom%domain%hi(1)-geom%domain%lo(1)+1
    write(cnx,fmt='(i5.5)') nx
    u_name = "u_centerline_"//trim(adjustl(cnx))
    v_name = "v_centerline_"//trim(adjustl(cnx))

    call MPI_File_Open(MPI_COMM_WORLD, u_name, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, u_unit, ierr)
    call MPI_File_Open(MPI_COMM_WORLD, v_name, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, v_unit, ierr)

    dx = geom%dx(1)

    ! centerline indices: u(i_cl,:), v(:,j_cl)
    i_cl = (geom%domain%hi(1)+1)/2
    j_cl = (geom%domain%hi(2)+1)/2


    call amrex_mfiter_build(mfi, field)
    do while (mfi%next())
       bx = mfi%tilebox()

       ! process u-centerline
       if (bx%lo(1).le.i_cl.and.bx%hi(1).ge.i_cl) then
          lo = bx%lo(2)
          hi = bx%hi(2)
          n_ = hi-lo+1
          allocate(u(n_), pos(n_))

          dp => field%dataptr(mfi)
          u(:) = dp(i_cl,lo:hi,1,U_i)
          pos(:) = [ (dx*(0.5d0+real(j,amrex_real)),j=lo,hi) ]

          pos_off = 8*lo
          vel_off = pos_off+8*nx

          call MPI_File_Write_At(u_unit, pos_off, pos, n_, &
               MPI_DOUBLE_PRECISION, status, ierr)
          call MPI_File_Write_At(u_unit, vel_off, u, n_, &
               MPI_DOUBLE_PRECISION, status, ierr)

          deallocate(u, pos)
       end if

       ! process v-centerline
       if (bx%lo(2).le.j_cl.and.bx%hi(2).ge.j_cl) then
          lo = bx%lo(1)
          hi = bx%hi(1)
          n_ = hi-lo+1
          allocate(v(n_), pos(n_))

          dp => field%dataptr(mfi)
          v(:) = dp(lo:hi,j_cl,1,V_i)
          pos(:) = [ (dx*(0.5d0+real(i,amrex_real)),i=lo,hi) ]

          pos_off = 8*lo
          vel_off = pos_off+8*nx

          call MPI_File_Write_At(v_unit, pos_off, pos, n_, &
               MPI_DOUBLE_PRECISION, status, ierr)
          call MPI_File_Write_At(v_unit, vel_off, v, n_, &
               MPI_DOUBLE_PRECISION, status, ierr)

          deallocate(v, pos)
       end if
    end do
    call amrex_mfiter_destroy(mfi)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_File_Close(u_unit, ierr)
    call MPI_File_Close(v_unit, ierr)
  end subroutine io_write_centerlines

end module io_module
