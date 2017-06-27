module io_module
  use amrex_base_module
  implicit none

  private

  public :: init_io, check_io

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


  subroutine check_io(step, dt, time, P, U, V, geom)
    implicit none
    integer, intent(in) :: step
    real(amrex_real), intent(in) :: time, dt
    type(amrex_multifab), intent(in) :: P, U, V
    type(amrex_geometry), intent(in) :: geom


    if (dumpp(step, time, dt)) then
       call data_to_io(P, U, V)
       call write_io(geom, time)
    end if

    plot_time = plot_time+dt

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
  subroutine data_to_io(P, U, V)
    implicit none
    type(amrex_multifab), intent(in) :: P, U, V
    !-
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: io_dp, p_dp, &
         u_dp, v_dp
    integer :: i, j, c

    call amrex_mfiter_build(mfi, io)
    do while (mfi%next())
       bx = mfi%tilebox()
       io_dp => io%dataptr(mfi)

       P_dp => P%dataptr(mfi)
       U_dp => U%dataptr(mfi)
       V_dp => V%dataptr(mfi)

       do j=bx%lo(2),bx%hi(2)
          do i=bx%lo(1),bx%hi(1)
             io_dp(i,j,1,1) = P_dp(i,j,1,1)
             io_dp(i,j,1,2) = 0.5d0*(U_dp(i,j,1,1)+U_dp(i+1,j,1,1))
             io_dp(i,j,1,3) = 0.5d0*(V_dp(i,j,1,1)+V_dp(i,j+1,1,1))
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

    call amrex_string_build(varname(1), "P")
    call amrex_string_build(varname(2), "U")
    call amrex_string_build(varname(3), "V")

    call amrex_write_plotfile(name, nlevs, phi, varname, geom_, &
         time, stepno, rr)
  end subroutine write_io

end module io_module
