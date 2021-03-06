
module amrex_geometry_module

  use iso_c_binding
  use amrex_fort_module, only : ndims => amrex_spacedim, amrex_real
  use amrex_box_module

  implicit none

  private

  public :: amrex_pmask, amrex_problo, amrex_probhi
  public :: amrex_geometry_build, amrex_geometry_destroy, amrex_geometry_init_data

  logical, save :: amrex_pmask(3)  = .false.  
  real(amrex_real), save :: amrex_problo(3) = 0.0_amrex_real
  real(amrex_real), save :: amrex_probhi(3) = 1.0_amrex_real

  logical, save :: amrex_geometry_initialzied = .false.

!$omp threadprivate(amrex_geometry_initialzied,amrex_pmask,amrex_problo,amrex_probhi)

  type, public :: amrex_geometry
     logical          :: owner     = .false.
     type(c_ptr)      :: p         = c_null_ptr
     !
     real(amrex_real) :: dx(3)     = 0.0_amrex_real
     type(amrex_box)  :: domain
   contains
     generic :: assignment(=) => amrex_geometry_assign  ! shallow copy
     procedure, private :: amrex_geometry_assign
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_geometry_destroy
#endif
  end type amrex_geometry

  ! interfaces to c++ functions

  interface
     subroutine amrex_fi_new_geometry (geom,lo,hi) bind(c)
       import
       implicit none
       type(c_ptr) :: geom
       integer, intent(in) :: lo(3), hi(3)
     end subroutine amrex_fi_new_geometry

     subroutine amrex_fi_delete_geometry (geom) bind(c)
       import
       implicit none
       type(c_ptr), value :: geom
     end subroutine amrex_fi_delete_geometry

     subroutine amrex_fi_geometry_get_pmask (pmask) bind(c)
       import
       implicit none
       integer(c_int) :: pmask(3)
     end subroutine amrex_fi_geometry_get_pmask

     subroutine amrex_fi_geometry_get_probdomain (problo,probhi) bind(c)
       import
       implicit none
       real(amrex_real) :: problo(3), probhi(3)
     end subroutine amrex_fi_geometry_get_probdomain

     subroutine amrex_fi_geometry_get_intdomain (geom,lo,hi) bind(c)
       import
       implicit none
       type(c_ptr), value :: geom
       integer(c_int), intent(out) :: lo(3), hi(3)
     end subroutine amrex_fi_geometry_get_intdomain
  end interface

contains

  subroutine amrex_geometry_init ()
    integer :: imask(3)
    imask = 0
    call amrex_fi_geometry_get_pmask(imask)
    where (imask .eq. 1) amrex_pmask = .true.
    call amrex_fi_geometry_get_probdomain(amrex_problo, amrex_probhi)
  end subroutine amrex_geometry_init

  subroutine amrex_geometry_build (geom, domain)
    type(amrex_geometry) :: geom
    type(amrex_box), intent(in) :: domain
    geom%owner = .true.
    call amrex_fi_new_geometry(geom%p, domain%lo, domain%hi)
    call amrex_geometry_init_data(geom)
  end subroutine amrex_geometry_build

  subroutine amrex_geometry_init_data (geom)  ! geom%p must be valid!
    type(amrex_geometry), intent(inout) :: geom
    integer :: i, lo(3), hi(3)
    if (.not.amrex_geometry_initialzied) then
       call amrex_geometry_init()
       amrex_geometry_initialzied = .true.
    end if
    call amrex_fi_geometry_get_intdomain(geom%p, lo, hi)
    geom%domain = amrex_box(lo, hi)
    do i = 1, ndims
       geom%dx(i) = (amrex_probhi(i)-amrex_problo(i)) / dble(hi(i)-lo(i)+1)
    end do
  end subroutine amrex_geometry_init_data

  impure elemental subroutine amrex_geometry_destroy (geom)
    type(amrex_geometry), intent(inout) :: geom
    if (geom%owner) then
       if (c_associated(geom%p)) then
          call amrex_fi_delete_geometry(geom%p)
       end if
    end if
    geom%owner = .false.
    geom%p = c_null_ptr
  end subroutine amrex_geometry_destroy

  subroutine amrex_geometry_assign (dst, src)
    class(amrex_geometry), intent(inout) :: dst
    type (amrex_geometry), intent(in   ) :: src
    call amrex_geometry_destroy(dst)
    dst%owner  = .false.
    dst%p      = src%p
    dst%dx     = src%dx
    dst%domain = src%domain
  end subroutine amrex_geometry_assign

end module amrex_geometry_module
