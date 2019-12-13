!> \file modruraldata.f90
!! By Michael Koene, TU Delft, section Atmospheric Physics, 27 September 2019

module modruraldata

  implicit none
  save
  public :: applydamping

  logical :: lruralboundary = .false.        !< Switch to enable rural boundary representation
  logical :: ldefrural      = .false.        !< Switch to customize the rural boundary area
  logical :: lnoslip        = .false.        !< Switch to use a no slip condition for the walls
  logical :: lwallfunc      = .false.        !< Switch to use wallfunctions to describe wallshear

  integer :: imaxb,jmaxb ! Integers imax and jmax to prevent circular dependency
  
  !< Field for the wall shear stress
  real, allocatable    :: damping(:,:,:)
  !< Field for the immersed boundary height
  integer, allocatable :: bc_height(:,:)     !< Height of immersed boundary at grid pos x,y
  !< Pressurefield for the correction term
  real, allocatable    :: pres0(:,:,:)

contains
  subroutine applydamping(ekm,i,j,k)
    
    implicit none
    real, intent(inout) :: ekm
    integer, intent(in) :: i,j,k

    if ((.not. (lruralboundary)).or.(.not.(lwallfunc))) return
    ekm=ekm*damping(i-1,j-1,k)
    return
  end subroutine applydamping
end module modruraldata
