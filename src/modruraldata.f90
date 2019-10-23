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

  !< Field for the wall shear stress
  real, allocatable    :: damping(:,:,:)

contains
  subroutine applydamping(ekm,i,j,k)
    
    implicit none
    real, intent(inout) :: ekm
    integer, intent(in) :: i,j,k

    if ((.not. (lruralboundary)).or.(.not.(lwallfunc))) return
    !if (i>17) write(6,*) 'i = ', i
    !if (j>17) write(6,*) 'j = ', j
    !if (k>20) write(6,*) 'k = ', k
    ekm=ekm*damping(i-1,j-1,k)
    !write(6,*) 'wall damping applied'
    return
  end subroutine applydamping
end module modruraldata
