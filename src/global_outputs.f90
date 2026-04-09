module global_outputs

  !=========================
  ! Flags and global outputs
  ! Module added by A.JANIN 14.02.2026
  !=========================

  implicit none

  !=========================
  ! Output flags ofr optional arrays
  !=========================

  logical :: o_orient
  logical :: o_invariant
  logical :: o_failure
  logical :: o_ugrad

  !=========================
  ! Global output arrays
  !=========================
  real(8), allocatable :: xfgrid(:), yfgrid(:), zfgrid(:)
  real(8), allocatable :: goutarray_displ(:,:)
  real(8), allocatable :: goutarray_stress(:,:)
  real(8), allocatable :: goutarray_strain(:,:)
  real(8), allocatable :: goutarray_failure(:,:)
  real(8), allocatable :: goutarray_orient(:,:)
  real(8), allocatable :: goutarray_invariant(:,:)
  real(8), allocatable :: goutarray_ugrad(:,:,:)
  real(8), allocatable :: goutarray_elements(:,:)

  !=========================
  ! Output flags
  !=========================
  logical :: goutflag_converged


contains

  subroutine allocate_global_outputs(npts, ndis)
    implicit none
    integer, intent(in) :: npts, ndis

    ! global output arrays
    allocate(xfgrid(npts), yfgrid(npts), zfgrid(npts))
    allocate(goutarray_displ(npts,3))
    allocate(goutarray_stress(npts,6))
    allocate(goutarray_strain(npts,6))
    allocate(goutarray_failure(npts,6))
    allocate(goutarray_orient(npts,9))
    allocate(goutarray_invariant(npts,4))
    allocate(goutarray_ugrad(npts,3,3))
    allocate(goutarray_elements(ndis,12))

  end subroutine allocate_global_outputs


  subroutine deallocate_global_outputs()
    implicit none

    deallocate(xfgrid, yfgrid, zfgrid)
    deallocate(goutarray_displ)
    deallocate(goutarray_stress)
    deallocate(goutarray_strain)
    deallocate(goutarray_failure)
    deallocate(goutarray_orient)
    deallocate(goutarray_invariant)
    deallocate(goutarray_ugrad)
    deallocate(goutarray_elements)

  end subroutine deallocate_global_outputs

end module global_outputs
