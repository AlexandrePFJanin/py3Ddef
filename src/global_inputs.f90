module global_inputs

  !=========================
  ! Flags and global inputs
  ! Module added by A.JANIN 27.02.2026
  !=========================

  implicit none

  !=========================
  ! Global input arrays
  !=========================
  integer, allocatable :: i_kode(:)
  real*4, allocatable  :: i_ss(:)
  real*4, allocatable  :: i_ds(:)
  real*4, allocatable  :: i_ts(:)
  integer, allocatable :: i_fcode(:)
  real*4, allocatable  :: i_sdrop(:)
  real*4, allocatable  :: i_rhoLitho(:)
  real*4, allocatable  :: i_rhoFluid(:)
  real*4, allocatable  :: i_cohes(:)
  real*4, allocatable  :: i_disfric(:)

  !=========================
  ! Friction solver parameter
  !=========================
  integer*4 :: i_maxiter
  real*4    :: i_tolsolver


contains

  subroutine allocate_global_inputs(ndis)
    implicit none
    integer, intent(in) :: ndis

    allocate(i_kode(ndis))
    allocate(i_ss(ndis))
    allocate(i_ds(ndis))
    allocate(i_ts(ndis))
    allocate(i_fcode(ndis))
    allocate(i_sdrop(ndis))
    allocate(i_rhoLitho(ndis))
    allocate(i_rhoFluid(ndis))
    allocate(i_cohes(ndis))
    allocate(i_disfric(ndis))

  end subroutine allocate_global_inputs

  subroutine deallocate_global_inputs()
    implicit none

    deallocate(i_kode)
    deallocate(i_ss)
    deallocate(i_ds)
    deallocate(i_ts)
    deallocate(i_fcode)
    deallocate(i_sdrop)
    deallocate(i_rhoLitho)
    deallocate(i_rhoFluid)
    deallocate(i_cohes)
    deallocate(i_disfric)

  end subroutine deallocate_global_inputs

end module global_inputs
