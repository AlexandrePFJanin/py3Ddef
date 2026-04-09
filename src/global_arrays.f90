module global_arrays

  !=========================
  ! Dynamic memory allocation for global arrays
  ! Module added by A.JANIN 13.02.2026
  !=========================

  implicit none

  !=========================
  ! Dimensions runtime
  !=========================
  integer :: npts_glob
  integer :: ndis_glob

  !=========================
  ! Boundary conditions
  !=========================
  integer, allocatable :: ISPACE(:)
  real(4), allocatable :: SPACE(:)

  !=========================
  ! Influence matrix and derivates
  !=========================
  real(4), allocatable :: XMATRIX(:,:) ! The matrix: influence coeff + extra (last) column for b.c.s
  real(4), allocatable :: AMATRIX(:,:) ! Original coeff matrix only, never shrinked or reshaped
  integer :: NUM_Ds_SAVED ! keep track of NUM_Ds, index of column where the b.c.s are written
  logical :: debug
  real(4), allocatable :: DVEC(:,:)    ! Vectors of *displacement* on each element (along-strike,down-dip,normal-out): returned by SOLVE
  real(4), allocatable :: DVEC2(:,:)   ! Same as DVEC but ALL relative displacements have been set to 0. (usefull in friction, to not count twice the effect of KODE>=10)
  real(4), allocatable :: DVEC0(:,:)   ! Reference DVEC: vector DVEC at the previous iteration - used during the convergence search
  real(4), allocatable :: DVECI(:,:) ! Displacement built iterativelly when frictional elements are introduced

  !=========================
  ! Stresses on each dislocations
  !=========================
  real(4), allocatable :: SSTORED(:,:) ! not released (only when friction)
  real(4), allocatable :: SDRIVER(:,:) ! released for motion

  !=========================
  ! Transformation matrices
  !=========================
  real(4), allocatable :: UG2P(:,:,:)
  real(4), allocatable :: SG2P(:,:,:)

  !=========================
  ! Element descriptors
  !=========================
  real(4), allocatable :: ZCE(:) ! z coordinates (depths) of the central point of elements
  real(4), allocatable :: XO(:),YO(:),ZO(:)
  real(4), allocatable :: C(:),S(:),DIP(:)
  real(4), allocatable :: CDIP(:),SDIP(:)
  real(4), allocatable :: BWX1(:),BWX2(:)
  integer, allocatable :: NBX1(:),NBX2(:)

  !=========================
  ! Station coordinates
  !=========================
  real(8), allocatable :: xu(:),yu(:),zu(:)

  !=========================
  ! Frictional status of the elements
  !=========================
  logical :: any_frictional    ! True if there is any frictional element
  integer(4), allocatable :: element_fstatus(:)       ! globally, what is the frictional status of the element (if the elements slides once in the iteration process, will be set as "sliding" for ever): that's the definitive fstatus.
  integer(4), allocatable :: element_iter_fstatus(:)  ! frictional status of each element only from one iteration to another (so, one element can goes from sliding to locked to sliding again with the iterative process)
  !  -10: unknown
  !   -2: not frictional
  !   -1: frictional but not determined
  !    0: locked
  !    1: sliding

  !=========================
  ! Final solution status
  !=========================
  integer(4), allocatable:: codeStatus(:)

contains

  subroutine allocate_global_arrays(npts,ndis)
    implicit none
    integer, intent(in) :: npts, ndis

    npts_glob = npts
    ndis_glob = ndis

    ! Reserved space for b.c. (boundary condition) codes, and b.c.s
    allocate(ISPACE(ndis))
    allocate(SPACE(3*ndis))

    ! Influence coefficient matrix; last column contains the 
    ! b.c.s and then the solution (relative displacements)
    allocate(XMATRIX(3*ndis,3*ndis+1))

    ! Original coeff matrix
    allocate(AMATRIX(3*ndis,3*ndis))

    ! Displacement vector on each element
    allocate(DVEC(ndis,3))
    allocate(DVEC2(ndis,3))
    allocate(DVEC0(ndis,3))
    allocate(DVECI(ndis,3))
    DVEC  = 0.0    ! Init the displacement on all elements to 0.0
    DVEC2 = 0.0    ! Init at 0
    DVEC0 = 0.0    ! Init at 0, no displacement
    DVECI = 0.0    ! Init at 0

    ! Stresses on each dislocation
    allocate(SSTORED(ndis,3))
    allocate(SDRIVER(ndis,3))
    SSTORED = 0.0  ! Initialize all elements to 0.0
    SDRIVER = 0.0  ! Initialize all elements to 0.0

    ! MAX_PLN = ndis in py3Ddef: no subelements

    ! Displacement and stress transformation matrices from
    ! global to planar (in plane) coordinates
    allocate(UG2P(3,3,ndis)) ! displacement transformation matrix for each plane 
    allocate(SG2P(3,6,ndis)) ! stress transformation matrix for each plane 

    ! Element descriptor parameters; for each plane
    allocate(ZCE(ndis))
    allocate(XO(ndis),YO(ndis),ZO(ndis)) ! X,Y,Z reference point in global coordinates
    allocate(C(ndis),S(ndis),DIP(ndis))  ! cosine & sine of the strike (cw wrt N)
    allocate(CDIP(ndis),SDIP(ndis))   ! dip (wrt horizontal), cosine & sine of the dip
    allocate(BWX1(ndis),BWX2(ndis))   ! sub-element widths in the strike & dip directions
    allocate(NBX1(ndis),NBX2(ndis))   ! number of sub-elements in the strike & dip directions

    ! user coordinates
    allocate(xu(npts),yu(npts),zu(npts))

    ! frictional status of elements
    allocate(element_fstatus(ndis))
    allocate(element_iter_fstatus(ndis))
    element_fstatus = -10        ! init on "unknown"
    element_iter_fstatus = -10   ! init on "unknown"

    ! final solution status when existing
    allocate(codeStatus(2)) ! codeStatus(1) = converged?, 1=True ,0=No
                            ! codeStatus(2) = final Niter

  end subroutine allocate_global_arrays


  subroutine deallocate_global_arrays()
    implicit none

    deallocate(ISPACE)
    deallocate(SPACE)

    deallocate(XMATRIX)
    deallocate(AMATRIX)
    deallocate(DVEC)
    deallocate(DVEC2)
    deallocate(DVEC0)
    deallocate(DVECI)

    deallocate(SSTORED)
    deallocate(SDRIVER)

    deallocate(UG2P)
    deallocate(SG2P)

    deallocate(ZCE)
    deallocate(XO,YO,ZO)
    deallocate(C,S,DIP)
    deallocate(CDIP,SDIP)
    deallocate(BWX1,BWX2)
    deallocate(NBX1,NBX2)

    deallocate(xu,yu,zu)

    deallocate(element_fstatus)
    deallocate(element_iter_fstatus)

    deallocate(codeStatus)

  end subroutine deallocate_global_arrays


  subroutine reset_xmatrix()
    implicit none

    ! -------------------------------------------------------
    ! Reset XMATRIX from AMATRIX
    ! -------------------------------------------------------

    XMATRIX(:, 1:3*ndis_glob) = AMATRIX
    XMATRIX(:, 3*ndis_glob+1) = 0.
  
  end subroutine reset_xmatrix


  subroutine update_DVEC()
    implicit none
    
    integer :: K, IBC

    ! -------------------------------------------------------
    ! Update the vector DVEC: To be put just before each RETURN in SOLVE
    ! -------------------------------------------------------
    DO K = 1, ndis_glob
        IBC = (K-1)*3
        DVEC(K,1) = XMATRIX(IBC+1,NUM_Ds_SAVED)
        DVEC(K,2) = XMATRIX(IBC+2,NUM_Ds_SAVED)
        DVEC(K,3) = XMATRIX(IBC+3,NUM_Ds_SAVED)
    ENDDO

    ! debug mode
    if (debug) then
      write(*,*) ''
      write(*,*) ' [DEBUG] (updated DVEC)'
      DO K = 1, ndis_glob
        IBC = (K-1)*3
        write(*,*) 'DVEC:',K,DVEC(K,1),DVEC(K,2),DVEC(K,3)
      ENDDO
    endif

  end subroutine update_DVEC


  subroutine update_DVEC2(i_kode)
    implicit none
    
    integer :: K, IBC
    integer :: i_kode(:)

    ! -------------------------------------------------------
    ! Update the vector DVEC2
    ! -------------------------------------------------------
    DO K = 1, ndis_glob
        IBC = (K-1)*3
        DVEC2(K,1) = DVEC(K,1)
        DVEC2(K,2) = DVEC(K,2)
        DVEC2(K,3) = DVEC(K,3)
        IF (i_kode(K).EQ.10) THEN
          DVEC2(K,1) = 0.
          DVEC2(K,2) = 0.
          DVEC2(K,3) = 0.
        ELSEIF (i_kode(K).EQ.11.OR.i_kode(K).EQ.12) THEN
          DVEC2(K,3) = 0.
        ELSEIF (i_kode(K).EQ.13) THEN
          DVEC2(K,2) = 0.
          DVEC2(K,3) = 0.
        ELSEIF (i_kode(K).EQ.14) THEN
          DVEC2(K,1) = 0.
          DVEC2(K,3) = 0.
        ELSEIF (i_kode(K).EQ.15) THEN
          DVEC2(K,1) = 0.
          DVEC2(K,2) = 0.
        ENDIF
    ENDDO

  end subroutine update_DVEC2


  subroutine update_XMATRIX_DISPL_DVECI()
    implicit none
    
    integer :: K, IBC

    ! -------------------------------------------------------
    ! Update the column of XMATRIX corresponding to the displacements
    ! with the elements of DVECI
    ! -------------------------------------------------------
    DO K = 1, ndis_glob
        IBC = (K-1)*3
        XMATRIX(IBC+1,NUM_Ds_SAVED) = DVECI(K,1)
        XMATRIX(IBC+2,NUM_Ds_SAVED) = DVECI(K,2)
        XMATRIX(IBC+3,NUM_Ds_SAVED) = DVECI(K,3)
    ENDDO

  end subroutine update_XMATRIX_DISPL_DVECI


  subroutine sum_DVECI_DVEC2()
    implicit none
    integer :: K
    ! -------------------------------------------------------
    ! Sum the vectors DVECI and DVEC2
    ! -------------------------------------------------------
    DO K = 1, ndis_glob
        DVECI(K,1) = DVECI(K,1) + DVEC2(K,1)
        DVECI(K,2) = DVECI(K,2) + DVEC2(K,2)
        DVECI(K,3) = DVECI(K,3) + DVEC2(K,3)
    ENDDO
  end subroutine sum_DVECI_DVEC2


  subroutine print_DVECI()
    implicit none
    integer :: K
    ! -------------------------------------------------------
    ! print DVECI in the terminal
    ! -------------------------------------------------------
    WRITE(*,*) '       ID,     Ds (cm),     Dd (cm),    Dn (cm)'
    DO K = 1, ndis_glob
        WRITE(*,*) K, DVECI(K,1), DVECI(K,2), DVECI(K,3)
    ENDDO
  end subroutine print_DVECI


  subroutine write_xmatrix(fname)
    implicit none

    character(7) :: fname
    integer :: K, L

    ! -------------------------------------------------------
    ! Dump XMATRIX to ASCII file for inspection
    ! -------------------------------------------------------
    open(unit=20, file=fname, status='unknown')
    WRITE(20,*) ' --- XMATRIX ---'
    WRITE(20,*) ''
    WRITE(20,*) 'SHAPE OF XMATRIX', SHAPE(XMATRIX)
    WRITE(20,*) 'LBOUND:', LBOUND(XMATRIX,1), LBOUND(XMATRIX,2)
    WRITE(20,*) 'UBOUND:', UBOUND(XMATRIX,1), UBOUND(XMATRIX,2)
    WRITE(20,*) ''
    WRITE(20,*) 'INDEX OF COLUMN CONTAINNIG B.C.S:'
    WRITE(20,*) NUM_Ds_SAVED
    WRITE(20,*) ''
    WRITE(20,*) 'BOUNDARY CONDITION VECTOR:'
    DO K = LBOUND(XMATRIX,1), UBOUND(XMATRIX,1)
      WRITE(20, *) K, XMATRIX(K,NUM_Ds_SAVED)
    ENDDO
    WRITE(20,*) ''
    WRITE(20,*) 'XMATRIX:'
    DO K = LBOUND(XMATRIX,1), UBOUND(XMATRIX,1)
      DO L = LBOUND(XMATRIX,2), UBOUND(XMATRIX,2)
        WRITE(20,*) K, L, XMATRIX(K,L)
      ENDDO
    ENDDO
    close(20)
  
  end subroutine write_xmatrix


  subroutine write_amatrix(fname)
    implicit none

    character(7) :: fname
    integer :: K, L

    ! -------------------------------------------------------
    ! Dump AMATRIX to ASCII file for inspection
    ! -------------------------------------------------------
    open(unit=20, file=fname, status='unknown')
    WRITE(20,*) ' --- AMATRIX ---'
    WRITE(20,*) ''
    WRITE(20,*) 'SHAPE OF AMATRIX', SHAPE(AMATRIX)
    WRITE(20,*) 'LBOUND:', LBOUND(AMATRIX,1), LBOUND(AMATRIX,2)
    WRITE(20,*) 'UBOUND:', UBOUND(AMATRIX,1), UBOUND(AMATRIX,2)
    WRITE(20,*) ''
    WRITE(20,*) 'AMATRIX:'
    DO K = LBOUND(AMATRIX,1), UBOUND(AMATRIX,1)
      DO L = LBOUND(AMATRIX,2), UBOUND(AMATRIX,2)
        WRITE(20,*) K, L, AMATRIX(K,L)
      ENDDO
    ENDDO
    close(20)
  
  end subroutine write_amatrix

end module global_arrays
