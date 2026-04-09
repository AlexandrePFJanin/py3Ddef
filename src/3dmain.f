	SUBROUTINE COMPUTE3DDEF(xg,yg,zg,xdis,ydis,zdis,
     &              length_dis,width_dis,strike_dis,dip_dis,
     &              input_kode,input_ss,input_ds,input_ts,
     &              input_fcode,input_sdrop,input_rhoLitho,
     &              input_rhoFluid,input_cohes,input_disfric,
     &              input_V,input_E,input_friction,
     &              oflag_orient,oflag_invariant,
     &              oflag_failure,oflag_ugrad,
     &              oflag_debug,
     &              maxiter,tolsolver,
     &              input_BG_flag,input_BG,
     &              npts,ndis,
     &              inlout_displ,inlout_stress,inlout_strain,
     &              inlout_orient,inlout_failure,
     &              inlout_invariant,inlout_ugrad,
     &              inlout_elements,inlout_fstatus,
     &              solutionStatus)
C******************************************************************************
C Program to calculate the deformation field using a 3-D
C boundary element algorithm.
C
C
C                        *** Some definitions ***
C   Plane = Element = Fault or fault segment with constant strike and dip.
C   Sub-element = each plane can be subdivided into smaller sub-elements.
C                 For each plane these have the same dimensions equal to
C                 the total length/number of sub-elements along strike and
C                 the total width/number of sub-elements along the dip direction.
C   Global coordinates = Coordinate system for all planes.  X is E, Y is N, and
C 	          +Z is vertical pointing up wrt the free surface.
C   Reference point,Xo,Yo,Zo = For each plane determines the relative position
C		  of all planes in the global coordinate system. 
C		  For each plane the user specifies it as the corner
C		  of the plane closest to the surface such that 
C                 looking away from the reference
C		  point the fault strikes with the hanging wall on the right.
C		  Internally Xo,Yo,Zo is redefined as the same end of the
C		  fault plane but at the lower corner.
C   Local coordinates = Coordinate system for each plane as defined by Okada;
C                 local XL axis is along strike,YL axis perpendicular to strike,
C      		  ZL axis is vertical. The system is right handed such that
C	          looking in the strike direction (+XL) the hanging-wall is on 
C	          the right and the +YL axis points into the footwall.
C                 The reference 
C                 point,Xo,Yo,Zo (specified in global coordinates) as specified 
C                 by the user is the corner of the plane closest to the surface 
C                 such that looking away from the reference point the fault 
C                 strikes in the +XL direction with the hanging wall on the right.
C 	Planar (in-plane) coordinates = Coordinate system for each plane such 
C                 that XP is along the strike direction (=XL), YP is along the 
C                 dip direction and is positive going up-dip from the internal 
C                 Xo,Yo,Zo. +ZP is
C		  normal to the plane with positive defined so that the system
C		  is right handed.
C
C
C
C                      *** Boundary Condition Codes ***
C Boundary conditions currently allowed:     
C	KODE	strike (shear)	dip (shear)		normal
C	 1	displ.		displ.			displ.
C	 2	stress		stress			stress
C	 3	stress		displ.			stress
C	 4	displ.		stress			stress
C	 5	stress		stress			displ.
C	10	relative	relative		relative 
C		displ.		displ.			displ.
C	11	angle 		stress (along	        relative 
C		                angle)			displ.
C	12	stress		stress			relative
C		                                        displ.
C	13	stress		relative		relative 
C	                	displ.			displ.
C	14	relative	stress			relative 
C		displ.			                displ.
C  Note that 11 and 12 are the same; for KODE=11 the user specifies
C  the direction of maximum shear as an angle measured from the
C  strike direction which is then resolved into shear stress conditions
C  along the strike and dip directions.
C
C
C                             *** Directions ***
C Fault (element) orientations are specified as follows:
C   1. When looking in the strike direction the hanging wall is on the right.
C   2. The strike angle is defined clockwise from north.
C   3. Dip is defined from the horizontal, from 0 to 90 degrees.
C
C Fault (element) relative displacement conventions are as follows:
C Disloc Component	Sign		Faulting Mode
C	   Strike		>0			left-lateral
C	    			<0			right-lateral
C	    Dip		    >0			thrust
C	  	 		    <0			normal
C	  Tensile (=normal to plane)
C                   >0			opening
C                   <0			closing
C
C Displacement boundary conditions refer to displacement of the footwall side
C of the element.
C
C                             *** Units ***
C  Specify input element (fault) dimensions in km, displacements in cm.
C  Output displacements with have units of cm, strains will be dimensionless.
C  Stresses (input and output) have units of whatever the Young's modulus
C  is given in.
C
C
C
C                             *** Misc. ***
C Stresses are written as a vector internally as
C  Vector Index		Stress Component
C  		1				Sxx
C		2				Sxy
C		3				Sxz
C		4				Syy
C		5				Syz
C		6				Szz
C
C  The maximum shear stress (magnitude and orientation) is given by tmax, tmaxo,
C  and tmaxa, each determined in subroutine TOPLANEa. tmaxo are components along
C  x and y for QUIVER in MATLAB, and tmaxa is the pitch angle (or rake) of tmax.
C
C****************************************************************************
C     load module for global arrays: Added A.JANIN 13.02.2026
      use global_arrays
      use global_outputs
      use global_inputs

C Material constants
C     Lame parameter (lambda)
      REAL*4 DMULT
C     1-(Vp/Vs)**2 = (lambda+mu)/(lambda+2*mu)
      REAL*4 ALPHA
C     rigidity
      REAL*4 XMU
C     2*rigidity
      REAL*4 XMU2
      COMMON/CONSTANTS/ALPHA,XMU,XMU2,DMULT

C Displacement and stress tensors, temporary storage for 
C influence coefficient matrices' calculations
C      matrix for displacements
      REAL*4 DSPL 
C      matrix for stresses
      REAL*4 STR
      COMMON/TEMPS/DSPL(3,3),STR(6,3)

C Displacement and stress tensors, storage for 
C influence coefficient matrices' calculations
      REAL*4 UCOEF 
      REAL*4 SCOEF
      COMMON/COEFS/ UCOEF(3,3),SCOEF(6,3)
C ** Added to main program on Aug. 29, 2000
C  Displacement gradient tensor values
      REAL*4 DGRAD
      COMMON/DGRADS/DGRAD(3,3,3)
		
C Common blocks s for Okada's routines
      REAL*8 ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,
     &       SDSD,CDCD,SDCD,S2D,C2D  
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,
     &           SDSD,CDCD,SDCD,S2D,C2D  
      REAL*8 XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,  
     &       X11,Y11,X32,Y32,EY,EZ,FY,FZ,GY,GZ,HY,HZ 
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,  
     &       X11,Y11,X32,Y32,EY,EZ,FY,FZ,GY,GZ,HY,HZ                                
		
C  Common blocks for output routines
      REAL*4 STRESS,EWE,dgrten,FRICTION,tmax,tmaxo,tmaxa
      COMMON/MIKE/STRESS(6),EWE(3),dgrten(3,3),V,E,FRICTION,
     &  tmax,tmaxo(2),tmaxa

      REAL*4 SSTRIKE,CSTRIKE,SDYP,CDYP
      COMMON/TOPLN/SSTRIKE,CSTRIKE,SDYP,CDYP

C  Flag and common block for background deformation field
      CHARACTER*4 BFLAG
      REAL*4 BSTRESS
      COMMON/BKGRND/BSTRESS(9),BFLAG

C Added by Alex. JANIN 16.06.24, define the input/output arguments of COMPUTE3DDEF
      integer npts, ndis
      real*8, intent(in)  :: xg(npts), yg(npts), zg(npts)
      real*8, intent(out) :: inlout_displ(npts,3)
      real*8, intent(out) :: inlout_stress(npts,6)
      real*8, intent(out) :: inlout_strain(npts,6)
      real*8, intent(out) :: inlout_failure(npts,6)
      real*8, intent(out) :: inlout_orient(npts,9)
      real*8, intent(out) :: inlout_invariant(npts,4)
      real*8, intent(out) :: inlout_ugrad(npts,3,3)
      real*8, intent(out) :: inlout_elements(ndis,12)
      integer*4, intent(out) :: inlout_fstatus(ndis)
      integer*4, intent(out) :: solutionStatus(2)
      real*4, intent(in)  :: xdis(ndis), ydis(ndis), zdis(ndis)
      real*4, intent(in)  :: length_dis(ndis), width_dis(ndis)
      real*4, intent(in)  :: strike_dis(ndis), dip_dis(ndis)
      integer, intent(in) :: input_kode(ndis)
      real*4, intent(in)  :: input_sdrop(ndis)
      real*4, intent(in)  :: input_rhoLitho(ndis)
      real*4, intent(in)  :: input_rhoFluid(ndis)
      real*4, intent(in)  :: input_cohes(ndis)
      real*4, intent(in)  :: input_disfric(ndis)
      integer, intent(in) :: input_fcode(ndis)
      integer accepted_kod(12)
      real*4, intent(in)  :: input_ss(ndis), input_ds(ndis)
      real*4, intent(in)  :: input_ts(ndis)
      real*4, intent(in)  :: input_V, input_E, input_friction
      integer*4, intent(in) :: maxiter
      real*4, intent(in)    :: tolsolver
      character*4, optional :: input_BG_flag
      real*4, optional      :: input_BG(9)

      logical :: oflag_orient
      logical :: oflag_invariant
      logical :: oflag_failure
      logical :: oflag_ugrad
      logical :: oflag_debug

      INTEGER*4 :: MAX_ELEM

C   scaling from degrees to radians
      DATA TORAD/.0174532925199/

C   ------- SUBROUTINE START ---------

      MAX_ELEM = ndis

      accepted_kod = [1,2,3,4,5,6,10,11,12,13,14,15]

      ! Load output flags in the module
      o_orient    = oflag_orient
      o_invariant = oflag_invariant
      o_failure   = oflag_failure
      o_ugrad     = oflag_ugrad

      ! Dynamic allocation of global array & output sizes
      CALL allocate_global_arrays(npts,ndis)
      CALL allocate_global_outputs(npts,ndis)
      CALL allocate_global_inputs(ndis)
      debug = oflag_debug ! export to global module

      ! Transfer grid data
      xfgrid = xg
      yfgrid = yg
      zfgrid = zg

      ! Transfer input b.c.s and frictional parameters
      i_kode    = input_kode
      i_ss      = input_ss
      i_ds      = input_ds
      i_ts      = input_ts
      i_fcode   = input_fcode
      i_sdrop   = input_sdrop
      i_rhoLitho = input_rhoLitho
      i_rhoFluid = input_rhoFluid
      i_cohes   = input_cohes
      i_disfric = input_disfric

      ! Transfer friction solver parameter
      i_maxiter   = maxiter
      i_tolsolver = tolsolver

      PRINT*, ''
      PRINT*, ' ------ COMPUTE3DDEF ------------------------'
      print*, ''

      NPLANE = ndis
      V = input_V
      E = input_E
      FRICTION = input_friction

      if (present(input_BG_flag)) then
          if (input_BG_flag.eq.'STRE'.or. 
     &        input_BG_flag.eq.'stre') then
              BFLAG = 'STRE'
          elseif (input_BG_flag.eq.'STRA'.or. 
     &            input_BG_flag.eq.'stra') then
              BFLAG = 'STRA'
          elseif (input_BG_flag.eq.'DISP'.or. 
     &            input_BG_flag.eq.'disp') then
              BFLAG = 'DISP'
          else
              BFLAG = 'NONE'
          endif
      else
          BFLAG = 'N'
      endif

      print*, ' MATERIAL PROPERTIES:'
      write(*,100) NPLANE,V,E,FRICTION
 100  FORMAT('   -> PLANES=',I2,' POISSIONS RATIO=',
     1      F5.3,' YOUNGS MOD: ',G10.4,' COEF. FRICTION=',F5.3)

C    - Calculate needed material constants
C      1-(Vs/Vp)**2 = (lambda+mu)/(lambda+2*mu)
      ALPHA=.5/(1.-V)
C     rigidity
      XMU=E*.5/(1.+V)
      XMU2=XMU*2.
      DMULT=V/(1.-2.*V)

C       Added by JANIN
      PRINT*, ''
      PRINT*, ' SIZE OF THE COMPUTATION GRID:'
      PRINT*, '   -> ',npts

C    - Check input KODE: Added by A.JANIN
      DO I=1,ndis
          if (any(input_kode(I)==accepted_kod)) then
              CONTINUE
          else
              print*, ''
              print*, ' ================= ERROR ================='
              print*, ' >>  ERROR IN KODE OF THE ELEMENT: ', I
              print*, '     AVAILABLE VALUES ARE:'
              print*, accepted_kod
              print*, ' >>  PROGRAM INTERRUPTED'
              print*, ' ========================================='
              STOP
          endif
      ENDDO

C    - Read in background deformation if desired
      IF(BFLAG.EQ.'STRE'.OR.BFLAG.EQ.'stre') then
C         read in background stresses
          DO J=1,6
              BSTRESS(J) = input_BG(J)
          END DO
          DO J=7,9
              BSTRESS(J) = 0.
          END DO
      ELSEIF(BFLAG.EQ.'STRA'.OR.BFLAG.EQ.'stra') then
C         read in background strains
          DO J=1,6
              BSTRESS(J) = input_BG(J)
          END DO
          DO J=7,9
              BSTRESS(J) = 0.
          END DO
C         convert to stresses &  read in background strains, rotations
          CALL FROM_STRAIN
      ELSEIF(BFLAG.EQ.'DISP'.OR.BFLAG.EQ.'disp') then
C         read in background displacement gradients
          DO J=1,9
              BSTRESS(J) = input_BG(J)
          END DO
          print*, ' INIT BACKGROUND DEFORMATION: ', BSTRESS
C         convert to stresses and rotations
          CALL FROM_DISPL
      ELSE
          DO J=1,9
              BSTRESS(J) = 0.
          END DO
      ENDIF

      print*, ''
      print*, ' BACKGROUND DEFORMATION: ', BFLAG
      print*, '   -> ',BSTRESS
C *** Read in the element structure of each plane ***      

      print*, ''
      print*, ' PLANE DETAILS:'
      write(*,300)
 300  FORMAT(
     1      '  PLANE-ID   ORG:XO      YO      ZO    ',
     2      ' WD:STK     DIP    #SUB-EL:STK DIP  ',
     3      'STK    DIP')

C     counter for max number of sub-elements along strikes
      MAXB1=0
C     counter for max number of sub-elements along dips
      MAXB2=0
      
      ! Global detection of any frictional element
      any_frictional = .FALSE.  ! Init

      ! Iterate on each plane
      DO  N=1,NPLANE

          ! Added A.JANIN 12.03.26: frictional status
          IF (i_fcode(N) .GT. 0 .AND. input_kode(N) .EQ. 2) THEN
              ! Frictional element
              any_frictional = .TRUE.
              element_fstatus(N) = -1 ! frictional but not determined
          ELSEIF (i_fcode(N) .GT. 0 .AND. input_kode(N) .NE. 2) THEN
              ! Raise an error and stop the program
              print*, ''
              print*, ' ================= ERROR ================='
              print*, ' >>  ERROR IN KODE OF THE ELEMENT: ', N
              print*, ' THE ELEMENT IS DEFINED AS FRICTIONAL'
              print*, ' (FCODE>0) AND THEREFORE IS ONLY'
              print*, ' COMPATIBLE WITH KODE=2 AND B.C.S'
              print*, ' REFLECTING PRE-EXISITING STRESSES.'
              print*, ' >>  PROGRAM INTERRUPTED'
              print*, ' ========================================='
              STOP
          ELSE
              element_fstatus(N) = -2 ! not frictional
          ENDIF

          X = xdis(N)
          Y = ydis(N)
          Z = zdis(N)
          W_STK = length_dis(N)
          W_DIP = width_dis(N)
C         NBX1 and NBX2 set fixed to 1 for the wrapping by A.JANIN
          NBX1(N) = 1
          NBX2(N) = 1
          STRIKE  = strike_dis(N)
          DIP(N)  = dip_dis(N)
C         Convert to sub-element widths internally
          BWX1(N)=W_STK/NBX1(N)
          BWX2(N)=W_DIP/NBX2(N)
C         Calculate sines and cosines of the strike and dip	
          DP=DIP(N)*TORAD
          CDIP(N)=COS(DP)
          SDIP(N)=SIN(DP)
          IF(DIP(N).EQ.180) SDIP(N)=0
          IF(DIP(N).EQ.180) CDIP(N)=-1
          STK=STRIKE*TORAD
          C(N)=COS(STK)
          S(N)=SIN(STK)
C         Internally the origin is the lower left corner of the plane
          TMP=W_DIP*CDIP(N)
          XO(N)=X+TMP*C(N)
          YO(N)=Y-TMP*S(N)
C         reference Z point is positive
          ZO(N)=-(Z-W_DIP*SDIP(N))
	
          write(*,200) N,X,Y,Z,W_STK,W_DIP,
     1            NBX1(N),NBX2(N),STRIKE,DIP(N)
 200      FORMAT(1x,I4,7x,3(f7.2,1x),3x,2g7.2,9X,I2,
     1                 2x,I2,3x,F6.2,1X,F6.2)

C         Find maximum number of sub-elements in each direction for 
C         array dimensioning flexibility
          MAXB1=MAX0(NBX1(N),MAXB1)
          MAXB2=MAX0(NBX2(N),MAXB2)

C         Calculate the transformation matrices to convert the 
C         influence coefficient matrices from global to in-plane coordinates
          CALL MK_TRANS(N,S(N),C(N),SDIP(N),CDIP(N))

          ! Compute the Z coordinate of the center of the patch
          ZCE(N) = Z - 0.5*SDIP(N)*W_DIP
      END DO

C     Check that the dimensions are not exceeded 
      ITOTAL=MAXB1*MAXB2*NPLANE

C   - Finish the remaining calculations in a subroutine to
C     allow flexibility in how planes and sub-elements are specified.
      CALL DO_ALL(MAXB1,MAXB2,NPLANE,ISPACE,SPACE,
     &            input_kode,input_ss,input_ds,input_ts,
     &            npts,ndis)

      PRINT*, '--------------------------------------------'
      PRINT*, ''

      ! Transfer data from global output arrays to output fields 
      inlout_displ     = goutarray_displ
      inlout_stress    = goutarray_stress
      inlout_strain    = goutarray_strain
      inlout_failure   = goutarray_failure
      inlout_orient    = goutarray_orient
      inlout_invariant = goutarray_invariant
      inlout_ugrad     = goutarray_ugrad
      inlout_elements  = goutarray_elements
      inlout_fstatus   = element_fstatus
      solutionStatus   = codeStatus

C     deallocate allocated memory to avoid memory leak in ipython terminal
      call deallocate_global_arrays()
      call deallocate_global_outputs()
      call deallocate_global_inputs()
      
      RETURN
      END

		 

      SUBROUTINE DO_ALL(MAXB1,MAXB2,NPLANE,KODE,BC,                        
     &              input_kode,input_ss,input_ds,input_ts,
     &              npts,ndis)

C******************************************************************************
C  This routine controls all the calculating now that almost everything 
C  is read in.  It is separated from the main program so that the numbers
C  of planes and sub-elements is flexible (array dimensions of KODE
C  defined via variables passed from the main routine).
C******************************************************************************

      use global_arrays
      use global_outputs
      use global_inputs

      INTEGER*4 FIXED, NUME, NUM_Ds
      INTEGER*4 FIXED_INIT, NUME_INIT, NUM_Ds_INIT
C   - FIXED =-1 if all elements have fixed rel. displ.
C           = 0 if no elements have fixed rel. displ.
C           = 1 if some but not all elements have fixed rel. displ.

      INTEGER*4 MAXB1, MAXB2, NPLANE
      INTEGER*4 NITER
      LOGICAL CONVERGED
      REAL*4 tol

C     boundary condition (b.c.) code
      INTEGER*4 KODE
C     b.c. in strike,dip,normal (tensile) directions
      REAL*4 BC
      DIMENSION KODE(MAXB1,MAXB2,NPLANE),BC(3,MAXB1,MAXB2,NPLANE)

C     Added by A.JANIN
      integer npts, ndis
	  integer, intent(in) :: input_kode(ndis)
	  real*4, intent(in)  :: input_ss(ndis), input_ds(ndis)
	  real*4, intent(in)  :: input_ts(ndis)

C     initially assume all elements have fixed rel. displ.
      FIXED = -1

      IF (any_frictional) THEN
          PRINT*, ''
          PRINT*, ' PLANE-ID,  KODE,  taus_i,',
     1 '  taud_i,  sigman_i,  FCODE,  SDROP,  C,  mu,  ',
     1 'rhoLitho,  rhoFluid' 
      ENDIF

      DO 10 N=1,NPLANE
			
C       For each element along strike
		DO 10 I=1,NBX1(N)

			DO 11 J=1,NBX2(N)
	         
				 KODE(I,J,N) = input_kode(N)
				 BC(1,I,J,N) = input_ss(N)
				 BC(2,I,J,N) = input_ds(N)
				 BC(3,I,J,N) = input_ts(N)

C       		some or no elements have unconstrained rel. displ.
				IF (KODE(I,J,N).NE.10) THEN
					FIXED = 0
				ENDIF

                IF (any_frictional) THEN
                    WRITE(*,100) N,KODE(I,J,N),
     &              BC(1,I,J,N),BC(2,I,J,N),BC(3,I,J,N),
     &              i_fcode(N),i_sdrop(N),i_cohes(N),i_disfric(N),
     &              i_rhoLitho(N),i_rhoFluid(N)

100                 FORMAT(I4,1X,I3,1X,3(G12.4),1X,I3,1X,5(G12.4))
                ENDIF

11				IF(KODE(I,J,N).EQ.11)
     &	 		    CALL RESOLVE(BC(1,I,J,N),BC(2,I,J,N))

10    	CONTINUE
		print*, ''
      
        IF(FIXED.GE.0) THEN
C       	some or no elements have fixed rel. displ. but not all
          	DO 12 N=1,NPLANE
	  			DO 12 I=1,NBX1(N)
	  				DO 12 J=1,NBX2(N)
      	  				IF(KODE(I,J,N).GE.10) FIXED=1 
12        	CONTINUE
        ENDIF
		
C     *** Calculate influence coefficient matrix
      IF (FIXED.GE.0) THEN
          CALL CALCUL(NUME,MAXB1,MAXB2,NPLANE,KODE)
          ! save the original coefficent matrix in AMATRIX
          AMATRIX = XMATRIX(:, 1:3*NDIS)
      ENDIF
	  
      ! Before any iterations, solve the initial b.c.s. parameters:
      ! will be restored and reinjected in the SOLVE routine.
      NUME_INIT   = NUME
      NUM_Ds_INIT = NUM_Ds
      FIXED_INIT  = FIXED
      NUM_Ds_SAVED= 0

      ! equalize the frictional status
      element_iter_fstatus = element_fstatus

C     *** Solve for shear and normal displacement discontinuities.
      ! Added A.JANIN 20.03.2026 Iterative scheme for frictional elements
      IF (any_frictional) THEN
          CONVERGED = .FALSE.
          goutflag_converged = .FALSE.
          NITER = 1
          DO WHILE (.NOT.CONVERGED .AND. NITER.LE.i_maxiter)
              WRITE(*,*) ''
              WRITE(*,'(A,I0,A)') '=== ITERATION ', NITER, ' ==='
              CALL SOLVE(NUME,NUM_Ds,FIXED,MAXB1,MAXB2,NPLANE,KODE,BC)
              ! Superposition of the displacements
              IF (NITER.EQ.1) THEN
                  DVECI = DVEC           ! init with every displacement
              ELSE
                  CALL update_DVEC2(KODE(1,1,:))
                  CALL sum_DVECI_DVEC2() ! avoid repeating fixed relative displacement already included
              ENDIF
              if (debug) then
                  ! -> Visualize the current state of the displacement vector
                  WRITE(*,*) ''
                  WRITE(*,*) ' [DEBUG] (Current displacement vector)'
                  CALL print_DVECI()
              endif
              ! Test the convergence: compute the maximum change in DVECI, have to be below tolsolver
              tol = MAXVAL(ABS(DVECI - DVEC0)) ! tolerance
              IF (tol.LE.i_tolsolver) THEN
                  CONVERGED = .TRUE.
                  goutflag_converged = .TRUE.
                  WRITE(*,*) ''
                  WRITE(*,*) 'ITERATIVE SOLVER REACHED CONVERGENCE!'
                  WRITE(*,*) '  -> FINAL NITER =',NITER
              ENDIF
              ! Decision
              IF (.NOT.CONVERGED .AND. NITER.LT.i_maxiter) THEN
                  ! Reset
                  NUME   = NUME_INIT
                  NUM_Ds = NUM_Ds_INIT
                  FIXED  = FIXED_INIT
                  NUM_Ds_SAVED = 0
                  ! Iterate
              ENDIF
              NITER  = NITER + 1
              DVEC0  = DVECI
          ENDDO
          ! Update XMATRIX with the displacement from DVECI
          CALL update_XMATRIX_DISPL_DVECI()
          NITER = NITER-1
          ! Fill the exiting code status vector
          codeStatus(2) = NITER ! final Niter
          IF (.NOT.CONVERGED) THEN
              WRITE(*,*) ''
              WRITE(*,*) '************* WARNING ***************'
              WRITE(*,*) 'N-ITER MAX REACHED BEFORE CONVERGENCE'
              WRITE(*,*) '*************************************'
              codeStatus(1) = 0 ! converged?, 1=True ,0=No
          ELSE
              codeStatus(1) = 1 ! converged?, 1=True ,0=No
          ENDIF
          WRITE(*,*)
          WRITE(*,*) '================='
      ELSE
          ! No frictional element, no iterative scheme
          CALL SOLVE(NUME,NUM_Ds,FIXED,MAXB1,MAXB2,NPLANE,KODE,BC)
          goutflag_converged = .TRUE.
          ! Fill the exiting code status vector
          codeStatus(1) = 1 ! converged?, 1=True ,0=No
          codeStatus(2) = 1 ! final Niter
      ENDIF

      PRINT*, ''
      PRINT*, 'SOLVING ON THE GRID'
      PRINT*, ''

      CALL GRID(NUM_Ds,NPLANE,npts)

      if (debug) then
          ! -> Export the coeff matrix AMATRIX
          !    before the inversion (with stress b.c.s)
          call write_amatrix('amatrix')
      endif

      PRINT*, 'COMPUTE 3D DEFORMATIONS: DONE!'

      RETURN
      END


	SUBROUTINE FROM_DISPL
C******************************************************************************
C  convert background displacement gradients to stresses and rotations    
C  Upon input BSTRESS(J),J=1,9 
C   = dUx/dx,dUx/dy,dUx/dz,dUy/dx,dUy/dy,dUy/dz,dUz/dx,dUz/dy,dUz/dz
C  Upon output BSTRESS(J),J=1,6
C   = Sxx,Sxy,Sxz,Syy,Syz,Szz
C******************************************************************************
	CHARACTER*4 BFLAG
	REAL*4 BSTRESS
	COMMON/BKGRND/BSTRESS(9),BFLAG
	  
	REAL*4 ALPHA,XMU,XMU2,DMULT
      	COMMON/CONSTANTS/ALPHA,XMU,XMU2,DMULT

C	( dUy/dx - dUx/dy)/2
	ROT_Z=.5*(BSTRESS(4)-BSTRESS(2))
C	( dUx/dz - dUz/dx)/2
	ROT_Y=.5*(BSTRESS(3)-BSTRESS(7))
C	( dUz/dy - dUy/dz)/2
	ROT_X=.5*(BSTRESS(8)-BSTRESS(6))

C       Convert to stresses - see Jaeger and Cook, pg. 110
      	DILAT=(BSTRESS(1)+BSTRESS(5)+BSTRESS(9))*DMULT

      	BSTRESS(1)=XMU2*(DILAT+BSTRESS(1))
      	BSTRESS(2)=XMU*(BSTRESS(2)+BSTRESS(4))
      	BSTRESS(3)=XMU*(BSTRESS(3)+BSTRESS(7))
      	BSTRESS(4)=XMU2*(DILAT+BSTRESS(5)) 
      	BSTRESS(5)=XMU*(BSTRESS(6)+BSTRESS(8))
      	BSTRESS(6)=XMU2*(DILAT+BSTRESS(9))

      	BSTRESS(7)=ROT_X
      	BSTRESS(8)=ROT_Y
      	BSTRESS(9)=ROT_Z
	  
	RETURN
	END



	SUBROUTINE FROM_STRAIN
C******************************************************************************
C  convert background strains to stresses      
C  Upon input BSTRESS(J),J=1,6
C   = Exx,Exy,Exz,Eyy,Eyz,Ezz
C  Upon output BSTRESS(J),J=1,6
C   = Sxx,Sxy,Sxz,Syy,Syz,Szz
C******************************************************************************
	CHARACTER*4 BFLAG
	REAL*4 BSTRESS
	COMMON/BKGRND/BSTRESS(9),BFLAG
	  
	REAL*4 ALPHA,XMU,XMU2,DMULT
      	COMMON/CONSTANTS/ALPHA,XMU,XMU2,DMULT

C	scaling from degrees to radians	  
      	DATA TORAD/.017453293/

      	DILAT=(BSTRESS(1)+BSTRESS(4)+BSTRESS(6))*DMULT
      
      	BSTRESS(1)=XMU2*(DILAT+BSTRESS(1))
      	BSTRESS(2)=XMU2*BSTRESS(2)
      	BSTRESS(3)=XMU2*BSTRESS(3)
      	BSTRESS(4)=XMU2*(DILAT+BSTRESS(4)) 
      	BSTRESS(5)=XMU2*BSTRESS(5)
      	BSTRESS(6)=XMU2*(DILAT+BSTRESS(6))
	  
	DO 2 J=7,9
2	    BSTRESS(J)=TORAD*BSTRESS(J)

	RETURN
	END

	 
	SUBROUTINE MK_TRANS(NS,SSTK,CSTK,SDIP2,CDIP2)
C******************************************************************************
C--- 	Create transformation matrices (from global to local in-plane coordinates)
C	 for plane NS.
C******************************************************************************
	
	use global_arrays
	  
C	index of the current plane
	INTEGER*4 NS
C	sine & cosine of the strike of the NSth plane
	REAL*4 SSTK,CSTK
C	sine & cosine of the dip of the NSth plane	
	REAL*4 SDIP2,CDIP2
	
C - 	UG2P(3,6,NS) transforms displacements from global coords to obtain
C   	displacements in NS in-plane coords.
C				in-plane component	global component
C				strike direction	strike direction 
      UG2P(1,1,NS)=SSTK
C				strike direction	dip direction
      UG2P(1,2,NS)=CSTK
C				strike direction	normal direction
      UG2P(1,3,NS)=0.
C				dip direction		strike direction
      UG2P(2,1,NS)=-CSTK*CDIP2
C				dip direction		dip direction
      UG2P(2,2,NS)=SSTK*CDIP2
C				dip direction		normal direction
      UG2P(2,3,NS)=SDIP2
C				normal direction	strike direction
      UG2P(3,1,NS)=CSTK*SDIP2
C				normal direction	dip direction	
      UG2P(3,2,NS)=-SSTK*SDIP2
C				normal direction	normal direction
      UG2P(3,3,NS)=CDIP2


C - SG2P(3,6,NS) transforms stresses from global coords to obtain stresses
C   in NS in-plane coords.
C						in-plane comp. global component
C							Sxz		Sxx
      SG2P(1,1,NS)=UG2P(3,1,NS)*UG2P(1,1,NS)
C							Sxz		Sxy	
      SG2P(1,2,NS)=UG2P(3,1,NS)*UG2P(1,2,NS)+
     &             UG2P(3,2,NS)*UG2P(1,1,NS)
C	    						Sxz		Sxz
      SG2P(1,3,NS)=UG2P(3,3,NS)*UG2P(1,1,NS)
C							Sxz		Syy
      SG2P(1,4,NS)=UG2P(3,2,NS)*UG2P(1,2,NS)
C	    						Sxz		Syz
      SG2P(1,5,NS)=UG2P(3,3,NS)*UG2P(1,2,NS)
C	                				Sxz		Szz
      SG2P(1,6,NS)=0.

C							Syz		Sxx
      SG2P(2,1,NS)=UG2P(2,1,NS)*UG2P(3,1,NS)
C							Syz		Sxy
      SG2P(2,2,NS)=UG2P(2,1,NS)*UG2P(3,2,NS)+
     &                 UG2P(2,2,NS)*UG2P(3,1,NS)
C							Syz		Sxz	
      SG2P(2,3,NS)=UG2P(2,1,NS)*UG2P(3,3,NS)+
     &                 UG2P(2,3,NS)*UG2P(3,1,NS)
C							Syz		Syy
      SG2P(2,4,NS)=UG2P(2,2,NS)*UG2P(3,2,NS)
C							Syz		Syz
      SG2P(2,5,NS)=UG2P(2,2,NS)*UG2P(3,3,NS)+
     &                 UG2P(2,3,NS)*UG2P(3,2,NS)
C							Syz		Szz
      SG2P(2,6,NS)=UG2P(2,3,NS)*UG2P(3,3,NS)


C							Szz		Szz
      SG2P(3,1,NS)=UG2P(3,1,NS)*UG2P(3,1,NS)
C							Szz		Szz
      SG2P(3,2,NS)=2.*UG2P(3,1,NS)*UG2P(3,2,NS)
C							Szz		Szz
      SG2P(3,3,NS)=2.*UG2P(3,1,NS)*UG2P(3,3,NS)
C							Szz		Szz
      SG2P(3,4,NS)=UG2P(3,2,NS)*UG2P(3,2,NS)
C							Szz		Szz
      SG2P(3,5,NS)=2.*UG2P(3,2,NS)*UG2P(3,3,NS)
C							Szz		Szz
      SG2P(3,6,NS)=UG2P(3,3,NS)*UG2P(3,3,NS)

      RETURN
      END


      SUBROUTINE RESOLVE(BC1,BC2)
C******************************************************************************
C  Subroutine to resolve a shear stress (BC2) in direction BC1 (angle
C  from strike) into shear stresses along in-plane strike and dip directions.
C******************************************************************************
          REAL*4 BC1,BC2,ANGLE,TORAD
	  DATA TORAD/.017453293/
      
          ANGLE=BC1*TORAD
	  BC1=BC2*COS(ANGLE)
	  BC2=BC2*SIN(ANGLE)
	  RETURN
	  END


	SUBROUTINE CALCUL(NUMEBC,MAXB1,MAXB2,NPLANE,KODE)
C******************************************************************************
C Subroutine to calculate influence coefs. These are placed in the matrix
C  XMATRIX(NUMEBC,NUMEBC).
C******************************************************************************  
	
	use global_arrays

	IMPLICIT NONE

	INTEGER*4 KODE, MAXB1, MAXB2, NPLANE
	
	DIMENSION KODE(MAXB1,MAXB2,NPLANE)

	INTEGER*4 NP, NS
	INTEGER*4 JP, IP, IP2, JS, IS

	REAL*8 BX1P, BX2P
	REAL*8 XNP, ZNP, YNP
	REAL*8 X, Y, Z
	REAL*8 AL1, AL2, AW1, AW2
    
	INTEGER*4 NUMEBC, NUMED, IGRAD, IRET

	WRITE(*,'(A)')
     &  ' BEGIN CALCULATION OF INFLUENCE COEFFICIENT MATRIX'    
	  
C   initialize number of boundary conditions
	NUMEBC=0
      
C For each plane NP	  
C     *** loop over all planes ***
      DO 450 NP=1,NPLANE
C   - retrieve sub-element widths, number of sub-elements
      BX1P=BWX1(NP)*.5
      BX2P=BWX2(NP)*.5
		  
C       *** loop over IP,JP sub-elements of plane NP ***
        DO 445 JP=1,NBX1(NP)
        DO 445 IP=1,NBX2(NP)

C    *** Find center of this sub-element in global coordinates ***        
C      - center of sub-element IP,JP in local NP coords	 
C	 center in local strike direction
	 XNP=BX1P*(2*JP-1)
	 IP2=2*IP-1
	 ZNP=-ZO(NP)+BX2P*IP2*SDIP(NP)
C	 center in local dip direction
	 YNP=BX2P*IP2*CDIP(NP)
      	 
C      - center of sub-element IP,JP in global coords
	 Z=ZNP
	 X=S(NP)*XNP-C(NP)*YNP+XO(NP)
	 Y=C(NP)*XNP+S(NP)*YNP+YO(NP)
		 
C	 initialize number of unknown displacement discont.s
	 NUMED=0

C        *** loop over NS planes ***	
         DO 170 NS=1,NPLANE
					       
C 	 - Define center of sub-element IP,JP,NP in local NS coords
	 XNP=S(NS)*(X-XO(NS))+C(NS)*(Y-YO(NS))
	 YNP=-C(NS)*(X-XO(NS))+S(NS)*(Y-YO(NS))
	 ZNP=Z
C        *** Offset center point for absolute displacement boundary
C        conditions to avoid problems with the Greens functions being
c        not being single-valued right on the element.  Offset is
C        on the minus side.  Added August 1999.
         if(KODE(JP,IP,NP).eq.1.or
     &    .KODE(JP,IP,NP).eq.3.or
     &    .KODE(JP,IP,NP).eq.4.or
     &    .KODE(JP,IP,NP).eq.5.or
     &    .KODE(JP,IP,NP).eq.6) YNP=YNP+1.e-6 ! from 1e-4 to 1e-6, A.JANIN 13/03/2026
      		
C	 *** loop over IS,JS sub-elements of plane NS ***
	 DO 115 JS=1,NBX1(NS)
C        - sub-element length range along strike
     	 AL1=FLOAT(JS-1)*BWX1(NS)
         AL2=AL1+BWX1(NS)
	       
  	 DO 115 IS=1,NBX2(NS)
C        - sub-element width range along dip
    	 AW1=FLOAT(IS-1)*BWX2(NS)
	     AW2=AW1+BWX2(NS)
					  
C	 - Calculate 6 stress and 3 displacement influence coefs in NS coords
C	 Put these in arrays STR(6,3),DSPL(3,3)
         IGRAD=0
 	 CALL DCD3D(XNP,YNP,ZNP,ZO(NS),DIP(NS),
     &    AL1,AL2,AW1,AW2,IGRAD,IRET)
	 IF(IRET.NE.0) THEN
		WRITE(*,'(A)') 
     &          ' ERROR CALCULATING INFLUENCE COEFFICIENT '
		STOP
	 ENDIF
		   
C        - Transform 3 (normal) stress and 3 displacement influence coeffs to
C          planar NP coordinates.  Put these in arrays STR(3,3),DSPL(3,3)
 	 CALL COEF2NP(NS,NP)
	
C	 - Create influence coef. matrix appropriate to b.c.s 
C          of sub-element IP,JP,NP
 	 CALL MK_MATRIX(KODE(JP,IP,NP),NUMED,NUMEBC)

C	 increment the displ. discont.s count
	 NUMED=NUMED+3
		   				
115      CONTINUE
C	 *** end of loop over blocks of plane NS ***
 	
170	 CONTINUE
C         *** end of loop over NS planes ***	
        
C	increment the boundary condition count
	NUMEBC=NUMEBC+3
			  
445     CONTINUE 
C       *** end of loop over blocks of plane NP ***

450    CONTINUE
C      *** end of loop over all planes ***

       return
       end



      SUBROUTINE SOLVE(NUME,NUM_Ds,FIXED,MAXB1,MAXB2,NPLANE,KODE,BC)
      ! **************************************************************************
      ! Subroutine to solve for the unknown relative displacement discontinuities.
      ! Rewritten in Fortran 90+ by A.JANIN 25.02.2026
      ! **************************************************************************
      use global_arrays
      use global_inputs
      INTEGER*4 FIXED
  
      INTEGER*4 KODE, MAXB1, MAXB2, NPLANE, SIDE
      REAL*4 BC
      DIMENSION KODE(MAXB1,MAXB2,NPLANE),BC(3,MAXB1,MAXB2,NPLANE)

      CHARACTER*4 BFLAG
      REAL*4 BSTRESS
      COMMON/BKGRND/BSTRESS(9),BFLAG

      DIMENSION BSTR(3)
      
      INTEGER(4) NUME, IBC, NP, I, J, K, JP, IP, L

      ! Save all initial boundary conditions - copy into 
      ! last column of XMATRIX
      IF(FIXED.GE.0) THEN
          ! index of column containing boundary conditions
          NUMEP1=NUME+1
      ELSE
          NUMEP1=1
          NUM_Ds=1
          NUM_Ds_SAVED = NUM_Ds ! keep track of NUM_Ds outside
      ENDIF

      ! *** reset XMATRIX from AMATRIX ***
      CALL reset_xmatrix()

      IBC = 0

      ! *** loop over all planes ***
      DO NP = 1, NPLANE
  
          ! *** correct for background stresses ***
          IF (BFLAG(1:1) .NE. 'N') THEN

              ! Transform background stresses to in-plane NP coords
              DO I = 1, 3
                  BSTR(I) = 0.0
                  DO K = 1, 6
                      BSTR(I) = BSTR(I) + SG2P(I,K,NP) * BSTRESS(K)
                  ENDDO
              ENDDO

          ELSE

              DO I = 1, 3
                  BSTR(I) = 0.0
              ENDDO

          ENDIF

          ! *** Garantee KODE=2 for frictional elements ***
          ! (KODE of frictional elements change after each iteration)
          IF (i_fcode(NP) .GT. 0) THEN
              ! Frictional element
              KODE(1,1,NP) = 2      ! Restore KODE = 2
              SIDE = -1             ! BC = external forces = initial (pre-existing/stored) stresses + background
              ! Restore the B.C.S (inital/pre-existing stresses)
              BC(1,1,1,NP) = i_ss(NP)
              BC(2,1,1,NP) = i_ds(NP)
              BC(3,1,1,NP) = i_ts(NP)
          ELSE
              SIDE = 1              ! substract the background stresses
          ENDIF


          ! *** loop over IP,JP sub-elements of plane NP ***
          DO JP = 1, NBX1(NP)
              DO IP = 1, NBX2(NP)

                  IF (KODE(JP,IP,NP) .EQ. 2) THEN
                      ! 3 stress b.c.s
                      ! Modified A.JANIN 14/03/2026: added SIDE for frictional elements
                      DO J = 1, 3
                          BC(J,JP,IP,NP) = BC(J,JP,IP,NP) -
     &                                     SIDE*BSTR(J)
                      ENDDO

                  ELSEIF (KODE(JP,IP,NP) .EQ. 3) THEN
                      ! strike-shear & normal stress, dip displ. b.c.s
                      BC(1,JP,IP,NP) = BC(1,JP,IP,NP) - BSTR(1)
                      BC(3,JP,IP,NP) = BC(3,JP,IP,NP) - BSTR(3)

                  ELSEIF (KODE(JP,IP,NP) .EQ. 4) THEN
                      ! dip-shear & normal stress, strike displ. b.c.s
                      BC(2,JP,IP,NP) = BC(2,JP,IP,NP) - BSTR(2)
                      BC(3,JP,IP,NP) = BC(3,JP,IP,NP) - BSTR(3)

                  ELSEIF (KODE(JP,IP,NP) .EQ. 5) THEN
                      ! strike & dip-shear stresses, normal displ. b.c.s
                      BC(2,JP,IP,NP) = BC(2,JP,IP,NP) - BSTR(2)
                      BC(1,JP,IP,NP) = BC(1,JP,IP,NP) - BSTR(1)

                  ELSEIF (KODE(JP,IP,NP) .EQ. 6) THEN
                      ! strike & dip-shear displ, normal stress. b.c.s added by A.JANIN
                      BC(3,JP,IP,NP) = BC(3,JP,IP,NP) - BSTR(3)

                  ELSEIF (KODE(JP,IP,NP) .EQ. 11 .OR. 
     &                    KODE(JP,IP,NP) .EQ. 12) THEN
                      ! strike & dip stress, normal displ. fixed b.c.s
                      BC(1,JP,IP,NP) = BC(1,JP,IP,NP) - BSTR(1)
                      BC(2,JP,IP,NP) = BC(2,JP,IP,NP) - BSTR(2)

                  ELSEIF (KODE(JP,IP,NP) .EQ. 13) THEN
                      BC(1,JP,IP,NP) = BC(1,JP,IP,NP) - BSTR(1)

                  ELSEIF (KODE(JP,IP,NP) .EQ. 14) THEN
                      BC(2,JP,IP,NP) = BC(2,JP,IP,NP) - BSTR(2)

                  ELSEIF (KODE(JP,IP,NP) .EQ. 15) THEN
                      ! added by A.JANIN
                      BC(3,JP,IP,NP) = BC(3,JP,IP,NP) - BSTR(3)

                  ENDIF

                  DO J = 1, 3
                      IBC = IBC + 1
                      ! Save all init b.c.s into last column of XMATRIX      
                      XMATRIX(IBC,NUMEP1) = BC(J,JP,IP,NP)

                  ENDDO

              ENDDO
          ENDDO

      ENDDO

      ! Update stress condition on each dislocation (in local coordinate system of the plane)
      ! A.JANIN 25.02.2026
      DO L=1, NPLANE
          IBC = 3*(L-1)
          IF (i_fcode(L) .EQ. 0) THEN
              ! NOT a frictional element:
              SDRIVER(L,1) = XMATRIX(IBC+1,NUMEP1) ! tau_strike(L)
              SDRIVER(L,2) = XMATRIX(IBC+2,NUMEP1) ! tau_dip(L)
              SDRIVER(L,3) = XMATRIX(IBC+3,NUMEP1) ! sigma_n(L)
          ! For frictional element SDRIVER is managed in FRICTIONALDISLOC directly
          ENDIF
      ENDDO

      ! If all boundary conditions are fixed rel. displ. return
      IF (FIXED.LT.0) THEN
          CALL update_DVEC()
          RETURN
      ENDIF

      ! SHRINK: (1) update XMATRIX coefficients if some dislocations are a source of motion
      !         (2) remove columns and rows from influence coefficient matrix that correspond
      !             to element/components with fixed displacement discontinuities
      
      IF (FIXED.GT.0) THEN
          ! If frictional patches, call FRICTIONALDISLOC inside SHRINK, after the
          ! update XMATRIX coefficients if some dislocations are a source of motion
          CALL SHRINK(NUME,MAXB1,MAXB2,NPLANE,KODE,0)
          CALL GET_FIXED(KODE, FIXED, MAXB1, MAXB2, NPLANE)
          ! If all boundary conditions are fixed rel. displ. return: no inversion needed, just return
          IF (FIXED.LT.0) THEN
              ! index of column containing boundary conditions
              NUM_Ds=NUME+1
              NUM_Ds_SAVED = NUM_Ds ! keep track of NUM_Ds outside
              CALL update_DVEC()
              RETURN
          ELSE
          ENDIF
      ELSE
          CALL FRICTIONALDISLOC(NUME,NPLANE,KODE)
          CALL GET_FIXED(KODE, FIXED, MAXB1, MAXB2, NPLANE)
          ! If all boundary conditions are fixed rel. displ. return: no inversion needed, just return
          IF (FIXED.LT.0) THEN
              ! index of column containing boundary conditions
              NUM_Ds=NUME+1
              NUM_Ds_SAVED = NUM_Ds ! keep track of NUM_Ds outside
              CALL update_DVEC()
              RETURN
          ENDIF
          IF (FIXED.GT.0) THEN
              ! FIXED is modified (.GE. 1): call SHRINK (step 2 only) to adapt XMATRIX
              CALL SHRINK(NUME,MAXB1,MAXB2,NPLANE,KODE,2)
          ENDIF
      ENDIF

      CALL GET_FIXED(KODE, FIXED, MAXB1, MAXB2, NPLANE) ! normally, FIXED wil not change here
      ! index of column containing boundary conditions
      NUM_Ds=NUME+1
      NUM_Ds_SAVED = NUM_Ds ! keep track of NUM_Ds outside

      ! Invert XMATRIX to solve for unknown relative displacement discontinuities
      WRITE(*,'(/,A,2(I5,A))') ' STARTING INVERSION OF ',NUME,
     1      ' BY ',NUME,' MATRIX - PLEASE BE PATIENT!'
      CALL INVERT(NUME)
      ! Before INVERT(), the last row of XMATRIX contained the vector stresses S
      ! (last column) at the center of each disloc and the coefficents matrix A
      ! INVERT() do: U = A^{-1} \dot S to compute the displacement U from S and A
      ! and replace S by U in XMATRIX.

      PRINT*, ''
      WRITE(*,'(A)') ' INVERSION DONE!'

      ! no fixed rel. displ. boundary conditions (else, we need to restore XMATRIX before returning)
      IF (FIXED .EQ. 0) THEN
          ! index of column containing boundary conditions
          NUM_Ds=NUME+1
          NUM_Ds_SAVED = NUM_Ds ! keep track of NUM_Ds outside
          CALL update_DVEC()
          RETURN
      ENDIF

      ! Add in any fixed relative displacements that were specified as 
      ! boundary conditions (and removed in SHRINK).

      WRITE(*,*)
      WRITE(*,*) 'RESTORE ANY FIXED RELATIVE DISPL. REMOVED BY SHRINK'

      IBC = 1

      ! *** loop over all planes ***
      DO NP = 1, NPLANE
		  
          ! *** loop over IP,JP sub-elements of plane NP ***
          DO JP = 1, NBX1(NP)
              DO IP = 1, NBX2(NP)

                  IF (KODE(JP,IP,NP).GT.7) THEN
                  ! There are fixed relative displacement boundary conditions on a non-frictional element
                  ! Condition changed by A.JANIN:
                  !    - from .GT.6 to .GT.7 because new b.c. added
                  !    - exclude the frictional dislocations: displacements treated separately

                      IF (KODE(JP,IP,NP) .EQ. 10) THEN
                      ! three fixed displacements - make room for 3 displacements
                      ! start shift at element/component IBC

                          DO J = IBC, NUME
                              JJ  = NUME - J + IBC
                              !   shift down 3
                              JJP = JJ + 3
                              XMATRIX(JJP,NUM_Ds) = XMATRIX(JJ,NUM_Ds)
                          ENDDO

                          IF (i_fcode(NP) .EQ. 0) THEN
                          ! condition (i_fcode(NP) .EQ. 0) added by A.JANIN 04.03.2026 to shift down
                          ! see above) but not udpate frictional patches, done separetly, below.
                              DO J = 1, 3
                                  JJ = IBC + J - 1
                                  XMATRIX(JJ,NUM_Ds) = BC(J,JP,IP,NP)
                              ENDDO
                          ELSEIF (i_fcode(NP) .GT. 0) THEN
                          ! for patches that are frictional and end up having KODE = 10: LOCKED patches
                              DO J = 1, 3
                                  JJ = IBC + J - 1
                                  XMATRIX(JJ,NUM_Ds) = 0.       ! Ds = Dd = Dn = 0
                              ENDDO
                          ENDIF

                          NUME = NUME + 3

                      ELSEIF (KODE(JP,IP,NP) .EQ. 11 .OR.
     &                        KODE(JP,IP,NP) .EQ. 12) THEN
                      ! fixed normal displacement only - make room for 1 displacement
                      ! start shift at element/component IBC+2

                          DO J = IBC+2, NUME
                              JJ  = NUME - J + IBC + 2
                              !   shift down 1
                              JJP = JJ + 1
                              XMATRIX(JJP,NUM_Ds) = XMATRIX(JJ,NUM_Ds)
                          ENDDO

                          IF (i_fcode(NP) .EQ. 0) THEN
                          ! condition (i_fcode(NP) .EQ. 0) added by A.JANIN 04.03.2026 to shift down
                          ! see above) but not udpate frictional patches, done separetly, below.
                              XMATRIX(IBC+2,NUM_Ds) = BC(3,JP,IP,NP)
                          ELSEIF (i_fcode(NP) .GT. 0) THEN
                            ! for patches that are frictional and end up having KODE = 12: SLIDING patches
                              XMATRIX(IBC+2,NUM_Ds) = 0.         ! Dn = 0
                          ENDIF

                          NUME = NUME + 1

                      ELSEIF (KODE(JP,IP,NP) .EQ. 13) THEN
                      ! fixed dip and normal displacement - make room for 2 displacements
                      ! start shift at element/component IBC+1

                          DO J = IBC+1, NUME
                              JJ  = NUME - J + IBC + 1
                              !   shift down 2
                              JJP = JJ + 2
                              XMATRIX(JJP,NUM_Ds) = XMATRIX(JJ,NUM_Ds)
                          ENDDO

                          XMATRIX(IBC+1,NUM_Ds) = BC(2,JP,IP,NP)
                          XMATRIX(IBC+2,NUM_Ds) = BC(3,JP,IP,NP)
                          NUME = NUME + 2

                      ELSEIF (KODE(JP,IP,NP) .EQ. 14) THEN
                      ! fixed strike and normal displacement - make room for 2 displacements
                      ! start shift at element/component IBC+1

                          DO J = IBC+1, NUME
                              JJ  = NUME - J + IBC + 1
                              !   shift down 2
                              JJP = JJ + 2
                              XMATRIX(JJP,NUM_Ds) = XMATRIX(JJ,NUM_Ds)
                          ENDDO

                          XMATRIX(IBC+1,NUM_Ds) = XMATRIX(IBC,NUM_Ds)
                          XMATRIX(IBC,NUM_Ds)   = BC(1,JP,IP,NP)
                          XMATRIX(IBC+2,NUM_Ds) = BC(3,JP,IP,NP)
                          NUME = NUME + 2

                      ELSEIF (KODE(JP,IP,NP) .EQ. 15) THEN
                      ! fixed strike and dip displacement Added by A.JANIN
                      ! start shift at element/component IBC+1

                          DO J = IBC+1, NUME
                              JJ  = NUME - J + IBC + 1
                              ! shift down 2
                              JJP = JJ + 2
                              XMATRIX(JJP,NUM_Ds) = XMATRIX(JJ,NUM_Ds)
                          ENDDO

                          XMATRIX(IBC+2,NUM_Ds) = XMATRIX(IBC,NUM_Ds)
                          XMATRIX(IBC,NUM_Ds)   = BC(1,JP,IP,NP)
                          XMATRIX(IBC+1,NUM_Ds) = BC(2,JP,IP,NP)
                          NUME = NUME + 2

                      ENDIF

                  ENDIF

                  ! index of next element
                  IBC = IBC + 3

              ENDDO
          ENDDO

      ENDDO
      CALL update_DVEC()
      RETURN
      END


      subroutine GET_FIXED(KODE,FIXED,MAXB1,MAXB2,NPLANE)
      ! *************************************************
      ! Compute the value of FIXED according to KODE
      !
      ! FIXED = -1 if all elements have fixed rel. displ.
      !       =  0 if no elements have fixed rel. displ.
      !       =  1 if some but not all elements have fixed rel. displ.
      ! *************************************************

      use global_arrays
      implicit none
      integer, intent(out) :: FIXED
      integer :: FIXED_I
      integer :: KODE,MAXB1,MAXB2,NPLANE
      DIMENSION KODE(MAXB1,MAXB2,NPLANE)

      FIXED_I = FIXED   ! save for verbose output
      IF (MINVAL(KODE) == 10 .AND. MAXVAL(KODE) == 10) THEN
          FIXED = -1
      ELSEIF (ANY(KODE > 10)) THEN
          FIXED = 1
      ELSE
          FIXED = 0
      ENDIF

      !IF (debug) THEN
      !    WRITE(*,*) ''
      !    WRITE(*,*) ' [DEBUG] (determine fixed)'
      !    WRITE(*,*) 'KODE:', KODE
      !    WRITE(*,*) 'FIXED (IN)', FIXED_I
      !    WRITE(*,*) 'FIXED (OUT)', FIXED
      !ENDIF
  
      end subroutine GET_FIXED


      SUBROUTINE FRICTIONALDISLOC(N, NPLANE, KODE)
      ! *************************************************
      ! XMATRIX(IBC+1,NUMEP1) ! tau_strike(L) ext loading
      ! XMATRIX(IBC+2,NUMEP1) ! tau_dip(L) ext loading
      ! XMATRIX(IBC+3,NUMEP1) ! sigma_n(L) ext loading
      ! *************************************************
      use global_arrays
      use global_inputs

      IMPLICIT NONE

      integer(4) :: K, L, M, N, IBC, JBC
      integer(4) :: NUMEP1, NPLANE
      integer(4) :: CORRSIGN_DISP, CORRSIGN_GRFL
      integer(4), intent(inout) :: KODE(1,1,NPLANE)
      real(4) :: DSE, TAU, TAUC, NORM
      real(4) :: TAUS_EXT, TAUD_EXT, SIGN_EXT
      real(4) :: TAUS_DPL, TAUD_DPL, SIGN_DPL
      real(4) :: SIGN_GRAVFLD
      real(4) :: TAUS, TAUD, SIGN

      WRITE(*,*)
      WRITE(*,*) 'RESOLVE FRICTION ON FRICTIONAL ELEMENTS:'

      ! Skip if no frictional patches
      IF (.NOT. any_frictional) THEN
          WRITE(*,*) '  -> NO FRICTIONAL ELEMENT: SKIP'
          RETURN
      ENDIF

      ! Get the displacement vector (corrected from the imposed relative displacement)
      CALL update_DVEC2(i_kode)

      WRITE(*,*) 'PLANE-ID,  (TAU_S,  TAU_D,  SIGMA_N),     ',
     &             '      TAU          TAUC,       MOTION'

      NUMEP1 = N+1

      ! *** loop over all planes ***
      DO L=1, NPLANE

          ! Index on XMATRIX of the first element (-1) of the triplet of b.c.s. for the element L
          IBC = 3*(L-1)

          IF (i_fcode(L) .GT. 0) THEN
          ! Frictional element: no need of a condition on KODE because KODE = 2 is garantee for frictional patches

              ! External stresses
              IF (element_iter_fstatus(L).EQ.-1) THEN
                  ! The element is frictional but it's slip condition are unknown: tau_ext = tau_BG + tau_others
                  TAUS_EXT = XMATRIX(IBC+1,NUMEP1)
                  TAUD_EXT = XMATRIX(IBC+2,NUMEP1)
                  SIGN_EXT = XMATRIX(IBC+3,NUMEP1)
              ELSEIF (element_iter_fstatus(L).EQ.0 .OR.
     &                element_iter_fstatus(L).EQ.1) THEN
                  ! The element is frictional and slip has been already computed and stresses balanced with the stored stresses
                  TAUS_EXT = SSTORED(L,1)
                  TAUD_EXT = SSTORED(L,2)
                  SIGN_EXT = SSTORED(L,3)
              ELSE
                  WRITE(*,*) ' *** ERROR ***'
                  WRITE(*,*) 'INVALID VALUE', element_iter_fstatus(L),
     &                       'FOR THE VARIABLE element_iter_fstatus'
                  WRITE(*,*) '-- PROGRAM STOPS'
                  STOP
              ENDIF

              ! Evaluate the gravity and fluid effects on sigma_n
              SIGN_GRAVFLD = (i_rhoLitho(L)-i_rhoFluid(L))*
     &                       9.80665*ZCE(L)*1000.

              ! Compute the induced stresses due to the resolved slip on the other dislocations
              TAUS_DPL = 0
              TAUD_DPL = 0
              SIGN_DPL = 0
              DO K=1, NPLANE
                  JBC = 3*(K-1)
                  DO M=1, 3
                      ! Use DVEC2 to avoid counting twice driving element (KODE>=10)
                      TAUS_DPL =TAUS_DPL+AMATRIX(IBC+1,JBC+M)*DVEC2(K,M)
                      TAUD_DPL =TAUD_DPL+AMATRIX(IBC+2,JBC+M)*DVEC2(K,M)
                      SIGN_DPL =SIGN_DPL+AMATRIX(IBC+3,JBC+M)*DVEC2(K,M)
                  ENDDO
              ENDDO

              ! prepare the slip (if need to be added)
              IF (element_iter_fstatus(L).EQ.1) THEN
                  CORRSIGN_DISP = 0   ! 
                  CORRSIGN_GRFL = 0   ! already counted
              ELSEIF (element_iter_fstatus(L).EQ.0) THEN
                  CORRSIGN_DISP = 1   ! need to be added
                  CORRSIGN_GRFL = 0   ! already counted
              ELSEIF (element_iter_fstatus(L).EQ.-1) THEN
                  CORRSIGN_DISP = 0   ! unknown
                  CORRSIGN_GRFL = -1  ! need to be added (-1 because ZCE <0)
              ENDIF

              IF (debug) THEN
                  WRITE(*,*) ' [DEBUG] (contribution of Gravity-Fluid)'
                  WRITE(*,*) '       ->', CORRSIGN_GRFL, SIGN_GRAVFLD
              ENDIF

              IF (debug) THEN
200               FORMAT(A9,1X,I2,1X,ES12.4,1X,ES12.4,1X,ES12.4)
                  WRITE(*,*) ' [DEBUG] (contribution of slips)'
                  WRITE(*,200) '       ->', CORRSIGN_DISP, TAUS_DPL,
     &                         TAUD_DPL, SIGN_DPL
              ENDIF

              ! Build the final stresses, take into account (or not) the slip
              TAUS = TAUS_EXT + CORRSIGN_DISP*TAUS_DPL
              TAUD = TAUD_EXT + CORRSIGN_DISP*TAUD_DPL
              SIGN = SIGN_EXT + CORRSIGN_DISP*SIGN_DPL + 
     &                          CORRSIGN_GRFL*SIGN_GRAVFLD

              ! Use the final stresses evaluate the failure condition 
              TAU  = SQRT(TAUS**2 + TAUD**2)
              TAUC = i_disfric(L) * ABS(SIGN) + i_cohes(L)
              DSE = TAU - TAUC

100           FORMAT(I4,1X,A1,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,A1,
     1               1X,ES12.4,1X,ES12.4,3X,A6)

              IF (DSE .LE. 0) THEN

                  WRITE(*,100) L, '(', TAUS, TAUD, SIGN,
     &                ')', TAU, TAUC, 'LOCKED'

                  ! frictionally locked:
                  element_iter_fstatus(L) = 0
                  IF (element_fstatus(L).LE.0) THEN
                      ! update only if unknown or already locked, not if the element slided during a previous iteration
                      element_fstatus(L) = 0
                  ENDIF
                  ! stored stress
                  SSTORED(L,1) = TAUS
                  SSTORED(L,2) = TAUD
                  SSTORED(L,3) = SIGN
                  ! switch to KODE 10 (0,0,0): locked patch
                  KODE(1,1,L) = 10
                  XMATRIX(IBC+1,NUMEP1) = 0.   ! Ds
                  XMATRIX(IBC+2,NUMEP1) = 0.   ! Dd
                  XMATRIX(IBC+3,NUMEP1) = 0.   ! Dn
                  ! update SDRIVER (convention: add 0 if locked)
                  SDRIVER(L,1) = SDRIVER(L,1) + 0.
                  SDRIVER(L,2) = SDRIVER(L,2) + 0.
                  SDRIVER(L,3) = SDRIVER(L,3) + 0.

              ELSE
                  
                  WRITE(*,100) L, '(', TAUS, TAUD, SIGN,
     &                ')', TAU, TAUC, 'SLIDES'

                  ! sliding:
                  element_iter_fstatus(L) = 1
                  element_fstatus(L) = 1        ! update anyway to "sliding" if slides once
                  ! switch to KODE 12 (tau_s', tau_d', 0): no tensile opening
                  KODE(1,1,L) = 12
                  XMATRIX(IBC+3,NUMEP1) = 0.  ! Dn

                  ! define the amplitude of the stress drop: NORM
                  IF (i_fcode(L) .EQ. 1) THEN
                      NORM = DSE * i_sdrop(L)
                  ELSEIF (i_fcode(L) .EQ. 2) THEN
                      NORM = TAU * i_sdrop(L)
                  ELSEIF (i_fcode(L) .EQ. 3) THEN
                      NORM = DSE+TAUC*i_sdrop(L)
                  ELSEIF (i_fcode(L) .EQ. 4) THEN
                      NORM = DSE+i_sdrop(L)
                  ELSEIF (i_fcode(L) .EQ. 5) THEN
                      NORM = i_sdrop(L)
                  ENDIF
                  NORM = MIN(NORM, TAU) ! avoid dropping more than TAU

                  ! apply the stress drop
                  XMATRIX(IBC+1,NUMEP1) = (-1.)*TAUS/TAU*NORM ! tau_s
                  XMATRIX(IBC+2,NUMEP1) = (-1.)*TAUD/TAU*NORM ! tau_d
                  XMATRIX(IBC+3,NUMEP1) = 0.0 ! sigma_n cancelled

                  ! update SDRIVER (if slides at this iteration, add the driving stress to the one of the previous iterations)
                  SDRIVER(L,1) = SDRIVER(L,1) + XMATRIX(IBC+1,NUMEP1)
                  SDRIVER(L,2) = SDRIVER(L,2) + XMATRIX(IBC+2,NUMEP1)
                  SDRIVER(L,3) = SDRIVER(L,3) + 0
                  ! stored stress
                  SSTORED(L,1) = TAUS + XMATRIX(IBC+1,NUMEP1) ! tau_s_final
                  SSTORED(L,2) = TAUD + XMATRIX(IBC+2,NUMEP1) ! tau_d_final
                  SSTORED(L,3) = SIGN + 0                     ! sig_n_final
              ENDIF
          ENDIF
      ENDDO

      END


      SUBROUTINE INVERT(N)
      ! ******************************************************************************
      ! Invert NxN matrix.  Routine overwrites boundary condition vector.
      ! A = influence coefficient matrix (XMATRIX)
      ! A(1,N+1) = boundary condition vector (XMATRIX(1,N+1))
      ! X = temporary solution vector
      ! Rewritten in Fortran 90+ by A.JANIN 24.02.2026, +add debug option
      ! ******************************************************************************
	  use global_arrays
	  IMPLICIT NONE

	  real(4), allocatable :: X(:)
	  integer(4) :: NB, N
  	  integer(4) :: I, J, L, JJ
	  real(4) :: XM, SUM

	  allocate(X(3*ndis_glob))

      if (debug) then
          ! -> Export the coeff matrix XMATRIX
          !    before the inversion (with stress b.c.s)
          call write_xmatrix('imatrix')
      endif
      
      NB = N - 1

      ! Forward elimination (Triangularization)
	  ! This transforms the matrix into an upper triangular matrix.
	  ! standard Gaussian elimination without pivoting
      do J = 1, NB
          L = J + 1
          do JJ = L, N
              XM = XMATRIX(JJ, J) / XMATRIX(J, J)

              do I = J + 1, N
                  XMATRIX(JJ, I) = XMATRIX(JJ, I) - XMATRIX(J, I) * XM
              end do

              XMATRIX(JJ, N+1) = XMATRIX(JJ, N+1) - XMATRIX(J, N+1) * XM
          end do
      end do

      ! Back substitution
	  ! After elimination, the matrix is upper triangular:
      X(N) = XMATRIX(N, N+1) / XMATRIX(N, N)

      do J = 1, NB
          JJ = N - J
          L = JJ + 1
          SUM = 0.0

          do I = L, N
              SUM = SUM + XMATRIX(JJ, I) * X(I)
          end do

          X(JJ) = (XMATRIX(JJ, N+1) - SUM) / XMATRIX(JJ, JJ)
      end do

      ! Copy solution back into last column
      do I = 1, N
          XMATRIX(I, N+1) = X(I)
      end do

      if (debug) then
          ! -> Export the coeff matrix XMATRIX
          !    after the inversion (with resolved displ b.c.s)
          call write_xmatrix('omatrix')
      endif

      ! deallocate allocated memory to avoid memory leak in ipython terminal
	  deallocate(X)

      RETURN
      END
	  


      SUBROUTINE SHRINK(NUME,MAXB1,MAXB2,NPLANE,KODE,INTRUCTION)
      !******************************************************************************
      ! Subroutine to 
      !    (1) update XMATRIX coefficients if some dislocations are a source of motion
      !        (i.e. when kode>=10). Will take the condition on the driving dislocation
      !        and transform it on stress condition on all the other dislocations!
      !        and store in updating the last row of XMATRIX
      !    (2) remove columns and rows from influence coefficient matrix,
      !        XMATRIX(NUME,NUME), that correspond to element/components with fixed
      !        displacement discontinuities
      ! 
      !   NUME is the total number of element/components in XMATRIX - this will
      !        be changed to reflect the reduced size of the matrix such that
      !        the returned XMATRIX(NUME,NUME) contains no zero rows.
      !
      ! Rewritten in Fortran ~90+ by A.JANIN 25.02.2026
      !******************************************************************************

      USE global_arrays
      use global_inputs
      INTEGER(4) KODE, SIDE
	  INTEGER(4) NUME, NUMEP1, L, IBC
      INTEGER(4) INTRUCTION
      DIMENSION KODE(MAXB1,MAXB2,NPLANE)

      WRITE(*,'(/,A)') ' SHRINKING INFLUENCE COEF. MATRIX'

      ! index of column containing boundary conditions
      NUMEP1 = NUME + 1

      IF (INTRUCTION.EQ.0 .OR. INTRUCTION.EQ.1) THEN

      ! *** Correct boundary conditions for fixed relative displacements;
      ! This is equivalent to subtracting the column of XMATRIX corresponding
      ! to the element/component with fixed relative displacement
      ! from the boundary condition vector.
      !
      ! Recall that the row corresponding to the element/component
      ! with fixed relative displacement boundary conditions will
      ! contain all zeros.
      ! Thus even though the fixed elements are included in the loops,
      ! 'correcting' them does not effect the value of that fixed relative
      ! displacement boundary condition.

      IBC = 0

      ! *** loop over all planes ***
      DO NP = 1, NPLANE

            ! *** loop over IP,JP sub-elements of plane NP ***
            DO JP = 1, NBX1(NP)
                  DO IP = 1, NBX2(NP)

                        IF (KODE(JP,IP,NP) .GT. 7) THEN
                        ! condition changed by A.JANIN from .GT.6 to .GT.7 because new b.c. added

                              IF (KODE(JP,IP,NP) .EQ. 10) THEN
                                    ! correct all other boundary conditions for all three fixed displacements
                                    DO J = 1, 3
                                          ! index of element/component with fixed relative displacement
                                          IBCJ = J + IBC
                                          ! correct all boundary conditions
                                          DO K = 1, NUME
                                              IF (i_fcode(NP).GT.0) THEN
                                                SIDE = -1
                                              ELSE
                                                SIDE = 1
                                              ENDIF
        XMATRIX(K,NUMEP1)=XMATRIX(K,NUMEP1)-
     &           SIDE*XMATRIX(K,IBCJ)*XMATRIX(IBCJ,NUMEP1)
                                          ENDDO
                                    ENDDO

                              ELSEIF (KODE(JP,IP,NP) .EQ. 11 .OR. 
     &                                KODE(JP,IP,NP) .EQ. 12) THEN
                                    ! correct all other boundary conditions for fixed normal displacement only
                                    IBC3 = IBC + 3
                                    DO K = 1, NUME
        XMATRIX(K,NUMEP1)=XMATRIX(K,NUMEP1)-
     &             XMATRIX(K,IBC3)*XMATRIX(IBC3,NUMEP1)
                                    ENDDO

                              ELSEIF (KODE(JP,IP,NP) .EQ. 13) THEN
                                    ! correct all other boundary conditions for fixed dip and normal displacement
                                    IBC2 = IBC + 2
                                    IBC3 = IBC + 3
                                    DO K = 1, NUME
        XMATRIX(K,NUMEP1)=XMATRIX(K,NUMEP1)-
     &                        XMATRIX(K,IBC2)*XMATRIX(IBC2,NUMEP1)
        XMATRIX(K,NUMEP1)=XMATRIX(K,NUMEP1)-
     &                        XMATRIX(K,IBC3)*XMATRIX(IBC3,NUMEP1)
                                    ENDDO

                              ELSEIF (KODE(JP,IP,NP) .EQ. 14) THEN
                                    ! correct all other boundary conditions for fixed strike and normal displacement
                                    IBC1 = IBC + 1
                                    IBC3 = IBC + 3
                                    DO K = 1, NUME
        XMATRIX(K,NUMEP1)=XMATRIX(K,NUMEP1)-
     &                XMATRIX(K,IBC1)*XMATRIX(IBC1,NUMEP1)
        XMATRIX(K,NUMEP1)=XMATRIX(K,NUMEP1)-
     &                XMATRIX(K,IBC3)*XMATRIX(IBC3,NUMEP1)
                                       ENDDO

                              ELSEIF (KODE(JP,IP,NP) .EQ. 15) THEN
                                    ! Added by A.JANIN correct all other boundary conditions for fixed strike and dip displacement
                                    IBC1 = IBC + 1
                                    IBC2 = IBC + 2
                                    DO K = 1, NUME
        XMATRIX(K,NUMEP1)=XMATRIX(K,NUMEP1)-
     &                XMATRIX(K,IBC1)*XMATRIX(IBC1,NUMEP1)
        XMATRIX(K,NUMEP1)=XMATRIX(K,NUMEP1)-
     &                XMATRIX(K,IBC2)*XMATRIX(IBC2,NUMEP1)
                                    ENDDO
                              ENDIF

                        ENDIF

                        ! index of next sub-element
                        IBC = IBC + 3

                  ENDDO
            ENDDO
      ENDDO

      ! Update stress condition on each dislocation (in local coordinate system of the plane)
      ! A.JANIN 25.02.2026
      DO L=1, NPLANE
          IBC = 3*(L-1)
          IF (KODE(1,1,L) .LT. 7) THEN
          ! no fixed relative displacements
              SDRIVER(L,1) = XMATRIX(IBC+1,NUMEP1) ! tau_strike(L)
              SDRIVER(L,2) = XMATRIX(IBC+2,NUMEP1) ! tau_dip(L)
              SDRIVER(L,3) = XMATRIX(IBC+3,NUMEP1) ! sigma_n(L)
          ELSEIF (KODE(1,1,L) .EQ. 10) THEN
          ! fixed relative displacements imposed (assumed no stress)
              SDRIVER(L,1) = 0.
              SDRIVER(L,2) = 0.
              SDRIVER(L,3) = 0.
          ELSEIF (KODE(1,1,L) .EQ. 11 .OR. KODE(1,1,L) .EQ. 12) THEN
              SDRIVER(L,1) = XMATRIX(IBC+1,NUMEP1)
              SDRIVER(L,2) = XMATRIX(IBC+2,NUMEP1)
              SDRIVER(L,3) = 0.
          ELSEIF (KODE(1,1,L) .EQ. 13) THEN
              SDRIVER(L,1) = XMATRIX(IBC+1,NUMEP1)
              SDRIVER(L,2) = 0.
              SDRIVER(L,3) = 0.
          ELSEIF (KODE(1,1,L) .EQ. 14) THEN
              SDRIVER(L,1) = 0.
              SDRIVER(L,2) = XMATRIX(IBC+2,NUMEP1)
              SDRIVER(L,3) = 0.
          ELSEIF (KODE(1,1,L) .EQ. 15) THEN
              SDRIVER(L,1) = 0.
              SDRIVER(L,2) = 0.
              SDRIVER(L,3) = XMATRIX(IBC+3,NUMEP1)
		  ENDIF
      ENDDO


      ! CALL FOR FRICTION
      CALL FRICTIONALDISLOC(NUME, NPLANE, KODE)

      ! --- If all boundary conditions are fixed rel. displ. return
      IF (MINVAL(KODE) == 10 .AND. MAXVAL(KODE) == 10) THEN
          RETURN
      ENDIF

      ENDIF

      ! *** Remove rows and columns of XMATRIX corresponding to elements/components
      !     with fixed relative displacement

      IF (INTRUCTION.EQ.0 .OR. INTRUCTION.EQ.2) THEN
      IBC = 1

      ! *** loop over all planes ***
      DO NP = 1, NPLANE

            ! *** loop over IP,JP blocks of plane NP ***
            DO JP = 1, NBX1(NP)
                  DO IP = 1, NBX2(NP)

                        IF (KODE(JP,IP,NP) .LE. 7) THEN
                              ! nothing is fixed
                              ! index of next element/comp to work on
                              IBC = IBC + 3

                        ELSEIF (KODE(JP,IP,NP) .EQ. 10) THEN
                              ! fixed strike,dip,normal displ.; row/column index is for 1st comp.
                              ! this element - remove row & column IBC 3 times
                              CALL ROWCOL(IBC,NUME,NUMEP1)
                              CALL ROWCOL(IBC,NUME,NUMEP1)
                              CALL ROWCOL(IBC,NUME,NUMEP1)

                        ELSEIF (KODE(JP,IP,NP) .EQ. 11 .OR.
     &                          KODE(JP,IP,NP) .EQ. 12) THEN
                              ! fixed normal displacement only; row/column index is for 3rd comp.
                              IBC = IBC + 2
                              CALL ROWCOL(IBC,NUME,NUMEP1)

                        ELSEIF (KODE(JP,IP,NP) .EQ. 13) THEN
                              ! fixed dip and normal displacement
                              ! row/column index is for 2nd comp.
                              IBC = IBC + 1
                              CALL ROWCOL(IBC,NUME,NUMEP1)
                              CALL ROWCOL(IBC,NUME,NUMEP1)

                        ELSEIF (KODE(JP,IP,NP) .EQ. 14) THEN
                              ! fixed strike and normal displacement
                              CALL ROWCOL(IBC,NUME,NUMEP1)
                              IBC = IBC + 1
                              CALL ROWCOL(IBC,NUME,NUMEP1)

                        ELSEIF (KODE(JP,IP,NP) .EQ. 15) THEN
                              ! Added by A.JANIN: fixed strike and dip displacement
                              CALL ROWCOL(IBC,NUME,NUMEP1)
                              CALL ROWCOL(IBC,NUME,NUMEP1)
                              IBC = IBC + 1
                        ENDIF

                  ENDDO
            ENDDO
      ENDDO

      ENDIF

      RETURN
      END


      SUBROUTINE DCD3D(X,Y,Z,DEPTH,DIP,                        
     &              AL1,AL2,AW1,AW2,IGRAD,IRET) 
C******************************************************************************
C Subroutine to calculate influence coefficients of a plane with reference depth
C DEPTH, planar dimensions (AL2-AL1)x(AW2-AW1), dip=DIP, at a point X,Y,Z in
C the local coordinate system of the plane.
C******************************************************************************
      IMPLICIT REAL*8 (A-H,O-W)

C ** Added next line Aug. 29, 2000
      REAL*4 DEPTH,DIP
	 
      REAL*8 X,Y,Z
 
      REAL*4 ALPHA,XMU,XMU2,DMULT
      COMMON/CONSTANTS/ALPHA,XMU,XMU2,DMULT
      
      REAL*4 DSPL,STR
      COMMON/TEMPS/DSPL(3,3),STR(6,3)
      
      REAL*8 UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ
      DIMENSION UXX(3),UYX(3),UZX(3),UXY(3),
     &          UYY(3),UZY(3),UXZ(3),UYZ(3),UZZ(3)

      REAL*4 DGRAD
      COMMON/DGRADS/DGRAD(3,3,3)

      REAL*4 STRESS,EWE,dgrten,V,E,FRICTION,tmax,tmaxo,tmaxa
      COMMON/MIKE/STRESS(6),EWE(3),dgrten(3,3),V,E,FRICTION,
     &    tmax,tmaxo(2),tmaxa
	  
C     flag for displacement gradient calculation
      INTEGER*4 IGRAD
	   
C     - Calculate the influence coefficients for each component
C     of displacement (1=strike,2=dip,3=normal or tensile)
      CALL DCD3(ALPHA,X,Y,Z,DEPTH,DIP,AL1,AL2,AW1,AW2,
     &  DSPL(1,1),DSPL(1,2),DSPL(1,3),
     &  DSPL(2,1),DSPL(2,2),DSPL(2,3),
     &  DSPL(3,1),DSPL(3,2),DSPL(3,3),
     &  UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
      IF(IRET.NE.0) RETURN

      XMU5=XMU*1.E-5
      XMU25=XMU2*1.E-5
	  
C   - Store displacement gradients and 
C     convert from displacement gradients to stresses      
      DO 1 J=1,3

      DILAT=(UXX(J)+UYY(J)+UZZ(J))*DMULT
      
      DUXX=DILAT+UXX(J)
      IF(ABS(DUXX).LT.1.E-6) THEN
  	STR(1,J)=0. 
      ELSE
      	STR(1,J)=XMU25*DUXX
      ENDIF
      DUYY=DILAT+UYY(J)
      IF(ABS(DUYY).LT.1.E-6) THEN
  	 STR(4,J)=0.
      ELSE
         STR(4,J)=XMU25*DUYY
      ENDIF
      DUZZ=DILAT+UZZ(J)
      IF(ABS(DUZZ).LT.1.E-6) THEN
   	 STR(6,J)=0. 
      ELSE
         STR(6,J)=XMU25*DUZZ
      ENDIF
      DUXY=UXY(J)+UYX(J)
      IF(ABS(DUXY).LT.1.E-6) THEN
   	 STR(2,J)=0.
      ELSE
     	 STR(2,J)=XMU5*DUXY
      ENDIF
      DUXZ=UXZ(J)+UZX(J)
      IF(ABS(DUXZ).LT.1.E-6) THEN
  	 STR(3,J)=0.
      ELSE
         STR(3,J)=XMU5*DUXZ
      ENDIF
      DUYZ=UZY(J)+UYZ(J)
      IF(ABS(DUYZ).LT.1.E-6) THEN
  	 STR(5,J)=0.
      ELSE
      	 STR(5,J)=XMU5*DUYZ
      ENDIF
1     CONTINUE

      IF(IGRAD.LE.0) RETURN

      DO 2 J=1,3
 	  DGRAD(1,1,J)=UXX(J)*1.E-5
	  DGRAD(1,2,J)=UXY(J)*1.E-5
	  DGRAD(1,3,J)=UXZ(J)*1.E-5

	  DGRAD(2,1,J)=UYX(J)*1.E-5
	  DGRAD(2,2,J)=UYY(J)*1.E-5
	  DGRAD(2,3,J)=UYZ(J)*1.E-5

	  DGRAD(3,1,J)=UZX(J)*1.E-5
	  DGRAD(3,2,J)=UZY(J)*1.E-5
2	  DGRAD(3,3,J)=UZZ(J)*1.E-5
 
      RETURN
      END

       SUBROUTINE COEF2NP(NS,NP)
C******************************************************************************
C Transform coef matrices in STR,DSPL from local NS coords to planar NP coords 
C******************************************************************************
      use global_arrays
      
      REAL*4 UCOEF,SCOEF
      COMMON/COEFS/ UCOEF(3,3),SCOEF(6,3)

      REAL*4 DSPL,STR
      COMMON/TEMPS/DSPL(3,3),STR(6,3)

C---- Transform STR,DSPL to global coords, store in SCOEF(6,3),UCOEF(3,3)      
      IGRAD=0
      CALL COEF2GLOB(NS,IGRAD)

C---- Transform normal components of SCOEF(6,3),and UCOEF(3,3) to planar NP
C     coordinates, store back in STR(3,3), DSPL(3,3)
      
      DO 10 I=1,3
       DO 10 J=1,3
       
        DSPL(I,J)=0.
	STR(I,J)=0.

	DO 2 K=1,6       
2	STR(I,J)=SG2P(I,K,NP)*SCOEF(K,J)+STR(I,J)

        DO 3 K=1,3
3       DSPL(I,J)=UG2P(I,K,NP)*UCOEF(K,J)+DSPL(I,J)

10    CONTINUE	
      RETURN
      END


       SUBROUTINE COEF2GLOB(NS,IGRAD)
C******************************************************************************
C transform coef matrices from local NS coords to global coords 
C******************************************************************************
      use global_arrays

      REAL*4 UCOEF,SCOEF
      COMMON/COEFS/ UCOEF(3,3),SCOEF(6,3)

      REAL*4 DSPL,STR
      COMMON/TEMPS/DSPL(3,3),STR(6,3)
      
      REAL*4 CSQ,SSQ,SC,SC2,C2MS2

      REAL*4 DGRAD
      COMMON/DGRADS/DGRAD(3,3,3)
	  
C     flag for displacement gradients
      INTEGER*4 IGRAD
  
C     cosine(strike)**2
      CSQ=C(NS)*C(NS)
C     sine(strike)**2	
      SSQ=S(NS)*S(NS)
C     cosine(strike)*sine(strike)
      SC=S(NS)*C(NS)
C     2*cosine(strike)*sine(strike)	
      SC2=2.*SC	
C     cosine(strike)**2-sine(strike)**2	
      C2MS2=CSQ-SSQ

C    - Rotate stress tensor and displacements to global coordinates      
      DO 1 J=1,3
       SCOEF(1,J)=SSQ*STR(1,J)  -SC2*STR(2,J)+CSQ*STR(4,J)
       SCOEF(2,J)= SC*(STR(1,J)-STR(4,J)) -C2MS2*STR(2,J) 
       SCOEF(3,J)=S(NS)*STR(3,J)-C(NS)*STR(5,J)
       SCOEF(4,J)=CSQ*STR(1,J)  +SC2*STR(2,J)+SSQ*STR(4,J)
       SCOEF(5,J)=C(NS)*STR(3,J)+S(NS)*STR(5,J)
       SCOEF(6,J)=STR(6,J)

       UCOEF(1,J)=S(NS)*DSPL(1,J)-C(NS)*DSPL(2,J)
       UCOEF(2,J)=C(NS)*DSPL(1,J)+S(NS)*DSPL(2,J)
1      UCOEF(3,J)=DSPL(3,J)

       IF(IGRAD.LE.0) RETURN
	   
C    - Rotate displacement gradien tensor to global coordinates
       DO 2 J=1,3
       tmp12=SC*(DGRAD(1,2,J)+DGRAD(2,1,J))
       dxx=SSQ*DGRAD(1,1,J)  -tmp12 +CSQ*DGRAD(2,2,J)
       tmpxy=SC*(DGRAD(1,1,J)-DGRAD(2,2,J))
       dxy= tmpxy +SSQ*DGRAD(1,2,J) -CSQ*DGRAD(2,1,J)
       dyx= tmpxy +SSQ*DGRAD(2,1,J) -CSQ*DGRAD(1,2,J)
       dxz=S(NS)*DGRAD(1,3,J)-C(NS)*DGRAD(2,3,J)
       dzx=S(NS)*DGRAD(3,1,J)-C(NS)*DGRAD(3,2,J)
       dyy=CSQ*DGRAD(1,1,J)  +tmp12 +SSQ*DGRAD(2,2,J)
       dyz=C(NS)*DGRAD(1,3,J)+S(NS)*DGRAD(2,3,J)
       dzy=C(NS)*DGRAD(3,1,J)+S(NS)*DGRAD(3,2,J)
       dzz=DGRAD(3,3,J)

       DGRAD(1,1,J)=dxx
       DGRAD(1,2,J)=dxy
       DGRAD(1,3,J)=dxz

       DGRAD(2,1,J)=dyx
       DGRAD(2,2,J)=dyy
       DGRAD(2,3,J)=dyz

       DGRAD(3,1,J)=dzx
       DGRAD(3,2,J)=dzy
2      DGRAD(3,3,J)=dzz

       RETURN
       END

       SUBROUTINE MK_MATRIX(KODNP,NUMED,NUMEBC)
C******************************************************************************
C  Subroutine to make up the influence coefficient matrix for the given
C  boundary conditions.
C******************************************************************************
      use global_arrays
	  
      INTEGER*4 KODNP,NUMED,NUMEBC
	  	  
      REAL*4 DSPL,STR
      COMMON/TEMPS/DSPL(3,3),STR(6,3)

      IF(KODNP.EQ.1) THEN
C	3 displacement b.c.s
	DO 1 J=1,3
	DO 1 I=1,3
1	XMATRIX(NUMEBC+I,NUMED+J)=DSPL(I,J)
      ELSEIF(KODNP.EQ.2) THEN
C	3 stress b.c.s
	DO 2 J=1,3
	DO 2 I=1,3
2	XMATRIX(NUMEBC+I,NUMED+J)=STR(I,J)
      ELSEIF(KODNP.EQ.3) THEN
C	strike-shear & normal stress, dip displ. b.c.s
	DO 3 I=1,3
	NUMEDI=NUMED+I
	XMATRIX(NUMEBC+1,NUMEDI)=STR(1,I)
	XMATRIX(NUMEBC+2,NUMEDI)=DSPL(2,I)
3	XMATRIX(NUMEBC+3,NUMEDI)=STR(3,I)
      ELSEIF(KODNP.EQ.4) THEN
C 	dip-shear & normal stress, strike displ. b.c.s
	DO 4 I=1,3
	NUMEDI=NUMED+I
	XMATRIX(NUMEBC+1,NUMEDI)=DSPL(1,I)
	XMATRIX(NUMEBC+2,NUMEDI)=STR(2,I)
4	XMATRIX(NUMEBC+3,NUMEDI)=STR(3,I)
      ELSEIF(KODNP.EQ.5) THEN
C 	strike & dip-shear stresses, normal displ. b.c.s
	DO 9 I=1,3
	NUMEDI=NUMED+I
	XMATRIX(NUMEBC+1,NUMEDI)=STR(1,I)
	XMATRIX(NUMEBC+2,NUMEDI)=STR(2,I)
9	XMATRIX(NUMEBC+3,NUMEDI)=DSPL(3,I)
      ELSEIF(KODNP.EQ.6) THEN
C 	strike & dip-shear displ, normal stress. b.c.s, added by A.JANIN
	DO I=1,3
	NUMEDI=NUMED+I
	XMATRIX(NUMEBC+1,NUMEDI)=DSPL(1,I)
	XMATRIX(NUMEBC+2,NUMEDI)=DSPL(2,I)
	XMATRIX(NUMEBC+3,NUMEDI)=STR(3,I)
	ENDDO
      ELSEIF(KODNP.EQ.10) THEN
C	strike, dip, normal displ. fixed b.c.s
	DO 7 J=1,3
        DO 7 I=1,3
7	XMATRIX(NUMEBC+I,NUMED+J)=0.
      ELSEIF(KODNP.EQ.11.OR.KODNP.EQ.12) THEN
C     strike & dip stress, normal displ. fixed b.c.s
	DO 5 I=1,3
	NUMEDI=NUMED+I
	XMATRIX(NUMEBC+1,NUMEDI)=STR(1,I)
 	XMATRIX(NUMEBC+2,NUMEDI)=STR(2,I)
5	XMATRIX(NUMEBC+3,NUMEDI)=0.
      ELSEIF(KODNP.EQ.13) THEN
	NUMEB1=NUMEBC+1
	XMATRIX(NUMEB1,NUMED+1)=STR(1,1)	
C	strike stress, normal & dip displ. fixed b.c.s
	XMATRIX(NUMEB1,NUMED+2)=STR(1,2)
	XMATRIX(NUMEB1,NUMED+3)=STR(1,3)
	DO 6 I=2,3
	NUMEBI=NUMEBC+I
	DO 6 J=1,3
6	XMATRIX(NUMEBI,NUMED+J)=0.
      ELSEIF(KODNP.EQ.14) THEN
	NUMEB2=NUMEBC+2
	XMATRIX(NUMEB2,NUMED+1)=STR(2,1)
C	dip stress, normal & strike displ. fixed b.c.s
	XMATRIX(NUMEB2,NUMED+2)=STR(2,2)
	XMATRIX(NUMEB2,NUMED+3)=STR(2,3)
	K=1
	DO 8 I=1,2
	NUMEBI=NUMEBC+K
	DO 8 J=1,3
 	XMATRIX(NUMEBI,NUMED+J)=0.
8	K=K+2
      ELSEIF(KODNP.EQ.15) THEN
C	Added by A.JANIN: strike and dip displ., normal stress. fixed b.c.s
	NUMEB3=NUMEBC+3
	XMATRIX(NUMEB3,NUMED+1)=STR(3,1)
	XMATRIX(NUMEB3,NUMED+2)=STR(3,2)
	XMATRIX(NUMEB3,NUMED+3)=STR(3,3)
	DO 10 I=1,2
	NUMEBI=NUMEBC+I
	DO 10 J=1,3
 	XMATRIX(NUMEBI,NUMED+J)=0.
10	K=K+2
	ENDIF
      RETURN
      END


	  SUBROUTINE ROWCOL(ICOL,NUME,NUMEP1)
C******************************************************************************
C  Subroutine to remove rows, columns from influence coef. matrix corresponding
C  to elements with fixed displacement discontinuities.  Corresponding 
C  b.c. (containing the fixed relative displacement discontinuities) are removed
C  from the b.c. vector.
C******************************************************************************
      use global_arrays

      INTEGER*4 ICOL,NUME,NUMEP1

C *** Remove a column
C     loop over columns
      DO 1 I=ICOL,NUMEP1   
          II=I+1    
C	  shift column I+1 to column I			     
	  DO 1 J=1,NUME
1	  XMATRIX(J,I)=XMATRIX(J,II)

C *** Remove a row
C     loop over rows
      DO 4 J=ICOL,NUME
	  JJ=J+1
C         shift row J+1 to row J;using NUME instead of  
C         NUMEP1 accounts for the column just removed
	  DO 4 I=1,NUME
4	  XMATRIX(J,I)=XMATRIX(JJ,I)

C *** Reduce the array dimension by 1
      NUME=NUME-1
      NUMEP1=NUMEP1-1
	  
      RETURN
      END


      SUBROUTINE CALCPO(X,Y,Z,NPLANE,NUM_Ds,IRET)
      ! **************************************************************************
      ! Subroutine to calculate deformation at a point X,Y,Z in global coordinates
      ! Rewritten in Fortran 90+ by A.JANIN 27.02.2026
      ! **************************************************************************
      
      use global_arrays
	  
      INTEGER*4     NPLANE
      REAL*8        X,Y

      REAL*4 UCOEF,SCOEF
      COMMON/COEFS/ UCOEF(3,3),SCOEF(6,3)
      
      REAL*8 XNP,YNP,Z,AL1,AL2,AW1,AW2

      REAL*4 DGRAD
      COMMON/DGRADS/DGRAD(3,3,3)

      REAL*4 STRESS,EWE,dgrten,V,E,FRICTION,tmax,tmaxo,tmaxa
      COMMON/MIKE/STRESS(6),EWE(3),dgrten(3,3),V,E,FRICTION,
     &   tmax,tmaxo(2),tmaxa

      CHARACTER*4 BFLAG
      REAL*4 BSTRESS
      COMMON/BKGRND/BSTRESS(9),BFLAG

      ! Flag for displacement gradient calculation
      INTEGER*4 IGRAD

      IGRAD=1

      ! Solution is the sum of deformation due to all elements; this is stored
      ! in arrays DSPL(3,1) (displacements) and STR(6,1) (stresses)

      ! Initialize solution vectors
      DO I=1,3
         DO J=1,3
            dgrten(I,J)=0.
         ENDDO
         EWE(I)=0.
      ENDDO

      DO I=1,6
         STRESS(I)=0.
      ENDDO

      ! For each point X,Y,Z

      NUMD3=0

      ! *** loop over NS planes ***
      DO NS=1,NPLANE

         ! Define point in local NS coords
         XNP=S(NS)*(X-XO(NS))+C(NS)*(Y-YO(NS))
         YNP=-C(NS)*(X-XO(NS))+S(NS)*(Y-YO(NS))

         ! *** loop over IS,JS sub-elements of plane NS ***
         DO JS=1,NBX1(NS)

            ! sub-element length range along dip
            AL1=FLOAT(JS-1)*BWX1(NS)
            AL2=AL1+BWX1(NS)

            DO IS=1,NBX2(NS)

               ! sub-element width range along strike
               AW1=FLOAT(IS-1)*BWX2(NS)
               AW2=AW1+BWX2(NS)

               ! Calculate 6 stress and 3 displacement influence coefs in NS coords
               ! Put these in arrays STR(6,3),DSPL(3,3)

               CALL DCD3D(XNP,YNP,Z,ZO(NS),DIP(NS),
     &                    AL1,AL2,AW1,AW2,IGRAD,IRET)

               IF(IRET.NE.0) THEN
                  WRITE(*,*)
     &            ' Singular Point in Calcpo; <x,y,z>=',XNP,YNP,Z

                  DO I=1,6
                     STRESS(I)=999999.
                  ENDDO

                  DO I=1,3
                     EWE(I)=999999.
                  ENDDO

                  DO I=1,3
                     DO J=1,3
                        dgrten(J,I)=999999.
                     ENDDO
                  ENDDO

                  RETURN
               ENDIF

               ! Transform STR,DSPL,DGRAD to global coords, store in
               ! SCOEF(6,3),UCOEF(3,3),DGRAD
               CALL COEF2GLOB(NS,IGRAD)

               ! Weight influence coefs by solved displacement
               ! discontinuities (Ds) and add to solution
               DO J=1,3

                  DO I=1,6
                     STRESS(I)=STRESS(I)
     &                    +SCOEF(I,J)*XMATRIX(NUMD3+J,NUM_Ds)
                  ENDDO

                  DO I=1,3
                     EWE(I)=EWE(I)
     &                  +UCOEF(I,J)*XMATRIX(NUMD3+J,NUM_Ds)
                  ENDDO

                  DO I=1,3
                     DO K=1,3
                        dgrten(K,I)=dgrten(K,I)
     &                       +DGRAD(K,I,J)
     &                       *XMATRIX(NUMD3+J,NUM_Ds)
                     ENDDO
                  ENDDO

               ENDDO

               ! increment displacement discontinuity count
               NUMD3=NUMD3+3

            ENDDO
         ENDDO
         ! *** end of loop over blocks of plane NS ***

      ENDDO
      ! *** end of loop over NS planes ***

      ! *** add in effects of background deformation
      IF(BFLAG.NE.'N') THEN

         DO I=1,6
            STRESS(I)=STRESS(I)+BSTRESS(I)
         ENDDO

         sum=(1.0+V)/E

         dgrten(1,2)=dgrten(1,2)
     &        +(BSTRESS(1)-V*(BSTRESS(4)+BSTRESS(6)))/E
         dgrten(1,2)=dgrten(1,2)
     &        +sum*BSTRESS(2)-BSTRESS(9)

         dgrten(1,3)=dgrten(1,3)
     &        +sum*BSTRESS(3)+BSTRESS(8)

         dgrten(2,1)=dgrten(2,1)
     &        +sum*BSTRESS(2)+BSTRESS(9)

         dgrten(2,2)=dgrten(2,2)
     &        +(BSTRESS(4)-V*(BSTRESS(1)+BSTRESS(6)))/E

         dgrten(2,3)=dgrten(2,3)
     &        +sum*BSTRESS(5)-BSTRESS(7)

         dgrten(3,1)=dgrten(3,1)
     &        +sum*BSTRESS(3)-BSTRESS(8)

         dgrten(3,2)=dgrten(3,2)
     &        +sum*BSTRESS(5)+BSTRESS(7)

         dgrten(3,3)=dgrten(3,3)
     &        +(BSTRESS(6)-V*(BSTRESS(1)+BSTRESS(4)))/E

         ! *** account for rigid body rotations
         IF(BFLAG.EQ.'D') THEN

            UX=EWE(3)*BSTRESS(8)-EWE(2)*BSTRESS(9)
            UY=EWE(1)*BSTRESS(9)-EWE(3)*BSTRESS(7)

            EWE(3)=EWE(3)
     &           +(EWE(2)*BSTRESS(7)-EWE(1)*BSTRESS(8))

            EWE(1)=UX+EWE(1)
            EWE(2)=UY+EWE(2)

         ENDIF

      ENDIF

      RETURN
      END

   	SUBROUTINE TOPLANE(SSTK,CSTK,SDIP,CDIP)
C******************************************************************************
C--- Create transformation matrices (from global to in-plane coordinates)
C	 for plane NS.
C******************************************************************************	  
	  
	REAL*4 ROT,U,TSTR,TE,S
        DIMENSION ROT(3,3),U(3),TSTR(3,3),TE(3,3),S(3,3)

	REAL*4 STRESS,EWE,dgrten,pr,ym,FRICTION,tmax,tmaxo,tmaxa
	COMMON/MIKE/STRESS(6),EWE(3),dgrten(3,3),
     &              pr,ym,FRICTION,tmax,tmaxo(2),tmaxa

C ** Added dum to REAL*4 declaration
        REAL*4 E,dum
        COMMON/TEMPS/E(3,3),dum(6,3)
		  
       S(1,1)=STRESS(1)  
       S(2,1)=STRESS(2)  
       S(3,1)=STRESS(3)  
       S(2,2)=STRESS(4)  
       S(2,3)=STRESS(5)  
       S(3,3)=STRESS(6)  
       S(1,2)=S(2,1)  
       S(3,2)=S(2,3)  
       S(1,3)=S(3,1)  
C				in-plane component	global component
C				strike direction	strike direction 
      ROT(1,1)=SSTK
C				strike direction	dip direction
      ROT(1,2)=CSTK
C				strike direction	normal direction
      ROT(1,3)=0.
C				dip direction		strike direction
      ROT(2,1)=-CSTK*CDIP
C				dip direction		dip direction
      ROT(2,2)=SSTK*CDIP
C				dip direction		normal direction
      ROT(2,3)=SDIP
C				normal direction	strike direction
      ROT(3,1)=CSTK*SDIP
C				normal direction	dip direction	
      ROT(3,2)=-SSTK*SDIP
C				normal direction	normal direction
      ROT(3,3)=CDIP

      DO 1 K=1,3
        DO 1 I=1,3
	TSTR(I,K)=0.
	TE(I,K)=0.
	DO 1 J=1,3
	TE(I,K)=TE(I,K)+E(I,J)*ROT(K,J)
1       TSTR(I,K)=TSTR(I,K)+S(I,J)*ROT(K,J)

      DO 2 K=1,3  
        DO 2 I=1,3  
	S(I,K)=0.
	E(I,K)=0.
	DO 2 J=1,3
	E(I,K)=E(I,K)+TE(J,K)*ROT(I,J)
2       S(I,K)=S(I,K)+TSTR(J,K)*ROT(I,J)

	  
      STRESS(1)=S(1,1)  
      STRESS(2)=S(2,1)  
      STRESS(3)=S(3,1)  
      STRESS(4)=S(2,2)  
      STRESS(5)=S(2,3)  
      STRESS(6)=S(3,3)  
C
C  Code added 2/15/94 to calcultate the max shear stress across the plane
C  of interest. tmax is the magnitude, tmaxo is a vector whose x,y components give
C  the orientation of tmax (this is intended for the MATLAB function QUIVER), and
C  tmaxa is the rake or pitch angle of tmax.
C  See manual for derivation of these results.
C
	tmaxa=atan(stress(5)/stress(3))
	tmaxo(1)=stress(3)
	tmaxo(2)=stress(5)
	tmax=stress(3)*cos(tmaxa)+stress(5)*sin(tmaxa)
C
C
        DO 6 K=1,3
        DO 6 I=1,3
	TSTR(I,K)=0.
	DO 6 J=1,3
6       TSTR(I,K)=TSTR(I,K)+dgrten(I,J)*ROT(K,J)

        DO 7 K=1,3  
        DO 7 I=1,3  
	dgrten(I,K)=0.
	DO 7 J=1,3	
7       dgrten(I,K)=dgrten(I,K)+TSTR(J,K)*ROT(I,J)

	DO 4 K=1,3
	U(K)=0.  
        DO 4 I=1,3
4	U(K)=U(K)+EWE(I)*ROT(K,I)
	    
	DO 5 I=1,3
5	EWE(I)=U(I)
	  
      RETURN
      END


      subroutine invariants(volchg,critic,octshr,work)
*
*********************************************************************
*	Subroutine invariants calculates the stress difference between a
*	an assumed failure stress, f(z), and the stress state following
*	deformation.  It is also necessary to assume some reference
*	stress state in the crust.  This is tricky, as any reference state
*	is likely to be both temporally and spatailly variable in a region 
*	of active faulting.  For simplicity, we will assume that the average
*	state is lithostatic.  Bear in mind, however, that this subroutine is
*	trying to predict where aftershocks happen, and so dealing with an 
*	"average" stress state is not the best thing to do.  
*	
*********************************************************************
*
C      include 'sizes.inc'
	  
C  ** Added declaration of tmax,tmaxo(2),tmaxa as real*4 on Aug. 29,2000
      real*4 stress,ewe,dgrten,pr,ym,FRICTION,tmax,tmaxo,tmaxa
      common/mike/stress(6),ewe(3),dgrten(3,3),pr,ym,FRICTION,
     &    tmax,tmaxo(2),tmaxa
  
      real*4 volchg,critic,octshr,work
      real*4 s(3,3),p(3),v(3,3)
      real*4 dev1,dev3,phi
	  	  
      real*4 e,dum
      common/temps/e(3,3),dum(6,3)

      rho=2650.
c     g=9.78
C     angle of internal friction in radians
C     corresponding to a coefficient of friction=0.6
      phi=0.540
					
c     meanst=(abs(z)*1000.*g*rho)
      meanst=0.
	  
      s(1,1)=stress(1)+meanst
	  s(1,2)=stress(2)
	  s(1,3)=stress(3)
	  s(2,1)=stress(2)
	  s(2,2)=stress(4)+meanst
	  s(2,3)=stress(5)
	  s(3,1)=stress(3)
	  s(3,2)=stress(5)
	  s(3,3)=stress(6)+meanst
	  
	  call jacobi(s,p,v)
	  
	  dev1=max(p(1),p(2),p(3))
	  dev3=min(p(1),p(2),p(3))
	  
	  critic=(dev1+dev3)/2.

 	  volchg=e(1,1)+e(2,2)+e(3,3)
 
C	  Calculate stress invariants I and II*
 	  sinv1=stress(1)*stress(4)*stress(6)
 	  temp=stress(6)*(stress(4)+stress(1))+stress(1)*stress(4)
 	  tmp=stress(5)*stress(5)+stress(3)*stress(3)+stress(2)*stress(2)
 	  sinv2=-1.0*temp+tmp
 
C	  calculate the octahedral shear stress
 	  octshr=0.4714*sqrt(sinv1*sinv1+3.0*(sinv2))
 
C	  calculate the total strain energy (sum of energy due to distortion
C	  and volume change).  See Ramsay, 1967, p.288.
 	  work=sinv1*sinv1-2.0*(sinv1+pr)*sinv2

	  return
	  end
	  
	  Subroutine Third_Results(ex,ey,ez,px,tx,py,ty,pz,tz)
********************************************************************************
*
*	Third_Results finds the principal strains and their orientations with respect
*	to a geographic reference frame.
*
********************************************************************************
*
	  real*4 eval(3)
	  real*4 ex,ey,ez,VX(3),VY(3),VZ(3),TX,PX,TY,PY,TZ,PZ

	  real*4 e,dum
	  common/temps/e(3,3),dum(6,3)

          REAL*4 d,fing,evec
	  COMMON/DGRADS/d(3,3),fing(3,3),evec(3,3)

C         Converts strain tensor to deformation	gradient tensor.  
	  DO 400 I=1,3
              DO 410 J=1,3
C	      Note: this is not the same as the displacement gradient tensor.
	      D(I,J)=E(I,J)
	      IF(I.EQ.J) D(I,J)=E(I,J)+1.0
 410	      CONTINUE
 400	  CONTINUE

C         Converts the deformation gradient tensor to Finger's tensor.
	  CALL FINGER(d,fing)

C         Obtains eigenvalues and eigenvectors from Finger's tensor, 
C         corresponding to squares of the principal strains and their 
C         orientations.
	  CALL JACOBI(FING,EVAL,EVEC)

C	  Gets strains and sorts them by size.
	  CALL SORT(EVAL,EVEC,EX,EY,EZ,VX,VY,VZ)

 	  ex=sqrt(ex)
	  ey=sqrt(ey)
	  ez=sqrt(ez)

C	  Calculates geographic orientation of principal strains.
	  CALL TREND(VX,VY,VZ,TX,PX,TY,PY,TZ,PZ)

	  return
	  end


	  subroutine STRAIN (e,delV)
*************************************************************************
*
* 	STRAIN takes the output from CALCPO and calculates the strain tensor 
*	components and the volume change (where 0 = no change)
*
*   The stress tensor components are stored according to:
*         Vector 	       Tensor	      Component
*	   Index               Index
*  		1		1,1	        Sxx
*		2		1,2 2,1		Sxy
*		3		1,3 3,1		Sxz
*		4		2,2		Syy
*		5		2,3		Syz
*		6		3,3		Szz
*************************************************************************
*
*
	  real*4 sum,e(3,3)
	  real*8, intent(out) :: delV
*
	  REAL*4 S,EWE,dgrten,pr,ym,FRICTION,tmax,tmaxo,tmaxa
	  COMMON/MIKE/S(6),EWE(3),dgrten(3,3),pr,ym,FRICTION,
     &    tmax,tmaxo(2),tmaxa
*
	  sum=(1.0+PR)/YM
*
	  e(1,1)=(s(1)-PR*(s(4)+s(6)))/YM
	  e(2,2)=(s(4)-PR*(s(1)+s(6)))/YM
	  e(3,3)=(s(6)-PR*(s(1)+s(4)))/YM
	  e(1,2)=sum*s(2)
	  e(2,3)=sum*s(5)
	  e(3,1)=sum*s(3)
	  e(2,1)=e(1,2)
	  e(3,2)=e(2,3)
	  e(1,3)=e(3,1)
*
	  delV=e(1,1)+e(2,2)+e(3,3)
*
	  return
	  end



	  SUBROUTINE TREND(VX,VY,VZ,TX,PX,TY,PY,TZ,PZ)
**************************************************************************
*
*	Subroutine TREND finds the plunge and trend of the principal 
*	strains.  It does it in an extraordinarily difficult and 
*	cumbersome way, not even taking advantage og their orthogonal
*	nature.  Never mind.  Labels in the range 1000-1090.
*
**************************************************************************
 	REAL*4 VX(1,3),VY(1,3),VZ(1,3)
	REAL*4 TX,PX,TY,PY,TZ,PZ,INC
	REAL*4 L,T,CON
 	INC=1.5707963
	CON=180.0/3.1415927

	TX=0.
	PX=0.
	TY=0.
	PY=0.
	TZ=0.
	PZ=0.
*
*	FIRST THE X AXIS
*
*	This IF statement prevents the null result of a vector
*	being perfectly vertical by forcing it to plunge steeply 
*	to the north.
*
	IF(VX(1,1).EQ.0.0.AND.VX(1,2).EQ.0.0)THEN
		PX=INC
		TX=0.0
		GO TO 1005
	ENDIF
*
*
	IF(VX(1,3).GE.0.0.AND.VX(1,1).EQ.0.0.AND.VX(1,2).GT.0.0)THEN
		TX=INC
		GO TO 1000
	ELSEIF(VX(1,3).GT.0.0.AND.VX(1,1).EQ.0.0.AND.VX(1,2).LT.0.0)THEN
		TX=3.0*INC
		GO TO 1000
	ELSEIF(VX(1,3).LT.0.0.AND.VX(1,1).EQ.0.0.AND.VX(1,2).GT.0.0)THEN
		TX=3.0*INC
		GO TO 1000
	ELSEIF(VX(1,3).LE.0.0.AND.VX(1,1).EQ.0.0.AND.VX(1,2).LT.0.0)THEN
		TX=INC
		GO TO 1000
	ENDIF
*
	IF(VX(1,1).EQ.0.0) WRITE(*,*)'VX IS IT'
	T=ATAN(VX(1,2)/VX(1,1))
*
	IF(VX(1,3).GT.0.0.AND.VX(1,1).GT.0.0.AND.VX(1,2).GT.0.0)THEN
		TX=2.0*INC-T
	ELSEIF(VX(1,3).GT.0.0.AND.VX(1,1).LT.0.0.AND.VX(1,2).GT.0.0)THEN
		TX=-1.0*T
	ELSEIF(VX(1,3).GT.0.0.AND.VX(1,1).LT.0.0.AND.VX(1,2).LE.0.0)THEN
		TX=4.0*INC-T
	ELSEIF(VX(1,3).GT.0.0.AND.VX(1,1).GT.0.0.AND.VX(1,2).LE.0.0)THEN
		TX=2.0*INC-T
*
	ELSEIF(VX(1,3).LE.0.0.AND.VX(1,1).GT.0.0.AND.VX(1,2).GE.0.0)THEN
		TX=4.0*INC-T
	ELSEIF(VX(1,3).LE.0.0.AND.VX(1,1).LT.0.0.AND.VX(1,2).GE.0.0)THEN
		TX=2.0*INC-T
	ELSEIF(VX(1,3).LE.0.0.AND.VX(1,1).LT.0.0.AND.VX(1,2).LT.0.0)THEN
		TX=2.0*INC-T
	ELSEIF(VX(1,3).LE.0.0.AND.VX(1,1).GT.0.0.AND.VX(1,2).LT.0.0)THEN
		TX=-1.0*T
	ENDIF
*
*
 1000	CONTINUE
*
	  L=SQRT(VX(1,1)*VX(1,1)+VX(1,2)*VX(1,2))
	  PX=ATAN(VX(1,3)/L)
*
 1005	TX=TX*CON
		IF(TX.GE.359.9) TX=0.0
	  PX=ABS(PX*CON)
*
*	NOW THE Y AXIS          ******************************************
*
*	This IF statement prevents the null result of a vector
*	being perfectly vertical by forcing it to plunge steeply 
*	to the north.
*
	IF(VY(1,1).EQ.0.0.AND.VY(1,2).EQ.0.0)THEN
		PY=INC
		TY=0.0
		GO TO 1015
	ENDIF
*
*
	IF(VY(1,3).GE.0.0.AND.VY(1,1).EQ.0.0.AND.VY(1,2).GT.0.0)THEN
		TY=INC
		GO TO 1010
	ELSEIF(VY(1,3).GT.0.0.AND.VY(1,1).EQ.0.0.AND.VY(1,2).LT.0.0)THEN
		TY=3.0*INC
		GO TO 1010
	ELSEIF(VY(1,3).LT.0.0.AND.VY(1,1).EQ.0.0.AND.VY(1,2).GT.0.0)THEN
		TY=3.0*INC
		GO TO 1010
	ELSEIF(VY(1,3).LE.0.0.AND.VY(1,1).EQ.0.0.AND.VY(1,2).LT.0.0)THEN
		TY=INC
		GO TO 1010
	ENDIF
*
	IF(VY(1,1).EQ.0.0) WRITE(*,*)'VY IS IT'
	TT=ATAN(VY(1,2)/VY(1,1))
*
	IF(VY(1,3).GT.0.0.AND.VY(1,1).GT.0.0.AND.VY(1,2).GT.0.0)THEN
		TY=2.0*INC-TT
	ELSEIF(VY(1,3).GT.0.0.AND.VY(1,1).LT.0.0.AND.VY(1,2).GT.0.0)THEN
		TY=-1.0*TT
	ELSEIF(VY(1,3).GT.0.0.AND.VY(1,1).LT.0.0.AND.VY(1,2).LE.0.0)THEN
		TY=4.0*INC-TT
	ELSEIF(VY(1,3).GT.0.0.AND.VY(1,1).GT.0.0.AND.VY(1,2).LE.0.0)THEN
		TY=2.0*INC-TT
*
	ELSEIF(VY(1,3).LE.0.0.AND.VY(1,1).GT.0.0.AND.VY(1,2).GE.0.0)THEN
		TY=4.0*INC-TT
	ELSEIF(VY(1,3).LE.0.0.AND.VY(1,1).LT.0.0.AND.VY(1,2).GE.0.0)THEN
		TY=2.0*INC-TT
	ELSEIF(VY(1,3).LE.0.0.AND.VY(1,1).LT.0.0.AND.VY(1,2).LT.0.0)THEN
		TY=2.0*INC-TT
	ELSEIF(VY(1,3).LE.0.0.AND.VY(1,1).GT.0.0.AND.VY(1,2).LT.0.0)THEN
		TY=-1.0*TT
	ENDIF
*
*
 1010	CONTINUE
*
	L=SQRT(VY(1,1)*VY(1,1)+VY(1,2)*VY(1,2))
	PY=ATAN(VY(1,3)/L)
*
 1015	TY=TY*CON
		IF(TY.GE.359.9) TY=0.0
    	PY=ABS(PY*CON)
*
*	AND NOW THE Z AXIS (THE HARD WAY!)    ************************************
*
*	This IF statement prevents the null result of a vector
*	being perfectly vertical by forcing it to plunge steeply 
*	to the north.
*
	IF(VZ(1,1).EQ.0.0.AND.VZ(1,2).EQ.0.0)THEN
		PZ=INC
		TZ=0.0
		GO TO 1025
	ENDIF
*
*
	IF(VZ(1,3).GE.0.0.AND.VZ(1,1).EQ.0.0.AND.VZ(1,2).GT.0.0)THEN
		TZ=INC
		GO TO 1020
	ELSEIF(VZ(1,3).GT.0.0.AND.VZ(1,1).EQ.0.0.AND.VZ(1,2).LT.0.0)THEN
		TZ=3.0*INC
		GO TO 1020
	ELSEIF(VZ(1,3).LT.0.0.AND.VZ(1,1).EQ.0.0.AND.VZ(1,2).GT.0.0)THEN
		TZ=3.0*INC
		GO TO 1020
	ELSEIF(VZ(1,3).LE.0.0.AND.VZ(1,1).EQ.0.0.AND.VZ(1,2).LT.0.0)THEN
		TZ=INC
		GO TO 1020
	ENDIF
*
	IF(VZ(1,1).EQ.0.0) WRITE(*,*)'VZ IS IT'
	TTT=ATAN(VZ(1,2)/VZ(1,1))
*
	IF(VZ(1,3).GT.0.0.AND.VZ(1,1).GT.0.0.AND.VZ(1,2).GT.0.0)THEN
		TZ=2.0*INC-TTT
	ELSEIF(VZ(1,3).GT.0.0.AND.VZ(1,1).LT.0.0.AND.VZ(1,2).GT.0.0)THEN
		TZ=-1.0*TTT
	ELSEIF(VZ(1,3).GT.0.0.AND.VZ(1,1).LT.0.0.AND.VZ(1,2).LE.0.0)THEN
		TZ=4.0*INC-TTT
	ELSEIF(VZ(1,3).GT.0.0.AND.VZ(1,1).GT.0.0.AND.VZ(1,2).LE.0.0)THEN
		TZ=2.0*INC-TTT
*
	ELSEIF(VZ(1,3).LE.0.0.AND.VZ(1,1).GT.0.0.AND.VZ(1,2).GE.0.0)THEN
		TZ=4.0*INC-TTT
	ELSEIF(VZ(1,3).LE.0.0.AND.VZ(1,1).LT.0.0.AND.VZ(1,2).GE.0.0)THEN
		TZ=2.0*INC-TTT
	ELSEIF(VZ(1,3).LE.0.0.AND.VZ(1,1).LT.0.0.AND.VZ(1,2).LT.0.0)THEN
		TZ=2.0*INC-TTT
	ELSEIF(VZ(1,3).LE.0.0.AND.VZ(1,1).GT.0.0.AND.VZ(1,2).LT.0.0)THEN
		TZ=-1.0*TTT
	ENDIF
*
*
 1020	CONTINUE
*
	L=SQRT(VZ(1,1)*VZ(1,1)+VZ(1,2)*VZ(1,2))
	PZ=ATAN(VZ(1,3)/L)
*
 1025	TZ=TZ*CON
		IF(TZ.GE.359.9) TZ=0.0
		PZ=ABS(PZ*CON)
*
*
*	Adjustment (10/19/92) for the fact that the global Y axis is geographic
*	north and the X is east.
 
		tx=tx+90.
		if(tx.gt.360.)tx=tx-360.
		ty=ty+90.
		if(ty.gt.360.)ty=ty-360.
		tz=tz+90.
		if(tz.gt.360.)tz=tz-360.
	  return
	  END
	
	
	  SUBROUTINE SORT(E,EVEC,X,Y,Z,VX,VY,VZ)
**************************************************************************
*
*	Subroutine SORT sorts the eigenvalues by size, X being the largest.
*	Labels in range 800-890. It also copies the appropriate eigenvectors
*	from EVEC to VX, VY, and VZ.
*
**************************************************************************
	REAL*4 E(3),EVEC(3,3),VX(3),VY(3),VZ(3)
	REAL*4 X,Y,Z
*
	IF(E(1).GE.E(2).AND.E(1).GE.E(3)) THEN
		X=E(1)
		GO TO 800
	ELSE IF(E(2).GE.E(1).AND.E(2).GE.E(3)) THEN
		X=E(2)
		GO TO 810
	ELSE 
		X=E(3)
	END IF
	IF(E(2).GE.E(1)) THEN
		Y=E(2)
		Z=E(1)
	ELSE IF(E(1).GE.E(2)) THEN
		Y=E(1)
		Z=E(2)
	END IF
	GO TO 820
 800	IF(E(2).GE.E(3)) THEN
 		Y=E(2)
		Z=E(3)
	ELSE
		Y=E(3)
		Z=E(2)
	END IF
	GO TO 820
 810	IF(E(1).GE.E(3)) THEN
 		Y=E(1)
		Z=E(3)
	ELSE
		Y=E(3)
		Z=E(1)
	END IF
 820	CONTINUE 
 		IF(E(1).EQ.X)THEN
			DO 825 I=1,3
825			VX(I)=EVEC(I,1)
 		ELSEIF(E(2).EQ.X)THEN
			DO 830 I=1,3
830			VX(I)=EVEC(I,2)
 		ELSE
			DO 835 I=1,3
835			VX(I)=EVEC(I,3)
 		ENDIF
*
		IF(E(1).EQ.Y)THEN
			DO 840 I=1,3
840			VY(I)=EVEC(I,1)
 		ELSEIF(E(2).EQ.Y)THEN
			DO 845 I=1,3
845			VY(I)=EVEC(I,2)
 		ELSE
			DO 850 I=1,3
850			VY(I)=EVEC(I,3)
 		ENDIF
*
		IF(E(1).EQ.Z)THEN
			DO 855 I=1,3
855			VZ(I)=EVEC(I,1)
 		ELSEIF(E(2).EQ.Z)THEN
			DO 860 I=1,3
860			VZ(I)=EVEC(I,2)
 		ELSE
			DO 865 I=1,3
865			VZ(I)=EVEC(I,3)
 		ENDIF
*
*	The following lines redefine the eigenvector of the minimum principal axis
*	so that the three axes define a right-handed coordinate system, by using the
*	cross-product of the other two principal axes.
*
	VZ(1)=VX(2)*VY(3)-VX(3)*VY(2)
	VZ(2)=VX(3)*VY(1)-VX(1)*VY(3)
	VZ(3)=VX(1)*VY(2)-VX(2)*VY(1)

	END
	
	
	
	SUBROUTINE FINGER (D, F)
******************************************************************************
	REAL*4 D(3,3), F(3,3)
*
	F(1,1)=D(1,1)*D(1,1)+D(1,2)*D(1,2)+D(1,3)*D(1,3)
	F(2,2)=D(2,2)*D(2,2)+D(2,1)*D(2,1)+D(2,3)*D(2,3)
	F(3,3)=D(3,3)*D(3,3)+D(3,2)*D(3,2)+D(3,1)*D(3,1)
	F(1,2)=D(2,1)*D(1,1)+D(2,2)*D(1,2)+D(2,3)*D(1,3)
	F(2,1)=F(1,2)
	F(1,3)=D(3,1)*D(1,1)+D(3,2)*D(1,2)+D(3,3)*D(1,3)
	F(3,1)=F(1,3)
	F(2,3)=D(3,1)*D(2,1)+D(3,2)*D(2,2)+D(3,3)*D(2,3)
	F(3,2)=F(2,3)
*
 	END


	SUBROUTINE JACOBI(M,D,V)
******************************************************************************
*
*	Subroutine JACOBI finds the eigenvalues and eigenvectors of a symmetric
*	3x3 matrix.  The matrix is input as M and transferred to A, since the 
*	contents of A are destroyed. D is a vector containing the unsorted 
*	eigenvalues, and V is a 3x3 matrix whose columns contain the eigenvectors
*	corresponding to the eigenvalues in D.  The routine has been (very slightly)
*	modified from that in NUMERICAL RECIPES, pp 346.
*
*******************************************************************************
	PARAMETER (NMAX=100)
	REAL*4 A(3,3),M(3,3),D(3),V(3,3),B(NMAX),Z(NMAX)
	DO 700 I=1,3
		DO 710 J=1,3
		A(I,J)=M(I,J)
 710	CONTINUE
 700	CONTINUE
*
	DO 730 IP=1,3
		DO 720 IQ=1,3
		V(IP,IQ)=0.0
 720	CONTINUE
 	V(IP,IP)=1.0
 730	CONTINUE
 	DO 740 IP=1,3
		B(IP)=A(IP,IP)
		D(IP)=B(IP)
		Z(IP)=0.0
 740	CONTINUE
 	NROT=0
	DO 775 I=1,50
		SM=0.0
		DO 760 IP=1,2
			DO 750 IQ=IP+1,3
			SM=SM+ABS(A(IP,IQ))
 750	CONTINUE
 760	CONTINUE
 	IF(SM.EQ.0.0)RETURN
	IF(I.LT.4)THEN
		TRESH=0.2*SM/9.
	ELSE
		TRESH=0.0
	ENDIF
	DO 798 IP=1,2
		DO 796 IQ=IP+1,3
		G=100.0*ABS(A(IP,IQ))
		IF((I.GT.4).AND.(ABS(D(IP))+G.EQ.ABS(D(IP)))
     >	.AND.(ABS(D(IQ))+G.EQ.ABS(D(IQ))))THEN
     		A(IP,IQ)=0.0
		ELSEIF(ABS(A(IP,IQ)).GT.TRESH)THEN
		H=D(IQ)-D(IP)
		IF(ABS(H)+G.EQ.ABS(H))THEN
		T=A(IP,IQ)/H
		ELSE
		THETA=0.5*H/A(IP,IQ)
		T=1.0/(ABS(THETA)+SQRT(1.0+THETA**2))
		IF(THETA.LT.0.0)T=-T
		ENDIF
		C=1.0/SQRT(1.0+T**2)
		S=T*C
		TAU=S/(1.0+C)
		H=T*A(IP,IQ)
		Z(IP)=Z(IP)-H
		Z(IQ)=Z(IQ)+H
		D(IP)=D(IP)-H
		D(IQ)=D(IQ)+H
		A(IP,IQ)=0.0
		DO 770 J=1,IP-1
			G=A(J,IP)
			H=A(J,IQ)
			A(J,IP)=G-S*(H+G*TAU)
			A(J,IQ)=H+S*(G-H*TAU)
 770	CONTINUE
 		DO 780 J=IP+1,IQ-1
			G=A(IP,J)
			H=A(J,IQ)
			A(IP,J)=G-S*(H+G*TAU)
			A(J,IQ)=H+S*(G-H*TAU)
 780	CONTINUE
 		DO 790 J=IQ+1,3
			G=A(IP,J)
			H=A(IQ,J)
			A(IP,J)=G-S*(H+G*TAU)
			A(IQ,J)=H+S*(G-H*TAU)
 790	CONTINUE
 		DO 792 J=1,3
			G=V(J,IP)
			H=V(J,IQ)
			V(J,IP)=G-S*(H+G*TAU)
			V(J,IQ)=H+S*(G-H*TAU)
 792	CONTINUE
 	NROT=NROT+1
	ENDIF
 796	CONTINUE
 798	CONTINUE
 	DO 785 IP=1,3
		B(IP)=B(IP)+Z(IP)
		D(IP)=B(IP)
		Z(IP)=0.0
 785 	CONTINUE
 775	CONTINUE
* 	WRITE(*,*)'50 ITERATIONS SHOULD NEVER HAPPEN'
	RETURN
	END

	subroutine failure_planes(sp1,d1,rake1,sp2,d2,rake2)
****************************************************************************8
*
*	Subroutine to find the potential failure planes assuming a Coulomb
*	failure criterion in which the two optimal planes of failure are 
*	oriented at +/- 0.5invtan(1/mew), where mew is the coefficient of 
*	internal friction.  Each failure plane contains the intermediate 
*	principal stress axis.  The Coulomb criterion is probably not
*	appropriate for complex strain fields, but there are no other
*	well-known or well-established criteria, so Coulomb will suffice
*	for now.
*	OUTPUT:
*	  sp1 = strike of plane 1
*	  sp2 = strike of plane 2
*	  d1  = dip of plane 1
*	  d2  = dip of plane 2
*	  rake1 = rake of plane 1
*	  rake2 = rake of plane 2
*
*****************************************************************************
*
	real*4 rot(3,3),s(3,3),vx(3),vy(3),vz(3),p(3),v(3,3)
	real*4 sfv1(3),sfv2(3),sv1(3),sv2(3),pn1(3),pn2(3)
	real*4 p1(3),p2(3)
	real*4 sig1,sig2,sig3,sp1,d1,sp2,d2,x,y,z
        real*4 rake1,rake2
* 
	real*4 stress,dum,dum2,vee,eee,friction 
	real*4 rot1(3,3),rot2(3,3),reig(3,3)
	real*4 rprin1(3,3),rprin2(3,3)
	real*4 s1inp(3),s2inp(3)
	common/mike/stress(6),dum(3),dum2(3,3),vee,eee,friction,
     *               junk1,junk2,junk3,junk4
*
*	beta is the angle between the failure plane normal and sigma1
*
	beta=1.5707963-0.5*atan(1./friction)

	todeg=57.295779513
	pi = 3.141592654
	pi2 = pi/2.

	s(1,1)=stress(1)
	s(1,2)=stress(2)
	s(1,3)=stress(3)
	s(2,1)=s(1,2)
	s(2,2)=stress(4)
	s(2,3)=stress(5)
	s(3,1)=s(1,3)
	s(3,2)=s(2,3)
	s(3,3)=stress(6)
C
C	Find Eigenvectors and Eigenvalues
C
	call jacobi(s,p,v)
C
C	The  eigenvectors will be sorted and will also be output
C	in a RHS.
C
	call sort(p,v,sig1,sig2,sig3,vx,vy,vz)
*
*	Make the vectors for failure_plane_normals    
*	with respect to the principal axes.
*
	pn1(1)=cos(beta)
	pn1(2)=0.0
	pn1(3)=sin(beta)
*
	pn2(1)=cos(-beta)
	pn2(2)=0.0
	pn2(3)=sin(-beta)
C
C	Put together rotation matrix
C
	rot(1,1) = vx(2)
	rot(2,1) = -vx(1)
	rot(3,1) = vx(3)
	
	rot(1,2) = vy(2)
	rot(2,2) = -vy(1)
	rot(3,2) = vy(3)

	rot(1,3) = vz(2)
	rot(2,3) = -vz(1)
	rot(3,3) = vz(3)
*
*	Now rotate these into the geopgraphic axes, using the 
*	rotation matrix just formed from the eigenvectors.       		
*
	do 20 i=1,3
	   p1(i)=rot(i,1)*pn1(1)+rot(i,2)*pn1(2)+rot(i,3)*pn1(3)
	   p2(i)=rot(i,1)*pn2(1)+rot(i,2)*pn2(2)+rot(i,3)*pn2(3)
20	continue
C
C	Find the strike and dip of both failure planes
C	   sp1: strike of plane 1 in radians
C	   sp2: strike of plane 2 in radians
C	   d1: dip of plane 1 in radians
C	   d2: dip of plane 2 in radians
C
C	Note:
C	   x == north
C	   y == west
C	   z == vertical
C
	x = p1(1)
	y = p1(2)
	z = p1(3)
*
	angle=abs(atan(x/(-y)))
	d1=acos(abs(z))

	if(z.le.0.0)then
		if(y.ge.0.0.and.x.ge.0.0)sp1=angle
		if(y.ge.0.0.and.x.lt.0.0)sp1=2*pi-angle
		if(y.lt.0.0.and.x.le.0.0)sp1=pi+angle
		if(y.lt.0.0.and.x.gt.0.0)sp1=pi-angle
	else
		if(y.ge.0.0.and.x.ge.0.0)sp1=pi+angle
		if(y.ge.0.0.and.x.lt.0.0)sp1=pi-angle
		if(y.lt.0.0.and.x.le.0.0)sp1=angle
		if(y.lt.0.0.and.x.gt.0.0)sp1=2*pi-angle
	endif
*
	x=p2(1)
	y=p2(2)
	z=p2(3)
*
	angle=abs(atan(x/(-y)))
	d2=acos(abs(z))

	if(z.le.0.0)then
		if(y.ge.0.0.and.x.ge.0.0)sp2=angle
		if(y.ge.0.0.and.x.lt.0.0)sp2=2*pi-angle
		if(y.lt.0.0.and.x.le.0.0)sp2=pi+angle
		if(y.lt.0.0.and.x.gt.0.0)sp2=pi-angle
	else
		if(y.ge.0.0.and.x.ge.0.0)sp2=pi+angle
		if(y.ge.0.0.and.x.lt.0.0)sp2=pi-angle
		if(y.lt.0.0.and.x.le.0.0)sp2=angle
		if(y.lt.0.0.and.x.gt.0.0)sp2=2*pi-angle
	endif
C
C	Create rotation matrix
C
	rot1(1,1) = cos(sp1-pi2) 
	rot1(1,2) = cos(sp1) 
	rot1(1,3) =  0.0
	rot1(2,1) = cos(pi-sp1)*cos(d1)
	rot1(2,2) = cos(sp1-pi2)*cos(d1)
	rot1(2,3) = sin(d1) 
	rot1(3,1) = cos(pi-sp1)*cos(pi2+d1) 
	rot1(3,2) = cos(sp1-pi2)*cos(pi2+d1) 
	rot1(3,3) = cos(d1) 
	
	rot2(1,1) = cos(sp2-pi2) 
	rot2(1,2) = cos(sp2) 
	rot2(1,3) =  0.0
	rot2(2,1) = cos(pi-sp2)*cos(d2)
	rot2(2,2) = cos(sp2-pi2)*cos(d2)
	rot2(2,3) = sin(d2) 
	rot2(3,1) = cos(pi-sp2)*cos(pi2+d2) 
	rot2(3,2) = cos(sp2-pi2)*cos(pi2+d2) 
	rot2(3,3) = cos(d2) 
C
C	Now rotate to in-failure plane
C
	s1inp(3)=rot1(3,1)*vx(1)+rot1(3,2)*vx(2)+rot1(3,3)*vx(3)
	s2inp(3)=rot2(3,1)*vx(1)+rot2(3,2)*vx(2)+rot2(3,3)*vx(3)
C
C	Rotation matrix from global to principal axis of the stress tensor
C
	reig(1,1) = vx(1)
	reig(2,1) = vx(2)
	reig(3,1) = vx(3)
	reig(1,2) = vy(1)
	reig(2,2) = vy(2)
	reig(3,2) = vy(3)
	reig(1,3) = vz(1)
	reig(2,3) = vz(2)
	reig(3,3) = vz(3)
C
C	Is Sigma 1 in the Foot wall or the hanging wall?
C	To determine:
C	   Foot wall if s1inp(3) > 0.0
C	   Hanging wal if s1inp(3) >= 0.0
C
	if (abs(s1inp(3)).lt.1e-6) then
	   write(77,*)'ERROR: s1inp(3) is zero!!!!'
	   write(77,*)'     ****** STOP NOW!'
	   stop
	endif
	if (abs(s2inp(3)).lt.1e-6) then
	   write(77,*)'ERROR: s2inp(3) is zero!!!!'
	   write(77,*)'     ****** STOP NOW!'
	   stop
	endif
	if (s1inp(3).lt.0.0) then
CC	   write(77,*)'For Plane 1 sigma1 is in the footwall.'
	   sv1(1) = sin(beta)
	   sv1(2) = 0.0
	   sv1(3) = cos(pi-beta)
	else
CC	   write(77,*)'For Plane 1 sigma1 is in the hanging-wall.'
	   sv1(1) = cos(pi2 + beta)
	   sv1(2) = 0.0
	   sv1(3) = cos(beta)
	endif
	if (s2inp(3).lt.0.0) then
CC	   write(77,*)'For Plane 2 sigma1 is in the foot wall.'
	   sv2(1) = sin(beta)
	   sv2(2) = 0.0
	   sv2(3) = cos(beta)
	else
CC	   write(77,*)'For Plane 2 sigma1 is in the hanging wall.'
	   sv2(1) = cos(pi2 + beta)
	   sv2(2) = 0.0
	   sv2(3) = cos(pi - beta)
	endif
C
C	Compute Rotation matrix needed to go from in-failure plane
C	to principal axis 
C	   rprin1: principal axis number 1
C	   rprin2: printcipal axis number 2
C	Note:
C	  (in-failure):(Geograph) * (Geographic):(Principal axis)
C	is the same as
C	  (in-failure):(global) * (global):(Principal axis)
C
	ndim = 3
	call matmult(rot1,reig,rprin1,ndim)
	call matmult(rot2,reig,rprin2,ndim)
C
C	Use rotation matrix to rotate into the in-failure plane
C
	do 70 i=1,3
	   sfv1(i)=rprin1(i,1)*sv1(1)+rprin1(i,2)*sv1(2)
     *           +rprin1(i,3)*sv1(3)
	   sfv2(i)=rprin2(i,1)*sv2(1)+rprin2(i,2)*sv2(2)
     *           +rprin2(i,3)*sv2(3)
70	continue
C
C	Compute Rake
C
	alpha1 = atan2(sfv1(2),sfv1(1))
	alpha2 = atan2(sfv2(2),sfv2(1))

	if (sfv1(2).gt.0.0) then
	   rake1 = abs(alpha1)
	else
	   rake1 = -abs(alpha1)
	endif
C
	if (sfv2(2).gt.0.0) then
	   rake2 = abs(alpha2)
	else
	   rake2 = -abs(alpha2)
	endif
C
C	Convert everything to Degrees
C
	rake1 = rake1*todeg
	rake2 = rake2*todeg

	sp1 = sp1*todeg
	sp2 = sp2*todeg

	d1 = d1*todeg
	d2 = d2*todeg

	end

	subroutine matmult(a,b,c,ndim)
C *********************************************************************
C
C	Matrix multiplication of two square matrices of dimension
C	ndim.  Matrix a and b are input....the result is passed
C	back in matrix c
C
C *********************************************************************

	real*4 a(ndim,ndim),b(ndim,ndim),c(ndim,ndim)

	do 20 i=1,ndim
	   do 10 j=1,ndim
	      c(i,j) = a(i,1)*b(1,j) + a(i,2)*b(2,j) + a(i,3)*b(3,j)
10	   continue
20	continue

	return
	end
