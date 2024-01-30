	  SUBROUTINE GRID(NUM_Ds,NPLANE,xg,yg,zg,
     &              inlout_displ,inlout_stress,inlout_strain,
     &              inlout_orient,inlout_failure,inlout_elements,
     &              inlout_ugrad,
     &              npts)
C******************************************************************************
C   This code uses 3-d dislocation solutions for the stress tensor, displacement
C   gradients and the displacement vectors as a function of position. 
C
C	The input file is assumed to be still open.  If the program is being run
C	in a non-interactive mode, the user specifies all the options describing
C	the grid parameters, etc. to the bottom of the input file.  If it is run 
C	interactively, then the user must enter grid parameters interactively 
C   (there are lots of questions to answer!).
C
C   All routines that control output to X,Y,Z,Data formatted files are 
C   contained in this file.  These output files are read by routines in
C	the file OUTPUT_VECTOR.F and re-formatted in a vector format (e.g.
C   as required by plotting programs such as MATLAB); this re-formatting
C   is of course optional.
C******************************************************************************

	  INTEGER*4 NUM_Ds,NPLANE

C  Added by A.JANIN
      integer npts
	  real*8 xg(npts), yg(npts), zg(npts)
	  real*8 inlout_displ(npts,3)
	  real*8 inlout_stress(npts,6)
	  real*8 inlout_strain(npts,6)
	  real*8 inlout_failure(npts,6)
	  real*8 inlout_orient(npts,9)
	  real*8 inlout_ugrad(npts,3,3)
	  real*8 inlout_elements(NPLANE,6)

	  call userco(NPLANE,NUM_Ds,xg,yg,zg,
     &              inlout_displ,inlout_stress,inlout_strain,
     &              inlout_orient,inlout_failure,inlout_elements,
     &              inlout_ugrad,
     &              npts)

	return
	end
	  
	  subroutine final_write(NPLANE,NUM_Ds,inlout_elements)
C****************************************************************************
C  Write out element displacements in X,Y,Z,Data format if desired.  Read
C  X,Y,Z,Data formatted files are write new ones in vector format if
C  desired.
C
C  Close all files.
C****************************************************************************
	  include 'sizes.inc'
	  
	  integer*4 NPLANE
	  
	  logical o_elem,o_stress,o_strain,o_dgrad,o_rbr,
     &  o_orient,o_disp,o_inv,xyz_out,o_fail,o_proj
	  character*24 sufx
	  common/outputs/sufx,o_elem,o_stress,o_strain,o_dgrad,
     &  o_rbr,o_orient,o_disp,o_inv,xyz_out,o_fail,o_proj

	  real*8 inlout_elements(NPLANE,6)
 
	  if(o_elem) then
	      call XYZ_elems(NPLANE,NUM_Ds,inlout_elements)
	  endif


	  return
	  end


        subroutine get_fopts
C****************************************************************************
C This routine reads the input file options and sets output file options
C accordingly.  The input file should look like
C       * local coordinate output (global otherwise) - only for planar grids
C       global
C  	* Output information (each line results in a separate output file)
C	stress
C	strain
C	gradients of displacement
C	rigid body rotations
C	orientation of principal strain axes
C	displacements
C       failure plane 
C	invarients of the strain field
C	element relative displacements
C       * Output file suffix
C       dat1
C Remember that lines beginning with * are only descriptions and aren't used
C by the program (they can say anything).  The list of files containing various
C characteristics of the calculated deformation field can be in any order. 
C The word on each line sets an appropriate flag so that an output file 
C containing the parameters described by the word will be generated.  
C The program only interprets the first letter (the first 4 letters for
C 'stress' and 'strain') so that their is some flexibility in what words
C are used.  For example, you must write 'gradients of displacement' or
C simply 'gradients' to generate an output file containing the components
C of the displacement gradient tensor; use of 'displacement gradients'
C would be incorrect (the key word must start with a 'g') and result in
C an output file containing displacements.
C****************************************************************************
		  
	logical o_elem,o_stress,o_strain,o_dgrad,o_rbr,
     &  o_orient,o_disp,o_inv,xyz_out,o_fail,o_proj
	  character*24 sufx
	  common/outputs/sufx,o_elem,o_stress,o_strain,o_dgrad,
     &  o_rbr,o_orient,o_disp,o_inv,xyz_out,o_fail,o_proj
		
C       xyz_out  - obsolete argument now - A.JANIN
	    xyz_out  = .true.

		o_elem   = .true.
		o_disp   = .true.
		o_stress = .true.
		o_strain = .true.
		o_orient = .true.
		o_inv    = .true.
		o_rbr    = .true.
		o_dgrad  = .false.
		o_fail   = .true.
		o_proj   = .false.
        sufx     = 'out'

        write(*,'(a)') ' '
        return
	end

	  
	  subroutine XYZ_Results(x,y,z,NPLANE,NUM_Ds)
C******************************************************************************
C  Subroutine to calculate various parameters describing the deformation field.
C  These are output to ascii files. 
C
C  In this routine seven files are written to; at each point X,Y,Z in global 
C  coordinates the following is output (the coordinates of each point are also
C  written to the file):
C    File 		         Contents
C    stressten.??? - Six independent elements of the stress tensor, 
C		     Sxx,Sxy,Sxz,Syy,Syz,Szz. If output on a plane is 
C                    requested these correspond to the in-plane coordinates 
C		     (e.g. Szz is the stress normal to the fault plane,
C                    Sxy is the traction in the strike direction,
C                        is the traction in the up-dip direction, etc.)
C    strainten.??? - Six independent elements of the strain tensor, Exx,
C		     Exx,Exy,Exz,Eyy,Eyz,Ezz, and the volumetric strain 
C		     (Exx+Eyy+Ezz,independent of the coordinate syst
C		     If output on a plane is requested these correspond to  
C                    the in-plane coordinates as in stressten.out.
C   displacements.??? - Displacements in the x,y,z directions. If output 
C	             on a plane is requested these correspond to thespond to the 
C		     in-plane coordinates; Ux is in the strike direction,directio
C                    Uy in the up-dip direction, Uz is the normal direction.l to the plane.
C    dgten.??? -     Displacement gradient tensor; dUx/dx,dUx/dy,dUx/dz,
C		     dUy/dx,dUy/dy,dUy/dz, dUz/dx,dUz/dy,dUz/dz.  As in,
C                    the above, these correspond to the in-plane coordinates
C                    if the user has selected the plane option.
C    strainorient.??? -	 Principal strains,their plunge and trend with respect
C		     to the global axes.  They are listed from max. to min.from
C    invariants.??? -Strain & stress invarients: the volume change, failure 
C	             stresses,octahedral shear,strain energy density.
C    rbr.??? -	     Rotations about the global x,y,and z axes
C    failure.??? -   Failure plane descriptions: strike,dip,slip vector
C                    direction cosines for each plane.
C  Remember that if you change the format of these files, you must also
C  change the routines in OUTPUT_VECTOR.f since they must read these files.
C******************************************************************************

	  logical o_elem,o_stress,o_strain,o_dgrad,o_rbr,
     &  o_orient,o_disp,o_inv,xyz_out,o_fail,o_proj
	  character*24 sufx
	  common/outputs/sufx,o_elem,o_stress,o_strain,o_dgrad,
     &  o_rbr,o_orient,o_disp,o_inv,xyz_out,o_fail,o_proj

	  REAL*8  x,y,z,delV
	  REAL*8 zd
	  	  
      REAL*4 E,dum
      COMMON/TEMPS/e(3,3),dum(6,3)

	  REAL*4 STRESS,EWE,dgrten,pr,ym,friction,tmax,tmaxo,tmaxa
	  COMMON/MIKE/STRESS(6),EWE(3),dgrten(3,3),
     &                pr,ym,friction,tmax,tmaxo(2),tmaxa

	  REAL*4 SSTRIKE,CSTRIKE,SDYP,CDYP
	  COMMON/TOPLN/SSTRIKE,CSTRIKE,SDYP,CDYP

C Added by A.JANIN - export these global output to userco
	  real*4 ex,ey,ez,px,tx,py,ty,pz,tz
	  real*4  volchg,critic,octshr,work
	  real*4  sp1,d1,rake1,sp2,d2,rake2
	  common/JANIN/ex,ey,ez,px,tx,py,ty,pz,tz,
     &             volchg,critic,octshr,work,
     &             sp1,d1,rake1,sp2,d2,rake2
	  
	  DATA TODEG/57.295779/

	  zd=z
	  
C         Gets stress and displacements
	  call CALCPO(x,y,zd,NPLANE,NUM_Ds,ierr)

C	  PRINT*, x,y,zd,ewe
 
C         Gets strain tensor and vol. change
	  call strain(e,delV)

C Gets principal strains and orientations and write to output file
	  if(o_orient) then
		call third_results(ex,ey,ez,px,tx,py,ty,pz,tz)
	  endif
 	
C     Rotate to in-plane coordinates 
      if(o_proj) CALL TOPLANE(SSTRIKE,CSTRIKE,SDYP,CDYP)
C
C  Code added 2/15/94 to get max shear stress and its orientation on the inspection 
C  plane.  See subroutine TOPLANE in 3dmain.f and manual for more info.
C
	  tmaxa=tmaxa*TODEG

	  if(o_inv) then
	    call invariants(volchg,critic,octshr,work)
	  endif
 	
	  if(o_fail) then
	    call failure_planes(sp1,d1,rake1,sp2,d2,rake2)
	  endif

      return
      end

	  subroutine XYZ_elems(NPLANE,NUM_Ds,inlout_elements)
C********************************************************************
C Subroutine to write out element displacements (relative) in
C X,Y,Z,Data format.  Each line of the output file will contain the
C coordinates of the center of the element in global and in-plane
C coordinate systems,and the relative
C displacements in the strike, dip, and tensile (normal)
C directions.
C********************************************************************

      INCLUDE 'sizes.inc'
      REAL*4 XMATRIX
      COMMON/SOLN/XMATRIX(MAX3_ELEM,MAX3_ELEM+1)
      REAL*4 UG2P,SG2P
      COMMON/TRANS/UG2P(3,3,MAX_PLN),SG2P(3,6,MAX_PLN)

      REAL*4 XO,YO,ZO,C,S,DIP,CDIP,SDIP,BWX1,BWX2
      INTEGER*4 NBX1,NBX2
      COMMON/DEFS/XO(MAX_PLN),YO(MAX_PLN),ZO(MAX_PLN),
     &		   C(MAX_PLN),S(MAX_PLN),DIP(MAX_PLN),
     &             CDIP(MAX_PLN),SDIP(MAX_PLN),
     &             BWX1(MAX_PLN),BWX2(MAX_PLN),
     &		   NBX1(MAX_PLN),NBX2(MAX_PLN)

	  real*8, intent(out) :: inlout_elements(NPLANE,6)
			
      IBC=0
      DO 45 NP=1,NPLANE
C----- sub-element half-widths -----
      BX1P=BWX1(NP)*.5
      BX2P=BWX2(NP)*.5

      DO 45 JP=1,NBX1(NP)
      DO 45 IP=1,NBX2(NP)

C----- Find center of this sub-element in global coordinates -----    
C      - center of sub-element IP,JP in local NP coords	 
C                center in local strike direction
	 	 XNP=BX1P*(2*JP-1)
	 	 ZNP=BX2P*(2*IP-1)*SDIP(NP)
C                center in local dip direction
	 	 YNP=BX2P*(2*IP-1)*CDIP(NP)
      	 
C      - center of sub-element IP,JP in global coords (Zo shifted)
	 	 Z=ZNP
	 	 X=S(NP)*XNP-C(NP)*YNP
	 	 Y=C(NP)*XNP+S(NP)*YNP

C      - center of sub-element IP,JP in in-plane coords
         XNP=UG2P(1,1,NP)*X+UG2P(1,2,NP)*Y+UG2P(1,3,NP)*Z
         YNP=UG2P(2,1,NP)*X+UG2P(2,2,NP)*Y+UG2P(2,3,NP)*Z
         ZNP=UG2P(3,1,NP)*X+UG2P(3,2,NP)*Y+UG2P(3,3,NP)*Z

	 	 X=X+XO(NP)
	 	 Y=Y+YO(NP)
	 	 Z=-ZO(NP)+Z
		 
		 inlout_elements(NP,1) = X
		 inlout_elements(NP,2) = Y
		 inlout_elements(NP,3) = Z
		 inlout_elements(NP,4) = XMATRIX(1+IBC,NUM_Ds)
		 inlout_elements(NP,5) = XMATRIX(2+IBC,NUM_Ds)
		 inlout_elements(NP,6) = XMATRIX(3+IBC,NUM_Ds)

 45	     IBC=IBC+3
	  
	  return
	  end

	  subroutine userco(NPLANE,NUM_Ds,xg,yg,zg,
     &              inlout_displ,inlout_stress,inlout_strain,
     &              inlout_orient,inlout_failure,inlout_elements,
     &              inlout_ugrad,
     &              npts)
C****************************************************************************
C  Subroutine to read parameters for output at user coordinates
c        found on the input file
C
C  Note that output parameters that are referenced to the global 
C  coordinate system include:
C		rigid body rotations
C		orientation of principle strain axes
C               failure plane descriptions
C  output parameters that are in the coordinate system of the gridded plane 
C  (X=strike,Y=dip,Z=normal directions) include:
C  		stress tensor elements
C  		strain tensor elements
C  		displacement gradients tensor elements
C		displacements
C  Element relative displacements correspond to the strike,dip,normal directions 
C  of each element's plane.
c c.meertens 24feb93
C***************************************************************************

	  include 'sizes.inc'
	  integer*4 NPLANE

C  Added by A.JANIN
	  integer npts
	  real*8 xg(npts), yg(npts), zg(npts)
	  real*8 inlout_displ(npts,3)
	  real*8 inlout_stress(npts,6)
	  real*8 inlout_strain(npts,6)
	  real*8 inlout_failure(npts,6)
	  real*8 inlout_orient(npts,9)
	  real*8 inlout_ugrad(npts,3,3)
	  real*8 inlout_elements(NPLANE,6)
	  
	  logical o_elem,o_stress,o_strain,o_dgrad,o_rbr,
     &  o_orient,o_disp,o_inv,xyz_out,o_fail,o_proj
	  character*24 sufx
	  common/outputs/sufx,o_elem,o_stress,o_strain,o_dgrad,
     &  o_rbr,o_orient,o_disp,o_inv,xyz_out,o_fail,o_proj

	  real*8 x,y,z

          real*4 stress,ewe,dgrten,pr,ym,friction,tmax,tmaxo,tmaxa
	  common/mike/stress(6),ewe(3),dgrten(3,3),pr,ym,friction,
     &                tmax,tmaxo(2),tmaxa
	  
	  real*4 e,dum
	  common/temps/e(3,3),dum(6,3)
c user coordinates
      real*8 xu,yu,zu
      common/ucoordr/xu(maxco),yu(maxco),zu(maxco)

C         scaling from degrees to radians
	  DATA TORAD/.017453293/

C Added by A.JANIN - global import from XYZ_Results
	  real*4 ex,ey,ez,px,tx,py,ty,pz,tz
	  real*4  volchg,critic,octshr,work
	  real*4  sp1,d1,rake1,sp2,d2,rake2
	  common/JANIN/ex,ey,ez,px,tx,py,ty,pz,tz,
     &             volchg,critic,octshr,work,
     &             sp1,d1,rake1,sp2,d2,rake2



C     Determine what information to write to output files
      call get_fopts

C---- Ensure that no coordinate projections are attempted
      o_proj=.false.

C---- Do calculations on the grid and write the results. -------	      	  
	  do j=1,npts
		x = xg(j)
		y = yg(j)
		z = zg(j)
        
		percent=float(j)*100./float(npts)
		
C		write(*,*)'working.....',percent,' percent done'
  		
		call XYZ_Results(x,y,z,NPLANE,NUM_Ds)
		inlout_displ(j,1) = ewe(1)
		inlout_displ(j,2) = ewe(2)
		inlout_displ(j,3) = ewe(3)
		inlout_stress(j,1) = stress(1)
		inlout_stress(j,2) = stress(2)
		inlout_stress(j,3) = stress(3)
		inlout_stress(j,4) = stress(4)
		inlout_stress(j,5) = stress(5)
		inlout_stress(j,6) = stress(6)
		inlout_strain(j,1) = e(1,1)
		inlout_strain(j,2) = e(1,2)
		inlout_strain(j,3) = e(1,3)
		inlout_strain(j,4) = e(2,2)
		inlout_strain(j,5) = e(2,3)
		inlout_strain(j,6) = e(3,3)
		inlout_ugrad(j,:,:) = dgrten
		inlout_failure(j,:) = [sp1,d1,rake1,sp2,d2,rake2]
		inlout_orient(j,:) = [ex,px,tx,ey,py,ty,ez,pz,tz]
	  end do 
 
      call final_write(NPLANE,NUM_Ds,inlout_elements)

      return
      end
