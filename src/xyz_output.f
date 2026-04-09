	  SUBROUTINE GRID(NUM_Ds,NPLANE,npts)
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

      use global_arrays

	  INTEGER*4 NUM_Ds,NPLANE,npts

      if (debug) then
          ! -> Export the final corrected XMATRIX
          !    (will be use to compute the solution on the grid
          !    and returned as element displacement vector)
          call write_xmatrix('zmatrix')
      endif

	  call userco(NPLANE,NUM_Ds,npts)

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

      use global_outputs

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
	  common/GOUTPUTS/ex,ey,ez,px,tx,py,ty,pz,tz,
     &             volchg,critic,octshr,work,
     &             sp1,d1,rake1,sp2,d2,rake2
	  
	  
	  DATA TODEG/57.295779/
	  zd=z
	  
C     Gets stress and displacements
	  call CALCPO(x,y,zd,NPLANE,NUM_Ds,ierr)
 
C     Gets strain tensor and vol. change
	  call strain(e,delV)

C     Gets principal strains and orientations and write to output file
	  if (o_orient) then
		  call third_results(ex,ey,ez,px,tx,py,ty,pz,tz)
	  endif
	  
C     Rotate to in-plane coordinates: REMOVED: A.JANIN 14.02.2026
C     if(o_proj) CALL TOPLANE(SSTRIKE,CSTRIKE,SDYP,CDYP)
C
C     Code added 2/15/94 to get max shear stress and its orientation on the inspection 
C     plane. See subroutine TOPLANE in 3dmain.f and manual for more info.
C
	  tmaxa=tmaxa*TODEG

C     Get Stress/Strain invairants
	  if (o_invariant) then
		  call invariants(volchg,critic,octshr,work)
      endif

C     Get Optimal failure planes
	  if (o_failure) then
      	call failure_planes(sp1,d1,rake1,sp2,d2,rake2)
	  endif

      return
      end


	  
	subroutine XYZ_elems(NPLANE,NUM_Ds)
C *******************************************************************
C Subroutine to write out element displacements (relative) in
C X,Y,Z,Data format.  Each line of the output file will contain the
C coordinates of the center of the element in global and in-plane
C coordinate systems,and the relative
C displacements in the strike, dip, and tensile (normal)
C directions.
C *******************************************************************

	! Function reworded A.JANIN 14.02.2026 to be more F90

	use global_arrays
	use global_outputs

	implicit none

	integer, intent(in) :: NPLANE, NUM_Ds

	integer :: NP, JP, IP
	integer :: IBC

	real(8) :: BX1P, BX2P
	real(8) :: XNP, YNP, ZNP
	real(8) :: X, Y, Z

	IBC = 0

	do NP = 1, NPLANE

		!----- sub-element half-widths -----
		BX1P = BWX1(NP)*0.5d0
		BX2P = BWX2(NP)*0.5d0

		do JP = 1, NBX1(NP)
			do IP = 1, NBX2(NP)

				! center in local strike direction
				XNP = BX1P*(2*JP-1)
				ZNP = BX2P*(2*IP-1)*SDIP(NP)

				! center in local dip direction
				YNP = BX2P*(2*IP-1)*CDIP(NP)

				! center in global coords
				Z = ZNP
				X = S(NP)*XNP - C(NP)*YNP
				Y = C(NP)*XNP + S(NP)*YNP

				! center in in-plane coords
				XNP = UG2P(1,1,NP)*X + UG2P(1,2,NP)*Y + UG2P(1,3,NP)*Z
				YNP = UG2P(2,1,NP)*X + UG2P(2,2,NP)*Y + UG2P(2,3,NP)*Z
				ZNP = UG2P(3,1,NP)*X + UG2P(3,2,NP)*Y + UG2P(3,3,NP)*Z

				X = X + XO(NP)
				Y = Y + YO(NP)
				Z = -ZO(NP) + Z

				goutarray_elements(NP,1) = X
				goutarray_elements(NP,2) = Y
				goutarray_elements(NP,3) = Z
				goutarray_elements(NP,4) = XMATRIX(1+IBC,NUM_Ds)
				goutarray_elements(NP,5) = XMATRIX(2+IBC,NUM_Ds)
				goutarray_elements(NP,6) = XMATRIX(3+IBC,NUM_Ds)

				goutarray_elements(NP,7) = SDRIVER(NP,1)
				goutarray_elements(NP,8) = SDRIVER(NP,2)
				goutarray_elements(NP,9) = SDRIVER(NP,3)

				goutarray_elements(NP,10) = SSTORED(NP,1)
				goutarray_elements(NP,11) = SSTORED(NP,2)
				goutarray_elements(NP,12) = SSTORED(NP,3)

				IBC = IBC + 3

			end do
		end do
	end do

	end subroutine XYZ_elems

	  subroutine userco(NPLANE,NUM_Ds,npts)
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

	  use global_arrays
	  use global_outputs

	  integer*4 NPLANE,NUM_Ds,npts

	  real*8 x,y,z

      real*4 stress,ewe,dgrten,pr,ym,friction,tmax,tmaxo,tmaxa
	  common/mike/stress(6),ewe(3),dgrten(3,3),pr,ym,friction,
     &                tmax,tmaxo(2),tmaxa
	  
	  real*4 e,dum
	  common/temps/e(3,3),dum(6,3)

C         scaling from degrees to radians
	  DATA TORAD/.017453293/

C Added by A.JANIN - global import from XYZ_Results
	  real*4 ex,ey,ez,px,tx,py,ty,pz,tz
	  real*4  volchg,critic,octshr,work
	  real*4  sp1,d1,rake1,sp2,d2,rake2
	  common/GOUTPUTS/ex,ey,ez,px,tx,py,ty,pz,tz,
     &             volchg,critic,octshr,work,
     &             sp1,d1,rake1,sp2,d2,rake2

C---- Do calculations on the grid and write the results. -------	      	  
	  do j=1,npts
		x = xfgrid(j)
		y = yfgrid(j)
		z = zfgrid(j)
        		  		
		call XYZ_Results(x,y,z,NPLANE,NUM_Ds)
		goutarray_displ(j,1) = ewe(1)
		goutarray_displ(j,2) = ewe(2)
		goutarray_displ(j,3) = ewe(3)
		goutarray_stress(j,1) = stress(1)
		goutarray_stress(j,2) = stress(2)
		goutarray_stress(j,3) = stress(3)
		goutarray_stress(j,4) = stress(4)
		goutarray_stress(j,5) = stress(5)
		goutarray_stress(j,6) = stress(6)
		goutarray_strain(j,1) = e(1,1)
		goutarray_strain(j,2) = e(1,2)
		goutarray_strain(j,3) = e(1,3)
		goutarray_strain(j,4) = e(2,2)
		goutarray_strain(j,5) = e(2,3)
		goutarray_strain(j,6) = e(3,3)

		if (o_ugrad) then
			goutarray_ugrad(j,:,:) = dgrten
		endif
		if (o_failure) then
			goutarray_failure(j,:) = [sp1,d1,rake1,sp2,d2,rake2]
		endif
		if (o_orient) then
			goutarray_orient(j,:) = [ex,px,tx,ey,py,ty,ez,pz,tz]
		endif
		if (o_invariant) then
			goutarray_invariant(j,:) = [volchg,critic,octshr,work]
		endif
	end do
 
C  simplified by A.JANIN 14.02.2026
	  call XYZ_elems(NPLANE,NUM_Ds)

      return
      end
