# -*- coding: utf-8 -*-
"""
@author: Alexandre JANIN
@aim:    Module containing misc computational routines
"""

# External dependencies:
import numpy as np


# ----------------- CLASS -----------------

class Stress:

    def __init__(self, array=None):
        self.xx = None
        self.xy = None
        self.xz = None
        self.yx = None
        self.yy = None
        self.yz = None
        self.zx = None
        self.zy = None
        self.zz = None

        if array is not None:
            self.load(array)

    def load(self, array):
        if not isinstance(array, np.ndarray):
            raise TypeError("The input stress array must be a numpy.ndarray")

        array = np.asarray(array)

        # Case 1: full 3x3 tensor
        if array.shape == (3, 3):
            self.xx = array[0, 0]
            self.xy = array[0, 1]
            self.xz = array[0, 2]
            self.yx = array[1, 0]
            self.yy = array[1, 1]
            self.yz = array[1, 2]
            self.zx = array[2, 0]
            self.zy = array[2, 1]
            self.zz = array[2, 2]

        # Case 2: flat full tensor (9,)
        elif array.shape == (9,):
            mat = array.reshape((3, 3))
            self.load(mat)

        # Case 3: flat symmetric tensor (6,)
        elif array.shape == (6,):
            sxx, sxy, sxz, syy, syz, szz = array

            self.xx = sxx
            self.xy = sxy
            self.xz = sxz

            self.yx = sxy
            self.yy = syy
            self.yz = syz

            self.zx = sxz
            self.zy = syz
            self.zz = szz

        else:
            raise ValueError(
                "Stress array must have shape (3,3), (9,) or (6,) "
                f"but got shape {array.shape}"
            )
        
    @property
    def sym6(self):
        return np.array([self.xx, self.xy, self.xz, self.yy, self.yz, self.zz], dtype=float)
    
    @property
    def flat(self):
        return np.array([self.xx, self.xy, self.xz, self.yx, self.yy, self.yz, self.zx, self.zy, self.zz], dtype=float)
    
    @property
    def mat(self):
        return np.array([[self.xx, self.xy, self.xz],
                [self.yx, self.yy, self.yz],
                [self.zx, self.zy, self.zz]], dtype=float)

    def get_inplane(self, strike, dip):
        return BSTR2P(self.sym6, strike, dip)
        





# ----------------- FUNCTIONS -----------------


def MK_TRANS(STK,DIP):
    '''
    Create transformation matrices from global to local in-plane coordinates for a given plane.
	  
    STK = strike of the plane (in deg)
	DIP = dip of the plane (in deg)
	UG2P,SG2P = transformation matrices
                UG2P(3,3) transforms displacements from global coords to obtain displacements in in-plane coords.
                SG2P(3,6) transforms stresses from global coords to obtain stresses in in-plane coords.
    '''
    
    SSTK = np.sin(STK*np.pi/180)
    CSTK = np.cos(STK*np.pi/180)
    SDIP = np.sin(DIP*np.pi/180)
    CDIP = np.cos(DIP*np.pi/180)
    UG2P = np.zeros((3,3))
    SG2P = np.zeros((3,6))

#   in-plane component  global component
#   -----------------   -------------------
#   strike direction    strike direction 
    UG2P[0,0] = SSTK
#   strike direction    dip direction
    UG2P[0,1] = CSTK
#   strike direction    normal direction
    UG2P[0,2] = 0.
#   dip direction       strike direction
    UG2P[1,0] = -CSTK*CDIP
#   dip direction       dip direction
    UG2P[1,1] = SSTK*CDIP
#   dip direction       normal direction
    UG2P[1,2] = SDIP
#   normal direction    strike direction
    UG2P[2,0] = CSTK*SDIP
#   normal direction    dip direction	
    UG2P[2,1] = -SSTK*SDIP
#   normal direction    normal direction
    UG2P[2,2] = CDIP

#   in-plane component  global component
#   -----------------   -------------------
#                        Sxz        Sxx
    SG2P[0,0] = UG2P[2,0]*UG2P[0,0]
#                        Sxz        Sxy	
    SG2P[0,1] = UG2P[2,0]*UG2P[0,1]+UG2P[2,1]*UG2P[0,0]
#                        Sxz        Sxz
    SG2P[0,2] = UG2P[2,2]*UG2P[0,0]
#                        Sxz        Syy
    SG2P[0,3] = UG2P[2,1]*UG2P[0,1]
#                        Sxz        Syz
    SG2P[0,4] = UG2P[2,2]*UG2P[0,1]
#                        Sxz        Szz
    SG2P[0,5] = 0.

#                        Syz        Sxx
    SG2P[1,0] = UG2P[1,0]*UG2P[2,0]
#                        Syz        Sxy
    SG2P[1,1] = UG2P[1,0]*UG2P[2,1]+UG2P[1,1]*UG2P[2,0]
#                        Syz        Sxz	
    SG2P[1,2] = UG2P[1,0]*UG2P[2,2]+UG2P[1,2]*UG2P[2,0]
#                        Syz        Syy
    SG2P[1,3] = UG2P[1,1]*UG2P[2,1]
#                        Syz        Syz
    SG2P[1,4] = UG2P[1,1]*UG2P[2,2]+UG2P[1,2]*UG2P[2,1]
#                        Syz        Szz
    SG2P[1,5] = UG2P[1,2]*UG2P[2,2]

#                        Szz        Sxx
    SG2P[2,0] = UG2P[2,0]*UG2P[2,0]
#                        Szz        Sxy
    SG2P[2,1] = 2.*UG2P[2,0]*UG2P[2,1]
#                        Szz        Sxz
    SG2P[2,2] = 2.*UG2P[2,0]*UG2P[2,2]
#                        Szz        Syy
    SG2P[2,3] = UG2P[2,1]*UG2P[2,1]
#                        Szz        Syz
    SG2P[2,4] = 2.*UG2P[2,1]*UG2P[2,2]
#                        Szz        Szz
    SG2P[2,5] = UG2P[2,2]*UG2P[2,2]

    return UG2P, SG2P


def BSTR2P(BSTRESS,STK,DIP):
    '''
    Project the background stress in the in-plane coordinate system
    '''
    UG2P, SG2P = MK_TRANS(STK,DIP)
    BSTRESS    = BSTRESS
    BSTR       = np.zeros(3)
    for i in range(3):
        for k in range(6):
            BSTR[i] += SG2P[i,k]*BSTRESS[k]
    return BSTR


def gradDispl2stress(gradDispl, E, nu):
    """
    Convert a displacement gradient tensor in stress
    using E the Young's modulus.

    Args:
        gradDispl (numpy.ndarray of shape (9,)):
                displacement gradiant as a flat matrix:
                (dUx/dx, dUx/dy, dUx/dz, dUy/dx, dUy/dy, dUy/dz, dUz/dx, dUz/dy, dUz/dz)
        E (scalar):
                Young's modulus: It's unit will define the unit of the output stress (e.g. Pa)
        nu (scalar):
                Poisson ratio

    Returns:
        stress (numpy.ndarray of shape (9,)):
                equivalent stress (same unit as E)
                (Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz)
                NOTE: will be symmetrical
    """
    error = False
    if not isinstance(gradDispl, np.ndarray):
        error = True
    else:
        if not gradDispl.shape == (9,):
            error = True
    if error:
        raise TypeError('The input displacement gradient must be a numpy.ndarray of shape (9,)')
    else:
        stress = np.zeros(9)
    # --- Material constants
    # rigidity
    XMU=E*.5/(1.+nu)
    # rigidity*2
    XMU2=XMU*2.
    # Lame parameter (lambda)
    DMULT=nu/(1.-2.*nu)
    # --- Convert to stresses - see Jaeger and Cook, pg. 110
    DILAT=(gradDispl[0]+gradDispl[4]+gradDispl[8])*DMULT
    stress[0]=XMU2*(DILAT+gradDispl[0])
    stress[1]=XMU*(gradDispl[1]+gradDispl[3])
    stress[2]=XMU*(gradDispl[2]+gradDispl[6])
    stress[3]=XMU*(gradDispl[1]+gradDispl[3])
    stress[4]=XMU2*(DILAT+gradDispl[4]) 
    stress[5]=XMU*(gradDispl[5]+gradDispl[7])
    stress[6]=XMU*(gradDispl[2]+gradDispl[6])
    stress[7]=XMU*(gradDispl[5]+gradDispl[7])
    stress[8]=XMU2*(DILAT+gradDispl[8])
    return Stress(stress)



def strain2stress(strain, E, nu):
    """
    Convert a background strain tensor in stress
    using E the Young's modulus.

    Args:
        strain (numpy.ndarray of shape (9,)):
                strain tensor as a flat matrix:
                (Exx, Exy, Exz, Eyx, Eyy, Eyz, Ezx, Ezy, Ezz)
                NOTE: The non-symmetrical part is ignored.
        E (scalar):
                Young's modulus: It's unit will define the unit of the output stress (e.g. Pa)
        nu (scalar):
                Poisson ratio

    Returns:
        stress (numpy.ndarray of shape (9,)):
                equivalent stress (same unit as E)
                (Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz)
                NOTE: will be symmetrical
    """
    error = False
    if not isinstance(strain, np.ndarray):
        error = True
    else:
        if not strain.shape == (9,):
            error = True
    if error:
        raise TypeError('The input strain tensor must be a numpy.ndarray of shape (9,)')
    else:
        stress = np.zeros(9)
     # --- Material constants
    # rigidity
    XMU=E*.5/(1.+nu)
    # rigidity*2
    XMU2=XMU*2.
    # Lame parameter (lambda)
    DMULT=nu/(1.-2.*nu)
    # --- compute
    DILAT=(strain[0]+strain[3]+strain[5])*DMULT

    stress[0]=XMU2*(DILAT+strain[0])    # sxx
    stress[1]=XMU2*strain[1]            # sxy
    stress[2]=XMU2*strain[2]            # sxz
    stress[3]=XMU2*strain[1]            # sxy
    stress[4]=XMU2*(DILAT+strain[3])    # syy
    stress[5]=XMU2*strain[4]            # syz
    stress[6]=XMU2*strain[2]            # sxz
    stress[7]=XMU2*strain[4]            # syz
    stress[8]=XMU2*(DILAT+strain[5])    # szz
    return Stress(stress)


def gradDispl2rotation(gradDispl):
    """
    Convert a displacement gradient tensor into
    rotation vector.

    Args:
        gradDispl (numpy.ndarray of shape (9,)):
                displacement gradiant as a flat matrix:
                (dUx/dx, dUx/dy, dUx/dz, dUy/dx, dUy/dy, dUy/dz, dUz/dx, dUz/dy, dUz/dz)

    Returns:
        rot (numpy.ndarray of shape (3,)):
                equivalent rotation (rotx, roty, rotz)
    """
    error = False
    if not isinstance(gradDispl, np.ndarray):
        error = True
    else:
        if not gradDispl.shape != (9,):
            error = True
    if error:
        raise TypeError('The input displacement gradient must be a numpy.ndarray of shape (9,)')
    # --- Convert to stresses - see Jaeger and Cook, pg. 110
    # ( dUy/dx - dUx/dy)/2
    ROT_Z=.5*(gradDispl[3]-gradDispl[1])
    # ( dUx/dz - dUz/dx)/2
    ROT_Y=.5*(gradDispl[2]-gradDispl[6])
    #( dUz/dy - dUy/dz)/2
    ROT_X=.5*(gradDispl[7]-gradDispl[5])
    return np.array([ROT_X, ROT_Y, ROT_Z])

