
# -*- coding: utf-8 -*-
"""
@author: Alexandre JANIN
@aim:    Main
"""

# External dependencies:
import numpy as np
import gc
from copy import deepcopy
import pickle

# Internal dependencies
from py3Ddef._all3Ddef import compute3ddef as computeSolution # export only the subroutine compute3ddef
from .generics import im
from .geotransform import displacement_sdt_to_xyz
from .geometry import UniformGrid, PatchCollection
from .viewer import patches2paraview, grid2paraview, plotFault2D, plotFault3D
from .compute import Stress, gradDispl2stress


# ----------------- ELEMENTARY DATA STRUCTURES -----------------


class Positions:

    def __init__(self, x=None, y=None, z=None):
        self.x = x
        self.y = y
        self.z = z
    
    def reshape(self, shape):
        self.x = self.x.reshape(shape)
        self.y = self.y.reshape(shape)
        self.z = self.z.reshape(shape)



class GridScalar:
    def __init__(self, a):
        self.solution = a
        self.v = a
    
    @property
    def shape(self):
        return self.v.shape
    
    @property
    def arrshape(self):
        return self.shape

    @property
    def array(self):
        return self.v

    def reshape(self, shape):
        self.v = self.v.reshape(shape)
    
    def rescale(self, k):
        self.v *= k



class GridVector:
    def __init__(self, a):
        self.solution = a
        self.x = a[:, 0]
        self.y = a[:, 1]
        self.z = a[:, 2]
    
    @property
    def shape(self):
        return self.x.shape
    
    @property
    def arrshape(self):
        return self.shape + (3, )
    
    @property
    def array(self):
        comps = [self.x, self.y, self.z]
        try:
            bshape = np.broadcast(*comps).shape
        except ValueError:
            raise ValueError("Vector components are not broadcastable together")
        # --- allocate tensor
        V = np.empty(bshape + (3,), dtype=self.x.dtype)
        # --- fill tensor
        V[..., 0] = self.x
        V[..., 1] = self.y
        V[..., 2] = self.z
        return V

    @property
    def norm(self):
        """L-2 norm"""
        a = self.array
        return np.sqrt(a[:,:,:,0]**2 + a[:,:,:,1]**2 + a[:,:,:,2]**2)
    
    def reshape(self, shape):
        self.x = self.x.reshape(shape)
        self.y = self.y.reshape(shape)
        self.z = self.z.reshape(shape)
    
    def rescale(self, k):
        self.x *= k
        self.y *= k
        self.z *= k
    


class GridTensor:
    def __init__(self, a):
        self.solution = a
        importFailure = False
        if len(a.shape) == 2:
            if a.shape[1] == 6:
                self.xx = a[:, 0]
                self.xy = a[:, 1]
                self.xz = a[:, 2]
                self.yx = self.xy.copy()
                self.yy = a[:, 3]
                self.yz = a[:, 4]
                self.zx = self.xz.copy()
                self.zy = self.yz.copy()
                self.zz = a[:, 5]
            elif a.shape[1] == 9:
                self.xx = a[:, 0]
                self.xy = a[:, 1]
                self.xz = a[:, 2]
                self.yx = a[:, 3]
                self.yy = a[:, 4]
                self.yz = a[:, 5]
                self.zx = a[:, 6]
                self.zy = a[:, 7]
                self.zz = a[:, 8]
            else:
                importFailure = True
        elif len(a.shape) == 3:
            if a.shape[1] == 3 and a.shape[2] == 3:
                self.xx = a[:, 0, 0]
                self.xy = a[:, 0, 1]
                self.xz = a[:, 0, 2]
                self.yx = a[:, 1, 0]
                self.yy = a[:, 1, 1]
                self.yz = a[:, 1, 2]
                self.zx = a[:, 2, 0]
                self.zy = a[:, 2, 1]
                self.zz = a[:, 2, 2]
            else:
                importFailure = True
        else:
            importFailure = True
        if importFailure:
            raise TypeError('Incorrect shape of the input. Should be either (:,6) or (:,3,3).')
        
    @property
    def shape(self):
        return self.xx.shape
    
    @property
    def arrshape(self):
        return self.shape + (3, 3)
    
    @property
    def array(self):
        comps = [self.xx, self.xy, self.xz, self.yx, self.yy, self.yz, self.zx, self.zy, self.zz]
        try:
            bshape = np.broadcast(*comps).shape
        except ValueError:
            raise ValueError("Tensor components are not broadcastable together")
        # --- allocate tensor
        T = np.empty(bshape + (3, 3), dtype=self.xx.dtype)
        # --- fill tensor
        T[..., 0, 0] = self.xx
        T[..., 0, 1] = self.xy
        T[..., 0, 2] = self.xz
        T[..., 1, 0] = self.yx
        T[..., 1, 1] = self.yy
        T[..., 1, 2] = self.yz
        T[..., 2, 0] = self.zx
        T[..., 2, 1] = self.zy
        T[..., 2, 2] = self.zz
        return T
    
    @property
    def norm(self):
        """Forbenius norm"""
        a = self.array
        return np.sqrt(np.einsum('...ij,...ij->...', a, a))
    
    def reshape(self, shape):
        self.xx = self.xx.reshape(shape)
        self.xy = self.xy.reshape(shape)
        self.xz = self.xz.reshape(shape)
        self.yx = self.yx.reshape(shape)
        self.yy = self.yy.reshape(shape)
        self.yz = self.yz.reshape(shape)
        self.zx = self.zx.reshape(shape)
        self.zy = self.zy.reshape(shape)
        self.zz = self.zz.reshape(shape)
    
    def rescale(self, k):
        self.xx *= k
        self.xy *= k
        self.xz *= k
        self.yx *= k
        self.yy *= k
        self.yz *= k
        self.zx *= k
        self.zy *= k
        self.zz *= k
    


class DislocVectorLocal():
    """
    Data structure for vectorial fields in the local dislocation
    frame (strike, dip, normal)
    """
    def __init__(self, a):
        self.solution = a
        self.s = a[:, 0]   # strike
        self.d = a[:, 1]   # dip
        self.n = a[:, 2]   # normal
    
    @property
    def shape(self):
        return self.n.shape
    
    @property
    def arrshape(self):
        return self.shape + (3, )
    
    @property
    def array(self):
        comps = [self.s, self.d, self.n]
        try:
            bshape = np.broadcast(*comps).shape
        except ValueError:
            raise ValueError("Vector components are not broadcastable together")
        # --- allocate tensor
        V = np.empty(bshape + (3,), dtype=self.n.dtype)
        # --- fill tensor
        V[..., 0] = self.s
        V[..., 1] = self.d
        V[..., 2] = self.n
        return V

    @property
    def norm(self):
        return np.sqrt(self.s**2 + self.d**2 + self.n**2)
    
    def reshape(self, shape):
        self.s = self.s.reshape(shape)
        self.d = self.d.reshape(shape)
        self.n = self.n.reshape(shape)
    
    def rescale(self, k):
        self.s *= k
        self.d *= k
        self.n *= k
    





# ----------------- DERIVED DATA STRUCTURE -----------------


class GridDisplacement(GridVector):
    """
    Displacements, Ux, Uy, Uz, in the x,y,z directions,
    respectively. If output on a plane in in-plane
    coordinates is requested, then these values may be
    output in in-plane coordinates; Ux is in the strike
    direction, Uy in the up-dip direction, Uz normal to the plane.
    """
    def __init__(self, a):
        super().__init__(a)



class GridStress(GridTensor):
    """
    Six independent elements of the stress tensor,
    in the order: Sxx, Sxy, Sxz, Syy, Syz, Szz, the
    maximum shear stress (tmax) on an inspection plane,
    the rake angle of tmax, and the components of tmax
    in the strike and dip directions. If output on a plane
    is requested and global output is not requested, then
    these correspond to the in-plane coordinates.
    (For example, Sss, Ssd, Ssn, Sdd, Sdn, Snn where 's'
    refers to the strike direction, 'd' to the dip direction,
    and 'n' to the normal direction. Thus, Szz is really Snn
    and corresponds to the stress that is normal (tensile) to
    the fault plane, Sxz corresponds to Ssn and is the shear
    stress (sometimes called a traction) in the strike direction.)
    Otherwise the tensor components are referred to the global
    coordinate system. The maximum shear stress is a useful set
    of results since the distribution of tmax can be compared to
    the slip direction. In stress inversion techniques, these two
    directions are assumed to be the same.
    """
    def __init__(self, a):
        super().__init__(a)



class GridStrain(GridTensor):
    """
    Six independent elements of the strain tensor, Exx, Exy, Exz, Eyy, Eyz, Ezz,
    and the volumetric strain (Exx+Eyy+Ezz, independent of the coordinate system).
    If output on a plane in in-plane coordinates is requested, then these values
    may be output in the in-plane coordinates (as described in GridStress).
    """
    def __init__(self, a):
        super().__init__(a)

        

class GridPrincipalStrainOrientation:
    """
    Principal strains, and their plunge and trend
    with respect to the global axes. The values
    corresponding to the maximum principal strain are
    output first, and the minimum principal strain
    values are output last.
    """
    def __init__(self, a):
        self.solution = a
        self.ex = a[:, 0]
        self.px = a[:, 1]
        self.tx = a[:, 2]
        self.ey = a[:, 3]
        self.py = a[:, 4]
        self.ty = a[:, 5]
        self.ez = a[:, 6]
        self.pz = a[:, 7]
        self.tz = a[:, 8]
    
    def reshape(self, shape):
        self.ex = self.ex.reshape(shape)
        self.px = self.px.reshape(shape)
        self.tx = self.tx.reshape(shape)
        self.ey = self.ey.reshape(shape)
        self.py = self.py.reshape(shape)
        self.ty = self.ty.reshape(shape)
        self.ez = self.ez.reshape(shape)
        self.pz = self.pz.reshape(shape)
        self.tz = self.tz.reshape(shape)
    
    @property
    def shape(self):
        return self.ex.shape
    
    @property
    def arrshape(self):
        return self.shape + (3, 3)
    
    @property
    def array(self):
        comps = [self.ex, self.px, self.tx, self.ey, self.py, self.ty, self.ez, self.pz, self.tz]
        try:
            bshape = np.broadcast(*comps).shape
        except ValueError:
            raise ValueError("Tensor components are not broadcastable together")
        # --- allocate tensor
        T = np.empty(bshape + (3, 3), dtype=self.ex.dtype)
        # --- fill tensor
        T[..., 0, 0] = self.ex
        T[..., 0, 1] = self.px
        T[..., 0, 2] = self.tx
        T[..., 1, 0] = self.ey
        T[..., 1, 1] = self.py
        T[..., 1, 2] = self.ty
        T[..., 2, 0] = self.ez
        T[..., 2, 1] = self.pz
        T[..., 2, 2] = self.tz
        return T



class GridOptimalFailurePlane:
    """
    For each of the two optimal failure planes,
    the strike and dip of the plane and the direction
    cosines of the slip vector are given. The strike is
    clockwise from north, the dip is from the horizontal,
    and the direction cosines are with respect to X(E), Y(N), Z(up).
    """
    def __init__(self, a):
        self.solution = a
        self.str1 = a[:, 0]
        self.dip1 = a[:, 1]
        self.rak1 = a[:, 2]
        self.str2 = a[:, 3]
        self.dip2 = a[:, 4]
        self.rak2 = a[:, 5]
    
    def reshape(self, shape):
        self.str1 = self.str1.reshape(shape)
        self.dip1 = self.dip1.reshape(shape)
        self.rak1 = self.rak1.reshape(shape)
        self.str2 = self.str2.reshape(shape)
        self.dip2 = self.dip2.reshape(shape)
        self.rak2 = self.rak2.reshape(shape)

    @property
    def shape(self):
        return self.str1.shape
    
    @property
    def arrshape(self):
        return self.shape + (6)
    
    @property
    def array(self):
        comps = [self.str1, self.dip1, self.rak1, self.str2, self.dip2, self.rak2]
        try:
            bshape = np.broadcast(*comps).shape
        except ValueError:
            raise ValueError("Tensor components are not broadcastable together")
        # --- allocate tensor
        T = np.empty(bshape + (6,), dtype=self.str1.dtype)
        # --- fill tensor
        T[..., 0] = self.str1
        T[..., 1] = self.dip1
        T[..., 2] = self.rak1
        T[..., 3] = self.str2
        T[..., 4] = self.dip2
        T[..., 5] = self.rak2
        return T

    def getFM1(self):
        """
        Get all the 1st focal mechanisms (str1, dip1, rak1) as a GridVector object.
        """
        return GridVector(self.solution[0:3])
    
    def getFM2(self):
        """
        Get all the 2nd focal mechanisms (str2, dip2, rak2) as a GridVector object.
        """
        return GridVector(self.solution[2:])



class GridDisplacementGradient(GridTensor):
    """
    Displacement gradient tensor, dUx/dx, dUx/dy, dUx/dz, dUy/dx, dUy/dy, dUy/dz, dUz/dx, dUz/dy, dUz/dz.
    As for the above, these may correspond to the in-plane coordinates if the user has selected the plane
    and in-plane (not global) options.
    """
    def __init__(self, a):
        super().__init__(a)



class GridStressStrainInvariants:
    """
    Strain and stress invariants including:
     - volume change (volumetric strain): volchg
     - failure stresses (critical failure stress): critic
     - octahedral shear stress: octshr=0.4714*sqrt(sinv1*sinv1+3.0*(sinv2))
     - strain energy density: work=sinv1*sinv1-2.0*(sinv1+pr)*sinv2
    
    Details:

    With e the strain matrix (3x3) then, e(1,1), e(2,2), e(3,3)
    are the diagonal term of e. We define the volume change as:
    volchg=e(1,1)+e(2,2)+e(3,3)

    We define p the list of the eigenvalues of the stress tensor and
    dev1 = max(p(1),p(2),p(3)) and dev3 = min(p(1),p(2),p(3)).
    Then, we define critic, the critical failure stress as:
    critic=(dev1+dev3)/2.

    We calculate stress invariants I and II* from the stress matrix:
 	  sinv1=stress(1)*stress(4)*stress(6)
 	  temp=stress(6)*(stress(4)+stress(1))+stress(1)*stress(4)
 	  tmp=stress(5)*stress(5)+stress(3)*stress(3)+stress(2)*stress(2)
 	  sinv2=-1.0*temp+tmp
    Where stress(1) = sxx, stress(2) = sxy, stress(3) = sxz
          stress(4) = syy, stress(5) = syz, stress(6) = szz
    Then, we define the octahedral shear stress as:
 	octshr=0.4714*sqrt(sinv1*sinv1+3.0*(sinv2))
    
    Then we define the the total strain energy (sum of energy
    due to distortion and volume change). See Ramsay, 1967, p.288. as
 	work=sinv1*sinv1-2.0*(sinv1+pr)*sinv2
    """
    def __init__(self, a):
        self.solution = a
        self.volchg = a[:, 0]
        self.critic = a[:, 1]
        self.octshr = a[:, 2]
        self.work   = a[:, 3]

    
    def reshape(self, shape):
        self.volchg = self.volchg.reshape(shape)
        self.critic = self.critic.reshape(shape)
        self.octshr = self.octshr.reshape(shape)
        self.work   = self.work.reshape(shape)
    
    @property
    def shape(self):
        return self.volchg.shape
    
    @property
    def arrshape(self):
        return self.shape + (4)
    
    @property
    def array(self):
        comps = [self.volchg, self.critic, self.octshr, self.work]
        try:
            bshape = np.broadcast(*comps).shape
        except ValueError:
            raise ValueError("Vector components are not broadcastable together")
        # --- allocate tensor
        T = np.empty(bshape + (4,), dtype=self.volchg.dtype)
        # --- fill tensor
        T[..., 0] = self.volchg
        T[..., 1] = self.critic
        T[..., 2] = self.octshr
        T[..., 3] = self.work
        return T



class ElementStrDispl:
    """    
    Stresses (driving and stored) and net displacement across the element ('relative displacement'
    using the vocabulary of 3D~def) i.e., motion of the footwall relative to the hangingwall.
    The footwall lies under the hangingwall, so that a hangingwall moves up and over the footwall
    across a reverse fault, and it slides down and to a level below the footwall across a normal fault.
    If you use the right-handed convention (which 3d~def and py3Ddef do), then the hangingwall lies
    to the right of the strike direction.
    """
    def __init__(self, e, fstatus=None):
        self.solution = e
        self.pos = Positions(e[:, 0], e[:, 1], e[:, 2])
        self.us = e[:, 3]   # strike-slip
        self.ud = e[:, 4]   # dip-slip
        self.un = e[:, 5]   # tensile (normal)-slip
        self.sdriver = DislocVectorLocal(e[:, 6:9])
        self.sstored = DislocVectorLocal(e[:, 9:])
        self.fstatus = fstatus
        self.nop = e.shape[0]
        self.x  = None
        self.y  = None
        self.z  = None
        self.moving_wall = 'footwall'   # Keep track of which wall is define as moving:
                                        # By defaults, the convention of 3D~def is that the displacement
                                        # across the element is natively given as the motion of the footwall
                                        # relative to the hangingwall.
    
    @property
    def array_sdn(self):
        output = np.zeros((self.nop,3))
        output[:, 0] = self.us
        output[:, 1] = self.ud
        output[:, 2] = self.un
        return output

    @property
    def shape(self):
        return self.us.shape

    @property
    def array_xyz(self):
        if self.x is not None and self.y is not None and self.z is not None:
            output = np.zeros((self.nop,3))
            output[:, 0] = self.x
            output[:, 1] = self.y
            output[:, 2] = self.z
        else:
            raise ValueError('You should first convert the displacement from sdn to xyz using the internal function convert2xyz.')
        return output
    
    @property
    def norm(self):
        return np.sqrt(self.us**2 + self.ud**2 + self.un**2)

    def convert2xyz(self, strike, dip):
        """
        Convert the Element displacement from the local reference frame (along-strike,
        along-dip, normal) to the cartesian reference frame (x,y,z).

        Fills the inner fields: self.x, self.y and self.z

        Args:
            strike, dip (scalar or np.ndarray): lists of strikes and dips of
                every patches. Can be a scalar if all the same or a np.ndarray
                if different.
                Strike in degree clockwise to the north   (same as py3Ddef/3D~def)
                Dip in degree (right) from the horizontal (same as py3Ddef/3D~def)
        """
        e_xyz = displacement_sdt_to_xyz(self.us, self.ud, self.un, strike, dip)
        self.x = e_xyz[:, 0]
        self.y = e_xyz[:, 1]
        self.z = e_xyz[:, 2]

    def reverse(self):
        """
        Reverse the reference fault wall.
        """
        self.us *= -1
        self.ud *= -1
        self.un *= -1
        if self.x is not None:
            self.x *= -1
            self.y *= -1
            self.z *= -1
        if self.moving_wall == 'footwall':
            self.moving_wall = 'hangingwall'
        elif self.moving_wall == 'hangingwall':
            self.moving_wall = 'footwall'
    
    def rescale(self, k):
        self.us *= k
        self.ud *= k
        self.un *= k
        if self.x is not None:
            self.x *= k
            self.y *= k
            self.z *= k

    def get_sdn_vector(self):
        """
        Return the (Strike,Dip,Normal) displacement as a GridVector object.
        """
        return GridVector(self.array_sdn)
    
    def get_xyz_vector(self):
        """
        Return the (x,y,z) displacement as a GridVector object.
        """
        return GridVector(self.array_xyz)

    def reshape(self, shape):
        self.pos.reshape(shape)
        self.us = self.us.reshape(shape)
        self.ud = self.ud.reshape(shape)
        self.un = self.un.reshape(shape)
        if self.x is not None:
            self.x = self.x.reshape(shape)
            self.y = self.y.reshape(shape)
            self.z = self.z.reshape(shape)


class BackgroundDeformation():
    """
    A data structure for the background deformation
    """

    def __init__(self, bg_flag=None, bg_field=None, verbose=True):
        # --- General
        self.pName = 'BackgroundDeform'
        self.verbose = verbose  # controls the verbose level
        # --- Background deformation for 3D~DEF
        self.flag  = None       # flag ('STRE', 'STRA', 'DISP', ...)
        self.field = None       # flat background deformation matrix (shape: (9,))
        self.type  = None       # type of the field: 'stress', 'strain', 'displ'
        self.load_backgroundDef(bg_flag, bg_field)
        self.stress = None
        # ---  Medium parameter used to compute self.stress (if needed)
        self.E  = None      # Young's modulus (it's unit determine the unit of self.stress if converted from displ grad or strain)
        self.nu = None      # Poisson's ratio


    def im(self, textMessage, error=False, warn=False, structure=False, end=False):
        """
        Internal message function.
        """
        im(textMessage, 'DeformationRun', self.verbose, error=error, warn=warn, structure=structure, end=end)


    def load_backgroundDef(self, bg_flag, bg_field):
        """
        Loads a background deformation.
        Fills the fields: self.bg_flag and self.bg_field

        Args:
            bg_flag (str): Flag describing the type of background
                    deformation (Defaults: None).
                    If None or 'NONE', no background deformation.
                    Then:
                        'STRE' or 'stre' -> stress background deformation
                        'STRA' or 'stra' -> strain background deformation
                        'DISP' or 'disp' -> displacement gradient background deformation
            bg_field (numpy.ndarray of shape (9,)): Background field described
                    as a flat 3x3 tensor: (xx, xy, xz, yx, yy, yz, zx, zy, zz)
                    (Defaults: None)
        """
        self.im('Load a background deformation')
        # Check the flag and the field
        if np.logical_or(bg_flag is None, bg_flag == 'NONE') and bg_field is not None:
            raise ValueError('Ambiguous BG field: the flag is None but not the field. Please check.')
        if np.logical_and(bg_flag is not None, bg_flag != 'NONE') and bg_field is None:
            raise ValueError('Ambiguous BG field: the background field is None but not the flag. Please check.')
        # Check the BG flag
        allow_gb_flag = ['STRE', 'stre', 'STRA', 'stra', 'DISP', 'disp', 'NONE']
        if bg_flag is not None and bg_flag not in allow_gb_flag:
            raise ValueError('bg_flag should be in '+str(allow_gb_flag))
        else:
            if bg_flag in ['STRE', 'stre']:
                bg_type = 'stress'
            elif bg_flag in ['STRA', 'stra']:
                bg_type = 'strain'
            elif bg_flag in ['DISP', 'disp']:
                bg_type = 'displ'
            elif bg_flag in ['NONE']:
                bg_type = None
            else:
                raise ValueError('bg_flag should be in '+str(allow_gb_flag))
        # Check the BG field
        if bg_field is not None and not isinstance(bg_field, np.ndarray):
            raise TypeError('The input background field argument bg_field should be a numpy.ndarray of shape (9,)')
        else:
            if isinstance(bg_field, np.ndarray):
                if bg_field.shape != (9,):
                    raise TypeError('The input background field argument bg_field should be a numpy.ndarray of shape (9,)')
                else:
                    pass
        # update
        self.flag  = bg_flag
        self.field = bg_field
        self.type  = bg_type
        self.im('  -> Background deformation loaded')


    def get(self):
        return {
            'bg_flag':  self.bg_flag,
            'bg_field': self.bg_field
        }

    def get_stress(self, E=None, nu=None):
        if self.type == 'stress':
            self.stress = Stress(self.field)
        elif self.type == 'displ':
            if E is None or nu is None:
                raise TypeError('You must enter both E and nu')
            self.stress = gradDispl2stress(self.field,E,nu)




# ----------------- MAIN -----------------


# ------ MAIN FUNCTION

def compute3Ddef(x, y, z, \
                 xd, yd, zd, length, width, strike, dip, \
                 kode, ss, ds, ts, \
                 nu, E, gmu, \
                 fcode=None, sdrop=None, C=None, mu=None, \
                 rhoLitho=None, rhoFluid=None, \
                 output_pstrainOri = False, output_invariants = False, \
                 output_fplanes = False, output_gradDispl = False, \
                 debug = False, maxiter=100, tol=0.1, \
                 bg_flag = None, bg_field = None):
    """
    Compute the 3D deformations both within and at the surface of an elastic
    half-space medium and on each planar dislocation using boundary element framework.

    Args:

        --- Geometry of the stations:

        x, y, z (np.ndarray of shape (N,)): Coordinates of points where the analytical
                solution will be computed [unit: km]
        
        --- Geometry of the dislocations:

        xd, yd, zd (np.ndarray of shape (M,)): Coordinates of the reference corner of
                each dislocation [unit: km]
        length, width, strike, dip (np.ndarray of shape (M,)): Geometry of each planar
                dislocation.
                Units:  length: km
                        width: km
                        strike: deg clockwise relative to North
                        dip: deg from horizontal
        
        --- Boundary condition on the dislocations:

        kode (numpy.ndarray, dype=int, shape (M,)): description of the boundary conditions for each
                dislocation. See the documentation of the package.
        ss, dt, ts (numpy.ndarray of shape (M,)): Components of the boundary conditions for each
                dislocations. See the documentation of the package
        
        --- Friction

        fcode (None, numpy.ndarray, dype=int, shape (M,), optional, default: None):
                Convention: A frictional patch slides (relative displacement) only
                            when  |tau| > mu * |sigma_n| + C
                            We write the excess of stress '\Delta\sigma_e' the stress 'above'
                            the Mohr-Coulomb envelope:
                            \Delta\sigma_e = |tau| - (mu * |sigma_n| + C)
                            and the stress at the envelope tau_c = mu * |sigma_n| + C
                            Below, when we say 'release': stress used to slip ('converted' in 
                            relative displacement on the dislocation).
                Friction code for each dislocation, controls the effect of the input argument 'sdrop':
                    fcode = 0: No friction, the patch is considered frictionless
                               The other frictional parameters (i.e. sdrop, C, mu, etc) have no effect.
                    fcode = 1: If failure, releases a fraction 'sdrop' of the excess of stress \Delta\sigma_e
                               (e.g. sdrop=0: no release, sdrop=1: back on the failure envelope)
                    fcode = 2: If failure, releases a fraction 'sdrop' of the total stress |tau|
                               (e.g. sdrop=0.5: release half of the total stress (may still be above
                               the failure envelope), sdrop=1: release all the stress: no more stress
                               on the dislocation (completely used to slip)).
                    fcode = 3: If failure, releases all the excess stress \Delta\sigma_e plus a fraction 'sdrop'
                               of the stress envelope.
                               Release: \Delta\sigma_e + sdrop * tau_c
                    fcode = 4: If failure, releases all the excess stress \Delta\sigma_e plus an absolute stress 'sdrop'
                               (unit of sdrop: same as E).
                               Release: \Delta\sigma_e + sdrop
                    fcode = 5: If failure, releases an absolute stress 'sdrop' (unit of sdrop: same as E).
                               Release: sdrop
                If None: fcode = 0 for each patch: all frictionless.
        sdrop (None, numpy.ndarray, dype=float, shape (M,), optional, default: None):
                Stress drop parameter for each patch
                When fcode in [1, 2, 3] -> no unit, is a fraction
                When fcode in [4, 4]    -> same unit as E
                -> Can be None only is fcode is None
        rhoLitho (None, numpy.ndarray, dype=float, shape (M,), optional, default: None):
                Average rock density in the vicinity of the each dislocations (unit: need to be consistent with E
                if E in [Pa] (prefered), then rhoLitho is in [kg.m^{-3}]).
                Warning: vertical variation of rhoLitho will be ignored.
                -> Can be None only is fcode is None
        rhoFluid (None, numpy.ndarray, dype=float, shape (M,), optional, default: None):
                Average fluid density for each dislocations (unit: need to be consistent with E
                if E in [Pa] (prefered), then rhoLitho is in [kg.m^{-3}]).
                -> Can be None only is fcode is None
        C (None, numpy.ndarray, dype=float, shape (M,), optional, default: 0):
                Cohesion of each patch
                -> Can be None only is fcode is None
        mu (None, numpy.ndarray, dype=float, shape (M,)):
                Local friction of each patch.
                -> Can be None only is fcode is None
        
        --- Parameters of the Medium:
        
        nu (scalar): Poisson's ratio [units: None]
        E (scalar): Young's modulus [units: will define the unit of
                the stress output: e.g. is E = 1 (GPa) then, the stress
                output will be in GPa.]
        gmu (scalar): (global) coefficient of internal friction [units: None]
        
        --- Flag for optional outputs:
        
        output_pstrainOri (bool, optional): Option controlling the calculation
                of the principale strain orientations (optional output).
                (Defaults: False)
        output_invariants (bool, optional): Option controlling the calculation
                of the invariants of the stress and strain fields (optional output).
                (Defaults: False)
        output_fplanes (bool, optional): Option controlling the calculation
                of the optimal failure planes (optional output).
                (Defaults: False)
        output_gradDispl (bool, optional): Option controlling the calculation
                of the displacement gradient field (optional output).
                (Defaults: False)
        debug (bool, optional): If set to True then, export in ASCII text file
                the coefficient matrices build by the program at diff steps:
                 - imatrix: coeff and b.c.s before the inversion (b.c.s. in stress,
                            without patches with fixed relative displ)
                 - omatrix: coeff and b.c.s after the inversion (b.c.s. in displacement,
                            without patches with fixed relative displ)
                 - fmatrix: final coeff and b.c.s (corrected to add patches with fixed
                            relative displacement): matrix used to compute the analytical
                            solutions on the grid.
                NOTE 1: The coefficient (Okada coeff) never change, they are computed once.
                NOTE 2: These 3 files will appear in the local directory where the function
                        is called.
        
        --- Iterative solver parameter (when at least one frictional element is set)

        maxiter (int, optional): maximum number of iteration for the iterative solver.
                (Defaults: maxiter=100)
        tol (float, optional): tolerance (in centimeters) to that define the convergence
                criterion for the iterative solver. The abs difference between the (i-1)th
                and (i)th displacement vectors on the dislocations have to be less that
                the tolerance 'tol'.
                (Defaults: tol=0.1)

        --- Background deformation: If not set, runs without background deformation.

        bg_flag (None, str, optional): Flag describing the type of background
            deformation (Defaults: None).
            If None or 'NONE', no background deformation.
            Then:
                'STRE' or 'stre' -> stress background deformation
                'STRA' or 'stra' -> strain background deformation
                'DISP' or 'disp' -> displacement gradient background deformation
        bg_field (None, numpy.ndarray of shape (9,)): Background field described
                as a flat 3x3 tensor: (xx, xy, xz, yx, yy, yz, zx, zy, zz)
                (Defaults: None)
    """
    # ==== Default frictional parameters:
    ndis = xd.shape[0]
    # --- sdrop
    if sdrop is None:
        if fcode is not None and np.any(fcode > 0):
            raise ValueError("You must enter a value for 'sdrop' as at least one element is frictional")
        else:
            sdrop = np.zeros(ndis)
    # --- cohesion
    if C is None:
        if fcode is not None and np.any(fcode > 0):
            raise ValueError("You must enter a value for 'C' as at least one element is frictional")
        else:
            C = np.zeros(ndis)
    else:
        if np.any(C < 0):
            raise ValueError('Cohesion is positive or null')
    # --- frition coeff
    if mu is None:
        if fcode is not None and np.any(fcode > 0):
            raise ValueError("You must enter a value for 'mu' as at least one element is frictional")
        else:
            mu = np.zeros(ndis)
    else:
        if np.any(mu < 0):
            raise ValueError('The friction coefficient on each element must be >= 0')
    # --- fcode
    if fcode is None:
        fcode = np.zeros(ndis)
    else:
        allowed_fcode = [0, 1, 2, 3, 4, 5]
        invalid = ~np.isin(fcode, allowed_fcode)
        if np.any(invalid):
            raise ValueError('Unknown value(s) for fcode')
    # --- rhoLitho
    if rhoLitho is None:
        if fcode is not None and np.any(fcode > 0):
            raise ValueError("You must enter a value for 'rhoLitho' as at least one element is frictional")
        else:
            rhoLitho = np.zeros(ndis)
    else:
        if np.any(rhoLitho < 0):
            raise ValueError('Rock density must be positive or null')
    # --- rhoFluid
    if rhoFluid is None:
        if fcode is not None and np.any(fcode > 0):
            raise ValueError("You must enter a value for 'rhoFluid' as at least one element is frictional")
        else:
            rhoFluid = np.zeros(ndis)
    else:
        if np.any(rhoFluid < 0):
            raise ValueError('Fluid density must be positive or null')
    # --- consistency of rhoLitho and rhoFluid
    if np.any(rhoFluid > rhoLitho):
        raise ValueError('The density of the fluid should always be lower or equal to the density of rocks')
    # --- Check iterative solver parameter
    if maxiter < 1:
        raise ValueError("The iterative solver parameter 'maxiter' must be at least 1.")
    
    # ==== Compute
    # Run the main function according to the value of the input BG field
    if bg_flag is None or bg_field is None:
        u,s,e,o,f,j,g,d,fstatus,runParam = computeSolution(x, y, z,\
                             xd, yd, zd, length, width, strike, dip,\
                             kode, ss, ds, ts,\
                             fcode, sdrop, rhoLitho, rhoFluid, C, mu, \
                             nu, E, gmu, \
                             output_pstrainOri, output_invariants,\
                             output_fplanes, output_gradDispl,\
                             debug, maxiter, tol)
    else:
        u,s,e,o,f,j,g,d,fstatus,runParam = computeSolution(x, y, z,\
                             xd, yd, zd, length, width, strike, dip,\
                             kode, ss, ds, ts,\
                             fcode, sdrop, rhoLitho, rhoFluid, C, mu,\
                             nu, E, gmu,\
                             output_pstrainOri, output_invariants,\
                             output_fplanes, output_gradDispl,\
                             debug, maxiter, tol,\
                             bg_flag, bg_field)

    # Prepare and format the outputs
    displ     = GridDisplacement(u)
    stress    = GridStress(s)
    strain    = GridStrain(e)
    dislocs   = ElementStrDispl(d, fstatus=fstatus)
    if output_pstrainOri:
        pstrainOri = GridPrincipalStrainOrientation(o)
    else:
        pstrainOri = None
        del o
    if output_fplanes:
        fplanes    = GridOptimalFailurePlane(f)
    else:
        fplanes    = None
        del f
    if output_invariants:
        invariants  = GridStressStrainInvariants(j)
    else:
        invariants  = None
        del j
    if output_gradDispl:
        gradDispl   = GridDisplacementGradient(g)
    else:
        gradDispl   = None
        del g
    
    # Garbage collector call
    gc.collect()

    return displ, stress, strain, pstrainOri, fplanes, invariants, gradDispl, dislocs, fstatus, runParam



# ------ MAIN CLASS

class DeformationRun:

    def __init__(self, grid=None, patches=None, nu=None, mu=None, E=None, bg=None, \
                 output_pstrainOri=False, output_invariants=False,\
                 output_fplanes=False, output_gradDispl=False,\
                 debug=False, maxiter=100, tol=0.5, auto_reshape=True, verbose=True):
        """
        Data structure for py3Ddef calculation.

        Description of fields:

            .verbose (bool): Option controlling the verbose output of the class
            .auto_reshape (bool): Reshape automatically the outputs computed on the
                    grid to have the same shape. (Default: True)

            ===== INPUT FIELDS =====

            --- Geometry of the problem

            .grid (py3Ddef.geometry.UniformGrid):
                    Grid object describing where the solution will be analytically computed [unit: km].
            .patches (py3Ddef.geometry.UniformGrid):
                    Set of dislocations in the elastic medium [unit: km].

            --- Medium Parameters:
        
            .nu (scalar): Poisson's ratio [units: None]
            .E (scalar): Young's modulus [units: will define the unit of
                    the stress output: e.g. is E = 1 (GPa) then, the stress
                    output will be in GPa.]
            .mu (scalar): coefficient of internal friction [units: None]
            
            -- Output flags:

            .output_pstrainOri (bool, optional): Option controlling the calculation
                    of the principale strain orientations (optional output).
                    (Defaults: False)
            .output_invariants (bool, optional): Option controlling the calculation
                    of the invariants of the stress and strain fields (optional output).
                    (Defaults: False)
            .output_fplanes (bool, optional): Option controlling the calculation
                    of the optimal failure planes (optional output).
                    (Defaults: False)
            .output_gradDispl (bool, optional): Option controlling the calculation
                    of the displacement gradient field (optional output).
                    (Defaults: False)
            .debug (bool, optional): If set to True then, export in ASCII text file
                    the coefficient matrices build by the program at diff steps:
                    - imatrix: coeff and b.c.s before the inversion (b.c.s. in stress,
                                without patches with fixed relative displ)
                    - omatrix: coeff and b.c.s after the inversion (b.c.s. in displacement,
                                without patches with fixed relative displ)
                    - fmatrix: final coeff and b.c.s (corrected to add patches with fixed
                                relative displacement): matrix used to compute the analytical
                                solutions on the grid.
                    NOTE 1: The coefficient (Okada coeff) never change, they are computed once.
                    NOTE 2: These 3 files will appear in the local directory where the function
                            .compute3Ddef() is called.
            
            -- Iterative solver parameter (when at least one frictional element is set)

            .maxiter (int, optional): maximum number of iteration for the iterative solver.
                    (Defaults: maxiter=100)
            .tol (float, optional): tolerance (in centimeters) to that define the convergence
                    criterion for the iterative solver. The abs difference between the (i-1)th
                    and (i)th displacement vectors on the dislocations have to be less that
                    the tolerance 'tol'.
                    (Defaults: tol=0.1)
                
            -- Background deformation:

            .bg (None or py3Ddef.BackgroundDeformation, optional):
                    Backbround deformation (Defaults: None).
            
            ===== OUTPUT FIELDS =====

            NOTE: For more details, see the documentation of the package of the
                  help function of the different object.

                  e.g. >> import py3Ddef
                       >> help(py3Ddef.GridDisplacement)

            --- Default outputs:

            .displ (py3Ddef.GridDisplacement): Displacement computed on the input grid [unit: cm]
            .stress (py3Ddef.GridStress): Stress field computed on the input grid [unit: same as .E]
            .strain (py3Ddef.GridStrain): Strain field computed on the input grid [unit: None]

            .dislocs (py3Ddef.ElementStrDispl): Slip and stresses on each dislocation [unit: cm]

            --- Optional outputs:

            .pstrainOri (py3Ddef.GridPrincipalStrainOrientation): Principal strain orientations on the grid
                    Will be computed by the function .compute3Ddef() only if output_pstrainOri is True. [unit: deg]
            .fplanes (py3Ddef.GridOptimalFailurePlane): Optimal failure planes on the grid
                    Will be computed by the function .compute3Ddef() only if output_fplanes is True. [unit: deg]
            .invariants (py3Ddef.GridStressStrainInvariants): Stress/Strain invariants on the grid
                    Will be computed by the function .compute3Ddef() only if output_invariants is True.
            .gradDispl (py3Ddef.GridDisplacementGradient): Displacement gradient field on the grid
                    Will be computed by the function .compute3Ddef() only if output_gradDispl is True.
        """
        self.verbose = verbose
        self.auto_reshape = auto_reshape
        self.im('Initialization of a DeformationRun instance')
        # Geometry
        self.grid    = None
        self.patches = None
        if grid is not None:
            self.load_grid(grid)
        else:
            self.grid = None
        if patches is not None:
            self.load_patches(patches)
        else:
            self.patches = None
        # Medium Parameters
        self.nu = nu
        self.E  = E
        self.mu = mu
        # Background deformation
        if bg is not None:
            self.load_backgroundDef(bg)
        else:
            self.bg = None
        # Output Flags
        self.output_pstrainOri = output_pstrainOri
        self.output_invariants = output_invariants
        self.output_fplanes    = output_fplanes
        self.output_gradDispl  = output_gradDispl
        self.debug = debug
        # Iterative solver
        self.maxiter = maxiter
        self.tol     = tol
        # Solution
        self.displ      = None
        self.stress     = None
        self.strain     = None
        self.pstrainOri = None
        self.fplanes    = None
        self.invariants = None
        self.gradDispl  = None
        self.dislocs    = None
        # run parameters
        self.runParam = None
        self.is_converged = None
        self.niter = None
        self.im('Instance initialized')
    

    def im(self, textMessage, error=False, warn=False, structure=False, end=False):
        """
        Internal message function.
        """
        im(textMessage, 'DeformationRun', self.verbose, error=error, warn=warn, structure=structure, end=end)


    def load_grid(self, grid):
        """
        Loads a computationnal grid object.
        Fills the field: self.grid

        Args:
            grid (py3Ddef.geometry.UniformGrid)
        """
        self.im('Load a grid object')
        if isinstance(grid, UniformGrid):
            self.grid     = grid
            self.gridtype = grid.type
        elif grid is None:
            self.grid     = None
            self.gridtype = None
        else:
            raise TypeError('The input grid should be either None or py3Ddef.geometry.UniformGrid')
        self.im('  -> Grid loaded')
    

    def load_patches(self, patches):
        """
        Loads a dislocation collection object.
        Fills the field: self.patches

        Args:
            patches (py3Ddef.geometry.PatchCollection)
        """
        self.im('Load a patch collection object')
        if isinstance(patches, PatchCollection):
            check_consistency_friction(patches)
            self.patches  = patches
        elif patches is None:
            self.patches  = patches
        else:
            raise TypeError('The input patch collection should be an instance of py3Ddef.geometry.PatchCollection')
        self.im('  -> Patch collection loaded')
    

    def load_backgroundDef(self,bg):
        """
        Loads a background deformation.

        Args:
            bg (None, py3Ddef.BackgroundDeformation): input background deformation
        """
        if bg is None or isinstance(bg, BackgroundDeformation):
            self.bg = bg
            if self.E is not None and self.nu is not None:
                self.bg.get_stress(self.E, self.nu)
        else:
            raise TypeError('Unknown type for the input background deformation. Must be either None or an instance of py3Ddef.BackgroundDeformation')


    def saveRun(self, fname, path='./'):
        """
        Save the current class instance in a pickle file.
        """
        if path[-1] != '/':
            path += '/'
        with open(path+fname, 'wb') as file:
            pickle.dump(self, file)

    
    @classmethod
    def loadRun(cls, fname, path='./'):
        """
        Initialized a class instance form a pickle file:
        
        Examples:
            solution1 = DeformationRun()
            [...]
            solution1.saveRun('./myRun.pickle')
            solution2 = DeformationRun.loadRun('./myRun.pickle')
        """
        if path[-1] != '/':
            path += '/'
        with open(path + fname, 'rb') as file:
            return pickle.load(file)


    @property
    def u(self):
        """Alias for the displacement field (.displ)"""
        return self.displ
    
    @property
    def s(self):
        """Alias for the stress field (.stress)"""
        return self.stress

    @property
    def e(self):
        """Alias for the strain field (.srrain)"""
        return self.strain

    @property
    def o(self):
        """Alias for the orientation of the principal strain components field (.pstrainOri)"""
        return self.pstrainOri
    
    @property
    def g(self):
        """Alias for the displacement gradient field (.gradDispl)"""
        return self.gradDispl
    
    @property
    def f(self):
        """Alias for the optimal failure plane field (.fplanes)"""
        return self.fplanes
    
    @property
    def j(self):
        """Alias for the stress/strain invariants (.invariants)"""
        return self.invariants
    
    @property
    def d(self):
        """Alias for the dislocation slip (.dislocs)"""
        return self.dislocs


    def compute3Ddef(self, grid=None, patches=None, nu=None, E=None, mu=None,\
                     output_pstrainOri=None, output_invariants=None, \
                     output_fplanes=None, output_gradDispl=None, \
                     debug=None, maxiter=None, tol=None, bg=None, auto_reshape=None):
        """
        Compute the 3D deformations both within and at the surface of an elastic
        half-space medium and on each planar dislocation using boundary element framework.

        Solves the stresses, strains, and displacements from a set of initial conditions
        (a population of dislocation, a set of points in the surrounding elastic medium). 
        
        Load the solution in the current class instance, filling the input (if set) 
        and output fields. See class description.

        Args:
            --- Others
                auto_reshape (bool): Reshape automatically the outputs computed on the
                        grid to have the same shape. (Default: True)
            
            --- Geometry of the problem: If not, take the corresponding fields in the current instance of DeformationRun.

                grid (py3Ddef.geometry.UniformGrid): grid describing where the solution will be
                        analytically computed in the medium (i.e. stations).
                        [units: km]
                patches (py3Ddef.geometry.PatchCollection): Object describing a set of
                        dislocations (discontinuities) in the elastic medium.
                        [units: km]
            
            --- Parameter of the medium: If not, take the corresponding fields in the current instance of DeformationRun.

                nu (scalar): Poisson's ratio [units: None]
                E (scalar): Young's modulus [units: will define the unit of
                        the stress output: e.g. is E = 1 (GPa) then, the stress
                        output will be in GPa.]
                mu (scalar): coefficient of internal friction [units: None]
            
            --- Output flags: If not, take the corresponding fields in the current instance of DeformationRun.

                output_pstrainOri (bool, optional): Option controlling the calculation
                        of the principale strain orientations (optional output).
                        (Defaults: False)
                output_invariants (bool, optional): Option controlling the calculation
                        of the invariants of the stress and strain fields (optional output).
                        (Defaults: False)
                output_fplanes (bool, optional): Option controlling the calculation
                        of the optimal failure planes (optional output).
                        (Defaults: False)
                output_gradDispl (bool, optional): Option controlling the calculation
                        of the displacement gradient field (optional output).
                        (Defaults: False)
                debug (bool, optional): If set to True then, export in ASCII text file
                        the coefficient matrices build by the program at diff steps:
                        - imatrix: coeff and b.c.s before the inversion (b.c.s. in stress,
                                    without patches with fixed relative displ)
                        - omatrix: coeff and b.c.s after the inversion (b.c.s. in displacement,
                                    without patches with fixed relative displ)
                        - fmatrix: final coeff and b.c.s (corrected to add patches with fixed
                                    relative displacement): matrix used to compute the analytical
                                    solutions on the grid.
                        NOTE 1: The coefficient (Okada coeff) never change, they are computed once.
                        NOTE 2: These 3 files will appear in the local directory where the function
                                is called.
            
            --- Iterative solver parameter (when at least one frictional element is set)

                maxiter (int, optional): maximum number of iteration for the iterative solver.
                        (Defaults: maxiter=100)
                tol (float, optional): tolerance (in centimeters) to that define the convergence
                        criterion for the iterative solver. The abs difference between the (i-1)th
                        and (i)th displacement vectors on the dislocations have to be less that
                        the tolerance 'tol'.
                        (Defaults: tol=0.1)
            
            --- Background deformation: If not, take the corresponding fields in the current instance of DeformationRun.

                bg (None or py3Ddef.BackgroundDeformation, optional):
                        Backbround deformation (Defaults: None).

        Returns:
            None
        """

        self.im('Compute the 3D deformation')

        # --- Geometry
        if grid is None and self.grid is None:
            raise TypeError('You must enter a grid object in input')        
        if patches is None and self.patches is None:
            raise TypeError('You must set a patch collection object in input')
        
        # Update the grid
        if grid is not None:
            self.load_grid(grid)
    
        # Update the dislocation collection
        if patches is not None:
            self.load_patches(patches)
        
        # --- Medium parameters

        # Poisson's ratio
        if nu is None and self.nu is None:
            raise ValueError("You should enter a Poisson's ratio nu")
        elif nu is None and self.nu is not None:
            pass
        else:
            self.nu = nu
        
        # Young's modulus
        if E is None and self.E is None:
            raise ValueError("You should enter a Young's modulus E")
        elif E is None and self.E is not None:
            pass
        else:
            self.E = E

        # coeff of internal friction
        if mu is None and self.mu is None:
            raise ValueError("You should enter a coefficient of internal friction mu")
        elif mu is None and self.mu is not None:
            pass
        else:
            self.mu = mu
        
        # --- Output flags

        # output_pstrainOri
        if output_pstrainOri is None and self.output_pstrainOri is None:
            self.output_pstrainOri = False # Default
        elif output_pstrainOri is None and self.output_pstrainOri is not None:
            pass
        else:
            self.output_pstrainOri = output_pstrainOri
        
        # output_invariants
        if output_invariants is None and self.output_invariants is None:
            self.output_invariants = False # Default
        elif output_invariants is None and self.output_invariants is not None:
            pass
        else:
            self.output_invariants = output_invariants
        
        # output_fplanes
        if output_fplanes is None and self.output_fplanes is None:
            self.output_fplanes = False # Default
        elif output_fplanes is None and self.output_fplanes is not None:
            pass
        else:
            self.output_fplanes = output_fplanes
        
        # output_gradDispl
        if output_gradDispl is None and self.output_gradDispl is None:
            self.output_gradDispl = False # Default
        elif output_gradDispl is None and self.output_gradDispl is not None:
            pass
        else:
            self.output_gradDispl = output_gradDispl
        
        # debug
        if debug is None and self.debug is None:
            self.debug = False # Default
        elif debug is None and self.debug is not None:
            pass
        else:
            self.debug = debug

        # --- Iterative solver
        if maxiter is not None:
            self.maxiter = maxiter
        if tol is not None:
            self.tol = tol
        
        # --- Background conditions
        if bg is not None:
            self.load_backgroundDef(bg)
        
        if self.bg is None:
            bg_flag  = None
            bg_field = None
        else:
            bg_flag  = self.bg.flag
            bg_field = self.bg.field

        # --- Call the function to compute the deformation
        self.displ, self.stress, self.strain, self.pstrainOri, \
        self.fplanes, self.invariants, self.gradDispl, self.dislocs, \
        self.fstatus, runParam = \
            compute3Ddef(*self.grid.get(), *self.patches.get(), \
                         self.nu, self.E, self.mu, \
                         **self.patches.getFriction(), \
                         output_pstrainOri = self.output_pstrainOri, \
                         output_invariants = self.output_invariants, \
                         output_fplanes = self.output_fplanes, \
                         output_gradDispl = self.output_gradDispl, \
                         debug = self.debug, maxiter=self.maxiter, tol=self.tol, \
                         bg_flag = bg_flag, \
                         bg_field = bg_field)
        # auto reshape
        if auto_reshape is None and self.auto_reshape is None:
            self.auto_reshape = False # Default
        elif auto_reshape is None and self.auto_reshape is not None:
            pass
        else:
            self.auto_reshape = auto_reshape
        if self.auto_reshape:
            self.reshapeGSolutions()
        # run parameters
        self.runParam = runParam
        if runParam[0] == 0:
            self.is_converged = False
        elif runParam[0] == 1:
            self.is_converged = True
        else:
            raise ValueError('Invalid value %s for runParam[0]')
        self.niter = runParam[1]


    def plotFault2D(self, **kwargs):
        """
        Fast 2D plot of fault patches projected onto a plane
        orthogonal to a given 3D vector.

        See the documentation of the function py3Ddef.viewer.plotFault2D
        """
        plotFault2D(self.patches, **kwargs)
    


    def plotFault3D(self, **kwargs):
        """
        Interactive display of a set of dislocation in 3D.

        See the documentation of the function py3Ddef.viewer.plotFault3D
        """
        plotFault3D(self.patches, **kwargs)



    def reshapeGSolutions(self, shape=None):
        """
        Reshape the solutions computed on the computational grid.
        The function uses the shape provided in input or the shape
        if the loaded computational grid is shape=None.

        Args:
            shape (None, tuple): new shape for all the fields computed 
                    on the computational grid. If None then, take the
                    shape of the grid.
                    (Defaults: None)
        """
        self.im('Reshape the solutions computed on the grid')
        if shape is None:
            shape = self.grid.shape
            self.im('  -> Use the shape of the grid')
        self.im('  Shape used: '+str(shape))
        # reshape all
        self.displ.reshape(shape)
        self.stress.reshape(shape)
        self.strain.reshape(shape)
        if self.pstrainOri is not None:
            self.pstrainOri.reshape(shape)
        if self.fplanes is not None:
            self.fplanes.reshape(shape)
        if self.invariants is not None:
            self.invariants.reshape(shape)
        if self.gradDispl is not None:
            self.gradDispl.reshape(shape)
        self.im('  -> Reshape done.')



    def patches2paraview(self, fname, path='./', fields=None, fieldnames=None, displacement=True, verbose=True):
        """
        Export the current patch collection structure and the linked field(s) in
        XDMF and HDF5 files for a display in Paraview.

        About the XDMF/HDF5 files:
            The function adds automatically three vector fields at patch centers
            in order to keep track of the patch-relative reference frame:
                - dir_strike
                - dir_dip
                - dir_normal

            Structure of the mesh:
                Grid 1 : FaultPatches (quad mesh)
                    - geometry
                    - user scalar/tensor fields (cell data)

                Grid 2 : PatchCenters (point cloud)
                    - centers of patches
                    - dir_strike, dir_dip, dir_normal vectors (node data)

            In Paraview:
                - Open the XDMF file with the "XDMF reader" option (not the "Xdmf3 Reader S" or "Xdmf3 Reader T")
                - Separate "FaultPatches" and "PatchCenters" with "Extract Block" (avoid issues with "Glyph" representation)


        Args:
            fname (str): Name of the output file (without extension).
            path (str, optional): Path to the export directory (defaults, path='./')
            fields (None, list): list of additional fields you want to attach to
                    to the grid and display in Paraview (default, fields=None).
                    If not None, fields should be a list. Each element of the list
                    can have one of the following accepted type:
                        - py3Ddef.ElementStrDispl
                    No shape constrains for the py3Ddef objects as long as they are
                    compatible with the shape of the grid (i.e. can be reshaped).
            fieldnames (None, list): list of names for the input argument 'fields'.
                    They should have the same length and the same type (default None).
                    If the name of a given field is not given (e.g. fieldnames[i]=None),
                    then, an automatic name will be given.
            displacement (bool): If true, auto export the displacement field (self.displ)
                    (Defaults: True)
            verbose (bool) Option controlling the verbose output of the function.
                    (Defaults: True)
        
        Returns:
            None
        """
        self.im('Exportation of the data and structure of the PatchCollection instance to XDMF+HDF5 files for Paraview')
        myfields = []
        myfnames = []
        # --- Check
        if path[-1] != '/':
            path += '/'
        # --- Proc the default fields
        if displacement:
            self.im('  -> Export the displacement field')
            myfields, myfnames = attach2dislocations(self.patches, myfields, myfnames, self.dislocs, name=None, verbose=verbose)
        # --- Additional fields
        if fields is None:
            pass
        else:
            if isinstance(fields, list) and isinstance(fieldnames, list) and len(fields) == len(fieldnames):
                for i in range(len(fields)):
                    self.im('  -> Export '+str(fieldnames[i]))
                    myfields, myfnames = attach2dislocations(self.patches, myfields, myfnames, fields[i], name=fieldnames[i], verbose=verbose)
            else:
                raise TypeError("Additional fields (arguments 'fields' and 'fieldnames') should be entered as lists, even for just one additional field.")
        # Check
        if len(myfields) == 0:
            self.im('   -> Empty field: Export point IDs')
            v = np.arange(self.nop).reshape(self.shape)
            myfields, myfnames = attach2dislocations(self.patches, myfields, myfnames, v, name='pointID', verbose=verbose)
        # Call the function
        patches2paraview(fname, self.patches, myfields, fieldnames=myfnames, path=path)
        self.im('Dislocation: Structure and fields exported in:')
        self.im('   - path:  '+str(path))
        self.im('   - files: '+str(fname)+ '  (.xdmf/.h5)')
    


    def grid2paraview(self, fname, path='./', fields=None, fieldnames=None,
                      displacement=True, stress=True, strain=True, 
                      pstrainOri=None, fplanes=None, invariants=None, gradDispl=None, 
                      verbose=True):
        """
        Export the current grid and the linked field(s) in XDMF and HDF5 files
        for a display in Paraview.

        Args:
            fname (str): Name of the output file (without extension).
            path (str, optional): Path to the export directory (defaults, path='./')
            fields (None, list): list of additional fields you want to attach to
                    to the grid and display in Paraview (default, fields=None).
                    If not None, fields should be a list. Each element of the list
                    can have one of the following accepted type:
                        - py3Ddef.GridDisplacement,
                        - py3Ddef.GridStress,
                        - py3Ddef.GridStrain,
                        - py3Ddef.GridDisplacementGradient,
                        - py3Ddef.GridPrincipalStrainOrientation,
                        - py3Ddef.GridOptimalFailurePlane,
                        - py3Ddef.GridStressStrainInvariants,
                        - py3Ddef.GridScalar,
                        - py3Ddef.GridVector,
                        - py3Ddef.GridTensor,
                        - numpy.ndarray
                    No shape constrains for the py3Ddef objects as long as they are
                    compatible with the shape of the grid (i.e. can be reshaped).
                    However numpy.ndarray should have a specific shape according
                    to the content of the field:
                        (N,)    -> scalar field
                        (N,3)   -> vector field
                        (N,6)   -> symmetric tensor
                        (N,9)   -> full tensor
                        (N,k)   -> generic attribute
            fieldnames (None, list): list of names for the input argument 'fields'.
                    They should have the same length and the same type (default None).
                    If the name of a given field is not given (e.g. fieldnames[i]=None),
                    then, an automatic name will be given.
            displacement (bool): If true, auto export the displacement field (self.displ)
                    (Defaults: True)
            stress (bool);  If true, auto export the stress field (self.stress)
                    (Defaults: True)
            strain (bool);  If true, auto export the strain field (self.strain)
                    (Defaults: True)
            pstrainOri (None bool);  If true, auto export the principle strain
                    orientation field (self.pstrainOri). If None, use the value of
                    the flag self.output_pstrainOri. (Defaults: None)
            fplanes (None bool);  If true, auto export the optimal failure
                    planes field (self.fplanes). If None, use the value of
                    the flag self.output_fplanes. (Defaults: None)
            invariants (None bool);  If true, auto export the stress/strain
                    invariants field (self.invariants). If None, use the value of
                    the flag self.output_invariants. (Defaults: None)
            gradDispl (None bool);  If true, auto export the displacement
                    gradient field (self.gradDispl). If None, use the value of
                    the flag self.output_gradDispl. (Defaults: None)
            verbose (bool) Option controlling the verbose output of the function.
                    (Defaults: True)
        
        Returns:
            None
        """
        self.im('Exportation of the grid instance to XDMF+HDF5 files for Paraview')
        myfields = []
        myfnames = []
        # --- Check
        if path[-1] != '/':
            path += '/'
        # --- Proc the default fields
        if displacement:
            self.im('  -> Export the displacement field')
            myfields, myfnames = attach2grid(self.grid, myfields, myfnames, self.displ, name=None, verbose=verbose)
        if stress:
            self.im('  -> Export the stress field')
            myfields, myfnames = attach2grid(self.grid, myfields, myfnames, self.stress, name=None, verbose=verbose)
        if strain:
            self.im('  -> Export the strain field')
            myfields, myfnames = attach2grid(self.grid, myfields, myfnames, self.strain, name=None, verbose=verbose)
        # --- Proc de optional fields
        if pstrainOri is None:
            pstrainOri = self.output_pstrainOri
        if fplanes is None:
            fplanes = self.output_fplanes
        if invariants is None:
            invariants = self.output_invariants
        if gradDispl is None:
            gradDispl = self.output_gradDispl
        if pstrainOri:
            self.im('  -> Export the principal strain orientation field')
            myfields, myfnames = attach2grid(self.grid, myfields, myfnames, self.pstrainOri, name=None, verbose=verbose)
        if fplanes:
            self.im('  -> Export the optimal failure plane field')
            myfields, myfnames = attach2grid(self.grid, myfields, myfnames, self.fplanes, name=None, verbose=verbose)
        if invariants:
            self.im('  -> Export the stress/strain invariant field')
            myfields, myfnames = attach2grid(self.grid, myfields, myfnames, self.invariants, name=None, verbose=verbose)
        if gradDispl:
            self.im('  -> Export the displacement gradient field')
            myfields, myfnames = attach2grid(self.grid, myfields, myfnames, self.gradDispl, name=None, verbose=verbose)
        # --- Additional fields
        if fields is None:
            pass
        else:
            if isinstance(fields, list) and isinstance(fieldnames, list) and len(fields) == len(fieldnames):
                for i in range(len(fields)):
                    self.im('  -> Export '+str(fieldnames[i]))
                    myfields, myfnames = attach2grid(self.grid, myfields, myfnames, fields[i], name=fieldnames[i], verbose=verbose)
            else:
                raise TypeError("Additional fields (arguments 'fields' and 'fieldnames') should be entered as lists, even for just one additional field.")
        # Check
        if len(myfields) == 0:
            raise ValueError('Empty fields. You should export at least one field on the grid')
        # Call the function
        grid2paraview(fname, self.grid, fields=myfields, fieldnames=myfnames, path=path, precision='Float')
        self.im('Grid: Structure and fields exported in:')
        self.im('   - path:  '+str(path))
        self.im('   - files: '+str(fname)+ '  (.xdmf/.h5)')

        



# ------ Auxilary functions


def attach2dislocations(patches, v, vnames, a, name=None, verbose=True):
    """
    Links a new field (a, name) to the PatchCollection and add
    it to the current list of linked fields (v, vnames).
    
    -> Prepare a field for the function patches2paraview.

    Args:
        patches (py3Ddef.geometry.PatchCollection): reference patch collection object
        v (list): reference list of field where you want to add another one
        vnames (list): list of reference field names
        a (py3Ddef.ElementStrDispl: new input field that will be added to 'v'.
                The shape have to be similar to the shape of the input patch collection (i.e. flat).
        name (None, str): name of the field that will appear in Paraview.
                If None (default), an automatic name will be given.
        verbose (bool): Option controlling the verbose output of the function.
                (Defaults: True)
    
    Returns:
        v (list): updated input argument 'v' with the new field (added as a numpy.ndarray)
        vnames (list): updated input argument 'vnames' with the new field name
    """
    # Check input
    if not isinstance(patches, PatchCollection):
        raise TypeError('The input patches should be an instance of py3Ddef.geometry.PatchCollection.')
    # Import
    if isinstance(a, ElementStrDispl):
        wall = a.moving_wall
        if name is None:
            name = wall+' displacement'
        im("Attach a displacement field (%s) to a PatchCollection instance"%name, 'Attach2Disloc', verbose=verbose)
        if a.shape == patches.shape:
            if a.x is None or a.y is None or a.z is None:
                a.convert2xyz(patches.strike, patches.dip)    # generate the x,y,z displacements
            else:
                # already done
                pass
            # take only the xyz array
            a = a.array_xyz
        else:
            raise TypeError('Error in the shape of the input ElementStrDispl.')
    else:
        im("Attach an additional field '%s' to a PatchCollection instance"%name, 'Attach2Disloc', verbose=verbose)
        raise TypeError('Unkown type for the input.')
    # Append:
    v.append(a)
    vnames.append(name)
    im("  -> Field '%s' attached"%name, 'Attach2Disloc', verbose=verbose)
    return v, vnames


def detach2dislocations(v, vnames, i):
    """
    Remove the field number i.
    
    Args:
        i (int): index of the field to remove
    """
    v.pop(i)
    vnames.pop(i)


def attach2grid(grid, v, vnames, a, name=None, verbose=True):
    """
    Links a new field (a, name) to a given grid object and add
    it to the current list of linked fields (v, vnames).
    
    -> Prepare a field for the function grid2paraview.

    Args:
        grid (py3Ddef.geometry.UniformGrid): reference grid object
        v (list): reference list of field where you want to add another one
        vnames (list): list of reference field names
        a (py3Ddef.GridDisplacement,
                py3Ddef.GridStress,
                py3Ddef.GridStrain,
                py3Ddef.GridDisplacementGradient,
                py3Ddef.GridPrincipalStrainOrientation,
                py3Ddef.GridOptimalFailurePlane,
                py3Ddef.GridStressStrainInvariants,
                py3Ddef.GridScalar,
                py3Ddef.GridVector,
                py3Ddef.GridTensor,
                numpy.ndarray): new input field that will be added to 'v'.
                No shape constrains for the py3Ddef objects as long as they are
                compatible with the shape of the grid (i.e. can be reshaped).
                However numpy.ndarray should have a specific shape according
                to the content of the field:
                    (N,)    -> scalar field
                    (N,3)   -> vector field
                    (N,6)   -> symmetric tensor
                    (N,9)   -> full tensor
                    (N,k)   -> generic attribute
        name (None, str): name of the field that will appear in Paraview.
                If None (default), an automatic name will be given.
        verbose (bool): Option controlling the verbose output of the function.
                (Defaults: True)
    
    Returns:
        v (list): updated input argument 'v' with the new field (added as a numpy.ndarray)
        vnames (list): updated input argument 'vnames' with the new field name
    """
    # Check input
    if not isinstance(grid, UniformGrid):
        raise TypeError('The input grid should be an instance of py3Ddef.geometry.UniformGrid.')
    # Import
    im("Attach a field to the grid instance", 'Attach2Grid', verbose=verbose)
    a = deepcopy(a)     # make sure to not modify the original field
    # ====== Check the type and the shape of input field
    im("  - attach a %s field of shape %s"%(str(type(a)), str(a.shape)), 'Attach2Grid', verbose=verbose)
    # --- py3Ddef high level objects
    if isinstance(a, GridDisplacement):
        im("  - detected: GridDisplacement", 'Attach2Grid', verbose=verbose)
        if a.shape != grid.shape:
            a.reshape(grid.shape)   # reshape
        a = a.array                 # convert to array
        if name is None:
            name = 'displacement'
    elif isinstance(a, GridStress) or isinstance(a, GridStrain) or isinstance(a, GridDisplacementGradient):
        if name is None:
            if isinstance(a, GridStress):
                im("  - detected: GridStress", 'Attach2Grid', verbose=verbose)
                name = 'stress'
            elif isinstance(a, GridStrain):
                im("  - detected: GridStrain", 'Attach2Grid', verbose=verbose)
                name = 'strain'
            elif isinstance(a, GridDisplacementGradient):
                im("  - detected: GridDisplacementGradient", 'Attach2Grid', verbose=verbose)
                name = 'displGradient'
        if a.shape != grid.shape:
            a.reshape(grid.shape)   # reshape
        a = a.array                 # convert to array
    elif isinstance(a, GridPrincipalStrainOrientation):
        im("  - detected: GridPrincipalStrainOrientation", 'Attach2Grid', verbose=verbose)
        if a.shape != grid.shape:
            a.reshape(grid.shape)   # reshape
        a = a.array                 # convert to array
        if name is None:
            name = 'prStrainOrient'
    elif isinstance(a, GridOptimalFailurePlane):
        im("  - detected: GridOptimalFailurePlane", 'Attach2Grid', verbose=verbose)
        if a.shape != grid.shape:
            a.reshape(grid.shape)   # reshape
        a = a.array                 # convert to array
        if name is None:
            name = 'optiFailPlane'
    elif isinstance(a, GridStressStrainInvariants):
        im("  - detected: GridStressStrainInvariants", 'Attach2Grid', verbose=verbose)
        if a.shape != grid.shape:
            a.reshape(grid.shape)   # reshape
        a = a.array                 # convert to array
        if name is None:
            name = 'invariants'
    # --- py3Ddef medium level objects
    elif isinstance(a, GridScalar):
        im("  - detected: GridScalar", 'Attach2Grid', verbose=verbose)
        if a.shape != grid.shape:
            a.reshape(grid.shape)   # reshape
        a = a.array                 # convert to array
        if name is None:
            name = add_nameExtension('scalar', vnames)
    elif isinstance(a, GridVector):
        im("  - detected: GridVector", 'Attach2Grid', verbose=verbose)
        if a.shape != grid.shape:
            a.reshape(grid.shape)   # reshape
        a = a.array                 # convert to array
        if name is None:
            name = add_nameExtension('vectorial', vnames)
    elif isinstance(a, GridTensor):
        im("  - detected: GridTensor", 'Attach2Grid', verbose=verbose)
        if a.shape != grid.shape:
            a.reshape(grid.shape)   # reshape
        a = a.array                 # convert to array
        if name is None:
            name = add_nameExtension('tensor', vnames)
    # --- numpy.ndarray
    elif isinstance(a, np.ndarray):
        original_name = name
        if len(a.shape) == 1:   # probably a scalar field
            im("  - suspected: GridScalar (N)", 'Attach2Grid', verbose=verbose)
            a = GridScalar(a)
            a.reshape(grid.shape)   # reshape
            a = a.array             # convert to array
            if name is None:
                name = add_nameExtension('scalar', vnames)
        elif len(a.shape) == 2:
            if a.shape[1] == 3: # probably a vectorial field
                im("  - suspected: GridVector (N,3)", 'Attach2Grid', verbose=verbose)
                a = GridVector(a)
                a.reshape(grid.shape)   # reshape
                a = a.array             # convert to array
                if name is None:
                    name = add_nameExtension('vectorial', vnames)
            elif a.shape[1] == 6: # probably a symmetric tensor
                im("  - suspected: GridTensor (N,6) symmetric", 'Attach2Grid', verbose=verbose)
                a = GridTensor(a)
                a.reshape(grid.shape)   # reshape
                a = a.array             # convert to array
                if name is None:
                    name = add_nameExtension('tensor', vnames)
            else:
                raise TypeError('Error in the shape of the input numpy.ndarray, incompatible with the grid.')
        elif len(a.shape) == 3:
            if a.shape == grid.shape: # probably scalar
                im("  - suspected: GridScalar (nx,ny,nz)", 'Attach2Grid', verbose=verbose)
                a = GridScalar(a)
                a.reshape(grid.shape)   # reshape
                a = a.array             # convert to array
                if name is None:
                    name = add_nameExtension('scalar', vnames)
            elif a.shape[1] == 3 and a.shape[2] == 3:
                im("  - suspected: GridTensor (N,3,3)", 'Attach2Grid', verbose=verbose)
                a = GridTensor(a)
                a.reshape(grid.shape)   # reshape
                a = a.array             # convert to array
                if name is None:
                    name = add_nameExtension('tensor', vnames)
            else:
                raise TypeError('Error in the shape of the input numpy.ndarray, incompatible with the grid.')
        else:
            raise TypeError('Error in the shape of the input numpy.ndarray, incompatible with the grid.')
        # add a name extension
        if original_name is None:
            name = add_nameExtension(name, vnames)
    else:
        raise TypeError('Unkown type for the input.')
    im("  - convert to numpy.ndarray of shape %s"%str(a.shape), 'Attach2Grid', verbose=verbose)
    im("  - name: %s"%name, 'Attach2Grid', verbose=verbose)
    # ====== Extract arrays and transpose shape for Paraview convention
    # transpose to have the same convention as Paraview (nz, ny, nx, ncomponent)
    im("  - transpose shape for Paraview (nz, ny, nx, ncomponent)", 'Attach2Grid', verbose=verbose)
    if len(a.shape) == 5: # first, if tensor, flatten it
        a = a.reshape(grid.shape + (9,))
    if len(a.shape) == 3: # scalar
        a = np.transpose(a, (2, 1, 0))
    elif len(a.shape) == 4: # vector or tensor
        a = np.transpose(a, (2, 1, 0, 3))
    else:
        raise ValueError('The input field %s does not have a shape compatible with the input grid'%name)
    # ====== Append:
    v.append(a)
    vnames.append(name)
    im("  -> Field '%s' attached"%name, 'Attach2Grid', verbose=verbose)
    return v, vnames


def detach2grid(v, vnames, i):
    """
    Remove the field number i.
    
    Args:
        i (int): index of the field to remove
    """
    v.pop(i)
    vnames.pop(i)


def add_nameExtension(name,reflist):
    """
    Add an extension to an input name to avoid
    duplicate in a reference list:

    e.g. name = 'displ'
         reflist = ['stress', 'displ']
         new_names = add_nameExtension(name,reflist)
         print(new_names)

         >> 'displ_1'
    """
    runCond = True
    ii = 1
    while runCond:
        nname = name+'_%s'%str(ii)
        if nname not in reflist:
            name = nname
            runCond = False
        else:
            ii += 1
    return name


def check_consistency_friction(patches):
    """
    A function that check the consistency between boundary conditions
    (kode) and friction on each patch of a PatchCollection object.

    Args:
        patches (py3Ddef.geometry.PatchCollection):
                checked patch collection.
    """
    for i in range(patches.nop):
        if patches.fcode[i] != 0 and patches.kode[i] != 2:
            raise ValueError('The frictional patch (id=%s, fcode=%s) must have 3 stress boundary condition (kode must be 2, but kode=%s)'%(str(patches.ids[i]), str(patches.fcode[i]), str(patches.kode[i])))