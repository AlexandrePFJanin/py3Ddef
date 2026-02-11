
# -*- coding: utf-8 -*-
"""
@author: Alexandre JANIN
@aim:    Main
"""

# External dependencies:
import numpy as np

# Internal dependencies
from py3Ddef._all3Ddef import compute3ddef as computeSolution # export only the subroutine compute3ddef
from .geotransform import displacement_sdt_to_xyz


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



class ElementDisplacement:
    """    
    Net displacement across the element ('relative displacement' using the vocabulary of 3D~def)
    i.e., motion of the footwall relative to the hangingwall.
    The footwall lies under the hangingwall, so that a hangingwall moves up and over the footwall
    across a reverse fault, and it slides down and to a level below the footwall across a normal fault.
    If you use the right-handed convention (which 3d~def and py3Ddef do), then the hangingwall lies
    to the right of the strike direction.
    """
    def __init__(self, e):
        self.solution = e
        self.pos = Positions(e[:, 0], e[:, 1], e[:, 2])
        self.ss = e[:, 3]
        self.ds = e[:, 4]
        self.ts = e[:, 5]
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
        output[:, 0] = self.ss
        output[:, 1] = self.ds
        output[:, 2] = self.ts
        return output

    @property
    def shape(self):
        return self.ss.shape

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
        e_xyz = displacement_sdt_to_xyz(self.ss, self.ds, self.ts, strike, dip)
        self.x = e_xyz[:, 0]
        self.y = e_xyz[:, 1]
        self.z = e_xyz[:, 2]

    def reverse(self):
        """
        Reverse the reference fault wall.
        """
        self.ss *= -1
        self.ds *= -1
        self.ts *= -1
        if self.x is not None:
            self.x *= -1
            self.y *= -1
            self.z *= -1
        if self.moving_wall == 'footwall':
            self.moving_wall = 'hangingwall'
        elif self.moving_wall == 'hangingwall':
            self.moving_wall = 'footwall'
    
    def rescale(self, k):
        self.ss *= k
        self.ds *= k
        self.ts *= k
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
        self.ss = self.ss.reshape(shape)
        self.ds = self.ds.reshape(shape)
        self.ts = self.ts.reshape(shape)
        if self.x is not None:
            self.x = self.x.reshape(shape)
            self.y = self.y.reshape(shape)
            self.z = self.z.reshape(shape)

    


# ----------------- MAIN -----------------


def compute3Ddef(x,y,z,
                 xd,yd,zd,length,width,strike,dip,\
                 kode,ss,ds,ts,\
                 nu,E,mu,\
                 bg_flag=None, bg_field=None):
    
    if bg_flag is None or bg_field is None:
        u,s,e,o,f,d,g = computeSolution(x,y,z,\
                             xd,yd,zd,length,width,strike,dip,\
                             kode,ss,ds,ts,\
                             nu,E,mu)
    else:
        u,s,e,o,f,d,g = computeSolution(x,y,z,\
                             xd,yd,zd,length,width,strike,dip,\
                             kode,ss,ds,ts,\
                             nu,E,mu,\
                             bg_flag, bg_field)

    return GridDisplacement(u),\
           GridStress(s),\
           GridStrain(e),\
           GridPrincipalStrainOrientation(o),\
           GridOptimalFailurePlane(f),\
           GridDisplacementGradient(g),\
           ElementDisplacement(d)
