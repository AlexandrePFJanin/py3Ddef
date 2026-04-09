# -*- coding: utf-8 -*-
"""
@author: Alexandre JANIN
@aim:    Geometrical data structures and functions
"""

# External dependencies:
import numpy as np

# Internal dependencies
from .generics import im
from .viewer import plotFault3D


# ----------------- FUNCTIONS -----------------


def discreteDislocation(x0, y0, z0, L, W, dip, strike, n_strike, n_dip,\
                        kode=10, ss=0, ds=0, ts=0, group=None, \
                        fcode=0, sdrop=0, rhoLitho=0, rhoFluid=0, C=0, mu=0.6,\
                        verbose=True):
    """
    Discretize a dipping elastic dislocation plane into subpatches with arbitrary strike.
    Uses the same convention as 3D~def.

    Args:
        
        --- Geometry:
        
        x0, y0, z0 (float):
            Coordinates of the top-left corner of the fault (km)
            x0, initial coordinate in the East direction
            y0, initial coordinate in the North direction
            z0, initial coordinate in the elevation (up) direction (<0: depth)
        L (float):
            Fault length along strike (km)
        W (float):
            Fault width along dip (km); different from vertical projection on z
        dip (float):
            Dip angle in degrees (0 = horizontal, 90 = vertical)
            Dip to the right relative to the strike direction.
        strike (float):
            Strike angle in degrees clockwise from North (y-axis)
        n_strike (int):
            Number of subpatches along strike
        n_dip (int):
            Number of subpatches along dip
        
        --- Boundary conditions:

        kode (int):
            Code of the b.c.
        ss (scalar):
            b.c. along strike
        ds (scalar):
            b.c. along dip
        ts (scalar):
            b.c. along normal
        
        --- Friction

        fcode (int):
            Friction code
        sdrop (scalar):
            Stress drop parameter
        rhoLitho (scalar):
            Average density of the rocks [unit: need to be consistent with the medium Young's modulus; denisty in kg.m^{-3} is the Young's modulus of the medium is in Pa (recommended)]
        rhoFluid (scalar):
            Average density of the fluid [unit: same as rhoLitho]
        C (scalar):
            Cohesion [unit: same as medium Young's modulus]
        mu (float):
            local friction
        
        --- Others:

        verbose (bool):
            Parameter controlling the verbose level in the terminal


    Returns:
        patches (instance of PatchCollection):
            Contains instances of Patch, themselfs defining all the sub
            element (patches) and containing:
                - x0, y0, z0 : top-left coordinates of the patch
                - strike : strike angle in degrees
                - dip : dip angle in degrees
                - length : along strike (km)
                - width : along dip (km)
                - All the boundary conditions
                - All the frictional conditions
    """
    pName = 'discretizeFault'

    im('Patch collection builder', pName, verbose)
    im('Geometry:', pName, verbose)
    im('    - x0 (East):  '+str(x0)+' (km)', pName, verbose)
    im('    - y0 (North): '+str(y0)+' (km)', pName, verbose)
    im('    - z0 (Up):    '+str(z0)+' (km)', pName, verbose)
    im('    - total length along strike: '+str(L)+' (km)', pName, verbose)
    im('    - total width along dip:     '+str(W)+' (km)', pName, verbose)
    im('    - number of subpatches along strike: '+str(n_strike), pName, verbose)
    im('    - number of subpatches along dip:    '+str(n_dip), pName, verbose)
    im('    - total number of subpatches:        '+str(int(n_strike*n_dip)), pName, verbose)
    im('Boundary Cond. applied to every subpatch',pName, verbose)
    if kode == 0 and ss==0 and ds==0 and ts == 0:
        txt = ' (default: No motion imposed, locked)'
    else:
        txt = ''
    im('    - kode: '+str(kode)+txt, pName, verbose)
    im('    - condition 1 (ss): '+str(ss)+txt, pName, verbose)
    im('    - condition 2 (ds): '+str(ss)+txt, pName, verbose)
    im('    - condition 2 (ts): '+str(ss)+txt, pName, verbose)
    im('Frictional Cond. applied to every subpatch', pName, verbose)
    im('    - friction code: '+str(fcode), pName, verbose)
    if fcode == 0:
        txt = '  (IGNORED)'
    else:
        txt = ''
    im('    - sdrop:    '+str(sdrop)+txt, pName, verbose)
    im('    - rhoLitho: '+str(rhoLitho)+txt, pName, verbose)
    im('    - rhoFluid: '+str(rhoFluid)+txt, pName, verbose)
    im('    - C:        '+str(C)+txt, pName, verbose)
    im('    - mu:       '+str(mu)+txt, pName, verbose)

    patches = []

    # Convert angles to radians
    dip_rad = np.deg2rad(dip)
    strike_rad = np.deg2rad(strike)

    # Patch dimensions
    dL = L / n_strike
    dW = W / n_dip

    # Strike unit vector in xy-plane
    strike_vector = np.array([np.sin(strike_rad), np.cos(strike_rad), 0.0])

    # Dip vector: perpendicular to strike
    dip_vector = np.array([
        +np.cos(strike_rad) * np.cos(dip_rad),
        -np.sin(strike_rad) * np.cos(dip_rad),
        -np.sin(dip_rad)   # depth is negative z
    ])

    for i in range(n_strike):
        for j in range(n_dip):
            # Top-left corner of patch
            x_patch = x0 + i * dL * strike_vector[0] + j * dW * dip_vector[0]
            y_patch = y0 + i * dL * strike_vector[1] + j * dW * dip_vector[1]
            z_patch = z0 + i * dL * strike_vector[2] + j * dW * dip_vector[2]

            patch = Patch(
                x0 = x_patch,
                y0 = y_patch,
                z0 = z_patch,
                strike = strike,
                dip = dip,
                L = dL,
                W = dW,
                kode = kode,
                ss = ss,
                ds = ds,
                ts = ts,
                fcode = fcode,
                sdrop = sdrop,
                rhoLitho = rhoLitho,
                rhoFluid = rhoFluid,
                C = C,
                mu = mu
            )

            patches.append(patch)

    patches = PatchCollection(patches, group=group)

    im('Collection built.', pName, verbose)
    return patches



def displacement_sdt_to_xyz(u_s, u_d, u_t, strike, dip):
    """
    Convert displacement(s) from (strike, dip, tensile)
    to (x, y, z) coordinates.

    Args:
        u_s, u_d, u_t (float or array-like, shape (N,)):
            Displacement components in strike, dip (downward), tensile directions
        strike, dip (float or array-like, shape (N,)):
            Strike (deg, clockwise from North) and dip (deg, right-dipping)

    Returns:
        u_xyz (ndarray, shape (3,) for scalar input or (N, 3) for array input):
            Output x,y,z displacement.
    """
    u_s = np.asarray(u_s)
    u_d = np.asarray(u_d)
    u_t = np.asarray(u_t)
    strike = np.asarray(strike)
    dip = np.asarray(dip)

    strike_rad = np.deg2rad(strike)
    dip_rad = np.deg2rad(dip)

    # Strike unit vector
    s = np.stack([
        np.sin(strike_rad),
        np.cos(strike_rad),
        np.zeros_like(strike_rad)
    ], axis=-1)

    # Dip unit vector (down, right-dipping)
    d = np.stack([
        np.cos(strike_rad) * np.cos(dip_rad),
       -np.sin(strike_rad) * np.cos(dip_rad),
       -np.sin(dip_rad)
    ], axis=-1)

    # Normal unit vector
    n = np.cross(s, d)
    n /= np.linalg.norm(n, axis=-1, keepdims=True)

    # Combine components
    u_xyz = (
        u_s[..., None] * s +
        u_d[..., None] * d +
        u_t[..., None] * n
    )

    return u_xyz





class UniformGrid:
    """
    Uniform Cartesian grid in x, y, z.
    """

    def __init__(self,
                 xmin, xmax, nx,
                 ymin, ymax, ny,
                 zmin, zmax, nz,
                 verbose = True):
        """
        Args:
            xmin, xmax (float)
                Grid extent in x (East)
                If xmin == xmax and nx = 1 then, only one layer of x
            ymin, ymax (float)
                Grid extent in y (North)
                If ymin == ymax and ny = 1 then, only one layer of y
            zmin, zmax (float)
                Grid extent in z (Up)
                If zmin == zmax and nz = 1 then, only one layer of z
            nx, ny, nz (int)
                Number of cells in x, y, z
        """

        self.xmin, self.xmax = xmin, xmax
        self.ymin, self.ymax = ymin, ymax
        self.zmin, self.zmax = zmin, zmax

        self.nx, self.ny, self.nz = int(nx), int(ny), int(nz)

        self.check()

        # Cell sizes and grid axis
        if nx > 1:
            self.dx = (xmax - xmin) / (nx-1)
            self.xaxis = np.linspace(xmin, xmax, nx)
        else:
            self.dx = 0
            if xmin == xmax:
                self.xaxis = np.array([xmin])
            else:
                raise ValueError('Inconsistent values: if nx = 1, xmin = xmax is expected.')
        if ny > 1:
            self.dy = (ymax - ymin) / (ny-1)
            self.yaxis = np.linspace(ymin, ymax, ny)
        else:
            self.dy = 0
            if ymin == ymax:
                self.yaxis = np.array([ymin])
            else:
                raise ValueError('Inconsistent values: if ny = 1, ymin = ymax is expected.')
        if nz > 1:
            self.dz = (zmax - zmin) / (nz-1)
            self.zaxis = np.linspace(zmin, zmax, nz)
        else:
            self.dz = 0
            if zmin == zmax:
                self.zaxis = np.array([zmin])
            else:
                raise ValueError('Inconsistent values: if nz = 1, zmin = zmax is expected.')

        # Meshed grid
        self.x, self.y, self.z = np.meshgrid(self.xaxis, self.yaxis, self.zaxis, indexing='ij')

        # others
        self.verbose = verbose
        self.im('Uniform grid initialized')
    

    def im(self, textMessage, error=False, warn=False, structure=False, end=False):
        """
        Internal message function.
        """
        im(textMessage, 'UniformGrid', self.verbose, error=error, warn=warn, structure=structure, end=end)


    @property
    def type(self):
        return type(self)

    @property
    def shape(self):
        """Grid shape (nx, ny, nz)."""
        return self.nx, self.ny, self.nz

    @property
    def size(self):
        """Total number of cells."""
        return self.nx * self.ny * self.nz

    @property
    def cell_volume(self):
        """Return the volume of a single cell."""
        return self.dx * self.dy * self.dz

    @property
    def extent(self):
        """Return grid extent as (xmin, xmax, ymin, ymax, zmin, zmax)."""
        return (self.xmin, self.xmax,
                self.ymin, self.ymax,
                self.zmin, self.zmax)

    def check(self):
        if not isinstance(self.nx, int) or not isinstance(self.ny, int) or not isinstance(self.nz, int) or \
            self.nx < 0 or self.ny < 0 or self.nz <0:
            raise TypeError('The grid dimensions (nx,ny,nz) have to be positive integers.')
        if self.xmin > self.xmax:
            raise ValueError('Condition xmin >= xmax not fullfilled')
        if self.ymin > self.ymax:
            raise ValueError('Condition ymin >= ymax not fullfilled')
        if self.zmin > self.zmax:
            raise ValueError('Condition zmin >= zmax not fullfilled')
    

    def get(self):
        return self.x.flatten(), self.y.flatten(), self.z.flatten()


    def points(self):
        """
        Return all cell center coordinates as an array.

        Returns
        -------
        pts : ndarray, shape (N, 3)
        """
        return np.column_stack([self.x.ravel(), self.y.ravel(), self.z.ravel()])




class Patch:

    def __init__(self, x0=None, y0=None, z0=None,\
                 L=None, W=None, dip=None, strike=None,\
                 kode=10, ss=0, ds=0, ts=0,\
                 fcode=0, sdrop=None, rhoLitho=None, rhoFluid=None, C=None, mu=None):
        # --- Geometry
        self.x0 = x0            # origin of the element in the x (East) direction [km]
        self.y0 = y0            # origin of the element in the y (North) direction [km]
        self.z0 = z0            # origin of the element in the z (Elevation) direction [km]
        self.L  = L             # length of the element in the strike direction [km]
        self.W  = W             # width of the element along the dip direction [km]
        self.dip = dip          # dip of the element [deg from horizontal]
        self.strike = strike    # strike azimuth [deg clockwise from North]
        # --- Boundary conditions
        self.kode = kode        # boundary condition code
        self.ss = ss            # b.c. along strike
        self.ds = ds            # b.c. along dip
        self.ts = ts            # b.c. along normal
        # --- Friction
        self.check_friction(fcode, sdrop, rhoLitho, rhoFluid, C, mu)
        self.fcode = fcode      # friction code
        self.sdrop = sdrop      # stress drop parameter
        self.rhoLitho = rhoLitho# average density of the rocks [unit: need to be consistent with the medium Young's modulus; denisty in kg.m^{-3} is the Young's modulus of the medium is in Pa (recommended)]
        self.rhoFluid = rhoFluid# average density of the fluid [unit: same as rhoLitho]
        self.C     = C          # cohesion [unit: same as medium Young's modulus]
        self.mu    = mu         # local friction coefficient
    
    @property
    def area(self):
        return self.L * self.W

    @property
    def length(self):
        return self.L
    
    @property
    def width(self):
        return self.W

    @property
    def center(self):
        """
        Geometric center of element (z<0 is depth, dip=90 vertical)
        """
        strike_rad = np.deg2rad(self.strike)
        dip_rad = np.deg2rad(self.dip)

        strike_vector = np.array([np.sin(strike_rad), np.cos(strike_rad), 0.0])
        dip_vector = np.array([
            +np.cos(strike_rad) * np.cos(dip_rad),
            -np.sin(strike_rad) * np.cos(dip_rad),
            -np.sin(dip_rad)   # depth is negative z
        ])

        center_vector = 0.5 * self.L * strike_vector + 0.5 * self.W * dip_vector
        x_center = self.x0 + center_vector[0]
        y_center = self.y0 + center_vector[1]
        z_center = self.z0 + center_vector[2]

        return x_center, y_center, z_center

    @property
    def corners(self):
        # Convert angles to radians
        strike_rad = np.deg2rad(self.strike)
        dip_rad = np.deg2rad(self.dip)
        # Strike vector in xy-plane
        strike_vector = np.array([np.sin(strike_rad), np.cos(strike_rad), 0.0])
        # Dip vector: perpendicular to strike
        dip_vector = np.array([
            +np.cos(strike_rad) * np.cos(dip_rad),
            -np.sin(strike_rad) * np.cos(dip_rad),
            -np.sin(dip_rad)   # depth is negative z
        ])
        corner1 = np.array([self.x0, self.y0, self.z0])  # top-left
        corner2 = corner1 + self.L * strike_vector       # along strike
        corner3 = corner2 + self.W * dip_vector          # along dip
        corner4 = corner1 + self.W * dip_vector          # along dip from corner1
        return corner1, corner2, corner3, corner4


    @property
    def xc(self):
        xc, _, _ = self.center
        return xc

    @property
    def yc(self):
        _, yc, _ = self.center
        return yc
    
    @property
    def zc(self):
        _, _, zc = self.center
        return zc

    def check_friction(self, fcode, sdrop, rhoLitho, rhoFluid, C, mu):
        if fcode not in [0, 1, 2, 3, 4, 5]:
            raise ValueError('Invalid value of fcode.')
        if fcode == 0:
            pass
        else:
            if sdrop is None or rhoLitho is None or C is None or mu is None or rhoFluid is None:
                raise ValueError("Partially defined frictional event detected\nYou must set sdrop, rhoLitho, rhoFluid, C and mu\nfor each frictional element.")
            if rhoLitho < 0 or rhoFluid < 0:
                raise ValueError('Densities (rocks and fluid) are defined positive')



class PatchCollection:

    # class instance counter
    _count = 0

    def __init__(self, listPatches=[], autoload=True, group=None, verbose=True):
        # --- Group
        PatchCollection._count += 1
        self.group_init = int(PatchCollection._count)   # original group identifier of the class instance
        self.group  = []        # list of group identifier
        # --- Verbose
        self.verbose= verbose   # verbose control
        # --- Collection
        self.patches= []        # list of loaded py3Ddef.Patch objects
        self.ids    = []        # indices of each patch
        self.nop    = 0         # number of loaded patch
        # --- Geometry
        self.x0     = []        # origin of each element in the x (East) direction [km]
        self.y0     = []        # origin of each element in the y (North) direction [km]
        self.z0     = []        # origin of the element in the z (Elevation) direction [km]
        self.W      = []        # length of the element in the strike direction [km]
        self.L      = []        # width of the element along the dip direction [km]
        self.strike = []        # strike azimuth [deg clockwise from North]
        self.dip    = []        # dip of the element [deg from horizontal]
        self.xc     = []        # center of each patch in the x direction
        self.yc     = []        # center of each patch in the y direction
        self.zc     = []        # center of each patch in the z direction
        # --- Boundary conditions
        self.kode   = []        # b.c. code
        self.ss     = []        # b.c. along strike
        self.ds     = []        # b.c. along dip
        self.ts     = []        # b.c. along normal
        # --- Friction
        self.fcode    = []      # friction code
        self.sdrop    = []      # stress drop parameter
        self.rhoLitho = []      # average density of rocks [unit: need to be consistent with the medium Young's modulus; denisty in kg.m^{-3} is the Young's modulus of the medium is in Pa (recommended)]
        self.rhoFluid = []      # average fluid density [unit: same as rhoLitho]
        self.C        = []      # cohesion [unit: same as medium Young's modulus]
        self.mu       = []      # local friction coefficient
        # init
        self.init(self.nop)
        if autoload and len(listPatches) != 0:
            self.add(listPatches, group=group)

    def init(self, n):
        self.patches  = [Patch()]*n
        self.ids      = np.zeros(n, dtype=np.int32)
        self.x0       = np.zeros(n, dtype=np.float64)
        self.y0       = np.zeros(n, dtype=np.float64)
        self.z0       = np.zeros(n, dtype=np.float64)
        self.W        = np.zeros(n, dtype=np.float64)
        self.L        = np.zeros(n, dtype=np.float64)
        self.strike   = np.zeros(n, dtype=np.float64)
        self.dip      = np.zeros(n, dtype=np.float64)
        self.kode     = np.zeros(n, dtype=np.int32)
        self.ss       = np.zeros(n, dtype=np.float64)
        self.ds       = np.zeros(n, dtype=np.float64)
        self.ts       = np.zeros(n, dtype=np.float64)
        self.xc       = np.zeros(n, dtype=np.float64)
        self.yc       = np.zeros(n, dtype=np.float64)
        self.zc       = np.zeros(n, dtype=np.float64)
        self.group    = np.ones(n, dtype=np.int32) * self.group_init
        self.nop      = n
        self.fcode    = np.zeros(n, dtype=np.int32)
        self.sdrop    = np.zeros(n, dtype=np.float64)
        self.rhoLitho = np.zeros(n, dtype=np.float64)
        self.rhoFluid = np.zeros(n, dtype=np.float64)
        self.C        = np.zeros(n, dtype=np.float64)
        self.mu       = np.zeros(n, dtype=np.float64)
    
    def add(self, listPatches=[], group=None):
        group0 = group
        if group is None:
            group = self.group_init
        if not isinstance(group, int):
            raise TypeError('group need to be an integer')
        if isinstance(listPatches, list):
            if len(listPatches) == 0:
                raise ValueError('You should provide a list of Patch objects')
            self.patches += listPatches
            n = len(listPatches)
            self.nop += n
            # init the fields that will be concatenate to the current instance
            x0       = np.zeros(n, dtype=np.float64)
            y0       = np.zeros(n, dtype=np.float64)
            z0       = np.zeros(n, dtype=np.float64)
            W        = np.zeros(n, dtype=np.float64)
            L        = np.zeros(n, dtype=np.float64)
            strike   = np.zeros(n, dtype=np.float64)
            dip      = np.zeros(n, dtype=np.float64)
            kode     = np.zeros(n, dtype=np.int32)
            ss       = np.zeros(n, dtype=np.float64)
            ds       = np.zeros(n, dtype=np.float64)
            ts       = np.zeros(n, dtype=np.float64)
            xc       = np.zeros(n, dtype=np.float64)
            yc       = np.zeros(n, dtype=np.float64)
            zc       = np.zeros(n, dtype=np.float64)
            fcode    = np.zeros(n, dtype=np.int32)
            sdrop    = np.zeros(n, dtype=np.float64)
            rhoLitho = np.zeros(n, dtype=np.float64)
            rhoFluid = np.zeros(n, dtype=np.float64)
            C        = np.zeros(n, dtype=np.float64)
            mu       = np.zeros(n, dtype=np.float64)
            group2add  = np.array([group]*n, dtype=np.int32)
            # iterative reading
            for i in range(n):
                x0[i]       = listPatches[i].x0
                y0[i]       = listPatches[i].y0
                z0[i]       = listPatches[i].z0
                L[i]        = listPatches[i].L
                W[i]        = listPatches[i].W
                strike[i]   = listPatches[i].strike
                dip[i]      = listPatches[i].dip
                kode[i]     = listPatches[i].kode
                ss[i]       = listPatches[i].ss
                ds[i]       = listPatches[i].ds
                ts[i]       = listPatches[i].ts
                xc[i], yc[i], zc[i] = listPatches[i].center
                fcode[i]    = listPatches[i].fcode
                sdrop[i]    = listPatches[i].sdrop
                rhoLitho[i] = listPatches[i].rhoLitho
                rhoFluid[i] = listPatches[i].rhoFluid
                C[i]        = listPatches[i].C
                mu[i]       = listPatches[i].mu
        elif isinstance(listPatches, PatchCollection):
            self.patches += listPatches.patches
            self.nop += listPatches.nop
            x0       = listPatches.x0
            y0       = listPatches.y0
            z0       = listPatches.z0
            W        = listPatches.W
            L        = listPatches.L
            strike   = listPatches.strike
            dip      = listPatches.dip
            kode     = listPatches.kode
            ss       = listPatches.ss
            ds       = listPatches.ds
            ts       = listPatches.ts
            xc       = listPatches.xc
            yc       = listPatches.yc
            zc       = listPatches.zc
            fcode    = listPatches.fcode
            sdrop    = listPatches.sdrop
            rhoLitho = listPatches.rhoLitho
            rhoFluid = listPatches.rhoFluid
            C        = listPatches.C
            mu       = listPatches.mu
            if group0 is None: # i.e. if no group was specified
                group2add  = listPatches.group # take the one preexisting
            else:
                if isinstance(group, int):
                    group2add = np.array([group]*n, dtype=np.int32)
                else:
                    raise TypeError('Trial to overwrite the group number: group need to be an integer')
        else:
            raise TypeError('The input Collection of Patch should be either a list of Patch objects of a PatchCollection instance.')
        # concatenate
        self.ids      = np.arange(self.nop)
        self.x0       = np.concatenate((self.x0, x0))
        self.y0       = np.concatenate((self.y0, y0))
        self.z0       = np.concatenate((self.z0, z0))
        self.L        = np.concatenate((self.L, L))
        self.W        = np.concatenate((self.W, W))
        self.strike   = np.concatenate((self.strike, strike))
        self.dip      = np.concatenate((self.dip, dip))
        self.kode     = np.concatenate((self.kode, kode))
        self.ss       = np.concatenate((self.ss, ss))
        self.ds       = np.concatenate((self.ds, ds))
        self.ts       = np.concatenate((self.ts, ts))
        self.xc       = np.concatenate((self.xc, xc))
        self.yc       = np.concatenate((self.yc, yc))
        self.zc       = np.concatenate((self.zc, zc))
        self.group    = np.concatenate((self.group, group2add), dtype=np.int32)
        self.fcode    = np.concatenate((self.fcode, fcode))
        self.sdrop    = np.concatenate((self.sdrop, sdrop))
        self.rhoLitho = np.concatenate((self.rhoLitho,rhoLitho))
        self.rhoFluid = np.concatenate((self.rhoFluid,rhoFluid))
        self.C        = np.concatenate((self.C,C))
        self.mu       = np.concatenate((self.mu,mu))

    @property
    def shape(self):
        return self.xc.shape
    
    def im(self, textMessage, error=False, warn=False, structure=False, end=False):
        """
        Internal message function.
        """
        im(textMessage, 'PatchCollection', self.verbose, error=error, warn=warn, structure=structure, end=end)


    def get(self):
        return self.x0, self.y0, self.z0, \
               self.L, self.W, \
               self.strike, self.dip,\
               self.kode, self.ss, self.ds, self.ts

    def getFriction(self):
        return {
            'fcode':    self.fcode,
            'sdrop':    self.sdrop,
            'rhoLitho': self.rhoLitho,
            'rhoFluid': self.rhoFluid,
            'C':        self.C,
            'mu':       self.mu
            }

    def plot3D(self, **kwargs):
        plotFault3D(self, **kwargs)








