import numpy as np
from py3Ddef import compute3Ddef
from py3Ddef.geometry import discreteDislocation, PatchCollection, UniformGrid


# ================ py3Ddef example ================

# ---- Example of a vertical dike with injection triggering motion on a nearby fault ----

# Context:  One vertical dike with one meter of tensile opening
#           One dipping fault close to the dike, dipping toward the dike
#           The fault is stress free in along-strike and along-dip directection, no tensile displacement.
#           No background deformation.

# This example is based on 'Ex04' and has the same context and the same deformation.
# The focus here is on the 3D visualization options provided by py3Ddef and Paraview.
# We just changed here the grid on wich the result is computed to increase the resolution
# close the dislocations.

# =================================================


# --- Medium parameters

nu = 0.25   # Poisson's ratio
E  = 7.e10  # Young's modulus (unit: Pa)
mu = 0.6    # coeff of internal friction

# --- Grid

# To generate easily a uniform grid, the recommended option is to use py3Ddef.geometry.UniformGrid
# Why? Numerous verification and visualization routines have been implemented to work specificaly with this object.

nx = 100    # Resolution increased
ny = 100
nz = 40     # 3D grid, the computation of the deformation will take more time than in 'Ex04'.
grid = UniformGrid(-20, 20, nx,\
                   -20, 20, ny,\
                    -10, 0, nz)     # bounds of the grid: closer to the dislocations

# When one of the dimension is 1 (e.g. nz = 1) and if the two bounds in this
# direction are the same (e.g. zmin = 0, zmax = 0), then UniformGrid will
# create a single layer in this direction at the value zmin.

# --- Dislocations geometry

patches = PatchCollection()

# Vertical dike
dike = discreteDislocation(x0=0, y0=0, z0=-10*1e-3, L=5, W=3.5, dip=90, strike=0, n_strike=10, n_dip=10, \
                           kode=10, ss=0, ds=0, ts=200)
# Normal fault
nfault = discreteDislocation(x0=-2, y0=0, z0=0, L=5, W=2, dip=60, strike=0, n_strike=10, n_dip=10, \
                             kode=12, ss=0, ds=0, ts=0)

patches.add(dike)
patches.add(nfault)

# --- Computation of 3D deformation

u,s,e,o,f,g,d = compute3Ddef(*grid.get(),\
                             *patches.get(),\
                             nu,E,mu)

# Unit of the output displacements: centimeters.
u.rescale(1/100)	# conversion from cm to m
d.rescale(1/100)	# conversion from cm to m

# --- 3D Visualization 1 - The deformation grid.

# Using py3Ddef.geometry.UniformGrid for creating your computation grid,
# you can easily link to it an output field of py3Ddef.
# This will then allows you to easily manage and visualize at the same time
# complexe 3D data.

# 1. -- Link fields to the grid

# You can directly link the output fields (you can link as many fields as wanted)

grid.link(u)        # link the displacement field
grid.link(s)        # link the stress tensor

# You can also create new fields and link them to the grid
# if they have the same dimension and shape as the grid.

# Example: create a point IDs field (i.e. a scalar field)
ptIDS = np.arange(grid.nx * grid.ny * grid.nz).reshape(grid.shape)
grid.link(ptIDS, name='ptID')   # link the field: you can set a name
                                # if the name is not set or if it's not obvious for py3Ddef,
                                # a name will be automatically affected to the new linked field


# Type of data you can link:
#   - numpy.ndarray
#   - py3Ddef.GridScalar
#   - py3Ddef.GridVector
#   - py3Ddef.GridTensor
#   - py3Ddef.GridDisplacement
#   - py3Ddef.GridStress
#   - py3Ddef.GridStrain
#   - py3Ddef.GridPrincipalStrainOrientation
#   - py3Ddef.GridOptimalFailurePlane
#   - py3Ddef.GridDisplacementGradient

# Shapes:
#   - Scalar fields:    (nx, ny, nz)
#   - Vectorial fields: (nx, ny, nz, 3)
#   - Tensorial fields: (nx, ny, nz, 3, 3), (nx, ny, nz, 6), (nx, ny, nz, 9)

# All the linked field will appears in 'grid.v' and their names in 'grid.vnames'

print('Fields linked to the grid:')
for i in range(len(grid.v)):
    print(' - (%s) %s:\t\t[%s]'%(str(i), grid.vnames[i], str(grid.v[i])))

# You can remove a linked field knowing it's position in 'grid.v'
grid.removeField(2) # Here we remove the 'ptID' field.

# 2. -- Export in XMDF+HDF5 for Paraview

# Export the grid and the linked fields to
# a file format readable by Paraview (XDMF+HDF5).

fname = 'Ex05_results-grid'
path  = './'

grid.export2paraview(fname, path=path) # based on py3Ddef.viewer.grid2paraview()

# Then open Paraview and load the file './results-grid.xdmf'
# to visualize in 3D the linked fields.


# --- 3D Visualization 2 - The input dislocation

# Similarly to the grid, you can link a displacement field
# to a PatchCollection object.

# 1. -- Link fields to the grid

patches.link(d)

# 2. -- Export in XMDF+HDF5 for Paraview

# Export the grid and the linked fields to
# a file format readable by Paraview (XDMF+HDF5).


fname = 'Ex05_result-dislocations'
path  = './'

patches.export2paraview(fname, path=path)  # based on py3Ddef.viewer.patches2paraview()

# Note: Also create three vectorial field, centered on the center of patches
#       and representing the three local direction for each patch: (along-strike,
#       along-dip, along-normal) in order to better read the displacement on the 
#       dislocation (as it's ever the motion of the hangingwall or the motion of
#       the footwall (default)).
#       With the right-handed convention (used by py3Ddef), the hangingwall lies
#       to the right of the strike direction.

# Then open Paraview and load the file './results-dislocations.xdmf'
# to visualize in 3D the linked fields.

# In Paraview:
#   - Open the XDMF file with the "XDMF reader" option (not the "Xdmf3 Reader S" or "Xdmf3 Reader T")
#   - Separate "FaultPatches" and "PatchCenters" with "Extract Block" (avoid issues with "Glyph" representation)
