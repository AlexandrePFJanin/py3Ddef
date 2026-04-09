import numpy as np
from py3Ddef import DeformationRun
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

# Define your solution model and solve the deformations
solution = DeformationRun(grid=grid, patches=patches, nu=nu, mu=mu, E=E)
solution.compute3Ddef()

# Unit of the output displacements: centimeters.
solution.displ.rescale(1/100)   # conversion from cm to m
solution.dislocs.rescale(1/100) # conversion from cm to m

# --- 3D Visualization

# Using py3Ddef.DeformationRun for solving your deformation,
# the outputs computed on the grid and on the dislocations
# are directly linked to them.
# This will then allows you to easily manage and visualize
# complexe 3D data using Paraview.


# 1. -- Export the dislocations and their linked outputs

# Export the dislocations and their linked outputs in XDMF/HDF5 files.

fname = 'Ex06_result-dislocations'  # Name of the files (without their extension)
path  = './'                        # Path to the export directory

solution.patches2paraview(fname, path=path)  # based on py3Ddef.viewer.patches2paraview()

# Note: Also creates three vectorial field, centered on the center of patches
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


# 2. -- Export the grid and its linked outputs

# Export the grid and its linked outputs in XDMF/HDF5 files.

fname = 'Ex06_results-grid'         # Name of the files (without their extension)
path  = './'                        # Path to the export directory

solution.grid2paraview(fname, path=path) # based on py3Ddef.viewer.grid2paraview()

# Then open Paraview and load the file './results-grid.xdmf'
# to visualize in 3D the linked fields.


# ---------------------
# NOTE, to go further:

# You can attach additional fields to either the dislocations or the grid
# for your display in Paraview e.g. a field you computed with your own
# function.
# For that check the documentation of the functions solution.grid2paraview
# and solution.patches2paraview. 
# (in a terminal >> help(solution.patches2paraview))
