import matplotlib.pyplot as plt
import numpy as np
from py3Ddef import DeformationRun
from py3Ddef.geometry import UniformGrid, discreteDislocation


# ================ py3Ddef example ================

# Exact same context as Ex01_Dike-induced-faulting.py,
# but now using the full interface of py3Ddef.

# Based on the work of Rubin, A. M., & Pollard, D. D. (1988)
# (Dike-induced faulting in rift zones of Iceland and Afar. Geology, 16(5), 413-417.)
# Figure 2.A.

# ---- Example of a vertical dike extension ----

# Context:  One meter of opening of a vertical dike from 0.4 to 4.4 km depth.
#           No background deformation.
#           We want to compute the surface displacement resulting from this
#           opening in the vertical and transverse directions of the dike.

# =================================================


# --- Medium parameters

nu = 0.25   # Poisson's ratio
E  = 30.e9  # Young's modulus (here in Pa)
mu = 0.6    # coeff of internal friction

# --- Grid

# To generate easily a uniform grid, the recommended option is to use py3Ddef.geometry.UniformGrid
# Why? Numerous verification and visualization routines have been implemented to work specificaly with this object.

nx = 400 + 1
ny = 1
nz = 1
grid = UniformGrid(-20, 20, nx,\
                     0,  0, ny,\
                     0,  0, nz)

# --- Dislocations geometry and boundary conditions

# You can use directly the function py3Ddef.geometry.discreteDislocation to efficiently and safely
# generate and manage the geometry and slip condition of your dislocations.
# The function will return a py3Ddef.geometry.PatchColletion object - will make your
# analysis of the results more efficient.

dike = discreteDislocation(x0=0, y0=-50, z0=-0.4, L=100, W=4, dip=90, strike=0, n_strike=1, n_dip=1, \
                           kode=10, ss=0, ds=0, ts=100)

# --- Computation of 3D deformation

solution = DeformationRun(grid=grid, patches=dike, nu=nu, mu=mu, E=E,\
                          output_invariants=True,\
                          output_pstrainOri=True,\
                          output_fplanes=True,\
                          output_gradDispl=True)

solution.compute3Ddef()

# Description of the outputs:

# 1. Outputs computed on the input grid:

# Default outputs (allways generated):
# solution.displ:       Displacement field: vectorial field (type: py3Ddef.GridDisplacement)
# solution.stress:      Stress field: symmetric tensor (type: py3Ddef.GridStress)
# solution.strain:      Strain field: symmetric tensor (type: py3Ddef.GridStrain)

# Additional outputs (depend on the flags):
# solution.pstrainOri:  Principal strain orientation: ex, px, tx, ey, py, ty, ez, pz, tz (type: py3Ddef.GridPrincipalStrainOrientation)
# solution.fplanes:     Optimal failure planes: strike1, dip1, rake1, strike2, dip2, rake2 (type: py3Ddef.GridOptimalFailurePlane)
# solution.invariants:  Stress/strain invariants: volumetric strain, critical failure stress, octahedral shear stress, strain energy density (type: py3Ddef.GridStressStrainInvariants)
# solution.gradDispl:   Displacement gradient field: asymmetric tensor (type: py3Ddef.GridDisplacementGradient)

# 2. Outputs computed on the dislocation patches:

# solution.dislocs:     Displacement on each patches given in (strike slip, dip slip, tensil sip) (type: py3Ddef.ElementDisplacement)
#                       Can be converted to (x-east, y-north, z-up) using the function d.convert2xyz() and the geometry of the dislocations.

# 3. Units of the different input: See the README of py3Ddef (similar to 3D~def.)

# Unit of the output displacements: centimeters.
solution.displ.rescale(1/100)   # conversion from cm to m
solution.dislocs.rescale(1/100) # conversion from cm to m

# --- Figure - displacement across the dike

fig = plt.figure(figsize=(10,3))
ax  = fig.add_subplot(111)
ax.set_title('Surface displacements resulting from the opening of a vertical dyke by %s cm'%str(solution.patches.ts[0]))
ax.hlines(0, np.amin(solution.grid.x), np.amax(solution.grid.x), alpha=0.5, colors='k', linestyles='--')
ln1 = ax.plot(solution.grid.x[:, 0, 0], solution.displ.z[:, 0, 0],label=r'$u_z(x,y=0)$')
ln3 = ax.plot(solution.grid.x[:, 0, 0], solution.displ.x[:, 0, 0],label=r'$u_x(x,y=0)$')
ax.set_xlabel('x-axis (km)')
ax.set_ylabel(r'$u_i(x)$'+ ' (m)')
ax.legend()
fig.tight_layout()
plt.show()
