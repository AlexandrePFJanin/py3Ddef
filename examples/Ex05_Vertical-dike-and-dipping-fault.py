import matplotlib.pyplot as plt
import numpy as np
from py3Ddef import DeformationRun
from py3Ddef.geometry import discreteDislocation, PatchCollection, UniformGrid
from py3Ddef.viewer import plotFault2D


# ================ py3Ddef example ================

# ---- Example of a vertical dike with injection triggering motion on a nearby frictionless dipping fault ----

# Context:  One vertical dike with one meter of tensile opening
#           One dipping fault close to the dike, dipping toward the dike
#           The fault is stress free in along-strike and along-dip directection, no tensile displacement.
#           No background deformation.

# =================================================


# --- Medium parameters

nu = 0.25   # Poisson's ratio
E  = 7.e10  # Young's modulus (unit: Pa)
mu = 0.6    # coeff of internal friction

# --- Grid

nx = 200
ny = 200
nz = 1
grid = UniformGrid(-50, 50, nx,\
                   -20, 20, ny,\
                      0, 0, nz)

# --- Dislocations geometry and boundary conditions

patches = PatchCollection()

# Vertical dike
dike = discreteDislocation(x0=0, y0=0, z0=-10*1e-3, L=5, W=3.5, dip=90, strike=0, n_strike=10, n_dip=10, \
                           kode=10, ss=0, ds=0, ts=200)
# Normal fault
nfault = discreteDislocation(x0=-2, y0=0, z0=0, L=5, W=2, dip=60, strike=0, n_strike=10, n_dip=10, \
                             kode=12, ss=0, ds=0, ts=0)

patches.add(dike)
patches.add(nfault)

# Note: You can quickly inspect the geometry of your PatchCollection object with:
inspect3D = True
if inspect3D:
    patches.plot3D(centers=False, refpts=False)

# --- Computation of 3D deformation

# Define your solution model and solve the deformations
solution = DeformationRun(grid=grid, patches=patches, nu=nu, mu=mu, E=E)
solution.compute3Ddef()

# Unit of the output displacements: centimeters.
solution.displ.rescale(1/100)   # conversion from cm to m
solution.dislocs.rescale(1/100) # conversion from cm to m

# --- Figure 1: map view

# Reshape specifically only the output displacement field with the shape of the grid
solution.displ.reshape(grid.shape)

# NOTE: You can use the following alias in the instance 'solution' to call the corresponding fields
#       .u for .displ
#       .s for .stress
#       .e for .strain
#       .o for .pstrainOri
#       .f for .fplanes
#       .j for .invariants
#       .g for .gradDispl
#       .d for .dislocs

uplot = solution.u.z[:,:,-1] # -1, at the surface

plot_surf = True
if plot_surf:
    fig = plt.figure(figsize=(10,4))
    ax  = fig.add_subplot(111)
    ax.set_title('Displacement field - map view')
    k = 1
    vmin = min(np.amin(uplot), -np.amax(uplot)) * k
    vmax = max(-np.amin(uplot), np.amax(uplot)) * k
    r = 5
    cmap = ax.pcolormesh(grid.x[:, :, 0], grid.y[:, :, 0], uplot, cmap=plt.cm.seismic, vmin=vmin, vmax=vmax)
    ax.quiver(grid.x[::r, ::r, 0], grid.y[::r, ::r, 0], solution.u.x[::r, ::r, 0], solution.u.y[::r, ::r, 0])
    fig.colorbar(cmap, ax=ax, shrink=0.5, label='uz (m)')
    ax.set_aspect('equal')
    ax.set_xlabel('x (i.e. east, in km)')
    ax.set_ylabel('y (i.e. west, in km)')
    fig.tight_layout()
    plt.show()

# --- Figure 2: Across the system

index_j = int(ny/2) # where along the y-axis

plot_section = True
if plot_section:
    fig = plt.figure(figsize=(7,6))
    ax1 = fig.add_subplot(211)
    ax1.plot(grid.x[:, index_j, -1], solution.u.z[:, index_j, -1])
    ax1.set_ylabel('uz (m)')
    ax1.set_xlabel('x (i.e. east, in km)')
    ax2 = fig.add_subplot(212)
    ax2.plot(grid.x[:, index_j, -1], solution.u.x[:, index_j, -1])
    ax2.set_ylabel('ux (m)')
    ax2.set_xlabel('x (i.e. east, in km)')
    fig.tight_layout()
    plt.show()


# --- Figure 3: Motion on the fault

# Call the function py3Ddef.viewer.plotFault2D
# It's the function called when you asked 'solution.plotFault2D'
# The call the same excepte that you need to give in input a PatchCollection object.

# NOTE: You can of course use directly solution.plotFault2D - just a demo.

plotFault2D(patches, v=solution.d.us, vec=solution.d, group=[nfault.group_init], cmap=plt.cm.seismic, vec_scale=1e1)
