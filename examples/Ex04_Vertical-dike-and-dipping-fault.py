import matplotlib.pyplot as plt
import numpy as np
from py3Ddef import compute3Ddef
from py3Ddef.geometry import discreteDislocation, PatchCollection, UniformGrid
from py3Ddef.viewer import plotFault2D


# ================ py3Ddef example ================

# ---- Example of a vertical dike with injection triggering motion on a nearby dipping fault ----

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

# To generate easily a uniform grid, the recommended option is to use py3Ddef.geometry.UniformGrid
# Why? Numerous verification and visualization routines have been implemented to work specificaly with this object.

nx = 50
ny = 60
nz = 1
grid = UniformGrid(-150, 150, nx,\
                   -150, 150, ny,\
                      0, 0, nz)

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

# Note: You can quickly inspect the geometry of your PatchCollection object with:
inspect3D = False
if inspect3D:
    patches.plot3D()

# --- Computation of 3D deformation

u,s,e,o,f,g,d = compute3Ddef(*grid.get(),\
                             *patches.get(),\
                             nu,E,mu)

# Unit of the output displacements: centimeters.
u.rescale(1/100)	# conversion from cm to m
d.rescale(1/100)	# conversion from cm to m

# --- Figure 1: map view

u.reshape(grid.shape)

uplot = u.z[:,:,-1] # -1, at the surface

plot_surf = True
if plot_surf:
    fig = plt.figure(figsize=(7,6))
    ax  = fig.add_subplot(111)
    ax.set_title('Displacement field - map view')
    k = 1
    vmin = min(np.amin(uplot), -np.amax(uplot)) * k
    vmax = max(-np.amin(uplot), np.amax(uplot)) * k
    cmap = ax.pcolormesh(grid.x[:, :, 0], grid.y[:, :, 0], uplot, cmap=plt.cm.seismic, vmin=vmin, vmax=vmax)
    ax.quiver(grid.x[:, :, 0], grid.y[:, :, 0], u.x[:, :, 0], u.y[:, :, 0])
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
    fig = plt.figure(figsize=(7,3))
    ax  = fig.add_subplot(111)
    ax.plot(grid.x[:, index_j, -1], u.z[:, index_j, -1])
    ax.set_ylabel('uz (m)')
    ax.set_xlabel('x (i.e. east, in km)')
    fig.tight_layout()
    plt.show()


# --- Figure 3: Motion on the fault

plotFault2D(patches, v=d.ds, vec=d, group=[nfault.group_init], cmap=plt.cm.seismic, vec_scale=1e1)
