import matplotlib.pyplot as plt
import numpy as np
from py3Ddef import compute3Ddef
from py3Ddef.geometry import discreteDislocation
from py3Ddef.viewer import plotFault2D


# ================ py3Ddef example ================

# Based on the Example 1 of 3D~def by Joan Gomberg  and Mike Elis

# ---- Example of a vertical strike-slip fault driven by dextral simple shear ----

# Context:  Vertical fault, stress free in along-strike and along-dip directection, no tensile displacement.
#           Background deformation: Imposed displacement gradient in dextral simple shear

# =================================================


# --- Medium parameters

nu = 0.25   # Poisson's ratio
E  = 7.e10  # Young's modulus (unit: Pa)
mu = 0.6    # coeff of internal friction

# --- Grid

# Define where to compute the solution outsite of the dislocation

nx = 50
ny = 60

x = np.linspace(-50, 150, nx, dtype=np.float64)    # in km
y = np.linspace(-50, 150, ny, dtype=np.float64)    # in km
z = np.array([-0], dtype=np.float64)               # in km

nz = len(z)

# Generate a mesh
x, y, z = np.meshgrid(x, y, z, indexing='ij')

# The input grid for py3Ddef should be flat
x = x.flatten()
y = y.flatten()
z = z.flatten()

# --- Background deformation (optional arguments)

# Background field:
# Name           flag     format
# strain:        'stra'   (Exx,    Exy,    Exz,    Eyy,    Eyz,    Ezz,    0,      0,      0)
# stress:        'stre'   (Sxx,    Sxy,    Sxz,    Syy,    Syz,    Szz,    0,      0,      0)
# displacement:  'disp'   (dUx/dx, dUx/dy, dUx/dz, dUy/dx, dUy/dy, dUy/dz, dUz/dx, dUz/dy, dUz/dz)
#    -> N.B. The three last zeros are mandatory for the stress and strain cond.

# Units of displacement: centimeters 

bg_flag  = 'disp'
bg_field = np.array([0, 1e-4, 0.0, 0, 0, 0.0, 0.0, 0.0, 0.0])

# --- Dislocations geometry

# You can use directly the function py3Ddef.geometry.discreteDislocation to efficiently and safely
# generate and manage the geometry and slip condition of your dislocations.
# The function will return a py3Ddef.geometry.PatchColletion object - will make your
# analysis of the results more efficient.

fault = discreteDislocation(x0=10, y0=10.1, z0=0, L=100, W=20, dip=90, strike=90, n_strike=10, n_dip=5, \
                            kode=12, ss=0, ds=0, ts=0)

# --- Computation of 3D deformation

# Conveniently, if you generated your dislocations using py3Ddef.geometry.discreteDislocation,
# or more generally, if you created a py3Ddef.geometry.PatchCollection object for your dislocations,
# you can use the .get() method to automatically provide py3Ddef.compte3Ddef with all the required arguments.
# If so, don't forget the unpack command '*'.

u,s,e,o,f,g,d = compute3Ddef(x,y,z,\
                             *fault.get(),\
                             nu,E,mu,\
                             bg_flag, bg_field)

# Unit of the output displacements: centimeters.
u.rescale(1/100)	# conversion from cm to m
d.rescale(1/100)	# conversion from cm to m

# --- Figure 1: map view

gshape = (nx, ny, nz)

x = x.reshape(gshape)
y = y.reshape(gshape)
z = z.reshape(gshape)

u.reshape(gshape)

uplot = u.z[:,:,0] # 0, at the surface

plot_surf = True
if plot_surf:
    fig = plt.figure(figsize=(7,6))
    ax  = fig.add_subplot(111)
    ax.set_title('Displacement field - map view')
    k = 1
    vmin = min(np.amin(uplot), -np.amax(uplot)) * k
    vmax = max(-np.amin(uplot), np.amax(uplot)) * k
    cmap = ax.pcolormesh(x[:, :, 0], y[:, :, 0], uplot, cmap=plt.cm.seismic, vmin=vmin, vmax=vmax)
    ax.quiver(x[:, :, 0], y[:, :, 0], u.x[:, :, 0], u.y[:, :, 0])
    fig.colorbar(cmap, ax=ax, shrink=0.5, label='uz (m)')
    ax.set_aspect('equal')
    ax.set_xlabel('x (i.e. east, in km)')
    ax.set_ylabel('y (i.e. west, in km)')
    fig.tight_layout()
    plt.show()


# --- Figure 2: Motion on the fault

# You can easily visualize the displacement of the resulting displacement on the dislocation
# using the function py3Ddef.viewer.plotFault2D. This function project the result in 2D
# on a plane normal to a given projection vector. By default, if not set, the vector will be normal
# to the 1st patch of the input PatchCollection (index 0). The reference patch can be changed.

# Below, we detail the different inputs of the function.

# PatchCollection object you want to represent
coll = fault

# Scalar field that will color each patches (can be None)
v = d.ss    # flat np.ndarray, whose the length correspond to the total number of patch

# Vector field that will be ploted at the center of the patches (can be None)
vec = d     # have to be a py3Ddef.ElementDisplacement object

# Plot the figure on the same plane as the patch number refID
normal = True   # switch
refID  = 0      # select the patch ID (will be important if all the dislocations are not in the same 2D plane)

# Otherwise, define your own projection vector
view_vec = None # Example of input: view_vec = (1, 1, 1)
                # This is the vector that will be use for the projection: projection on a plane normal to this input vector
                # Note: if not None, view_vec will overprint the instruction set on the arguments 'normal' and 'refID'.


# Call the function
plotFault2D(coll, v=v, vec=vec, normal=normal, refID=refID, view_vec=view_vec, \
            cmap=plt.cm.viridis, vec_color='white', cbar_title='Strike slip (m)')


