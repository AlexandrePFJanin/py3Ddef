import matplotlib.pyplot as plt
import numpy as np
from py3Ddef import DeformationRun, BackgroundDeformation
from py3Ddef.geometry import UniformGrid, discreteDislocation


# ================ py3Ddef example ================

# Based on the Example 1 of 3D~def by Joan Gomberg and Mike Elis

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
nz = 1
grid = UniformGrid(-50, 150, nx,\
                   -50, 150, ny,\
                     0,  0, nz)

# --- Background deformation (optional arguments)

# Background field:
# Name           flag     format
# strain:        'stra'   (Exx,    Exy,    Exz,    Eyy,    Eyz,    Ezz,    0,      0,      0)
# stress:        'stre'   (Sxx,    Sxy,    Sxz,    Syy,    Syz,    Szz,    0,      0,      0)
# displacement:  'disp'   (dUx/dx, dUx/dy, dUx/dz, dUy/dx, dUy/dy, dUy/dz, dUz/dx, dUz/dy, dUz/dz)
#    -> N.B. The three last zeros are mandatory for the stress and strain cond.

# Units of displacement: centimeters 

bg = BackgroundDeformation(bg_flag='disp',\
                           bg_field=np.array([0, 1e-4, 0.0, 0, 0, 0.0, 0.0, 0.0, 0.0]))

# --- Dislocations geometry and boundary conditions

# Define your fault geometry and boundarby conditions:

fault = discreteDislocation(x0=10, y0=10.1, z0=0, L=100, W=20, dip=90, strike=90, n_strike=10, n_dip=5, \
                            kode=12, ss=0, ds=0, ts=0)

# --- Computation of 3D deformation

# Define your solution model

solution = DeformationRun()

# Compute the deformation
# NOTE: You can enter all the input paramater either when defining the
#       the solution model or when calling the solver

solution.compute3Ddef(grid=grid, patches=fault, nu=nu, mu=mu, E=E, bg=bg)

# NOTE: As we didn't enter any specific values for the flag controlling the computation
# of the optional outputs (e.g. output_invariants), they are all skipped.

# Unit of the output displacements: centimeters.
solution.displ.rescale(1/100)   # conversion from cm to m
solution.dislocs.rescale(1/100) # conversion from cm to m

# Reshape all the outputs computed on the grid with the shape of the input grid
solution.reshapeGSolutions() # If the parameter 'auto_reshape' of the intance of DeformationRun is True (default set to True, then will be done automatically)

# --- Figure 1: map view

uplot = solution.displ.z[:,:,0] # 0, at the surface

plot_surf = True
if plot_surf:
    fig = plt.figure(figsize=(7,6))
    ax  = fig.add_subplot(111)
    ax.set_title('Displacement field - map view')
    k = 1
    vmin = min(np.amin(uplot), -np.amax(uplot)) * k
    vmax = max(-np.amin(uplot), np.amax(uplot)) * k
    cmap = ax.pcolormesh(grid.x[:, :, 0], grid.y[:, :, 0], uplot, cmap=plt.cm.seismic, vmin=vmin, vmax=vmax)
    ax.quiver(grid.x[:, :, 0], grid.y[:, :, 0], solution.displ.x[:, :, 0], solution.displ.y[:, :, 0])
    fig.colorbar(cmap, ax=ax, shrink=0.5, label='uz (m)')
    ax.set_aspect('equal')
    ax.set_xlabel('x (i.e. east, in km)')
    ax.set_ylabel('y (i.e. west, in km)')
    fig.tight_layout()
    plt.show()


# --- Figure 2: Motion on the fault

# You can easily visualize the displacement of the resulting displacement on the dislocation
# using the function solution.plotFault2D. This function project the result in 2D
# on a plane normal to a given projection vector. By default, if not set, the vector will be normal
# to the 1st patch of the input PatchCollection (index 0). The reference patch can be changed.

# NOTE: solution.plotFault2D(), directly calls py3Ddef.viewer.plotFault2D()



# Below, we detail the different inputs of the function.

# Scalar field that will color each patches (can be None)
v = solution.dislocs.norm  # norm of the displacemenet, flat np.ndarray, whose the length correspond to the total number of patch

# Vector field that will be ploted at the center of the patches (can be None)
vec = solution.dislocs     # have to be a py3Ddef.ElementDisplacement object

# Plot the figure on the same plane as the patch number refID
normal = True   # switch
refID  = 0      # select the patch ID (will be important if all the dislocations are not in the same 2D plane)

# Otherwise, define your own projection vector
view_vec = None # Example of input: view_vec = (1, 1, 1)
                # This is the vector that will be use for the projection: projection on a plane normal to this input vector
                # Note: if not None, view_vec will overprint the instruction set on the arguments 'normal' and 'refID'.


# Call the function
solution.plotFault2D(v=v, vec=vec, normal=normal, refID=refID, view_vec=view_vec, \
                     cmap=plt.cm.viridis, vec_color='k', cbar_title='|u| (m)')


