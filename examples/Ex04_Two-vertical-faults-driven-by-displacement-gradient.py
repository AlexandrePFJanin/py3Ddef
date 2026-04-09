import matplotlib.pyplot as plt
import numpy as np
from py3Ddef import DeformationRun, BackgroundDeformation
from py3Ddef.geometry import UniformGrid, PatchCollection, discreteDislocation


# ================ py3Ddef example ================

# Based on the Example 1 of 3D~def by Joan Gomberg  and Mike Elis
# Extension of the example 'Ex03' of py3Ddef with 2 faults

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

# --- Dislocations geometry

patches = PatchCollection()

fault1 = discreteDislocation(x0=-10, y0=10.1, z0=0, L=100, W=20, dip=90, strike=90, n_strike=10, n_dip=5, \
                             kode=12, ss=0, ds=0, ts=0)

fault2 = discreteDislocation(x0=30, y0=30.1, z0=0, L=100, W=20, dip=90, strike=90, n_strike=10, n_dip=5, \
                             kode=12, ss=0, ds=0, ts=0)

patches.add(fault1)
patches.add(fault2)

# NOTE: Instead of adding fault1 and fault2 to an empty PatchCollection object,
# you can also add fault2 directly to fault1 by doing:
# fault1.add(fault2)
# Nevertheless, creating an empty PatchCollection object make the script more readable.

# --- Computation of 3D deformation

# Define your solution model

solution = DeformationRun()

# Compute the deformation
# NOTE: You can enter all the input paramater either when defining the
#       the solution model or when calling the solver

solution.compute3Ddef(grid=grid, patches=patches, nu=nu, mu=mu, E=E, bg=bg)

# NOTE: As we didn't enter any specific values for the flag controlling the computation
# of the optional outputs (e.g. output_invariants), they are all skipped.

# Unit of the output displacements: centimeters.
solution.displ.rescale(1/100)   # conversion from cm to m
solution.dislocs.rescale(1/100) # conversion from cm to m

# Reshape all the outputs computed on the grid with the shape of the input grid
solution.reshapeGSolutions()

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

# See the example script 'Ex02' for more detail on solution.plotFault2D.
# New here: as you have in this example 2 faults not contain in the same 2D plane,
#           you may want to represent the displacement on them separetly, on two 2D
#           figures, each on a plane containing each fault individually.
#           The function plot2Dfault allows that very easily using the unique dislocation IDs.
#           Each dislocation as an unique ID, set during it's creation and a list of 'group' IDs
#           for each of the patches it contains. You can set using the argument 'group'
#           the group of patches you want to display.
#           To display only the patches of the PatchCollection named 'fault1' and not the one added
#           from the object 'fault2', give to group the information about the initial group ID
#           of the object 'fault1' at its creation (fault1.group_init).
#           The example below show an example.

solution.plotFault2D(group=fault1.group_init, title="Displacement field on 'Fault 1'", v=solution.dislocs.norm, vec=solution.dislocs, cmap=plt.cm.viridis, vec_color='k', cbar_title='|u| (m)')

solution.plotFault2D(group=fault2.group_init, title="Displacement field on 'Fault 2'", v=solution.dislocs.norm, vec=solution.dislocs, cmap=plt.cm.viridis, vec_color='k', cbar_title='|u| (m)')


# You can see the group ID of all the patches contained in 'patches':
# print(patches.group)
