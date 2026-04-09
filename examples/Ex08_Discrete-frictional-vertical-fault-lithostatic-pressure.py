import matplotlib.pyplot as plt
import numpy as np
from py3Ddef import DeformationRun, BackgroundDeformation
from py3Ddef.geometry import discreteDislocation, UniformGrid
from py3Ddef.viewer import plotFault2D


# ================ py3Ddef example ================

# ---- Example of a discretized frictional vertical fault driven by a background
#      stress field with lithostatic pressure ----

# Context:  Vertical fault discretized in 20x20 dislocations under a background stress
#           field promoting strike-slip displacement on the fault.
#           The lithostatic pressure is considered in this problem as acting on the
#           normal stress on the fault and acting against its sliding. The
#           listhostatic stresses cannot generate motion on their own.

# =================================================


# --- Medium parameters

nu = 0.25   # Poisson's ratio
E  = 7.e10  # Young's modulus (unit: Pa)
mu = 0.6    # coeff of internal friction

# --- Grid

# Minimal grid: we are focusing here only on the slip distribution on the dislocations
# not on the surrounding elastic medium.

nx = 1
ny = 1
nz = 1
grid = UniformGrid(1, 1, nx,\
                   1, 1, ny,\
                  -1,-1, nz)


# --- Background deformation

bg = BackgroundDeformation(bg_flag='stre',\
                        bg_field=np.array([0, 50e6, 0.0, 0, 0, 0.0, 0.0, 0.0, 0.0]))

# --- Dislocations geometry and boundary conditions

# fritional vertical strike-slip fault:
# ss = ds = ts = 0: no imposed stresses stored on the fault before the computation of the deformation
# fcode = 1 and sdrop = 1, frictional element and if failure, return to envelope
# C = 0, mu = 0.6, cohesion 0 Pa and coefficient of friction of 0.6
# rhoLitho = 3000 kg.m^{-3}: average rock density, used to compute the lithostatic pressure
# rhoFluid = 0 kg.m^{-3}: Density of the fluid: here no fluid
fault = discreteDislocation(x0=0, y0=0, z0=0, L=5, W=5, dip=90, strike=90, n_strike=20, n_dip=20, \
                            kode=2, ss=0, ds=0, ts=0, rhoLitho=3000, rhoFluid=0,\
                            fcode=1, sdrop=1, mu=0.6)


# --- Computation of 3D deformation

# Define your solution model and solve the deformations
maxiter  = 150  # you may want to increase the maximum number of iterations to reach convergence (unecessary here, just for the example, by default maxiter=100)
solution = DeformationRun(grid=grid, patches=fault, nu=nu, mu=mu, E=E, bg=bg)
solution.compute3Ddef(maxiter=maxiter)

# Unit of the output displacements: centimeters.
solution.displ.rescale(1/100)   # conversion from cm to m
solution.dislocs.rescale(1/100) # conversion from cm to m

# --- Figure: slip on the fault

# Make a mask 'patches2hatch' of the dislocations on the fault that are locked.
mfault = solution.patches.group == fault.group_init # isolate the fault elements (unecessary here as there is just the fault, but usefull if you want to isolate one specific fault or if you add a dike for instance).
locked = solution.d.fstatus[mfault] == 0     # fstatus contains the frictional status of each dislocations (it's an output of the run). Here we isolate the one corresponding to the fault.
    # fstatus = -10:    unknonw (can be frictional or not)
    # fstatus =  -2:    not frictional
    # fstatus =  -1:    frictional but not determined (don't know yet if sliding or locked)
    # fstatus =   0:    frictional and locked
    # fstatus =   1:    frictional and sliding

# Hatch the locked elements to highlight them: use the keyword argument 'patches2hatch'
plotFault2D(fault, v=solution.d.norm, vec=solution.d, group=[fault.group_init], cmap=plt.cm.magma, vec_scale=60,
            patches2hatch=locked, vec_color='C0', cbar_title='|u| (m)', vmin=0, vmax=None,vec_width=4e-3,
            center=False, figsize=(7,7))

