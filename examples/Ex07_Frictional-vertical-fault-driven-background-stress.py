import matplotlib.pyplot as plt
import numpy as np
from py3Ddef import DeformationRun, BackgroundDeformation
from py3Ddef.geometry import discreteDislocation, PatchCollection, UniformGrid


# ================ py3Ddef example ================

# ---- Example of a frictional vertical fault driven by a background stress field ----

# Context:  Comparison between a frictional and a frictionless dislocation
#           under the same background stress field promoting strike-slip displacement
#           on them.
#           For the frictional dislocation, the normal stress on the fault is imposed
#           at 10 MPa, cohesion is set to 0 Pa, and the stress drop is parametrise so
#           that if failure, the stress will return to the Mohr-Coulomb envelope.
#           Notice that the effect of fluids and lithostatic pressure is ignored here
#           to clearly visualise the effect of friction.

# =================================================


# --- Medium parameters

nu = 0.25   # Poisson's ratio
E  = 7.e10  # Young's modulus (unit: Pa)
mu = 0.6    # coeff of internal friction


# --- Background deformation

bg = BackgroundDeformation(bg_flag='stre',\
                           bg_field=np.array([0, 1e7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))


# --- Dislocations geometry and boundary conditions

# fritionless vertical strike-slip fault
flessfault = discreteDislocation(x0=0, y0=-100, z0=0, L=200, W=5, dip=90, strike=0, n_strike=1, n_dip=1, \
                                 kode=2, ss=0, ds=0, ts=0)

# fritional vertical strike-slip fault:
# ts = 1e7: imposed normal stress: 10 MPa
# fcode = 1 and sdrop = 1, frictional element and if failure, return to envelope
# C = 0, mu = 0.6, cohesion 0 Pa and coefficient of friction of 0.6
frictfault = discreteDislocation(x0=0, y0=-100, z0=0, L=200, W=5, dip=90, strike=0, n_strike=1, n_dip=1, \
                                 kode=2, ss=0, ds=0, ts=1e7, \
                                 fcode=1, sdrop=1, mu=0.6, C=0)

# --- Figure init
# Initialisation of the comparison figure

fig = plt.figure(figsize=(7,6))
fig.suptitle('Friction on a vertical strike-slip fault', fontsize=14)
ax1 = fig.add_subplot(211)
ax1.set_title('Surface')
ax1.set_ylabel('Stress (Pa)')
ax1.set_xlabel('x (i.e. east, in km)')
ax2 = fig.add_subplot(212)
ax2.set_title('Center of the fault')
ax2.set_ylabel('Stress (Pa)')
ax2.set_xlabel('x (i.e. east, in km)')


# --- Computation of 3D deformation

# iterate on the two faults (without/with friction)
for j in range(2):

    # Add the relevant fault to the patch collection
    patches = PatchCollection()

    if j == 0:
        # Frictionless
        patches.add(flessfault)
        colors = 'C0'
        suffix = ' fless'
    else:
        # Frictional
        patches.add(frictfault)
        colors = 'C1'
        suffix = ' frict'
    
    # --- Grid

    # 1-D grid, along x at a given depth.
    # One at the surface, one at the center of the dislocation
    nx = 1000
    ny = 1
    nz = 1

    # iterate on (surface/center of the fault)
    for i in range(2):

        if i == 0:
            # Surface
            grid = UniformGrid(-20, 20, nx,\
                                0,  0, ny,\
                                0,  0, nz)
        else:
            # Dislocation center
            grid = UniformGrid(-20, 20, nx,\
                                0,  0, ny,\
                            patches.zc[0], patches.zc[0], nz)
            

        # --- Computation of 3D deformation

        # Define your solution model and solve the deformations
        solution = DeformationRun(grid=grid, patches=patches, nu=nu, mu=mu, E=E)
        solution.compute3Ddef(bg=bg)

        # Unit of the output displacements: centimeters.
        solution.displ.rescale(1/100)   # conversion from cm to m
        solution.dislocs.rescale(1/100) # conversion from cm to m

        # --- Figure: Across the system

        if i == 0:
            axi = ax1
        else:
            axi = ax2

        axi.plot(grid.x[:, 0, 0], solution.stress.xy[:, 0, 0], color=colors, linestyle='-',  label=r'$\sigma_{xy}$'+suffix)
        axi.plot(grid.x[:, 0, 0], solution.stress.yz[:, 0, 0], color=colors, linestyle='--', label=r'$\sigma_{yz}$'+suffix)
        
        xlim0, xlim1 = np.amin(grid.x), np.amax(grid.x)
        if j == 1 and i == 1:
            ax2.scatter([0], [frictfault.mu[0]*abs(frictfault.ts[0]+-1*(frictfault.rhoLitho[0]-frictfault.rhoFluid[0])*9.81*frictfault.zc[0]*1000)+frictfault.C[0]], label=r'$\tau_c$', color=colors)
        axi.hlines([0], xlim0, xlim1, linestyles=':', linewidths=1, colors=['k'], alpha=0.5, zorder=-1)
        axi.set_xlim(xlim0, xlim1)

        if j == 1:
            ylim0, ylim1 = axi.get_ylim()
            axi.vlines(patches.x0, ylim0, ylim1, linestyles='--', alpha=0.75, colors=['grey'], label='fault')
            axi.set_ylim(ylim0, ylim1)
        
        axi.legend(fontsize=8, loc='lower right')

fig.tight_layout()
plt.show()


# ------------------
# What you will see:
# ------------------
#
# Far from the dislocation, the stress field is in both case the imposed background stress.
# For the frictionless fault, the stresses are zero at the center of the dislocation (by definition).
# However, notice that it is not the case at the surface.
# For the frictional fault: The failure condition was meet, i.e. the stress in the plane of the fault
# exceeded the normal stress times the friction plus the cohesion (here 0) and the dislocation slides.
# The envelope is at tau_c = 0.6*1e7= 6e6 Pa. The orange dot on the figure displays tau_c and we can 
# see that the stress at the center of the frictional patches sliding is excatly on the envelope.
