import matplotlib.pyplot as plt
import numpy as np
import py3Ddef



# ================ py3Ddef example ================

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

# Define where to compute the solution outsite of the dislocation

nx = 400 +1

x = np.linspace(-20,20,nx)                  # in km
y = np.zeros(x.shape[0],dtype=np.float64)   # in km
z = np.zeros(x.shape[0],dtype=np.float64)   # in km

# --- Dislocations geometry (i.e. Dike geometry)

xd = np.array([0.0],dtype=np.float64)       # in km
yd = np.array([-50.0],dtype=np.float64)     # in km
zd = np.array([-0.4],dtype=np.float64)      # in km

length = np.array([100],dtype=np.float64)   # in km
width  = np.array([4],dtype=np.float64)     # in km
strike = np.array([0],dtype=np.float64)     # in km
dip    = np.array([90],dtype=np.float64)    # in km

# --- Dislocations deformation

# Units of the input displacements: centimeters.

kode = np.array([10],dtype=np.int32)        # See the doc of 3D-def (Gomberg and Ellis, 1993, User manual)
ss   = np.array([0],dtype=np.float64)       # Strike slip (in centimeters)
ds   = np.array([0],dtype=np.float64)       # Dip slip (in centimeters)
ts   = np.array([100],dtype=np.float64)     # Tensile slip (in centimeters)


# --- Computation of 3D deformation

u,s,e,o,f,g,d = py3Ddef.compute3Ddef(x,y,z,\
                         xd,yd,zd,length,width,strike,dip,\
                         kode,ss,ds,ts,\
                         nu,E,mu)

# Description of the outputs:

# 1. Outputs computed on the input grid:

# u:    Displacement field: vectorial field (type: py3Ddef.GridDisplacement)
# s:    Stress field: symmetric tensor (type: py3Ddef.GridStress)
# o:    Principal strain orientation: ex, px, tx, ey, py, ty, ez, pz, tz (type: py3Ddef.GridPrincipalStrainOrientation)
# f:    Optimal failure planes: strike1, dip1, rake1, strike2, dip2, rake2 (type: py3Ddef.GridOptimalFailurePlane)
# e:    Strain field: symmetric tensor (type: py3Ddef.GridStrain)
# g:    Displacement gradient field: asymmetric tensor (type: py3Ddef.GridDisplacementGradient)

# 2. Outputs computed on the dislocation patches:

# d:    Displacement on each patches given in (strike slip, dip slip, tensil sip) (type: py3Ddef.ElementDisplacement)
#       Can be converted to (x-east, y-north, z-up) using the function d.convert2xyz() and the geometry of the dislocations.

# 3. Units of the different input: See the README of py3Ddef (similar to 3D~def.)

# Unit of the output displacements: centimeters.
u.rescale(1/100) # conversion from cm to m
d.rescale(1/100) # conversion from cm to m

# --- Figure - displacement across the dike

fig = plt.figure(figsize=(10,3))
ax  = fig.add_subplot(111)
ax.set_title('Surface displacements resulting from the opening of a vertical dyke by %s cm'%str(ts[0]))
axb = ax.twinx()
ax.hlines(0,np.amin(x),np.amax(x),alpha=0.5,colors='k',linestyles='--')
ln1 = ax.plot(x,u.z,label=r'$u_z(x,y=0)$')
ln3 = ax.plot(x,u.x,label=r'$u_x(x,y=0)$')
ax.set_xlabel('x-axis (km)')
ax.set_ylabel(r'$u_z(x)$'+ ' (m)')
axb.set_ylabel(r'$u_x(x)$'+ ' (m)')
ax.legend()
fig.tight_layout()
plt.show()
