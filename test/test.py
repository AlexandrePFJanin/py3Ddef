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
E  = 30.e9  # Young's modulus
mu = 0.6    # coeff of internal friction

# --- Grid

nx = 400 +1

x = np.linspace(-20,20,nx)                  # in km
y = np.zeros(x.shape[0],dtype=np.float64)   # in km
z = np.zeros(x.shape[0],dtype=np.float64)   # in km

# --- Background deformation (optional arguments)

# Background field:
# Name           flag     format
# strain:        'stra'   (Exx,    Exy,    Exz,    Eyy,    Eyz,    Ezz,    0,      0,      0)
# stress:        'stre'   (Sxx,    Sxy,    Sxz,    Syy,    Syz,    Szz,    0,      0,      0)
# displacement:  'disp'   (dUx/dx, dUx/dy, dUx/dz, dUy/dx, dUy/dy, dUy/dz, dUz/dx, dUz/dy, dUz/dz)
#    -> N.B. The three last zeros are mandatory for the stress and strain cond.

# Here below an example, not give to the main function:

bg_flag  = 'disp'
bg_field = np.array([-0.5e-4, 0.5e-4, 0.0, -0.5e-4, 0.5e-4, 0.0, 0.0, 0.0, 0.0])

# --- Dislocations geometry (i.e. Dike geometry)

xd = np.array([0.0],dtype=np.float64)       # in km
yd = np.array([-50.0],dtype=np.float64)     # in km
zd = np.array([-0.4],dtype=np.float64)      # in km

length = np.array([100],dtype=np.float64)   # in km
width  = np.array([4],dtype=np.float64)     # in km
strike = np.array([0],dtype=np.float64)     # in km
dip    = np.array([90],dtype=np.float64)    # in km

# --- Dislocations deformation

kode = np.array([10],dtype=np.int32)        # See the doc of 3D-def (Gomberg and Ellis, 1993, User manual)
ss   = np.array([0],dtype=np.float64)       # Strike slip (in meters)
ds   = np.array([0],dtype=np.float64)       # Dip slip (in meters)
ts   = np.array([1],dtype=np.float64)       # Tensile slip (in meters)


# --- Computation of 3D deformation

u,s,d,o,f,e,g = py3Ddef.compute3ddef(x,y,z,\
                         xd,yd,zd,length,width,strike,dip,\
                         kode,ss,ds,ts,\
                         nu,E,mu)


ux = u[:,0]
uy = u[:,1]
uz = u[:,2]

# --- Figure

fig = plt.figure(figsize=(10,3))
ax  = fig.add_subplot(111)
ax.set_title('Surface displacements resulting from the opening of a vertical dyke by %sm'%str(ts[0]))
axb = ax.twinx()
ax.hlines(0,np.amin(x),np.amax(x),alpha=0.5,colors='k',linestyles='--')
ln1 = ax.plot(x,uz,label=r'$u_z(x,y=0)$')
ln3 = ax.plot(x,ux,label=r'$u_x(x,y=0)$')
ax.set_xlabel('x-axis (km)')
ax.set_ylabel(r'$u_z(x)$'+ ' (m)')
axb.set_ylabel(r'$u_x(x)$'+ ' (m)')
ax.legend()
fig.tight_layout()
plt.show()


