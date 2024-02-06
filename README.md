# py3Ddef
Implementation in Python of the three-dimensional, boundary element modeling code 3D~def (http://www.ceri.memphis.edu/people/ellis/3ddef/)


## Citation
This is a python implementation of the solution Gomberg and Ellis (1993). Please cite:

Gomberg, J. S., & Ellis, M. (1993). 3D-DEF; a user's manual (a three-dimensional, boundary element modeling program) (No. 93-547). US Geological Survey,.

Gomberg, J., & Ellis, M. (1994). Topography and tectonics of the central New Madrid seismic zone: Results of numerical experiments using a three‐dimensional boundary element program. Journal of Geophysical Research: Solid Earth, 99(B10), 20299-20310.



## Install py3Ddef:

### Compilation
```
python setup.py build
```

### Linking

Link in a user module directory:
```
python setup.py install --user
```
On some OS, it is sometimes better to run :
```
python setup.py install --user --prefix=
```

## Use py3Ddef:

**3D-DEF** and **py3Ddef** are a three-dimensional boundary-element models allowing computing stresses, strains, and displacements within and on the surface of an elastic half-space (no bottom to the model). The strength of this code is that, unlike the standard implementation of Okada's (1992) solution, the patches representing discontinuities in the elastic medium can be either displacement or stress discontinuities, not just displacement.

### The main function

**py3Ddef** contains a main function `py3Ddef.compute3ddef` solving the stresses, strains, and displacements for the input boundary condictions on each patches.

```python
import py3Ddef
```

The following example shows you how to use it and get the following output fields computed on the input computation grid (defined by `x`, `y`, `y` in *km*): the displacement `u` (in *cm*), the stress `s` (unit of the Young's modulus `E`, *e.g. Pa*), the strain `d`, orientations of principal strains `o`, optimal failure planes `f` (in *deg*), the relative displacements on each patches `e` (in *cm*) and the displacement gradients `g`.

```python
u,s,d,o,f,e,g = py3Ddef.compute3ddef(x,y,z,\
                        xd,yd,zd,length,width,strike,dip,\
                        kode,ss,ds,ts,\
                        nu,E,mu)
```

In the overlying example, `xd,yd,zd,length,width,strike,dip` define the geometry of each patch (in *deg* for the strike and the dip and in *km* for the others), `kode,ss,ds,ts` define the type of boundary conditions (see next paragraph) applied on each patch (*i.e.* slip in the strike direction, stress in the dip direction, etc.), `nu` is the Poisson's ratio.

### Boundary Condition types by Components

Discontinuities are user-defined for each patch in the reference frame of patches (along the strike, in the dip direction, along the normal of the patch, i.e. tensile direction). The following table summarizes the different options. Note that several options have been added from the original version of 3D-DEF.

| Code | Strike Dir. | Dip Dir. | Normal Dir.|
| :---: | :---: | :---: | :---: |
| 1  | $$u_{s}^{-}$$ | $$u_{d}^{-}$$ | $$u_{n}^{-}$$ |
| 2  | $$\tau_{s}$$ | $$\tau_{d}$$ | $$\sigma_{n}$$ |
| 3  | $$\tau_{s}$$ | $$u_{d}^{-}$$ | $$\sigma_{n}$$ |
| 4  | $$u_{s}^{-}$$ | $$\tau_{d}$$ | $$\sigma_{n}$$ |
| 5  | $$\tau_{s}$$ | $$\tau_{d}$$ | $$u_{n}^{-}$$ |
| 6  | $$u_{s}^{-}$$ | $$u_{d}^{-}$$ | $$\sigma_{n}$$ |
| 10 | $$D_{s}$$ | $$D_{d}$$ | $$D_{n}$$ |
| 11 | $$\phi$$ | $$\tau(\phi)$$ | $$D_{n}$$ |
| 12 | $$\tau_{s}$$ | $$\tau_{d}$$ | $$D_{n}$$ |
| 13 | $$\tau_{s}$$ | $$D_{d}$$ | $$D_{n}$$ |
| 14 | $$D_{s}$$ | $$\tau_{d}$$ | $$D_{n}$$ |
| 15 | $$D_{s}$$ | $$D_{d}$$ | $$\sigma_{n}$$ |


where $u_{s}^{-}$, $u_{d}^{-}$ and $u_{n}^{-}$ represent the absolute displacements (the magnitude of the motion of the footwall) in the strike, dip and normal directions, respectively. $D_{s}$, $D_{d}$ and $D_{n}$ represent the relative displacement in the strike, dip and normal directions, respectively. $\tau_{s}$, $\tau_{d}$ and $\sigma_{n}$ are the stresses defined at the center of the patch in the strike, dip and normal directions, respectively. For *code=11*, the user specifies a stress magnitude $\tau(\phi)$ and the direction of shear $\phi$ as an angle measured from the strike direction (in *deg*). Then, the code resolves these conditions into shear stress conditions in the strike and dip directions.

Absolute displacement refers to the motion of the footwall side of the patch with respect to the global coordinates (surrounding medium).
Relative  displacement is the net dislocation or slip of one side of the patch with respect to the other side (similar to the fault offset observed after an earthquake).

Displacements, both relative and absolute are in *cm* and the stresses have the same unit as the Young's modulus.
