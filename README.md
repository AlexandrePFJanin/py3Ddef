# py3Ddef
**py3Ddef** is a python implementation of [3D~def](http://www.ceri.memphis.edu/people/ellis/3ddef/), a three-dimensional boundary-element framework used to compute stresses, strains, and displacements both within and at the surface of an elastic half-space and on each element (*i.e.* planar dislocation).

The strength of **3D~def** lies in the fact that, unlike classical dislocation models, boundary-element approaches treat slip (*i.e.*, the amplitudes of the dislocation components) as unknown variables. These quantities are solved for by minimizing the strain energy of the surrunding medium while satisfying boundary conditions set on either the stress or the displacement on each dislocation surface. By doing this, **3D~def** explicitly represents mechanical interactions between dislocations (faults, dikes, sills etc) and between the dislocation and the background deformation, whereas more simple dislocation models typically neglect these interactions and do not incorporate strain-energy minimization.

**py3Ddef** provides a high-level interface to **3D~def** with full integration into a Python environment. It also includes a library of functions and classes that streamline the creation of dislocation and grid geometries, the management of outputs, and the visualization of results, either directly in Python using **Matplotlib** or through [ParaView](https://www.paraview.org/). In addition, compared with the original version of **3D~def**, **py3Ddef** incorporates several additional boundary conditions.


#### Acknoledgements:
[3D~def](http://www.ceri.memphis.edu/people/ellis/3ddef/) is the product of the work of J. Gomberg and M. Ellis. 

If you use **py3Ddef** or **3D~def** in your research, we kindly ask that you cite the following references, acknowledging the work of J. Gomberg and M. Ellis:

 - Gomberg, J. S., & Ellis, M. (1993). 3D-DEF; a user's manual (a three-dimensional, boundary element modeling program) (No. 93-547). US Geological Survey,.
 - Gomberg, J., & Ellis, M. (1994). Topography and tectonics of the central New Madrid seismic zone: Results of numerical experiments using a three‐dimensional boundary element program. Journal of Geophysical Research: Solid Earth, 99(B10), 20299-20310.


## Table of contents:

- [Introduction](#py3ddef)
- [Install py3Ddef](#install-py3ddef)
    - [Download](#1-download)
    - [Compilation with f2py and Meson](#2-compilation-with-f2py-and-meson)
    - [Linking](#3-linking)
    - [Test](#4-test)
- [Uninstall](#uninstall-py3ddef)
- [Use py3Ddef](#use-py3ddef)
    - [How does it work?](#how-does-it-work)
        - [Introduction](#introduction)
        - [Solving the deformation](#solving-the-deformation)
        - [Geometry and coordinate systems](#geometries-and-coordinate-systems)
        - [Boundary conditions](#boundary-conditions)
        - [Background deformation](#background-deformation)
        - [Friction](#friction)
        - [Units](#units)
        - [Displacement and stress conventions](#displacement-and-stress-conventions)
    - [How to use py3Ddef?](#how-to-use-py3ddef)
        - [Importation and main function](#importation-and-main-function)
        - [Dislocations](#dislocations)
        - [Computational grid](#computation-grid)
        - [Visualisation](#visualisation)
    - [Examples](#examples)
- [Cite py3Ddef](#cite-py3ddef)
- [References](#references)


## Install py3Ddef:

In this section, we describe how to install **py3Ddef** from the source code available on this page. Please follow the steps below sequentially.

> **_NOTE:_** Although not detailed here, as a best practice we strongly recommend installing **py3Ddef** within a dedicated Python environment, for example using [conda](https://conda.io/projects/conda/en/latest/index.html).

### 1. Download

Download **py3Ddef** from this page and unzip the archive. Then, move into the main directory of the package:
```
cd py3Ddef-main/
```

### 2. Compilation with f2py and Meson

The source code of **py3Ddef** are written in fortran and need to be compile first. To do so, **py3Ddef** is build to be used with [**numpy f2py**](https://numpy.org/doc/stable/f2py/f2py-user.html) and [**Meson**](https://numpy.org/doc/stable/f2py/buildtools/meson.html).

First, place you in your environement where you want to install **py3Ddef**.

Then, make sure ninja and meson are installed:
```
meson --version
ninja --version
```

Or install them with [pip](https://pip.pypa.io/en/stable/):
```
pip install meson meson-python ninja
```

Then, compile the source code with [**numpy f2py**](https://numpy.org/doc/stable/f2py/f2py-user.html) and [**Meson**](https://numpy.org/doc/stable/f2py/buildtools/meson.html) using the following line in the project directory:

```bash
python -m numpy.f2py --backend meson -c \
    src/3dmain.f src/okada_sub.f src/xyz_output.f \
    -I "$(pwd)/src" \
    -m _all3Ddef
```

It will generate a `.so` file in the project directory. Then,

```bash
rm -f py3Ddef/_all3Ddef*.so
mv _all3Ddef*.so py3Ddef/
```

### 3. Linking

To link in a user module directory, use [pip](https://pip.pypa.io/en/stable/) and run 
```
python -m pip install -e . --no-build-isolation
```

### 4. Test

To test that the installation process went well, you can run one of the provided examples:
```
python examples/Ex01_Dike-induced-faulting.py
```
Note: If you installed **py3Ddef** in a specific environement, make sure to have it active before running python.

## Uninstall py3Ddef:

As you used [pip](https://pip.pypa.io/en/stable/) for the installation, use it to uninstall the package. In a terminal, run:
```
python -m pip uninstall py3Ddef
```

$\rightsquigarrow$ [Back to the Table of Content](#table-of-contents)

## Use py3Ddef:

### How does it work?

#### Introduction

**py3Ddef** uses the three-dimensional boundary-element model **3D~def** allowing computing stresses, strains, and displacements within and on the surface of an elastic half-space (no bottom to the model at depth) as well as on (across) the dislocation themselves. The strength of this model is that, unlike more simple implementations of Okada's (1992) solution, the dislocations representing discontinuities in the elastic medium can be either displacement or stress discontinuities, not just displacement. 

In **3D~def**, a **dislocation represents a surface of constant slip**. Displacements are discontinuous across the dislocations but stresses are continuous everywhere (across the dislocation and everywhere in the grided medium).

**3D~def** solve the displacement across the dislocations by minimizing the strain energy in the surrounding medium while satisfying the boundary conditions in either stress or displacement applied on each dislocation surface. Doing this, the model explicitly accounts for the interactions between dislocations and with the background deformation.

$\rightsquigarrow$ [Back to the Table of Content](#table-of-contents)

#### Solving the deformation

The deformation across the dislocations is resolved as a set of linear equations. For each dislocation number $i$, the displacement **across** the dislocation (strike slip $D^s_i$, dip slip $D^d_i$, tensile slip $D^n_i$) and stress (shear stresses in the along-strike direction $\tau^s_i$, in the along-dip direction $\tau^d_i$, and the stress normal to the dislocation $\sigma^n_i$) is specified at it's center.

> **_NOTE:_** The superscripts $^s$, $^d$ and $^n$ always refer to along-strike, along-did and along-normal directions, respectively.

The set of equations solved take the following form for $n$ dislocations:

```math
\left[ \begin{matrix}
\tau^s_1 \\
\tau^d_1 \\
\sigma^n_1 \\
\cdots\\
\tau^s_n \\
\tau^d_n \\
\sigma^n_n
\end{matrix} \right] = 
\left[ \begin{matrix}
A^{ss}_{11}, A^{sd}_{11}, A^{sn}_{11}, \cdots, A^{ss}_{1n}, A^{sd}_{1n}, A^{sn}_{1n}\\
A^{ds}_{11}, A^{dd}_{11}, A^{dn}_{11}, \cdots, A^{ds}_{1n}, A^{dd}_{1n}, A^{dn}_{1n}\\
A^{ns}_{11}, A^{nd}_{11}, A^{nn}_{11}, \cdots, A^{ns}_{1n}, A^{nd}_{1n}, A^{nn}_{1n}\\
\cdots \\
A^{ss}_{n1}, A^{sd}_{n1}, A^{sn}_{n1}, \cdots, A^{ss}_{nn}, A^{sd}_{nn}, A^{sn}_{nn}\\
A^{ds}_{n1}, A^{dd}_{n1}, A^{dn}_{n1}, \cdots, A^{ds}_{nn}, A^{dd}_{nn}, A^{dn}_{nn}\\
A^{ns}_{n1}, A^{nd}_{n1}, A^{nn}_{n1}, \cdots, A^{ns}_{nn}, A^{nd}_{nn}, A^{nn}_{nn}\\
\end{matrix} \right] 
\left[ \begin{matrix}
D^s_1 \\
D^d_1 \\
D^n_1 \\
\cdots\\
D^s_n \\
D^d_n \\
D^n_n
\end{matrix} \right]
```

The matrix $A$ contains the Green’s functions (computed using the Okada, 1992, routines). The displacement vector $D$ **across** each dislocation is obtained by inverting $A$. The resulting displacements are then used to compute analytically the deformation throughout the surrounding medium.

> **_NOTE:_** A numerical error may occur if a point of the computational grid (used to compute the solution in the surrounding medium) lies exactly on a dislocation.

$\rightsquigarrow$ [Back to the Table of Content](#table-of-contents)

#### Geometries and coordinate systems

**3D~def** and **py3Ddef** use the cartesian **x-(East), y-(North), z-(Up)** system to describe the coordinates and the **along-strike, along-dip, along-normal** to describe displacements and stresses in addition to the cartesian system. Both coordinates system verified the right-hand convention.

> **_NOTE:_**  When looking in the strike direction, the hanging wall is on the right.

Dislocations are represented as rectangular patches, each defined by a width (along dip), a length (along strike), an origin ($X_0, Y_0, Z_0$), a strike (azimuth measured clockwise from north), and a dip (angle measured to the right of the strike direction; 0° = horizontal, 90° = vertical).

> **_NOTE:_**  A main difference between **py3Ddef** and **3D~def** is that **py3Ddef** does not use *sub-elements* (*i.e.* sub-dislocations), but only *elements* (*i.e* dislocations), in order to provide a more homogeneous and lower-level definition of the dislocation geometry. If you want a result with more elements, implement them as dislocations.

**3D~def** and **py3Ddef** describe the displacements at a given dislocation as either **relative** displacements (dispalcements across the dislocation: the total slip in every direction) or **absolute** displacements (displacements of one side of the fault relative to its original position and not relative to the other side of the fault). Relative displacements if the $i^{th}$ dislocation are written $D_i^k$ whereas the absolute displacements are written $u_i^k$ (notation used here and by J. Gomberg and M. Ellis in the [official documentation of 3D~def](http://www.ceri.memphis.edu/people/ellis/3ddef/)). The relation between the absolute and relative displacement can be simply express by:

```math
D_i^s = u_{i_-}^s - u_{i_+}^s \\
D_i^d = u_{i_-}^d - u_{i_+}^d \\
D_i^n = u_{i_-}^n - u_{i_+}^n
```
The subscripts $-$ and $+$ refer to the absolute displacement of the **footwall** and **hangingwall**, respectively..

The distinction between absolute and relative displacements is crucial.

$\rightsquigarrow$ [Back to the Table of Content](#table-of-contents)

#### Boundary conditions

Discontinuities are user-defined for each dislocation in the reference frame of dislocation (along-strike, along-dip, along-normal). The following table summarizes the different options available. Note that several options have been added from the original version of **3D~def**. For readability, the superscipts below and subscripts are inverted. Superscipts describe the absolute motion of the walls ($^-$ for the footwall and $^+$ for the hangingwall) and the subscripts the directions.

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

The difference between absolute displacement and relative displacement is crucial here. **A confusion will lead to a wrong result.**

For *code=11*, the user specifies a stress magnitude $\tau(\phi)$ and the direction of shear $\phi$ as an angle measured from the strike direction (in *deg*). Then, the code resolves these conditions into shear stress conditions in the strike and dip directions.

$\rightsquigarrow$ [Back to the Table of Content](#table-of-contents)

#### Background deformation

**3D~def** allows the users to specify a background deformation field that will act as a driver (or as an *additional* driver) of deformation. The background deformation field can be specified as either a stress, strain, or displacement gradient tensor (internally these tensors are converted to a stress tensor).

The background stresses are subtracted from the stress boundary conditions before solving for the relative displacements on the dislocations. Because the background deformation is added to that associated with the dislocation displacements, the resultant stresses on the dislocations are thus guaranteed to equal the boundary condition stress value. 

> **_NOTE:_**  Use of a background stress or strain tensor will have no affect on the dislocation relative displacements if you use the boundary conditions codes 1 or 10 for your dislocation because these boundary conditions do not involve stresses.

If a displacement gradient tensor is specified, then the only affect will be to modify the output displacements by any rigid body rotations associated with the background displacement gradient field.

$\rightsquigarrow$ [Back to the Table of Content](#table-of-contents)

#### Friction

The coefficient of internal friction is only used to determine potential failure planes, assuming a Coulomb failure criterion.

Note that the stress boundary conditions placed on the central element via code 12 are ending conditions. That is, the stresses take the zero values at the end of the model run, and displacements that occur across the element do so in order to preserve the zero stress conditions. In other words, the element behaves like a frictionless fault. To simulate friction on the element, you would need to impose a non-zero level of shear stress opposite in sign to that imposed by the driving stresses. This is leads to complications of having to assume a priori slip directions, which we will discuss elsewhere.

$\rightsquigarrow$ [Back to the Table of Content](#table-of-contents)

#### Units

**py3Ddef** and **3D~def** use the same units. Units are defined in the table below:

| Features | Units |
| :---: | :---: |
| Dislocation geometry | km |
| Computational grid geometry | km |
| Dislocation displacements | cm |
| Stresses | same as the input Young's modulus |
| Strain | - |
| Angles | degree |

$\rightsquigarrow$ [Back to the Table of Content](#table-of-contents)

#### Displacement and stress conventions

Conventions on the displacement:

| Component | Sign | Displacement |
| :---: | :---: | :---: |
| Strike | $>0$ | left-lateral |
| Strike | $<0$ | right-lateral |
| Dip | $>0$ | thrust |
| Dip | $<0$ | normal |
| Normal (tensile) | $>0$ | opening |
| Normal (tensile) | $<0$ | closing |

Conventions on the stress:

| Sign | Syle |
| :---: | :---: |
| $<0$ | Compression |
| $>0$ | Extension |

$\rightsquigarrow$ [Back to the Table of Content](#table-of-contents)




### How to use py3Ddef?

In this section we will detail how to use the main features of **py3Ddef**. A complete set of commented examples using these features are provided in the directory [examples/](./examples/). Examples are detailed in the section [Examples](#examples).

#### Importation and main function

When installed and linked to your python environnement as shown above in the section [Install py3Ddef](#install-py3ddef), **py3Ddef** can be imported with:

```python
from py3Ddef
```

**py3Ddef** is structured to have a main function `py3Ddef.compute3Ddef` to solve the stresses, strains, and displacements from a set of initial conditions (a population of dislocation, a set of points in the surrounding elastic medium). 

```python
from py3Ddef import compute3Ddef
```

The function `py3Ddef.compute3Ddef` takes in input all the information about the geometry of the points where the solution should be computed and the geometry and boundary condition on each dislocation and returns the result as a tuple of python objects, one per type of output. The call of the function will look like:

```python
# Example without background deformation

u,s,e,o,f,g,d = compute3ddef(x,y,z,\
                        xd,yd,zd,length,width,strike,dip,\
                        kode,ss,ds,ts,\
                        nu,E,mu)
```

where the input variable `xd,yd,zd,length,width,strike,dip` define the geometry of each patch (in *deg* for the strike and the dip and in *km* for the others), `kode,ss,ds,ts` define the type of boundary conditions (see next paragraph) applied on each patch (*i.e.* slip in the strike direction, stress in the dip direction, etc.), `nu` is the Poisson's ratio. `x`, `y`, `z` (in *km*) are the coordiantes of the set of points where the analytical solution is computed
The outputs computed by the functions are: `u` the displacement field (in *cm*) on the input grid, `s` the stress (unit of the Young's modulus `E`, *e.g. Pa*), `e` the strain, `o` orientations of principal strains , `f` the optimal failure planes (in *deg*), `g` the displacement gradients and `d` the relative displacements ($D_i^k$ with our notations) on each dislocation (in *cm*).

Each of the output of `py3Ddef.compute3Ddef` has a specific type whose the properties are listed in the table below.

Results on the computational grid:

| Variable | Type | Numerical Type  | Fields | Comments |
| :---: | :---: | :---: | :---: | :---: |
| `u` | Vectorial | `py3Ddef.GridDisplacement` | `.x`, `.y`, `.z` |  |
| `s` | Tensorial (symmetric) | `py3Ddef.GridStress` | `.xx`, `.xy`, `.xz`, `.yx`, `.yy`, `.zz`, `.zx`, `.zy`, `.zz`|  |
| `e` | Tensorial (asymmetric) | `py3Ddef.GridStrain` | `.xx`, `.xy`, `.xz`, `.yx`, `.yy`, `.zz`, `.zx`, `.zy`, `.zz`|  |
| `o` | Tensorial (asymmetric) | `py3Ddef.GridPrincipalStrainOrientation` | `.ex`, `.px`, `.tx`, `.ey`, `.py`, `.ty`, `.ez`, `.pz`, `.tz`| magnitude of the max. principal strain (`e`), plunge with respect to horizontal of the max. principal strain (`p`), trend (clockwise with respect to North) of the max. principal strain (`t`)|
| `f` | Vectorial (double) | `py3Ddef.GridOptimalFailurePlane` | `.str1`, `.dip1`, `.rak1`, `.str2`, `.dip2`, `.rak2`| strike (`str`), dip (`dip`) and rake (`rak`) for the two possible failure planes labelled 1 and 2 |
| `g` | Tensorial (asymmetric) | `py3Ddef.GridDisplacementGradient` | `.xx`, `.xy`, `.xz`, `.yx`, `.yy`, `.zz`, `.zx`, `.zy`, `.zz`|  |
| `u` | Vectorial | `py3Ddef.GridDisplacement` | `.x`, `.y`, `.z` |  |

Results on the dislocations:

| Variable | Type | Numerical Type  | Fields | Comments |
| :---: | :---: | :---: | :---: | :---: |
| `d` | Vectorial | `py3Ddef.ElementDisplacement` | `.pos`, `.ss`, `.ds`, `.ts`, `x`, `.y`, `.z` | Position of the center of each dislocation ()`.pos`), net displacement across the dislocation in the along-strike, along-dip and along-normal reference frame (`.ss`, `.ds`, `.ts`) and in the xyz coordinate system (`x`, `.y`, `.z`, not computed by default: need to run the internal function `.convert2xyz`)|


> **_EXAMPLE:_** Check the example [Ex01_Dike-induced-faulting.py](./examples/Ex01_Dike-induced-faulting.py) to see a first basic example of how to use `py3Ddef.compute3Ddef`. The example is about a vertical dyke opening at depth, based on the work of Rubin & Pollard (1988)

$\rightsquigarrow$ [Back to the Table of Content](#table-of-contents)

#### Dislocations

**py3Ddef** offers a comprehensive framework for defining and managing dislocation geometries and boundary conditions. The `py3Ddef.geometry.Patch` object defines a unit dislocation, while `py3Ddef.geometry.PatchCollection` defines a collection of unit dislocations and provides tools for efficiently creating a set of dislocations joined (no gap between dislocations) in the same plane with the function `py3Ddef.geometry.discreteDislocation` returning an instance of `py3Ddef.geometry.PatchCollection`.

```python
from py3Ddef.geometry import discreteDislocation

# Creation of a set of 100 (n_strike=10, n_dip=10) dislocations representing a vertical (dip=90) north-south (strike=0) dike extending from the point (x=0,y=0,z=0) to 3.5 km depth (W), 5 km long (L) and opening by 1 meters (kode=10, ss=0, ds=0, ts=100).

dike = discreteDislocation(x0=0, y0=0, z0=0, L=5, W=3.5, dip=90, strike=0, n_strike=10, n_dip=10, \
                           kode=10, ss=0, ds=0, ts=100)
```

The boundary conditions on the set of dislocation is enter throught the input variables `kode` for the code (*i.e.* the index) of the boundary condition (see the section [Boundary conditions](#boundary-conditions)) and `ss`, `ds` and `ts` for each components (here, no strike slip, no dip slip and 1 meter of tensile opening).

Working with `py3Ddef.geometry.PatchCollection` objects also allows to easily pass to `py3Ddef.compute3Ddef` all the necessary information about the set of dislocations in the correct order and with the correct format. For that you can use `*dike.get()` (`*` unpack command) directly when calling `py3Ddef.compute3Ddef`.

> **_EXAMPLE:_** Check the example [Ex02_Vertical-fault-driven-displacement-gradient.py](./examples/Ex02_Vertical-fault-driven-displacement-gradient.py) to see an example of how to use `py3Ddef.geometry.discreteDislocation` to generate dislocation geometries and [Ex03_Two-vertical-faults-driven-by-displacement-gradient.py](./examples/Ex03_Two-vertical-faults-driven-by-displacement-gradient.py) to see how to combine several geometries with functions of `py3Ddef.geometry.PatchCollection` objects.

$\rightsquigarrow$ [Back to the Table of Content](#table-of-contents)


#### Computational grid

Similarly to dislocations, **py3Ddef** offers a comprehensive framework for defining and managing computational grid geometries for solving the deformation in the surrounding elastic medium through the class `py3Ddef.geometry.UniformGrid`.

```python
# Creation of a uniform grid of shape (50,60,1) extending from -150 km to 150 km in both the x (east) and y (north) direction and at just 1 depth z=0 km, the surface.

nx = 50
ny = 60
nz = 1
grid = UniformGrid(-150, 150, nx,\
                   -150, 150, ny,\
                      0, 0, nz)
```

Coordinates are then available at `grid.x`, `grid.y` and `grid.z`. Using `py3Ddef.geometry.UniformGrid` objects make then easier the management of the grid and the visualisation of fields computed on it.

Similarly to dislocation sets, `py3Ddef.geometry.UniformGrid` objects also allows to easily pass to `py3Ddef.compute3Ddef` all the necessary information about the grid geometry in the correct order and with the correct format. For that you can use `*grid.get()` (`*` unpack command) directly when calling `py3Ddef.compute3Ddef`.

> **_EXAMPLE:_** Check the example [Ex04_Vertical-dike-and-dipping-fault.py](./examples/Ex04_Vertical-dike-and-dipping-fault.py) to see an example of how to create and use a uniform grid with `py3Ddef.geometry.UniformGrid`.

$\rightsquigarrow$ [Back to the Table of Content](#table-of-contents)


#### Visualisation

The high level description of **py3Ddef** facilitates the management and visualization of inputs and outputs. This offers flexibility to plot the results with your favorite python visualisation package ([matplotlib](https://matplotlib.org/), [seaborn](https://seaborn.pydata.org/) ...). In addition, **py3Ddef** also offers direct options for visualizing the inputs and results of `py3Ddef.compute3Ddef`.

The function `py3Ddef.viewer.plotFault2D` allows to project a set of dislocations and optionally, a vectorial and a scalar field carried on the dislocations in a 2D plane. The projection direction can be set either as a normal vector to the projection plane (with the input argument `view_vec`) or a normal to a given dislocation.

> **_EXAMPLE:_** Check the example [Ex02_Vertical-fault-driven-displacement-gradient.py](./examples/Ex02_Vertical-fault-driven-displacement-gradient.py) to see an example of how to use `py3Ddef.viewer.plotFault2D`.

The internal function `.plot3D()` of a `py3Ddef.geometry.PatchCollection` object enable a quick and easy visualisation of the 3D geometry of a collection of dislocation (for instance, to rapidly check if the input geometry was correct, *e.g.* depth should be negative, positive values of z are for elevation.)

> **_EXAMPLE:_** Check the example [Ex04_Vertical-dike-and-dipping-fault.py](./examples/Ex04_Vertical-dike-and-dipping-fault.py) to see an example of how to use `.plot3D()`.

You can link output fields of `py3Ddef.compute3Ddef` to either the input grid object (`py3Ddef.geometry.UniformGrid`) or to the input dislocation collection object (`py3Ddef.geometry.PatchCollection`). This link will ensure the compatibility in the shapes and formats between the outputs fields and the host geometries and will make easy to export both in [ParaView](https://www.paraview.org/).

> **_EXAMPLE:_** Check the example [Ex05_Vertical-dike-and-dipping-fault-3D-visualisation.py](./examples/Ex05_Vertical-dike-and-dipping-fault-3D-visualisation.py) to see how to automatically export **py3Ddef** fields and geometries in a format (XDMF and HDF5) readable by **Paraview**. Note that this examples comes with the Paraview files already computed (.xdmf and .h5 files (`Ex05_results*`) available [here](./examples/)) so that you can load them directly and see if a visualisation of your data with Paraview can be usefull for you.


$\rightsquigarrow$ [Back to the Table of Content](#table-of-contents)


### Examples

The following table describes the examples provided with **py3Ddef**. They highlight different features of **py3Ddef** and were built to be explored sequentially, in the order shown.

| Examples | Descriptions |
| :---: | :---: |
| [Ex01_Dike-induced-faulting.py](./examples/Ex01_Dike-induced-faulting.py) | Dike-induced displacement field. From Rubin & Pollard (1988). Shows basic notions of how to use `py3Ddef.compute3Ddef` |
| [Ex02_Vertical-fault-driven-displacement-gradient.py](./examples/Ex02_Vertical-fault-driven-displacement-gradient.py) | Based on the Example 1 of **3D~def**: A vertical strike-slip fault driven by dextral simple shear. Shows how to drive deformation with background deformations and how to efficiently create a set of dislocations with `py3Ddef.geometry.discreteDislocation` . Shows how to use the visualisation function `py3Ddef.viewer.plotFault2D`. |
| [Ex03_Two-vertical-faults-driven-by-displacement-gradient.py](./examples/Ex03_Two-vertical-faults-driven-by-displacement-gradient.py) | Similar context as Ex02 but now with several dislocations. Shows how to combine sets of dislocations with the object `py3Ddef.geometry.PatchCollection` and how to visualise them using the notion of `group` of dislocation. |
| [Ex04_Vertical-dike-and-dipping-fault.py](./examples/Ex04_Vertical-dike-and-dipping-fault.py) | Original example of a vertical dike with injection triggering motion on a nearby dipping fault. Shows how to efficiently create a uniform grid for the computation of the solution on the surrounding elastic medium using `py3Ddef.geometry.UniformGrid`. Shows how to visualise in 3D the input dislocation geometry with the internal function `.plot3D` of a `py3Ddef.geometry.PatchCollection` object. |
| [Ex05_Vertical-dike-and-dipping-fault-3D-visualisation.py](./examples/Ex05_Vertical-dike-and-dipping-fault-3D-visualisation.py) | Similar context as Ex04. Shows how to easily visualise the results in [ParaView](https://www.paraview.org/) using the internal functions `.export2paraview` for the computational grid and the set of dislocations. Shows how to link output fields of `py3Ddef.compute3Ddef` to the input computational grid and the input dislocation set with the internal functions `.link` of the relevant objects.|

$\rightsquigarrow$ [Back to the Table of Content](#table-of-contents)

## Cite py3Ddef
[3D~def](http://www.ceri.memphis.edu/people/ellis/3ddef/) is the product of the work of J. Gomberg and M. Ellis. 

If you use **py3Ddef** or **3D~def** in your research, we kindly ask that you cite the following references, acknowledging the work of J. Gomberg and M. Ellis:

 - Gomberg, J. S., & Ellis, M. (1993). 3D-DEF; a user's manual (a three-dimensional, boundary element modeling program) (No. 93-547). US Geological Survey,.
 - Gomberg, J., & Ellis, M. (1994). Topography and tectonics of the central New Madrid seismic zone: Results of numerical experiments using a three‐dimensional boundary element program. Journal of Geophysical Research: Solid Earth, 99(B10), 20299-20310.


$\rightsquigarrow$ [Back to the Table of Content](#table-of-contents)

## References


 - Gomberg, J. S., & Ellis, M. (1993). 3D-DEF; a user's manual (a three-dimensional, boundary element modeling program) (No. 93-547). US Geological Survey,.
 - Gomberg, J., & Ellis, M. (1994). Topography and tectonics of the central New Madrid seismic zone: Results of numerical experiments using a three‐dimensional boundary element program. Journal of Geophysical Research: Solid Earth, 99(B10), 20299-20310.
 - Okada, Y. (1992). Internal deformation due to shear and tensile faults in a half-space. Bulletin of the seismological society of America, 82(2), 1018-1040.
 - Rubin, A. M., & Pollard, D. D. (1988). Dike-induced faulting in rift zones of Iceland and Afar. Geology, 16(5), 413-417.


$\rightsquigarrow$ [Back to the Table of Content](#table-of-contents)