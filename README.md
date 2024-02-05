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

### Boundary Condition types by Components

|Code|Strike|Dip|Normal|


| Code | Strike | Dip | Normal |
| :---: | :---: | :---: | :---: |
| 1 | $$u_{s}^{-}$$ | $$u_{d}^{-}$$ | $$u_{n}^{-}$$ |
