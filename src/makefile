#-------------------------------------------------------------
# MAKEFILE FOR THE SOURCES OF PY3DDEF
#-------------------------------------------------------------

# A.JANIN 15.02.2026

# NOTE: the compilation of source is integrated to the
# installation of py3Ddef with Meson + f2py (see README.txt).
# There is no need to compile the src files seperatly. However,
# if you want to modify the sources, you may want to test their
# compilation. This makefile is provided for that.

# Clean the src directory
clean :
	rm -f *.o *.mod

# Compile the source separetly
compile :
	gfortran -c global_inputs.f90
	gfortran -c global_arrays.f90
	gfortran -c global_outputs.f90
	gfortran -c 3dmain.f
	gfortran -c okada_sub.f
	gfortran -c xyz_output.f


