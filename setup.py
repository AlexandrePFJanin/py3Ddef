from os.path import join
import numpy as np

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    
    # Get the numpy direction
    nppath=np.__path__[0]
    npinclude = join(nppath, 'core/include/numpy')

    # Create a config object
    config = Configuration('py3Ddef', parent_package, top_path)
    
    # Where are the includes
    include_dirs = ['src/', npinclude]
    
    # Which are the sources
    sources = ['src/3dmain.f', 'src/okada_sub.f', 'src/xyz_output.f']

    # Additional flags
    CFLAGS = []

    # Create an extension
    config.add_extension('_all3Ddef', 
            sources=sources, 
            include_dirs=include_dirs,
            extra_compile_args=CFLAGS)

    # All done
    return config

# Main
if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration, version='1.0.0')