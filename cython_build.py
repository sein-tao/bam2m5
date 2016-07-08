from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import pysam

import sys
if len(sys.argv) -1 == 0:
    sys.argv.extend(['build_ext', '--inplace'])

extension = Extension("cbam2m5", ["cbam2m5.pyx",],
            include_dirs=pysam.get_include(),
            extra_compile_args=["-w",],
            # language='c++',
            )
setup(
        name = 'cbam2m5',
        ext_modules = cythonize(extension)
)

