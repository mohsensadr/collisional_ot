from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension(
        "collision_wrapper",
        ["collision_wrapper.pyx", "collision.c"],
        include_dirs=[numpy.get_include()],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]  # Suppress deprecated API warning
    )
]

setup(
    name="OT_Collision",
    ext_modules=cythonize(extensions)
)
