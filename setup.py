from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
import os

# Base path to the cython backend folder
cython_dir = os.path.join("src", "collisional_ot", "cython")

extensions = [
    Extension(
        "collisional_ot.cython.collision_wrapper",  # fully qualified module name
        [
            os.path.join(cython_dir, "collision_wrapper.pyx"),
            os.path.join(cython_dir, "collision.c"),
        ],
        include_dirs=[numpy.get_include()],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        extra_compile_args=["-O3"],
    )
]

setup(
    name="collisional_ot",
    version="0.1.0",
    description="Collision-based dynamics for optimal transport",
    packages=["collisional_ot"],
    package_dir={"": "src"},
    ext_modules=cythonize(extensions),
)

