from setuptools import setup
from Cython.Build import cythonize


setup(
    ext_modules = cythonize("hoenl_like.pyx"),
    compiler_directives={"language_level": "3"}
)
setup(
    ext_modules = cythonize("legendre_plynomials.pyx"),
    compiler_directives={"language_level": "3"}
)
setup(
    ext_modules = cythonize("constants_and_atomic_properties.pyx"),
    compiler_directives={"language_level": "3"}
)
setup(
    ext_modules = cythonize("matrix_coefficients_v2.pyx"),
    compiler_directives={"language_level": "3"}
)