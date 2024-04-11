from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup
import platform

__version__ = '0.1.0'

libraries = []
extra_objects = []
extra_compile_args = ["-O3", "-ffast-math", "-march=native"]
extra_link_args = []
if platform.system() == 'Darwin':
    extra_compile_args.insert(-1, "-Xpreprocessor -fopenmp")
    extra_link_args.insert(-1, "-lomp")
else:
    extra_compile_args.insert(-1, "-fopenmp")
    extra_link_args.insert(-1, "-fopenmp")

ext_modules = [
    Pybind11Extension(
        'snnlib',
        ['python_bindings/bindings.cpp'],
        cxx_std='2a',
        language='c++',
        libraries=libraries,
        extra_objects=extra_objects,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
    ),
]

setup(
    name='snnlib',
    version=__version__,
    author = "David C. Anastasiu",
    author_email = "danastasiu@scu.edu",
    description = "A library for nearest neighbor search and graph constuction in sparse data.",
    long_description="""SNNLib""",
    download_url= "https://github.com/davidanastasiu/snnlib",
    ext_modules=ext_modules,
    python_requires = ">=3.10",
    cmdclass={'build_ext': build_ext},
    zip_safe=False,
)