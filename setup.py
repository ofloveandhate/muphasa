import os
import platform
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from Cython.Build import cythonize
import numpy



# Compiler flags
extra_compile_args = ["-O3", "-ffast-math", "-D_hypot=hypot"]
extra_link_args = []

system = platform.system()
if system == "Windows":
    extra_compile_args.append("/std:c++14")
elif system == "Darwin":
    extra_compile_args.extend(["-std=c++11", "-mmacosx-version-min=10.9"])
    extra_link_args.extend(["-stdlib=libc++", "-mmacosx-version-min=10.9"])
else:
    extra_compile_args.append("-std=c++11")

ext_modules = [
    Extension(
        "mph.pyMPH",
        sources=[
            "src/python/mph/pyMPH.pyx",
            "src/cpp/main.cpp",
            "src/cpp/grade.cpp",
            "src/cpp/IO.cpp",
            "src/cpp/matrix.cpp",
            "src/cpp/utils.cpp",
            "src/cpp/presentation.cpp"
        ],
        include_dirs=[
            "src/python/mph",   # for .pxd files
            "src/cpp",      # for C++ headers
            numpy.get_include(),
            os.environ.get("BOOST_INC_DIR", "/opt/homebrew/include" if system == "Darwin" else "/usr/local/include")
        ],
        language="c++",
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args
    )
]

setup(
    ext_modules=cythonize(ext_modules, language_level=3),
    packages = find_packages(where="src/python"),
    package_dir = {"": "src/python"},
    # package_dir={"":"src/python"},
    zip_safe=False
)
