# 2022-December
# https://github.com/ludlows

from setuptools import find_packages, setup, Extension

with open("README.md", "r") as fh:
    long_description = fh.read()


class CythonExtension(Extension):
    def __init__(self, *args, **kwargs):
        self._include = []
        super().__init__(*args, **kwargs)

    @property
    def include_dirs(self):
        return self._include

    @include_dirs.setter
    def include_dirs(self, dirs):
        self._include = dirs


extensions = [
    CythonExtension(
        "libfftcy",
        ["libfft/libfftcy.pyx", "libfft/ludlows_libfft.h"],
        include_dirs=["libfft"],
        compiler_directives={'language_level': '3'},
        language="c")
]

setup(
    name="libfft",
    version="0.0.0",
    author="ludlows",
    description="library about fast fourier transform using the 2-base method",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ludlows/libfft",
    packages=find_packages(),
    package_data={'libfft': ["*.pyx", "*.h"]},
    ext_package='libfft',
    ext_modules=extensions,
    setup_requires=['setuptools>=18.0', 'cython'],
    classifiers=[
        "Programming Language :: Python",
        "Operating System :: OS Independent"
    ]
)

