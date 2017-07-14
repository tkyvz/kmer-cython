import os
import sys
from distutils.core import setup

try:
    from Cython.Build import cythonize
    from Cython.Distutils import Extension
except ImportError:
    print("""Cannot find Cython!
    Cython is required to correctly build kmer's pyx extensions.
    In most cases, running the following command should be sufficient:
        $ pip install Cython

    Exception: ImportError
    """)
    exit()

here = os.path.dirname(__file__)

ext_files = [
    'resources/dsk.pyx',
    'resources/bfcounter.pyx'
]

ext_modules = [
    Extension(
        'kmercounter',
        ext_files
    )
]

ext_modules = cythonize(ext_files)

setup(
    name='kmercounter',
    version='1.0.0',
    author='Utku Yavuz',
    author_email='utku.yavuz@yahoo.com.tr',
    url='https://github.com/tkyvz/kmer-cython',
    description='A K-mer counter based on BFCounter and DSK algorithms.',
    install_requires=[],
    ext_modules=ext_modules
)
