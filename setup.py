from distutils.core import setup
from setuptools import find_packages

with open('version') as f:
    version = f.read()

setup(
    name='ldscore',
    version=version,
    packages=find_packages(),
    url='',
    license='',
    author='knomics',
    author_email='',
    description='ldscore by Bulik rewritten for python 3.6'
)