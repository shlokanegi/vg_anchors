# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='sample',
    version='0.1.0',
    description='Assemble HiFi reads using vg',
    long_description=readme,
    author='Francesco Andreace',
    author_email='',
    url='https://github.com/frankandreace',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)
