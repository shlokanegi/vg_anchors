from setuptools import setup, find_packages, Extension
import pybind11
from setuptools.command.build_py import build_py
import os
import sys
import subprocess

cpp_args = ['-std=c++11', '-O3']

gtest_module = Extension(
    'assembler.gtest',
    sources=[
        'assembler/cpp/GTest.cpp',
        'assembler/cpp/bindings.cpp',
        'assembler/cpp/SHASTA_ASSERT.cpp'
    ],
    include_dirs=[
        pybind11.get_include(),
        'assembler/cpp'
    ],
    language='c++',
    extra_compile_args=cpp_args,
)

class BuildCommand(build_py):
    """Custom build command (kept minimal; no PyInstaller here)."""
    def run(self):
        super().run()

setup(
    name='vg-anchors',
    version='0.1.0',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'Click',
        'bdsg',
        'matplotlib',
        'seaborn',
        'pandas',
        'numpy',
        'flask',
        'biopython',
        'pybind11>=2.6'
    ],
    entry_points={
        'console_scripts': [
            'vg-anchors=assembler.cli:cli',
        ],
    },
    author="Shloka Negi",
    author_email="shnegi@ucsc.edu",
    description="Python tool to construct anchors from a pangenome using read alignments",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/shlokanegi/vg_anchors',
    setup_requires=['pybind11>=2.6'],
    ext_modules=[gtest_module],
    cmdclass={'build_py': BuildCommand}
)
