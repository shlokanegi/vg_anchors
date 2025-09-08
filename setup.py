from setuptools import setup, find_packages, Extension
import pybind11

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
        'assembler/cpp'  # Add the directory containing your C++ headers
    ],
    language='c++',
    extra_compile_args=cpp_args,
)

setup(
    name='assembler',
    version='0.1',
    packages=find_packages(),
    author="ShlokaNegi",
    author_email="shnegi@ucsc.edu",
    description="Python tool to construct anchors from a pangenome using read alignments",
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
    setup_requires=['pybind11>=2.6'],
    ext_modules=[gtest_module],
    entry_points='''
        [console_scripts]
        vg_anchor=assembler.cli:cli
    ''',
)
