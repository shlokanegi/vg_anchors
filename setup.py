from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import os
import sys
import subprocess

# Define the path to the bdsg library
bdsg_lib_path = os.path.join('libbdsg', 'lib')

# Ensure the library path is in the sys.path
if bdsg_lib_path not in sys.path:
    sys.path.insert(0, bdsg_lib_path)

class BuildCommand(build_py):
    """Custom build command to run PyInstaller."""
    def run(self):
        super().run()
        
        # Generate the default config.ini before building
        print("--- Generating default config.ini ---")
        subprocess.check_call([sys.executable, 'generate_config.py'])
        
        # Get the version from the setup() call
        version = self.distribution.get_version()
        executable_name = f"vg-anchors-{version}"
        
        print(f"--- Building executable: {executable_name} ---")
        
        # Build the PyInstaller command
        pyinstaller_command = [
            'pyinstaller',
            'assembler/cli.py',
            '--name', executable_name,
            '--onefile',
            '--add-data', f"{os.path.join('libbdsg', 'lib')}:lib",
            '--add-data', 'config.ini:.'
        ]
        
        subprocess.check_call(pyinstaller_command)

setup(
    name='vg-anchors',
    version='0.1.0',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'click',
        'matplotlib',
        'biopython',
        'pyinstaller',
    ],
    data_files=[
        ('lib', [os.path.join(bdsg_lib_path, f) for f in os.listdir(bdsg_lib_path) if os.path.isfile(os.path.join(bdsg_lib_path, f))])
    ],
    entry_points={
        'console_scripts': [
            'vg-anchors=assembler.cli:cli',
        ],
    },
    author='Shloka Negi',
    author_email='shnegi@ucsc.edu',
    description='A python tool to construct anchors from a pangenome using read alignments',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/shlokanegi/vg_anchors',
    cmdclass={
        'build_py': BuildCommand,
    }
)
