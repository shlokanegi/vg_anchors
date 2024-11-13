import os
import sys

# Add the lib directory to sys.path
lib_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'lib'))
if lib_path not in sys.path:
    sys.path.append(lib_path)

# Import bdsg here if you want it available throughout the package
import bdsg

# CLI
from .cli import cli