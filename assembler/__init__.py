import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))

#path to ../libbdsg/lib
lib_path = os.path.abspath(os.path.join(current_dir, '..', 'libbdsg', 'lib'))
if lib_path not in sys.path:
    sys.path.append(lib_path)

# Import bdsg 
import bdsg

# CLI
from .cli import cli