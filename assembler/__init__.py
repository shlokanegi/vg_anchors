import os
import ctypes

# Preload libhandlegraph.so with RTLD_GLOBAL to make symbols available for libbdsg.so
# This must happen BEFORE bdsg is imported anywhere in the package
_current_dir = os.path.dirname(os.path.abspath(__file__))
_lib_path = os.path.join(_current_dir, '..', 'libbdsg', 'libhandlegraph.so')
if os.path.exists(_lib_path):
    ctypes.CDLL(_lib_path, mode=ctypes.RTLD_GLOBAL)
