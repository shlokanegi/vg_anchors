import os
import sys

# Ensure bundled shared libraries are discoverable at runtime
this_dir = os.path.dirname(sys.executable if getattr(sys, 'frozen', False) else __file__)
lib_dir = os.path.join(this_dir, 'lib')
cur = os.environ.get('LD_LIBRARY_PATH', '')
paths = [this_dir, lib_dir]
if cur:
    paths.append(cur)
os.environ['LD_LIBRARY_PATH'] = ':'.join(paths)

