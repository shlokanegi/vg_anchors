#!/usr/bin/env python3
"""
Profile wrapper script for vg-anchors commands using line_profiler.

This script allows you to profile vg-anchors commands with line-by-line timing.
Usage:
    kernprof -l scripts/profile_vg_anchors.py <command> [arguments...]
    
After running, view the results with:
    python -m line_profiler profile_vg_anchors.py.lprof
"""

import sys
import os

# Add the parent directory to the path so we can import assembler
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from assembler.cli import cli

if __name__ == "__main__":
    # Click CLI automatically parses sys.argv
    cli()
