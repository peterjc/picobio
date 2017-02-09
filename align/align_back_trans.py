#!/usr/bin/env python
import sys


def sys_exit(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg)
    sys.exit(error_level)

sys_exit("""Python script for 'back-translating' a protein alignment.

This script was originally available from here:
https://github.com/peterjc/picobio/tree/master/align

It is now available from here instead, with an optional Galaxy wrapper:
https://github.com/peterjc/pico_galaxy/tree/master/tools/align_back_trans

The Galaxy tool is available from the Galaxy Tool Shed here:
http://toolshed.g2.bx.psu.edu/view/peterjc/align_back_trans
""")
