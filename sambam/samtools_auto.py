#!/usr/bin/env python
"""Dirty hack to allow mixing of samtools commands between versions.

It can be a downside that the samtools command line API is a single
binary which offers multiple (often independent) commands.

Right now, samtools 1.1 still lacks some functionality from 0.1.19,
for example "samtools index", "samtools depad" and "samtools rmdup"
are not yet fully functional. e.g.
 - https://github.com/samtools/samtools/issues/199
 - https://github.com/samtools/samtools/issues/291

Conversely, the "samtools bam2fq" from samtools 0.1.19 has had
several issues fixed.

This wrapper allows me to call "samtools" and route this to the
appropriate binary. In this case:

 - ``samtools`` (alone) will call samtools 1.1
 - ``samtools bam2fq [...]`` will call samtools 1.1
 - ``samtools depad [...]`` will call samtools 0.1.19
 - ``samtools rmdup [...]`` will call samtools 0.1.19
 - etc

Install this by putting the Python script (or a symlink to it) on
your ``$PATH`` as ``samtools``, for example under ``~/bin/``::

   $ cd ~/bin
   $ ln -s samtools_auto.py samtools

Also install binaries for samtools 0.1.19 and 1.1 and set their
paths below (variables ``samtools_old`` and ``samtools_new``).
"""
import os
import sys

samtools_old = "/mnt/galaxy/bin/samtools_0.1.19"
samtools_new = "/mnt/galaxy/bin/samtools_1.1"

def pick_binary():
    """Return new samtools unless known to be using a broken command.

    i.e. Avoid samtools commands with known regressions!
    """
    if len(sys.argv) == 1:
        return samtools_new
    elif sys.argv[1] in ["index", "depad", "rmdup"]:
        return samtools_old
    else:
        return samtools_new

#argv[0] is this python script
#Turn the argv list into a string, escaping as needed
def wrap(text):
    if " " in text and not text[0]=='"' and not text[-1]=='"':
        return '"%s"' % text
    else:
        return text

cmd = pick_binary() + " " + " ".join(wrap(arg) for arg in sys.argv[1:])

err = os.system(cmd)
if 0 < err < 128:
    sys.exit(err)
elif err:
    #Returning 512 gives 0 (odd)
    sys.exit(1)
