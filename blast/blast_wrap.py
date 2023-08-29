#!/usr/bin/env python
"""Simple Python wrapper for BLAST to locally cache databases.

This is intended for use within a computer cluster where there is
a central copy of the BLAST databases which is updated regularly,
but which is slow to access over the network. Therefore for speed,
and assuming sufficient disk space on each node, we want to copy
the BLAST databases to a local cache.

This can be setup as a scheduled task (e.g. with cron) for key
databases (see blast_sync.py), and/or uses via a wrapper script
to sync on demand (see blast_wrap.py).

Currently uses two hard coded settings for the network mounted
master copy of the databases (e.g. NFS, SAMBA), and the local
fast hard drive to use as the cache.

We're using /mnt/gfs/blast/galaxy as the master, and /var/blast/galaxy
as the local cache - see variables 'master' and 'local' below.

Intention is rather than this:

$ blastx -query=example.fasta -db=/mnt/gfs/blast/galaxy/nr ...

You do this:

$ blast_wrap.py blastx -query=example.fasta -db=/mnt/gfs/blast/galaxy/nr ...

This will cache the /mnt/gfs/blast/galaxy/nr.* files from the server
as /var/blast/galaxy/ncbi/nr.* and run this command for you:

$ blastx -query=example.fasta -db=/var/blast/galaxy/ncbi/nr ...

TODO: Work out the database path if not given explicitly (e.g. just nr)
but via the BLAST environment variable etc.
"""
from __future__ import print_function

# Old naming,
# master = "/mnt/shared/cluster/blast/galaxy"
# local = "/var/local/blast/galaxy"

master = "/mnt/shared/cluster/blast/galaxy"
local = "/mnt/scratch/local/blast/galaxy"

# e.g. db = "ncbi/nr"

import os
import sys
import time

# print("Given: %r" % sys.argv)

# argv[0] is this python script
# Turn the argv list into a string, escaping as needed


def wrap(text):
    if " " in text and not text[0] == '"' and not text[-1] == '"':
        return '"%s"' % text
    elif "|" in text:
        return '"%s"' % text
    else:
        return text


cmd = " ".join(wrap(arg) for arg in sys.argv[1:])

if master in cmd:
    # We have syncing to do!
    i = cmd.find(master + "/")
    db = cmd[i + len(master) + 1 :].split(None, 1)[0]
    if db.endswith('"'):
        db = db.rstrip('"')
    print("Synchronising database,")
    assert master + "/" + db in cmd
    start = time.time()
    sync = os.path.join(os.path.dirname(os.path.abspath(__file__)), "blast_sync.py")
    if not os.path.isfile(sync):
        sys.stderr.write("Syncing %s failed (missing %s script)\n" % (master, sync))
        sys.exit(1)
    err = os.system("%s %s %s %s > /dev/null" % (sync, master, local, db))
    taken = time.time() - start
    if 0 < err < 128:
        sys.stderr.write("Syncing %s failed (error code %i)\n" % (db, err))
        sys.exit(err)
    elif err:
        sys.stderr.write("Syncing %s failed (error code %i --> 1)\n" % (db, err))
        sys.exit(1)
    # Update the command
    cmd = cmd.replace(master + "/" + db, local + "/" + db)
    if taken > 100:
        print("%s done in %0.1fm" % (db, taken / 60.0))
    else:
        print("%s done in %is" % (db, int(taken)))

# Run the command
# print cmd
err = os.system(cmd)
if 0 < err < 128:
    sys.exit(err)
elif err:
    # Returning 512 gives 0 (odd)
    sys.exit(1)
