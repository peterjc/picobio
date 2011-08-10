#!/usr/bin/env python
"""Sync BLAST database(s)

Arguments:

(*) Master path, e.g. /data/blastdb
(*) Local path, e.g. /tmp/galaxy-blastdb
(*) Database name(s), e.g. ncbi/nr

Return codes:

0 - Worked, up to date
2 - Locked, aborted
3 - Failed

TODO - Try single call to rsync with wild card (exclude tar balls)
TODO - Different return code when already up to date
TODO - Automatically expand pal and nal alias files...
"""
#master = "/data/blastdb"
#local = "/tmp/galaxy-blastdb"
#db = "ncbi/nr"

import os
import sys
import time

if len(sys.argv) < 4:
    sys.stderr.write("Bad BLAST database sync arguments\n")
    sys.stderr.write("Expect source path, dest path, one or more DBs\n")
    sys.exit(3)

master = sys.argv[1]
local = sys.argv[2]
names = sys.argv[3:]

if not os.path.isdir(master):
    sys.stderr.write("Master directory %s not found\n" % master)
    sys.exit(3)

if not os.path.isdir(local):
    os.makedirs(local, 0777)

n_ext = [".nal", ".nhd", ".nhi", ".nhr", ".nin", ".nnd", ".nni",
         ".nog", ".nsd", ".nsi", ".nsq"]
p_ext = [".pal", ".phd", ".phi", ".phr", ".pin", ".pnd", ".pni",
         ".pog", ".ppd", ".ppi", ".psd", ".psi", ".psq"]

cmd = "rsync -v -rtz %s %s"

def file_list(master, local, db):
    d = os.path.split(db)[0]
    directory, name = os.path.split(os.path.join(master, db))
    for f in os.listdir(directory):
        if f.startswith(name + "."):
            ext = os.path.splitext(f)[-1]
            if ext in n_ext or ext in p_ext:
                yield os.path.join(master, d, f), os.path.join(local, d, f)


for db in names:
    print db
    lock = os.path.join(local, db + ".lock")
    if not os.path.isdir(os.path.split(lock)[0]):
        os.makedirs(os.path.split(lock)[0], 0777)
    #Wait a little, enough for other copies of the script
    #to finish the sync check if the files are current
    if os.path.isfile(lock):
        time.sleep(5)
    if os.path.isfile(lock):
        time.sleep(10)
    if os.path.isfile(lock):
        time.sleep(15)
    if os.path.isfile(lock):
        time.sleep(30)
    if os.path.isfile(lock):
        time.sleep(60)
    if os.path.isfile(lock):
        sys.stderr.write("BLAST Database already locked:\n")
        try:
            handle = open(lock)
            sys.stderr.write(handle.read())
            handle.close()
        except:
            pass
        sys.stderr.write("\nAborting sync\n")
        sys.exit(2)
    try:
        handle = open(lock, 'w')
        handle.write(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()))
        handle.close()
    except:
        sys.stderr.write("Could not create BLAST DB lock\n")
        sys.exit(2)

    start = time.time()
    for old, new in file_list(master, local, db):
        print "%s -> %s" % (old, new)
        err = os.system(cmd % (old, new))
        if err:
            sys.stderr("Return code %i from rsync:\n%s\n" % (err, cmd % (old, new)))
            os.remove(lock)
            sys.exit(3)
    taken = time.time() - start
    os.remove(lock)
    if taken > 100:
        print "%s done in %0.1fm" % (db, taken/60.0)
    else:
        print "%s done in %is" % (db, int(taken))
print "Done"
