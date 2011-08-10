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

TODO - Different return code when already up to date?
"""
#master = "/data/blastdb"
#local = "/tmp/galaxy-blastdb"
#db = "ncbi/nr"

import os
import sys
import time

try:
    from os.path import relpath
except ImportError:
    #Must be prior to Python 2.6
    #This implementation is based on that by James Gardner in
    #MIT licensed package BareNecessities
    import posixpath
    def relpath(path, start=posixpath.curdir):
        """Return a relative version of a path"""
        if not path:
            raise ValueError("no path specified")
        start_list = posixpath.abspath(start).split(posixpath.sep)
        path_list = posixpath.abspath(path).split(posixpath.sep)
        # Work out how much of the filepath is shared by start and path.
        i = len(posixpath.commonprefix([start_list, path_list]))
        rel_list = [posixpath.pardir] * (len(start_list)-i) + path_list[i:]
        if not rel_list:
            return curdir
        return posixpath.join(*rel_list)
    assert relpath("/data/blast/ncbi/nr.pal", "/data/blast") == "ncbi/nr.pal"

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

def sync_blast_alias_db(master, local, db, index):
    print "Syncing %s" % index
    handle = open(index)
    for line in handle:
        if line.startswith("DBLIST "):
            dbs = line[7:].strip().split()
            for d in dbs:
                d = os.path.join(os.path.split(index)[0], d)
                err = sync_blast_db(master, local, relpath(d, master))
                if err:
                    return err
    handle.close()
    return 0

def sync_blast_db(master, local, db):
    for index in [os.path.join(master, db + ".nal"),
                  os.path.join(master, db + ".pal")]:
        if os.path.isfile(index):
            err = sync_blast_alias_db(master, local, db, index)
            if err:
                return err
            #else continue to sync the alias file itself,

    cmd = "rsync -v -rtz --exclude=*.tar.gz --exclude=*.md5 %s %s"
    old = os.path.join(master, db + ".*")
    new = os.path.join(local, os.path.split(db)[0]) #Folder namer!
    #print "%s -> %s" % (old, new)
    print cmd % (old, new)
    err = os.system(cmd % (old, new))
    if err:
        sys.stderr.write("Return code %i from rsync:\n%s\n" % (err, cmd % (old, new)))
    return err

for db in names:
    print db
    lock = os.path.join(local, db + ".lock")
    if not os.path.isdir(os.path.split(lock)[0]):
        os.makedirs(os.path.split(lock)[0], 0777)
    #Wait a little, enough for other copies of the script
    #to finish the sync check if the files are current
    if os.path.isfile(lock):
        sys.stderr.write("BLAST Database already locked,\n")
        time.sleep(5)
    if os.path.isfile(lock):
        time.sleep(10)
    if os.path.isfile(lock):
        time.sleep(15)
    if os.path.isfile(lock):
        time.sleep(30)
    #if os.path.isfile(lock):
    #    time.sleep(60)
    if os.path.isfile(lock):
        try:
            handle = open(lock)
            sys.stderr.write(handle.read())
            handle.close()
        except:
            pass
        sys.stderr.write("Aborting sync\n")
        sys.exit(2)
    try:
        handle = open(lock, 'w')
        handle.write(time.strftime("%a, %d %b %Y %H:%M:%S +0000\n", time.gmtime()))
        handle.close()
    except:
        sys.stderr.write("Could not create BLAST DB lock\n")
        sys.exit(2)

    start = time.time()
    err = sync_blast_db(master, local, db)
    if err:
        os.remove(lock)
        sys.exit(3)
    taken = time.time() - start
    os.remove(lock)
    if taken > 100:
        print "%s done in %0.1fm" % (db, taken/60.0)
    else:
        print "%s done in %is" % (db, int(taken))
print "Done"
