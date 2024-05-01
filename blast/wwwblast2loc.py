#!/usr/bin/env python
# Short Python script to parse the blastwww setup files to extract a list of
# nucleotide and protein BLAST databases (with descriptions) and write them
# out as location files for use in Galaxy.
#
# Copyright 2010, Peter Cock.
#
# v001 - First version
# v002 - Use print as function
from __future__ import print_function

import os

# This gives us the list of databases and their type (nt vs aa):
blastrc = "/var/www/html/blast/blast.rc"
# This gives us a sensible order and their descriptions:
blastwww = "/var/www/html/blast/blast.html"

# BLAST DB path,
blastpath = "/data/blastdb"

# Output files
blast_nt = "blastdb.loc"
blast_aa = "blastdb_p.loc"


def load_blast_db_list(filename):
    nt = set()
    aa = set()
    handle = open(filename)
    for line in handle:
        if line.startswith("#") or not line.strip():
            continue
        elif line.startswith("NumCpuToUse"):
            continue
        elif line.startswith(("blastn ", "tblastn", "tblastx ")):
            nt.update(line.rstrip().split()[1:])
        elif line.startswith(("blastp ", "blastx")):
            aa.update(line.rstrip().split()[1:])
        else:
            raise ValueError(line)
    handle.close()
    return nt, aa


nt, aa = load_blast_db_list(blastrc)
# print(nt)
# print(aa)


def load_blast_db_descr(html_filename, nt, aa):
    nt_list = []
    aa_list = []
    handle = open(html_filename)
    for line in handle:
        line = line.strip()
        if not line.startswith("<option ") or not line.endswith("</option>"):
            continue
        i = line.find(">")
        assert i != -1, line
        descr = line[i + 1 : -9].strip()
        attrs = line[8:i]
        for db in nt:
            if '"%s"' % db in attrs:
                if descr.endswith(" (nt)"):
                    descr = descr[:-5].strip()
                # print("NT: %s -> %s" % (db, descr))
                nt_list.append((db, descr))
                nt.remove(db)
                break
        for db in aa:
            if '"%s"' % db in attrs:
                if descr.endswith(" (aa)"):
                    descr = descr[:-5].strip()
                # print("AA: %s -> %s" % (db, descr))
                aa_list.append((db, descr))
                aa.remove(db)
                break
    handle.close()
    return nt_list, aa_list


nt_list, aa_list = load_blast_db_descr(blastwww, nt, aa)
# print(nt_list)
# print(aa_list)


for filename, dbs in [(blast_nt, nt_list), (blast_aa, aa_list)]:
    handle = open(filename, "w")
    handle.write("#Automatically generated from\n")
    handle.write("#file %s\n" % blastrc)
    handle.write("#and %s\n" % blastwww)
    for db, descr in dbs:
        handle.write("\t".join([db, descr, os.path.join(blastpath, db)]) + "\n")
    handle.close()
print("Done, %i nt and %i aa BLAST databases" % (len(aa_list), len(nt_list)))
if aa:
    print("Missing protein databases:")
    for db in aa:
        print(db)
if nt:
    print("Missing nucleotide databases:")
    for db in nt:
        print(db)
