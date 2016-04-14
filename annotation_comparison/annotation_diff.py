#!/usr/bin/env python3
import sys
from Bio import SeqIO

MISSING_QUALIFIERS_TO_IGNORE = ["translation", "codon_start", "db_xref"]
QUALIFIERS_TO_IGNORE = ["inference", "note"]

def sniff(handle):
    offset = handle.tell()
    line = handle.readline()
    handle.seek(offset)

    if line.startswith("##gff-version"):
        raise NotImplementedError("TODO: GFF support")
    elif line.startswith("LOCUS "):
        return SeqIO.parse(handle, "gb")
    elif line.startswith("ID "):
        return SeqIO.parse(handle, "embl")
    else:
        sys.exit("Could not guess file type from first line:\n%s" % line)

def clean(value):
    if value is None:
        return None
    if isinstance(value, list):
        if len(value) == 1:
            value = value[0]  # unlist
        else:
            return [clean(v) for v in value]
    assert isinstance(value, str), value
    if value == 'Uncharacterised protein':
        return 'hypothetical protein'
    if "%2C" in value:
        value = value.replace("%2C", ",")
    return value

def diff_f(old, new):
    assert old.type == new.type
    assert str(old.location) == str(new.location)


    if "locustag" in old.qualifiers and "locustag" in new.qualifiers:
        if old.qualifiers["locustag"] == new.qualifiers["locustag"]:
            name = old.qualifiers["locustag"]

    keys = set(old.qualifiers).union(new.qualifiers).difference(QUALIFIERS_TO_IGNORE)
    for k in keys:
        if k in MISSING_QUALIFIERS_TO_IGNORE:
            if k not in old.qualifiers or k not in new.qualifiers:
                continue
        old_v = clean(old.qualifiers.get(k, None))
        new_v = clean(new.qualifiers.get(k, None))
        if old_v != new_v:
            print("\t".join([str(old.location), old.type, k, repr(old_v), repr(new_v)]))

# TODO: Proper command line API
old_filename, new_filename = sys.argv[1:]

old_handle = open(old_filename)
new_handle = open(new_filename)

old_iter = sniff(old_handle)
new_iter = sniff(new_handle)

for old, new in zip(old_iter, new_iter):
    print("# Comparing records %s vs %s" % (old.id, new.id))
    assert old.id == new.id or old.id == "XXX"
    assert len(old) == len(new)
    assert len(old.features) == len(new.features)
    for old_f, new_f in zip(old.features, new.features):
        diff_f(old_f, new_f)

print("# Done")
