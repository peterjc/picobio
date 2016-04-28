#!/usr/bin/env python3
import sys
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import UnknownSeq
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


def patch_gff(handle, diffs):
    """Quick hack to patch Bacterial GFF files from Prokka etc.

    Does NOT support multi-line features (i.e. splicing and
    multiple exons).
    """
    out_handle = sys.stdout

    line = handle.readline()
    assert line.startswith("##gff-version 3"), line
    out_handle.write(line)

    # print("Parsing GFF3")
    references = OrderedDict()
    for line in handle:
        if line.strip() == "##FASTA":
            out_handle.write(line)
            break
        elif line.startswith("#"):
            out_handle.write(line)
            continue
        elif line.count("\t") == 8:
            seqid, source, ftype, start, end, score, strand, phase, attributes = line.rstrip().split("\t")
            #assert seqid in references, seqid
            start = int(start)  # Leave this as one-based
            end = int(end)
            assert 0 <= start < end #< len(references[seqid])
            loc = "%i..%i" % (start, end)
            if strand == "-":
                loc = "complement(%s)" % loc
            elif strand not in "+.?":
                # "+" = Forward strand - do nothing
                # "." = Unstranded - do nothing to match INSDC
                # "?" = Stranded but missing - do nothing to match INSDC
                raise ValueError("Bad strand %r in line: %s" % (strand, line))
            diff_key = "%s:%s:%s" % (seqid, loc, ftype)
            if diff_key in diffs:
                a = attributes + ";"
                for key, old, new in diffs[diff_key]:
                    if key == "EC_number":
                        key = "eC_number"
                    if new is None:
                        new = ""  # Remove it!
                    else:
                        new = "%s=%s" % (key, new.replace(",", "%2C"))
                    if old is None:
                        a += new + ";"
                    else:
                        old = "%s=%s" % (key, old.replace(",", "%2C"))
                        assert old in a, (line, old, new)
                        a = a.replace(old, new)
                assert a.endswith(";")
                line = line.replace(attributes, a[:-1])
            assert line.count("\n") == 1 and line.endswith("\n"), repr(line)
            out_handle.write(line)
        else:
            raise NotImplementedError(line)
    # Deal with any FASTA block
    for line in handle:
        out_handle.write(line)

def load_diffs(handle):
    answer = dict()
    for line in handle:
        if line.startswith("#"):
            continue
        parts = line.strip("\n").split("\t")
        ref, loc, ftype, key, old, new = parts
        if old == "None":
            old = None
        elif old.startswith("'") and old.endswith("'"):
            old = old[1:-1]
        elif old.startswith('"') and old.endswith('"'):
            old = old[1:-1]
        else:
            raise NotImplementedError(old)
        if new == "None":
            new = None
        elif new.startswith("'") and new.endswith("'"):
            new = new[1:-1]
        elif new.startswith('"') and new.endswith('"'):
            new = new[1:-1]
        else:
            raise NotImplementedError(new)
        diff_key = "%s:%s:%s" % (ref, loc, ftype)
        if diff_key not in answer:
            answer[diff_key] = []
        answer[diff_key].append((key, old, new))
    return answer

def apply_diffs(handle, diffs):
    offset = handle.tell()
    line = handle.readline()
    handle.seek(offset)

    if line.startswith("##gff-version"):
        return patch_gff(handle, diffs)
    elif line.startswith("LOCUS "):
        raise NotImplementedError
    elif line.startswith("ID "):
        raise NotImplementedError
    else:
        sys.exit("Could not guess file type from first line:\n%s" % line)

# TODO: Proper command line API
try:
    diff_filename, old_filename = sys.argv[1:]
except ValueError:
    sys.exit("Want two arguments: patch/diff filename, and input filename\n")

if diff_filename == "-":
    diffs = load_diffs(sys.stdin)
else:
    with open(diff_filename) as handle:
        diffs = load_diffs(handle)
sys.stderr.write("Loaded %i differences to apply\n" % len(diffs))
#sys.stderr.write("%s\n" % list(diffs.keys())[:10])
#sys.stderr.write("%s\n" % list(diffs.values())[0])

old_handle = open(old_filename)
apply_diffs(old_handle, diffs)
sys.stderr.write("Done\n")
