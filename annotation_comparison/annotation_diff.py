#!/usr/bin/env python3
import sys
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import UnknownSeq
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

# TODO - expose this as a public API in Biopython:
from Bio.SeqIO.InsdcIO import _insdc_location_string as location_string

FEATURE_TYPE_TO_IGNORE = ["source"]
FEATURE_TYPE_WANTED = ["CDS"]  # Empty/None for any not in FEATURE_TYPE_TO_IGNORE
MISSING_QUALIFIERS_TO_IGNORE = ["translation", "codon_start", "db_xref", "ID", "transl_table"]
QUALIFIERS_TO_IGNORE = ["inference", "note"]

def parse_gff(handle):
    """Quick hack to parse Bacterial GFF files from Prokka etc.

    Does NOT support multi-line features (i.e. splicing and
    multiple exons). Will load EVERYTHING into memory!

    Iterator yielding SeqRecord objects, intended to fit into the
    Biopython SeqIO structure.
    """
    line = handle.readline()
    assert line.startswith("##gff-version 3"), line
    # print("Parsing GFF3")
    references = OrderedDict()
    for line in handle:
        # print(line)
        if line.startswith("##sequence-region "):
            _, name, start, end = line.split()
            assert start == "1"
            references[name] = SeqRecord(UnknownSeq(int(end)), id=name, name=name)
        elif line.strip() == "##FASTA":
            break
        elif line.startswith("#"):
            raise NotImplementedError(line)
        elif line.count("\t") == 8:
            seqid, source, ftype, start, end, score, strand, phase, attributes = line.split("\t")
            assert seqid in references, "Reference %r not declared with ##sequence-region line:\n%r" % (seqid, line)
            start = int(start) - 1
            end = int(end)
            assert 0 <= start < end < len(references[seqid])
            if ftype in FEATURE_TYPE_TO_IGNORE:
                continue
            if FEATURE_TYPE_WANTED and ftype not in FEATURE_TYPE_WANTED:
                continue
            if strand == "+":
                loc = FeatureLocation(start, end, +1)
            elif strand == "-":
                loc = FeatureLocation(start, end, -1)
            elif strand == ".":
                # Unstranded - should use zero but +1 to match EMBL/GB
                loc = FeatureLocation(start, end, +1)
            elif strand == "?":
                # Stranded by missing - should use None but +1 to match EMBL/GB
                loc = FeatureLocation(start, end, +1)
            else:
                raise ValueError("Bad strand %r in line:\n%r" % (strand, line))
            f = SeqFeature(loc, type=ftype)
            for part in attributes.strip().split(";"):
                if not part:
                    assert ";;" in line, line
                    sys.stderr.write("Warning - missing key=value or double semi-colon in line:\n%r\n" % line)
                    continue
                if "=" not in part:
                    sys.exit("Bad key=value entry %r in line:\n%r" % (part, line))
                key, value = part.split("=", 1)
                if key in MISSING_QUALIFIERS_TO_IGNORE:
                    continue
                if key == "eC_number":
                    key = "EC_number"
                value = value.replace("%2C", ",")
                try:
                    f.qualifiers[key].append(value)
                except KeyError:
                    f.qualifiers[key] = [value]
            references[seqid].features.append(f)
        else:
            raise NotImplementedError(line)
    # Deal with any FASTA block
    name = None
    seqs = []
    for line in handle:
        if line.startswith(">"):
            if name and seqs:
                seq = "".join(seqs)
                assert len(seq) == len(references[name]), \
                    "FASTA entry for %s was %i long, expected %i" % (name, len(seq), len(references[name]))
                references[name].seq = Seq(seq)
            name = line[1:].split(None, 1)[0]
            seqs = []
        elif name:
            seqs.append(line.strip())
        elif line.strip():
            raise NotImplementedError(line)
    if name and seqs:
        seq = "".join(seqs)
        assert len(seq) == len(references[name]), \
            "FASTA entry for %s was %i long, expected %i" % (name, len(seq), len(references[name]))
        references[name].seq = Seq(seq)
    # Return results
    for name, record in references.items():
        # print("%s length %i with %i features" % (name, len(record), len(record.seq)))
        yield record

def sniff(handle):
    offset = handle.tell()
    line = handle.readline()
    handle.seek(offset)

    if line.startswith("##gff-version"):
        return parse_gff(handle)
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

def diff_f(ref_name, ref_len, old, new):
    assert old.type == new.type
    assert str(old.location) == str(new.location), \
        "%s location %s vs %s" % (old.type, old.location, new.location)
    assert location_string(old.location, ref_len) == location_string(new.location, ref_len), \
        "%s location %s vs %s" % (old.type, location_string(old.location, ref_len), location_string(new.location, ref_len))

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
            if k == "locus_tag" and old_v.split("_", 1)[1] == new_v.split("_", 1)[1]:
                # Different prefix, ignore
                pass
            elif k == "product" and old_v.split() == new_v.split():
                # White space only, ignore
                pass
            else:
                print("\t".join([ref_name, location_string(old.location, ref_len), old.type, k, repr(old_v), repr(new_v)]))

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
    if FEATURE_TYPE_WANTED:
        old_fs = [f for f in old.features if f.type in FEATURE_TYPE_WANTED]
        new_fs = [f for f in new.features if f.type in FEATURE_TYPE_WANTED]
    else:
        old_fs = [f for f in old.features if f.type not in FEATURE_TYPE_TO_IGNORE]
        new_fs = [f for f in new.features if f.type not in FEATURE_TYPE_TO_IGNORE]

    assert len(old_fs) == len(new_fs), \
        "Have %i [%i] vs %i [%i] features, aborting" % (len(old_fs), len(old.features), len(new_fs), len(new.features))
    for old_f, new_f in zip(old_fs, new_fs):
        diff_f(old.id, len(old), old_f, new_f)

print("# Done")
