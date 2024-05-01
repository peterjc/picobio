#!/usr/bin/env python
"""Python script to insert gap features into EMBL files for ENA submission.

Given a FASTA assembly with runs of Ns in the sequence fed through Prokka
for annotation and then gff3_to_embl to make an EMBL file, it will fail
ENA validation:

$ java -jar embl-api-validator-1.1.146.jar BXF1.embl
...
ERROR: Sequence contains a stretch of 'n' characters between base {0} and {1} that is not represented with a "gap" feature (stretches of n greater than {2} gives a warning, greater than {3} gives an error). (26 occurrences) (SequenceToGapFeatureBasesCheck-1)
...

The validator reports a warning for any Ns without a gap feature, and an
error for runs of at least 10 Ns without a gap feature. Therefore seems
sensible to create gap features for any run of 10 or more N?

References
----------
Prokka: https://github.com/tseemann/prokka

GFF3 to ENA ready EMBL: https://github.com/sanger-pathogens/gff3toembl

ENA Validator: https://www.ebi.ac.uk/ena/software/flat-file-validator
and https://github.com/enasequence/sequencetools

This script: https://github.com/peterjc/picobio/tree/master/seq_manipulation

"""

import sys

if len(sys.argv) != 3:
    sys.exit("Expects two arguments: Input EMLB filename, output EMBL filename\n")
input_embl = sys.argv[1]
output_embl = sys.argv[2]
MIN_GAP = 10  # TODO: Could be a command line option?


try:
    from Bio import SeqIO
    from Bio._py3k import StringIO
    from Bio.SeqFeature import FeatureLocation
    from Bio.SeqFeature import SeqFeature
except ImportError:
    sys.exit("This script requires Biopython 1.69 or later")


def insert_feature(record, feature):
    pos = int(feature.location.start)
    i = 0
    for i, f in enumerate(record.features):
        if int(f.location.start) > pos:
            break
    record.features.insert(i, feature)


def insert_gaps(record):
    seq = str(record.seq).upper()
    sys.stderr.write(
        "Record %s (length %i bp) has %i N characters\n"
        % (record.id, len(seq), seq.count("N"))
    )
    gap = "N" * MIN_GAP
    try:
        i = seq.find(gap)
    except IndexError:
        sys.stderr.write(
            "No long gaps in record %s (length %i bp)\n" % (record.id, len(seq))
        )
        return record

    count = 0
    while i != -1:
        j = i + len(gap)
        while seq[j] == "N":
            j += 1
        sys.stderr.write(
            "Record %s (length %i bp) has run of %i N from %i to %i\n"
            % (record.id, len(seq), j - i, i + 1, j)
        )
        # WARNING - I suspect the validator is broken for features of one,
        # where I think the location is just X, rather than X..X instead?
        gap_feature = SeqFeature(
            FeatureLocation(i, j), type="gap", qualifiers={"estimated_length": j - i}
        )
        insert_feature(record, gap_feature)
        count += 1
        i = seq.find(gap, j)
    sys.stderr.write(
        "Added %i gap features to record %s (length %i bp)\n"
        % (count, record.id, len(seq))
    )
    return record


print("Adding gap features to %s making %s" % (input_embl, output_embl))

# This is the original short-and-sweet implementation, however right now
# as of November 2016 the Biopython EMBL round-trip is not close enough
# to avoid creating extra warnings from the ENA submission validator.
#
# fixed_records = (insert_gaps(r) for r in SeqIO.parse(input_embl, "embl"))
# count = SeqIO.write(fixed_records, output_embl, "embl")
#
# New version uses original header combined with Biopython's output of
# the feature table etc.


def get_header(embl_string):
    """Return everything up to but excluding the FH line.

    i.e. All the header lines prior to the features.
    """
    assert "\nFH" in embl_string
    answer = []
    for line in embl_string.split("\n"):
        if line.startswith("FH  "):
            break
        else:
            answer.append(line)
    return "\n".join(answer) + "\n"


def get_body(embl_string):
    """Return everything after and including the FH line.

    i.e. All the features and the sequence itself.
    """
    assert "\nFH" in embl_string
    answer = []
    for line in embl_string.split("\n"):
        if line.startswith("FH   "):
            answer = []
        answer.append(line)
    answer = "\n".join(answer)
    return answer.rstrip("\n") + "\n"


def break_up_embl_file(handle):
    record = []
    for line in handle:
        record.append(line)
        if line.startswith("//"):
            yield "".join(record)
            record = []
    if record:
        yield "".join(record)


count = 0
with open(input_embl) as in_handle:
    with open(output_embl, "w") as out_handle:
        for embl_string in break_up_embl_file(in_handle):
            assert embl_string.startswith("ID "), embl_string
            assert embl_string.endswith("//\n"), embl_string
            count += 1
            out_handle.write(get_header(embl_string))
            r = SeqIO.read(StringIO(embl_string), "embl")
            out_handle.write(get_body(insert_gaps(r).format("embl")))
print("Done, %i records" % count)
