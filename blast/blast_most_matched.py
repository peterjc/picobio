#!/usr/bin/env python
from __future__ import print_function

import sys

from collections import OrderedDict

from Bio import SeqIO

min_run = 50
min_bases = 100

contig_fasta = sys.argv[1]  # FASTA file used as BLAST query
contig_blast = sys.argv[2]  # Tabular 12 (std) columns + optional extras


def cull_runs(set_of_points, min_run):
    answer = set()
    start = None
    end = None
    for i in sorted(set_of_points):
        if start is None:
            # very first run
            start = i
            end = i
        elif i == end + 1:
            # Continues run
            end = i
        else:
            # End of run.
            if end - start + 1 >= min_run:
                answer.update(range(start, end + 1))
            start = i
            end = i
    # Final run,
    if start and end and end - start + 1 >= min_run:
        answer.update(range(start, end + 1))
    return answer


assert cull_runs(set([1, 2, 3, 5, 6, 10, 11, 12, 13]), 5) == set()
assert cull_runs(set([1, 2, 3, 5, 6, 10, 11, 12, 13]), 4) == set([10, 11, 12, 13])
assert cull_runs(set([1, 2, 3, 5, 6, 10, 11, 12, 13]), 3) == set(
    [1, 2, 3, 10, 11, 12, 13]
)
assert cull_runs(set([1, 2, 3, 5, 6, 10, 11, 12, 13]), 2) == set(
    [1, 2, 3, 5, 6, 10, 11, 12, 13]
)
assert cull_runs(set([1, 2, 3, 5, 6, 10, 11, 12, 13, 20]), 5) == set()
assert cull_runs(set([1, 2, 3, 5, 6, 10, 11, 12, 13, 20]), 4) == set([10, 11, 12, 13])
assert cull_runs(set([1, 2, 3, 5, 6, 10, 11, 12, 13, 20]), 3) == set(
    [1, 2, 3, 10, 11, 12, 13]
)
assert cull_runs(set([1, 2, 3, 5, 6, 10, 11, 12, 13, 20]), 2) == set(
    [1, 2, 3, 5, 6, 10, 11, 12, 13]
)

# key = contig id
# value = length of contig
contig_lengths = OrderedDict()

# key = contig id
# value = length of preceeding contigs
contig_starts = dict()

offset = 0
for record in SeqIO.parse(contig_fasta, "fasta"):
    length = len(record)
    contig_lengths[record.id] = length
    contig_starts[record.id] = offset
    offset += length

# key = subject id
# value = set of base coordinates (cummulative along all contigs)
contig_mapping = dict()

contig_species = dict()
for line in open(contig_blast):
    parts = line.rstrip("\n").split("\t")
    assert len(parts) > 12
    (
        qseqid,
        sseqid,
        pident,
        length,
        mismatch,
        gapopen,
        qstart,
        qend,
        sstart,
        send,
        evalue,
        bitscore,
    ) = parts[:12]
    if len(parts) >= 14:
        contig_species[sseqid] = parts[14]
    start = int(qstart) - 1
    end = int(qend)
    assert 0 <= start <= end <= contig_lengths[qseqid]
    offset = contig_starts[qseqid]
    if sseqid not in contig_mapping:
        contig_mapping[sseqid] = set()
    contig_mapping[sseqid].update(range(offset + start, offset + end))


def pop_most_mapped():
    global contig_mapping
    # Sort by most bases mapped to each subject
    contig_mapping_counts = sorted(((len(v), k) for k, v in contig_mapping.items()))
    most_mapped_count, most_mapped = contig_mapping_counts[-1]
    # Now remove all those bases from consideration!
    taken_bases = contig_mapping[most_mapped]
    del contig_mapping[most_mapped]
    for sseqid in list(contig_mapping):  # list as editing dict during loop
        contig_mapping[sseqid].difference_update(taken_bases)
        contig_mapping[sseqid] = cull_runs(contig_mapping[sseqid], min_run)
        if len(contig_mapping[sseqid]) < min_bases:
            # print("Culled %s" % sseqid)
            del contig_mapping[sseqid]
    return most_mapped, most_mapped_count


for sseqid in list(contig_mapping):  # list as editing dict during loop
    contig_mapping[sseqid] = cull_runs(contig_mapping[sseqid], min_run)
    if len(contig_mapping[sseqid]) < min_bases:
        del contig_mapping[sseqid]

# print("- Raw -")
# contig_mapping_counts = sorted(((len(v), k) for k, v in contig_mapping.items()), reverse=True)
# for sseqid, bases in contig_mapping_counts:
#    print(sseqid, bases)
print("- Culling 1st hit -")
while contig_mapping:
    sseqid, count = pop_most_mapped()
    print(sseqid, count, contig_species.get(sseqid, ""))
