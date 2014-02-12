#!/usr/bin/env python
import sys
from collections import OrderedDict
from Bio import SeqIO

contig_fasta = sys.argv[1]
contig_blast = sys.argv[2]

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
for line in open(contig_blast):
    parts = line.rstrip("\n").split("\t")
    assert len(parts) > 12
    qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = parts[:12]
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
    for sseqid in list(contig_mapping): #list as editing dict during loop
        contig_mapping[sseqid].difference_update(taken_bases)
        #TODO Remove isolated bases?
        if len(contig_mapping[sseqid]) < 100:
            #print("Culled %s" % sseqid)
            del contig_mapping[sseqid]
    return most_mapped, most_mapped_count

print "- Raw -"
contig_mapping_counts = sorted(((len(v), k) for k, v in contig_mapping.items()), reverse=True)
for sseqid, bases in contig_mapping_counts:
    print sseqid, bases
print "- Culling 1st hit -"
while contig_mapping:
    print pop_most_mapped()
