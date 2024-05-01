#!/usr/bin/env python
"""Python script for assembly contig de-duplication using BLASTN."""

import os
import shutil
import sys
import tempfile
from optparse import OptionParser

from Bio import SeqIO

usage = """Basic usage: dedup_assembly.py assembly.fasta dedup_output.fasta

This will sort the input assembly by contig length making a temporary
FASTA file and BLAST database.
"""

makeblastdb_binary = "makeblastdb"
blastn_binary = "blastn"


def sys_exit(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg.rstrip())
    sys.exit(error_level)


def run(cmd):
    print(cmd)
    return_code = os.system(cmd)
    if return_code:
        sys_exit("Return %i from: %s" % (return_code, cmd), return_code)


parser = OptionParser(usage=usage)
parser.add_option(
    "-m",
    "--min-contig-len",
    dest="min_len",
    type="int",
    help="Minimum contig length for FASTA output (default 1000)",
    default=1000,
)
parser.add_option(
    "-l",
    "--min-hit-len",
    dest="min_hit",
    type="int",
    help="Minimum BLAST hit length to consider (default 500)",
    default=500,
)
parser.add_option(
    "-p",
    "--min-perc-identity",
    dest="perc_identity",
    type="float",
    help="Minimum BLAST percentage identity to consider (default 95%)",
    default=95,
)
parser.add_option(
    "-c",
    "--min-cover",
    dest="min_cover",
    type="float",
    help="Minimum BLAST hit coverage to black-list (default 95%)",
    default=95,
)
(options, args) = parser.parse_args()

if len(args) != 2:
    sys_exit("Requires two arguments!\n\n" + usage)
assembly_fasta, output_fasta = args
min_len = int(options.min_len)
min_hit = int(options.min_hit)
min_cover = float(options.min_cover)
perc_identity = float(options.perc_identity)

if not os.path.isfile(assembly_fasta):
    sys_exit("Assembly FASTA file not found: %r" % assembly_fasta)

cols = "qseqid sseqid qlen slen length qstart qend"
c_query = 0
c_match = 1
c_qlen = 2
c_slen = 3
c_length = 4
c_qstart = 5
c_qend = 6


def prepare_sorted_fasta(assembly_fasta, sorted_fasta):
    # Sort & remove short contigs, being lazy and doing this in memory
    contigs = [r for r in SeqIO.parse(assembly_fasta, "fasta") if len(r) >= min_len]
    # Sort on length (longest first), tie break on identifier
    contigs.sort(key=lambda r: (-len(r), r.id))
    # Write sorted FASTA file
    count = SeqIO.write(contigs, sorted_fasta, "fasta")
    assert count == len(contigs)
    del contigs
    # Make BLAST database
    cmd = "%s -dbtype nucl -in %s" % (makeblastdb_binary, sorted_fasta)
    run(cmd)
    print("Prepared BLAST database of %i contigs passing minimum length" % count)


def prepare_blast(sorted_fasta, blast_file):
    cmd = '%s -db %s -query %s -out %s -perc_identity %f -outfmt="6 %s"' % (
        blastn_binary,
        sorted_fasta,
        sorted_fasta,
        blast_file,
        perc_identity,
        cols,
    )
    run(cmd)
    print("Ran self-BLAST")


def find_duplicates(blast_file):
    regions = {}
    lengths = {}
    for line in open(blast_file):
        fields = line.rstrip("\n").split("\t")
        query = fields[c_query]
        if query == fields[c_match]:
            # Ignore self hits
            continue
        qlen = int(fields[c_qlen])
        slen = int(fields[c_slen])
        if qlen > slen:
            # Ignore hits to smaller sequences
            continue
        elif qlen == slen and query > fields[c_match]:
            # If same length, using identifier as the tie break
            continue
        assert qlen <= slen
        if float(fields[c_length]) < min_hit:
            # Ignore short HSPs
            continue
        # Will next look at hit coverage etc
        qstart = int(fields[c_qstart])
        qend = int(fields[c_qend])
        assert qstart < qend
        try:
            regions[query].add((qstart, qend))
        except KeyError:
            regions[query] = {(qstart, qend)}
        lengths[query] = qlen
    for query in regions:
        regs = sorted(regions[query])
        qlen = lengths[query]
        if (1, qlen) in regs:
            # Short cut!
            # print("All of %s matched other longer contig(s)" % query)
            yield query
            continue
        # Now get effective total region covered by good hits
        i = 0
        q = 0
        while regs:
            qstart, qend = regs.pop(0)
            if i < qstart:
                # New region
                q += qend - qstart + 1
                i = qend
            elif i < qend:
                # Extending region
                q += qend - i
                i = qend
            else:
                # Already counted this region
                pass
        if len(regions[query]) == 1:
            # Sanity test
            qstart, qend = list(regions[query])[0]
            assert q == qend - qstart + 1
        # else:
        #    # For debugging:
        #    print qlen, sorted(regions[query])
        assert q <= qlen, "%0.2f of %s hit other longer contigs (%i of %i bp)" % (
            q * 100.0 / qlen,
            query,
            q,
            qlen,
        )
        if qlen * min_cover / 100 <= q:
            # print("%0.2f of %s hit other longer contigs" % (q * 100.0 / qlen, query))
            yield query


def dedup(assembly_fasta, blast_file, output_fasta):
    duplicates = set(find_duplicates(blast_file))
    print("Identified %i contigs to treat as duplicates" % len(duplicates))
    # Filter the original assembly FASTA file to retain its ordering.
    # Must repeat the minimum length filter here too...
    wanted = (
        r
        for r in SeqIO.parse(assembly_fasta, "fasta")
        if r.id not in duplicates and len(r) >= min_len
    )
    count = SeqIO.write(wanted, output_fasta, "fasta")
    print("Saved %i contigs to %s" % (count, output_fasta))


temp_dir = tempfile.mkdtemp(prefix="tmp_dedup_")
sorted_fasta = os.path.join(temp_dir, "pre_dedup_sorted.fasta")
blast_file = os.path.join(temp_dir, "pre_dedup_blast.tsv")

prepare_sorted_fasta(assembly_fasta, sorted_fasta)
prepare_blast(sorted_fasta, blast_file)

dedup(sorted_fasta, blast_file, output_fasta)

print("-" * 60)
count = 0
total = 0
for r in SeqIO.parse(assembly_fasta, "fasta"):
    count += 1
    total += len(r)
print("Input %i contigs, total length %i bp, in %s" % (count, total, assembly_fasta))
count = 0
total = 0
for r in SeqIO.parse(sorted_fasta, "fasta"):
    count += 1
    total += len(r)
print(
    "Min length gives %i contigs, total length %i bp, in %s"
    % (count, total, sorted_fasta)
)
count = 0
total = 0
for r in SeqIO.parse(output_fasta, "fasta"):
    count += 1
    total += len(r)
print("Output %i contigs, total length %i bp, in %s" % (count, total, output_fasta))
print("-" * 60)

shutil.rmtree(temp_dir)
