# Use case:
# - Have multiple annotated GenBank files
# - Aligned with Mauve, and orthologue file exported
#
# Want to copy the sister genome's gene identifiers
# into a reference GenBank file (as gene aliases, notes,
# etc) so they can be viewed/searched for easily.

from __future__ import print_function

from Bio import SeqIO

mauve_orthologues_file = "mauve_orthologues.txt"
reference_genbank_file = "reference.gbk"
reference_number_in_mauve = 0
output_genbank_file = "reference_with_aliases.gbk"

# Might be more than one contig
reference_records = list(SeqIO.parse(reference_genbank_file, "genbank"))
cds_dict = dict()
for r in reference_records:
    for f in r.features:
        if f.type == "CDS":
            name = f.qualifiers["gene"][0]
            key = "%i:%s:%i-%i" % (reference_number_in_mauve, name, f.location.start + 1, f.location.end)
            cds_dict[key] = f
#print(list(cds_dict.keys()))

for line in open(mauve_orthologues_file, "rU"):
    parts = sorted(line.strip().split("\t"))
    key = None
    #print(parts)
    for x in parts:
        if x in cds_dict:
        # if x.startswith("%i|" % reference_number_in_mauve):
            print("Using: %r" % parts)
            name = x.split(":")[1]
            names = [y.split(":")[1] for y in parts if y != x]
            cds_dict[x].qualifiers["name"] = [",".join([name] + names)]

SeqIO.write(reference_records, output_genbank_file, "genbank")
print("Wrote to %s" % output_genbank_file)
