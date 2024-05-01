#!/usr/bin/env python
"""Python script to rename locus tags in EMBL/GenBank files.

When submitting an annotated genome assembly to the ENA, you
are asked for your locus tag prefix, which must be globally
unique. It is possible to register the locus tag prefix in
advance, and I've found out the hard way why this is a good
idea.

If another project has already taken the prefix you wanted,
tough: You'll have to rename your features. Hopefully your
collaborators/consortium members will forgive you.

Script usage:

./rename_locustags OLD_PREFIX NEW_PREFIX old_file new_file

Include the trailing underscore on your prefix if used (this
allows the script to be used on input files where there is
no underscore after the prefix).

Using a minus sign for the input or output filenames will
select stdin or stdout respectively.

This will simply loop over the lines looking for feature
qualifiers using the old locus tag prefix, which it edits.
It is a quick-and-dirty approach (not using a full parser),
which assumes a locus tag is short enough not to cause a line
wrap, and that the prefix is followed by an underscore.

e.g. Given the following where locus tag prefix YP has already
been taken:

FT   gene            2334..3335
FT                   /locus_tag="YP_00001"
FT                   /gene="rbsR_1"
FT   CDS             2334..3335
FT                   /product="Ribose operon repressor"
FT                   /inference="ab initio prediction:Prodigal:2.60"
FT                   /db_xref="UniProtKB/Swiss-Prot:P0ACQ0"
FT                   /locus_tag="YP_00001"
FT                   /gene="rbsR_1"
FT                   /transl_table=11

Running the script with new old prefix ``YP_`` and new locus tag
prefix ``XYZYP_`` gives (note the trailing underscores):

FT   gene            2334..3335
FT                   /locus_tag="XYZYP_00001"
FT                   /gene="rbsR_1"
FT   CDS             2334..3335
FT                   /product="Ribose operon repressor"
FT                   /inference="ab initio prediction:Prodigal:2.60"
FT                   /db_xref="UniProtKB/Swiss-Prot:P0ACQ0"
FT                   /locus_tag="XYZYP_00001"
FT                   /gene="rbsR_1"
FT                   /transl_table=11

This example used EMBL format, but GenBank format is also
supported.
"""

import sys

# TODO - Proper API
if len(sys.argv) != 5:
    sys.exit(
        "Expects four arguments: Old prefix, new prefix, input filename, output filename"
    )
old_prefix, new_prefix, in_filename, out_filename = sys.argv[1:]

sys.stderr.write(
    'Replacing /locus_tag="%s..." with /locus_tag="%s..."\n' % (old_prefix, new_prefix)
)

if in_filename == "-":
    in_handle = sys.stdin
else:
    in_handle = open(in_filename)

if out_filename == "-":
    out_handle = sys.stdout
else:
    out_handle = open(out_filename, "w")

# Accept either EMBL or GenBank format
patterns = (
    'FT                   /locus_tag="%s' % old_prefix,
    '                     /locus_tag="%s' % old_prefix,
)
count = 0
locus_tags = set()
for line in in_handle:
    if line.startswith(patterns):
        count += 1
        line = line.replace(
            '/locus_tag="%s' % old_prefix, '/locus_tag="%s' % new_prefix
        )
        locus_tags.add(line[32:].strip())  # left quotes in place
    out_handle.write(line)

if in_filename != "-":
    in_handle.close()
if out_filename != "-":
    out_handle.close()

if not count:
    sys.stderr.write("Original locus tag prefix %s_... not found\n" % old_prefix)
    sys.exit(1)
sys.stderr.write(
    "Edited %s_... -> %s_... in %i lines.\n" % (old_prefix, new_prefix, count)
)
sys.stderr.write("Saw %i unique locus tags\n" % len(locus_tags))
