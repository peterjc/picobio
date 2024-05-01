#!/usr/bin/env python
"""Python script to convert HMMER3 table output into tab separated tables.

This handles both the per-sequence table and the per-domain table. The
script just takes two arguments, the input filename (or '-' for stdin),
and the output filename (or '-' for stdout).

This is v0.0.2 of the script.

Copyright 2012, Peter Cock, all rights reserved.

THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.
"""

from __future__ import print_function

import sys


def convert(input_handle, output_handle):
    """Convert HMMER space separated table into tab separated table.

    Can be given a per-sequence table,

    #                                                                       --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
    # target name        accession  query name                   accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
    #------------------- ----------         -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
    Amidinotransf        PF02274.12 Gpa_EST_02_04___C05_022|ORF2 -            3.3e-06   26.2   0.1   3.3e-06   26.2   0.0   1.0   1   0   0   1   1   1   1 Amidinotransferase

    Or, a per-domain table,

    #                                                                                    --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
    # target name        accession   tlen query name                   accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
    #------------------- ---------- -----         -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
    Amidinotransf        PF02274.12   281 Gpa_EST_02_04___C05_022|ORF2 -             77   3.3e-06   26.2   0.1   1   1   2.4e-10   3.3e-06   26.2   0.0   213   280    10    77     1    77 0.90 Amidinotransferase


    In both these examples, the columnation is clearly disturbed by the
    long query names.

    For these tables the HMMER3 manual states: "... is columnated neatly
    for human readability, but you should not write parsers that rely on
    this columnation; parse based on space-delimited fields instead."
    """
    h1 = input_handle.readline()
    assert h1.startswith("# ")
    h2 = input_handle.readline()
    assert h2.startswith("# target name")
    h3 = input_handle.readline()
    assert h3.startswith("#---")
    columns = len(h3.split())
    assert columns == 19 or columns == 23, columns

    # Hard code our expected header names (but allow for differences
    # in the spacing):
    if columns == 19:
        names = [
            "target name",
            "accession",
            "query name",
            "accession",
            "E-value",
            "score",
            "bias",
            "E-value",
            "score",
            "bias",
            "exp",
            "reg",
            "clu",
            "ov",
            "env",
            "dom",
            "rep",
            "inc",
            "description of target",
        ]
        assert " ".join(h2[2:-1].split()) == " ".join(names)
        # Now switch to longer names (including line one information):
        names = [
            "target name",
            "accession",
            "query name",
            "accession",
            "full sequence E-value",
            "full sequence score",
            "full sequence bias",
            "best 1 domain E-value",
            "best 1 domain score",
            "best 1 domain bias",
            "domain number estimation exp",
            "domain number estimation reg",
            "domain number estimation clu",
            "domain number estimation ov",
            "domain number estimation env",
            "domain number estimation dom",
            "domain number estimation rep",
            "domain number estimation inc",
            "description of target",
        ]
    else:
        names = [
            "target name",
            "accession",
            "tlen",
            "query name",
            "accession",
            "qlen",
            "E-value",
            "score",
            "bias",
            "#",
            "of",
            "c-Evalue",
            "i-Evalue",
            "score",
            "bias",
            "from",
            "to",
            "from",
            "to",
            "from",
            "to",
            "acc",
            "description of target",
        ]
        assert " ".join(h2[2:-1].split()) == " ".join(names)
        # Now switch to longer names (including line one information):
        names = [
            "target name",
            "accession",
            "tlen",
            "query name",
            "accession",
            "qlen",
            "full sequence E-value",
            "full sequence score",
            "full sequence bias",
            # The next two columns are for e.g. 1 of 3, 2 of 3, 3 of 3.
            "dom#",
            "ndom",
            "c-Evalue",
            "i-Evalue",
            "score",
            "bias",
            "hmm coord from",
            "hmm coord to",
            "ali coord from",
            "ali coord to",
            "env coord from",
            "env coord to",
            "acc",
            "description of target",
        ]
    assert len(names) == columns
    output_handle.write("#%s\n" % "\t".join(names))

    # Now the easy bit, tabify the data (using white space as instructed
    # in the HMMER3 manual, which says not to use the columnation).
    count = 0
    for line in input_handle:
        assert line[0] != "#"
        # There will be spaces in the last column (description of target)
        parts = line.rstrip("\n").split(None, columns - 1)
        assert len(parts) == columns, parts
        output_handle.write("\t".join(parts) + "\n")
        count += 1
    return count


if len(sys.argv) != 3:
    sys.exit("Expect two argument: input HMMER table, output tabular file")
in_file, out_file = sys.argv[1:]

if in_file == "-":
    inp = sys.stdin
else:
    inp = open(in_file, "rU")

if out_file == "-":
    out = sys.stdout
else:
    out = open(out_file, "w")

count = convert(inp, out)

if in_file != "-":
    inp.close()
if out_file != "-":
    out.close()
    print("Converted table with %i lines" % count)
