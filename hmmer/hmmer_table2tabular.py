#!/usr/bin/env python
"""Python script to convert HMMER3 table output into tab separated tables.

This handles both the per-sequence table and the per-domain table. The
script just takes two arguments, the input filename (or - for stdin),
and the output filename (or '-' for stdout).

TODO - Preserve the column headings

This is v0.0.1 of the script.

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
import os
import sys

def sys_exit(msg, err=1):
    sys.stderr.write(msg.rstrip()+"\n")
    sys.exit(err)

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

    """
    h1 = input_handle.readline()
    assert h1.startswith("# ")
    h2 = input_handle.readline()
    assert h2.startswith("# target name")
    h3 = input_handle.readline()
    assert h3.startswith("#---")
    columns = len(h3.split())
    assert columns == 19 or columns == 23, columns

    count = 0
    for line in input_handle:
        assert line[0] != "#"
        #There will be spaces in the last column (description of target)
        parts = line.rstrip("\n").split(None, columns-1)
        assert len(parts) == columns, parts
        output_handle.write("\t".join(parts) + "\n")
        count += 1
    return count

if len(sys.argv) != 3:
    sys_exit("Expect two argument: input HMMER table, output tabular file")
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
    print "Converted table with %i lines" % count
