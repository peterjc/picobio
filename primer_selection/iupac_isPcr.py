#!/usr/bin/env python
"""Generalise Jim Kent's isPcr to support IUPAC ambiguities by brute force.

As of v33 at least, ambiguous bases are rejected in the primers. So, this
script generalises the input file to record all the non-ambiguous
interpretations of the primer. Running isPcr will take several times longer,
and the output will probably need to be deduplicated.

The input and output are simple three-column TSV files with the name of
each primer pair, the forward primer sequence, and the reverse primer
sequence.
"""
import itertools
import sys

from Bio.Data.IUPACData import ambiguous_dna_values

expand_iupac = {
    # Treat I (inosine) like N
    "I": list(ambiguous_dna_values["N"].upper()),
    "i": list(ambiguous_dna_values["N"].lower()),
}
for base, expanded in ambiguous_dna_values.items():
    expand_iupac[base.upper()] = list(expanded.upper())
    expand_iupac[base.lower()] = list(expanded.lower())


def expand_iupac_bases(seq):
    """All possible unabmiguous sequences described with IUPAC ambiguities.

    e.g.

    >>> list(expand_iupac_bases("DAY"))
    ['AAC', 'AAT', 'GAC', 'GAT', 'TAC', 'TAT']
    """
    try:
        for alt in itertools.product(*[expand_iupac[base] for base in seq]):
            yield "".join(alt)
    except KeyError as err:
        sys.exit(f"ERROR - Problem with primer sequence {seq}, {err}")


before = after = 0
for line in sys.stdin:
    if line.startswith("#") or not line.strip():
        continue
    try:
        idn, fwd, rev = line.strip("\n").split("\t")[:3]
    except ValueError:
        t = line.count("\t")
        sys.exit(f"ERROR: Only {t} tabs in line: {line}")
    before += 1
    for fwd2 in expand_iupac_bases(fwd):
        for rev2 in expand_iupac_bases(rev):
            sys.stdout.write(f"{idn}\t{fwd2}\t{rev2}\n")
            after += 1
sys.stderr.write(f"Generalised {before} primer pairs into unabmiguous {after} pairs\n")
