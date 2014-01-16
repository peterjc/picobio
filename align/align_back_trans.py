#!/usr/bin/env python
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio import AlignIO

def stop_err(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg)
    sys.exit(error_level)

def _sequence_back_translate(aligned_protein_record, unaligned_nucleotide_record, gap=None):
    #TODO - Separate arguments for protein gap and nucleotide gap?
    if hasattr(aligned_protein_record.seq.alphabet, "gap_char"):
        if not gap:
            gap = aligned_protein_record.seq.alphabet.gap_char
        elif gap != aligned_protein_record.seq.alphabet.gap_char:
            raise ValueError("Gap %r does not match %r from %r aligned protein alphabet" \
                             % (gap, aligned_protein_record.seq.alphabet.gap_char,
                                aligned_protein_record.id))
    if not gap:
        raise ValueError("Please supply the protein alignment gap character")

    alpha = unaligned_nucleotide_record.seq.alphabet
    if hasattr(alpha, "gap_char"):
        gap_codon = alpha.gap_char * 3
        assert len(gap_codon) == 3
    else:
        from Bio.Alphabet import Gapped
        alpha = Gapped(alpha, gap)
        gap_codon = gap*3

    if len(aligned_protein_record.seq.ungap(gap))*3 != len(unaligned_nucleotide_record.seq):
        stop_err("Inconsistent lengths for %s, ungapped protein %i, "
                 "tripled %i vs ungapped nucleotide %i" %
                 (len(aligned_protein_record.seq.ungap(gap)),
                  len(aligned_protein_record.seq.ungap(gap))*3,
                  len(unaligned_nucleotide_record.seq)))

    seq = []
    nuc = str(unaligned_nucleotide_record.seq)
    for amino_acid in aligned_protein_record.seq:
        if amino_acid == gap:
            seq.append(gap_codon)
        else:
            seq.append(nuc[:3])
            nuc = nuc[3:]
    assert not nuc, "Nucleotide sequence for %r longer than protein %s" \
        % (unaligned_nucleotide_record.id, aligned_protein_record.id)

    aligned_nuc = unaligned_nucleotide_record[:] #copy for most annotation
    aligned_nuc.letter_annotation = {} #clear this
    aligned_nuc.seq = Seq("".join(seq), alpha) #replace this
    assert len(aligned_protein_record.seq) * 3 == len(aligned_nuc)
    return aligned_nuc

def alignment_back_translate(protein_alignment, nucleotide_records, key_function=None, gap=None):
    """Thread nucleotide sequences onto a protein alignment."""
    #TODO - Separate arguments for protein and nucleotide gap characters?
    if key_function is None:
        key_function = lambda x: x

    if hasattr(protein_alignment._alphabet, "gap_char"):
        if not gap:
            gap = protein_alignment._alphabet.gap_char
        elif gap != protein_alignment._alphabet.gap_char:
            raise ValueError("Gap %r does not match %r from protein alignment alphabet" \
                             % (gap, protein_alignment._alphabet.gap_char))

    aligned = []
    for protein in protein_alignment:
        try:
            nucleotide = nucleotide_records[key_function(protein.id)]
        except KeyError:
            raise ValueError("Could not find nucleotide sequence for protein %r" \
                             % protein.id)
        aligned.append(_sequence_back_translate(protein, nucleotide, gap))
    return MultipleSeqAlignment(aligned)


try:
    align_format, prot_align_file, nuc_fasta_file = sys.argv[1:]
except:
    stop_err("""This is a Python script for 'back-translating' a protein alignment,

It requires three arguments:
- alignment format (e.g. fasta, clustal),
- aligned protein file (in specified format),
- unaligned nucleotide file (in fasta format).

The nucleotide alignment is printed to stdout (in the specified format).

Example usage, capturing stdout to a file by redirection:

$ python align_back_trans.py fasta demo_prot_align.fasta demo_nucs.fasta > demo_nuc_align.fasta

This script is available with sample data here:
https://github.com/peterjc/picobio/tree/master/align
""")

prot_align = AlignIO.read(prot_align_file, align_format, alphabet=generic_protein)
nuc_dict = SeqIO.index(nuc_fasta_file, "fasta")
nuc_align = alignment_back_translate(prot_align, nuc_dict, gap="-")
AlignIO.write(nuc_align, sys.stdout, align_format)
