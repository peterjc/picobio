import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio import AlignIO

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
    sys.stderr.write("Three arguments: align_format, prot_align_file, nuc_fasta_file\n\nOutput to stdout.\n")
    sys.exit(1)

prot_align = AlignIO.read(prot_align_file, align_format, alphabet=generic_protein)
nuc_dict = SeqIO.index(nuc_fasta_file, "fasta")
nuc_align = alignment_back_translate(prot_align, nuc_dict, gap="-")
AlignIO.write(nuc_align, sys.stdout, align_format)
