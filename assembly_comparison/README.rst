Python scripts for visual assessment of (bacterial) assemblies.

These use Biopython and GenomeDiagram (calling ReportLab) to render
images of my assemblies, with the goal of a visual summary especially
where a reference assembly is available.


Dependencies
============

* Python, tested with Python 2.7, available from http://python.org
* Biopython, tested with Biopython 1.62, available from http://biopython.org
* ReportLab, tested with ReportLab 2.6, available from http://reportlab.com
* NCBI BLAST+, tested with BLAST 2.2.27+, available from
  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/


Sample Data
===========

As an example, we will use the first public assembly of the 2011 Shiga-toxin
producing *Escherichia coli* O104:H4 outbreak in Germany. This was part of the
open-source crowd-sourcing analysis described in Rohde et al. (2011) and here:
https://github.com/ehec-outbreak-crowdsourced/BGI-data-analysis/wiki

You can download this FASTA file with 3,057 sequences for either of these URLs,
for example using the ``wget`` command under Linux:

* http://static.xbase.ac.uk/files/results/nick/TY2482/TY2482.fasta.txt
* https://github.com/ehec-outbreak-crowdsourced/BGI-data-analysis/blob/master/strains/TY2482/seqProject/BGI/assemblies/NickLoman/TY2482.fasta.txt

This FASTA file ``TY2482.fasta.txt`` was the initial TY-2482 strain assembled
by Nick Loman from 5 runs of Ion Torrent data released by the BGI, using the
MIRA 3.2 assembler. It was initially released via his blog,
http://pathogenomics.bham.ac.uk/blog/2011/06/ehec-genome-assembly/

We will also need a complete reference genome, ideally one which has been
annotated. For example, although not a particularly close match, we could try
*E. coli* K-12 by downloading these two files:

* ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid57779/NC_000913.fna
* ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid57779/NC_000913.gbk

You will need to install the NCBI BLAST+ standalone tools, specifically we
will use ``makeblastdb`` and (from within the Python script) ``blastn``.
Now prepare a BLAST database from the reference FASTA file,

    $ makeblastdb -in NC_000913.fna -dbtype nucl

You can now run this script using this command:

    $ python assembly_comparison.py TY2482.fasta.txt  NC_000913.fna

This will call ``blastn`` to produce tabluar output, then produce a PDF diagram
comparing the TY-2482 assembly to the full circle of the reference *E. coli*
strain.

Note that with the default settings, these strains are not similar enough to
show matches along most of the reference - just small regions all the way
around it.


TODO
====

* Proper command line API including specification of output PDF filename
  and the tabular BLAST results.

* Control over sequence similarity thresholds.

* Try BLAT etc instead of BLASTN.

* Control over the colours?

* Galaxy wrapper?

* etc
