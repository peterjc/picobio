This is package is a Galaxy workflow for comparing three RXLR prediction
methods with a Venn Diagram, and creates a FASTA file of any proteins
passing all three methods.

See http://www.galaxyproject.org for information about the Galaxy Project.


Sample Data
===========

This workflow was developed and run on several Phytophthora species.
For example, try the "Phyca11" protein set for Phytophthora capsici:

http://genome.jgi-psf.org/Phyca11/download/Phyca11_filtered_proteins.fasta.gz

You can upload this directly into Galaxy via this URL. Galaxy will handle
removing the gzip compression to give you the FASTA protein file which
has 19,805 protein sequences. The expected results:

* 89 RXLRs using Whisson et al. (2007)
* 124 RXLRs using Win et al. (2007)
* 162 RXLRs using Bhattacharjee et al. (2006)

Of these, only 79 sequences pass all three of the RXLR prediction tools,
while 19643 have no RXLR matches at all.

.. image:: https://raw.github.com/peterjc/picobio/master/galaxy_workflows/rxlr_venn_workflow/Phyca11_example_output.png
   :height: 400px
   :width: 400px


Citation
========

If you use this workflow directly, or a derivative of it, in work leading
to a scientific publication, please cite:

Cock, P.J.A. and Pritchard, L. 2013. Galaxy as a platform for identifying
candidate pathogen effectors. Chapter 1 in "Plant-Pathogen Interactions:
Methods and Protocols (Second Edition)"; Methods in Molecular Biology.
Humana Press, Springer. In press.

Whisson, S.C., Boevink, C.V., Moleleki, L., et al. (2007)
A translocation signal for delivery of oomycete effector proteins into
host plant cells. Nature 450:115-118.
http://dx.doi.org/10.1038/nature06203

Win, J., Morgan, W., Bos, J., et al. (2007)
Adaptive evolution has targeted the C-terminal domain of the RXLR effectors
of plant pathogenic oomycetes. The Plant Cell 19:2349-2369.
http://dx.doi.org/10.1105/tpc.107.051037

Bhattacharjee, S., Luisa Hiller, N., Liolios, K., et al. (2006)
The malarial host-targeting signal is conserved in the Irish potato famine
pathogen. PLoS Pathogens 2(5):e50.
http://dx.doi.org/10.1371/journal.ppat.0020050


Availability
============

This workflow is available to download and/or install from the main
Galaxy Tool Shed:

http://toolshed.g2.bx.psu.edu/view/peterjc/rxlr_venn_workflow

Test releases (which should not normally be used) are on the Test Tool Shed:

http://testtoolshed.g2.bx.psu.edu/view/peterjc/rxlr_venn_workflow

Development is being done on github here:

https://github.com/peterjc/picobio/tree/master/galaxy_workflows/rxlr_venn_workflow


Dependencies
============

These dependencies should be resolved automatically via the Galaxy Tool Shed:

* http://toolshed.g2.bx.psu.edu/view/peterjc/tmhmm_and_signalp
* http://toolshed.g2.bx.psu.edu/view/peterjc/seq_filter_by_id
* http://toolshed.g2.bx.psu.edu/view/peterjc/venn_list

However, at the time of writing those Galaxy tools have their own dependencies
required for this workflow which require manual installation (SignalP v3.0,
HMMER v2.0, and the R/Bioconductor package limma).


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial release to Tool Shed (May, 2013)
        - Expanded README file to include example data
v0.0.2  - Added link to Tool Shed in the workflow annotation explaining there
          is a README file with sample data, and a requested citation.
======= ======================================================================


Developers
==========

This workflow is under source code control here:

https://github.com/peterjc/picobio/tree/master/galaxy_workflows/rxlr_venn_workflow

To prepare the tar-ball for uploading to the Tool Shed, I use this:

    $ tar -cf rxlr_venn_workflow.tar.gz README.rst repository_dependencies.xml rxlr_venn_workflow.ga

Check this,

    $ tar -tzf rxlr_venn_workflow.tar.gz
    README.rst
    repository_dependencies.xml
    rxlr_venn_workflow.ga
