This is package is a Galaxy workflow for comparing three RXLR prediction
methods with a Venn Diagram, and creates a FASTA file of any proteins
passing all three methods.

See http://www.galaxyproject.org for information about the Galaxy Project.


References
==========

Stephen C. Whisson, Petra C. Boevink, Lucy Moleleki, et al. (2007)
A translocation signal for delivery of oomycete effector proteins into
host plant cells. Nature 450:115-118.
http://dx.doi.org/10.1038/nature06203

Joe Win, William Morgan, Jorunn Bos, et al. (2007)
Adaptive evolution has targeted the C-terminal domain of the RXLR effectors
of plant pathogenic oomycetes. The Plant Cell 19:2349-2369.
http://dx.doi.org/10.1105/tpc.107.051037

Souvik Bhattacharjee, N. Luisa Hiller, Konstantinos Liolios, et al. (2006)
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


Developers
==========

This workflow is under source code control here:

https://github.com/peterjc/picobio/tree/master/galaxy_workflows/rxlr_venn_workflow

To prepare the tar-ball for uploading to the Tool Shed, I use this:

    $ tar -cf rxlr_venn_workflow.tar.gz README.txt repository_dependencies.xml rxlr_venn_workflow.ga

Check this,

    $ tar -tzf rxlr_venn_workflow.tar.gz
    README.txt
    repository_dependencies.xml
    rxlr_venn_workflow.ga
