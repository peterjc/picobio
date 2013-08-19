This is package is a Galaxy workflow for the identification of candidate
secreted proteins from a given protein FASTA file.

It runs SignalP v3.0 (Bendtsen et al. 2004) and selects only proteins with a
strong predicted signal peptide, and then runs TMHMM v2.0 (Krogh et al. 2001)
on those, and selects only proteins without a predicted trans-membrane helix.
This workflow was used in Kikuchi et al (2001), and is a simplification of
the candidate effector protocol described in Jones et al (2009).

See http://www.galaxyproject.org for information about the Galaxy Project.


Citation
========

If you use this workflow directly, or a derivative of it, in work leading
to a scientific publication, please cite:

Cock, P.J.A. and Pritchard, L. 2013. Galaxy as a platform for identifying
candidate pathogen effectors. Chapter 1 in "Plant-Pathogen Interactions:
Methods and Protocols (Second Edition)"; Methods in Molecular Biology.
Humana Press, Springer. In press.

Also consider citing:

Bendtsen, J.D., Nielsen, H., von Heijne, G., Brunak, S. (2004)
Improved prediction of signal peptides: SignalP 3.0. J Mol Biol 340: 783–95.
http://dx.doi.org/10.1016/j.jmb.2004.05.028

Krogh, A., Larsson, B., von Heijne, G., Sonnhammer, E. (2001)
Predicting transmembrane protein topology with a hidden Markov model:
application to complete genomes. J Mol Biol 305: 567- 580.
http://dx.doi.org/10.1006/jmbi.2000.4315


Additional References
=====================

Kikuchi, T., Cotton, J.A., Dalzell, J.J., Hasegawa. K., et al. (2011)
Genomic insights into the origin of parasitism in the emerging plant
pathogen Bursaphelenchus xylophilus. PLoS Pathog 7: e1002219.
http://dx.doi.org/10.1371/journal.ppat.1002219

Jones, J.T., Kumar, A., Pylypenko, L.A., Thirugnanasambandam, A., et al. (2009)
Identification and functional characterization of effectors in expressed
sequence tags from various life cycle stages of the potato cyst nematode
Globodera pallida. Mol Plant Pathol 10: 815–28.
http://dx.doi.org/10.1111/j.1364-3703.2009.00585.x


Availability
============

This workflow is available to download and/or install from the main
Galaxy Tool Shed:

http://toolshed.g2.bx.psu.edu/view/peterjc/secreted_protein_workflow

Test releases (which should not normally be used) are on the Test Tool Shed:

http://testtoolshed.g2.bx.psu.edu/view/peterjc/secreted_protein_workflow

Development is being done on github here:

https://github.com/peterjc/picobio/tree/master/galaxy_workflows/secreted_protein_workflow


Dependencies
============

These dependencies should be resolved automatically via the Galaxy Tool Shed:

* http://toolshed.g2.bx.psu.edu/view/peterjc/tmhmm_and_signalp
* http://toolshed.g2.bx.psu.edu/view/peterjc/seq_filter_by_id

However, at the time of writing those Galaxy tools have their own
dependencies required for this workflow which require manual
installation (SignalP v3.0 and TMHMM v2.0).


Developers
==========

This workflow is under source code control here:

https://github.com/peterjc/picobio/tree/master/galaxy_workflows/secreted_protein_workflow

To prepare the tar-ball for uploading to the Tool Shed, I use this:

    $ tar -cf secreted_protein_workflow.tar.gz README.rst repository_dependencies.xml secreted_protein_workflow.ga

Check this,

    $ tar -tzf secreted_protein_workflow.tar.gz 
    README.rst
    repository_dependencies.xml
    secreted_protein_workflow.ga
