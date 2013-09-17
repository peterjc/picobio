This is package is a Galaxy workflow for the identification of candidate
secreted proteins from a given protein FASTA file.

It runs SignalP v3.0 (Bendtsen et al. 2004) and selects only proteins with a
strong predicted signal peptide, and then runs TMHMM v2.0 (Krogh et al. 2001)
on those, and selects only proteins without a predicted trans-membrane helix.
This workflow was used in Kikuchi et al. (2011), and is a simplification of
the candidate effector protocol described in Jones et al. (2009).

As of 17 September 2013, development has moved from here:

* https://github.com/peterjc/picobio/tree/master/galaxy_workflows/secreted_protein_workflow

To here, along with the associated Galaxy tools:

* https://github.com/peterjc/pico_galaxy/tree/master/workflows/secreted_protein_workflow

This workflow is available to download and/or install from the main
Galaxy Tool Shed:

* http://toolshed.g2.bx.psu.edu/view/peterjc/secreted_protein_workflow
