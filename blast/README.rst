Misc. BLAST scripts.

Auto-caching of Databases
=========================

Files ``blast_sync.py`` and ``blast_wrap.py`` are used to
pre-cache our central BLAST databases onto a cluster node's
local hard drive (using ``rsync``).

This works by adding wrapper scripts like ``$HOME/bin/blastp``::

    $ more ~/bin/blastp
    #!/bin/bash
    #This bash script pretends to be an NCBI BLAST command line tool
    #acting as a proxy via a Python wrapper script to cache databases.
    #echo $@
    DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    $DIR/ncbi_blast/blast_wrap.py $DIR/ncbi_blast/blastp "$@"

This runs ``$HOME/ncbi_blast/blast_wrap.py`` which checks if a sync
is required via ``$HOME/ncbi_blast/blast_sync.py'', and then runs
the real NCBI BLAST+ binary named ``$HOME/bin/ncbi_blast/blastp``.


Converting wwwblast BLAST DB list to Galaxy loc files
=====================================================

We used to run a ``wwwblast`` server with a collection of
local BLAST databases, but transitioned to using BLAST+ via
Galaxy - see https://github.com/peterjc/galaxy_blast

The script ``wwwblast2loc.py`` was used during our transition
period to generate the Galaxy location files ``blastdb.loc``
and ``blastdb_p.loc`` from the ```wwwblast`` listing defined
in ``blast.rc`` and ``blast.html``.
