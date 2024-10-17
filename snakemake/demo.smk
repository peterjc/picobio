# Example usage at the command line:
#
# $ rm -rf *.md5; snakemake -q -s demo.smk -p $(for f in *.fna; do echo $f.md5; done); ls *.md5
#
# Here using a little bash loop to generate a listing of all the
# desired MD5 files based on the FASTA files present.
#
# Example usage from Python via the API (using same logic for targets):
#
# $ rm -rf *.md5; ./snakemake_progress_bar_demo.py ; ls *.md5
#
# The rule will sleep for between 1 and 10s, and then compute the MD5.
# However, 1 time in 20 it will fail instead.

rule fasta_checksum:
    input:
        "{genome}.fna"
    output:
        "{genome}.fna.md5"
    shell:
        #'X=$((1 + $RANDOM % 10)); if [ "$X" == "1" ]; then sleep 5; exit 1; else sleep $X; md5sum {input} > {output}; fi'
        'sleep $((1 + $RANDOM % 10)); if [ "$(($RANDOM % 20))" == "0" ]; then exit 1; else md5sum {input} > {output}; fi'
