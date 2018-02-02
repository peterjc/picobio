import os
import sys
try:
    from urllib.request import urlopen
except ImportError:
    from urllib import urlopen

project = "PRJEB2896"

# This was simple, but does not work anymore:
# fastq_url = "http://www.ebi.ac.uk/ena/data/view/reports/sra/fastq_files/internal/%s" % project
fastq_url = "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJEB2896&result=read_run&fields=study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp,fastq_galaxy,submitted_ftp,submitted_galaxy,sra_ftp,sra_galaxy,cram_index_ftp,cram_index_galaxy&download=txt"
# Or goto https://www.ebi.ac.uk/ena/data/view/PRJEB2896 and click on "TEXT" download link
fastq_file = "%s_fastq.tsv" % project


wanted = """ERS091755
ERS092427
ERS001595
ERS092081
ERS092426
ERS092348
ERS092525
ERS091953
ERS092579
ERS001598
ERS092349
ERS001809
ERS092350
ERS002001
ERS092351
ERS091952
ERS092525""".split()


def download_in_one(url, filename):
    print("Fetching %s" % url)
    n = urlopen(url)
    data = n.read()
    n.close()

    h = open(filename, "wb")
    h.write(data)
    h.close()
    print("Saved as %s" % filename)

if not os.path.isfile(fastq_file):
    download_in_one(fastq_url, fastq_file)

def process_fastq(project, fastq_filename):
    h = open(fastq_filename)
    line = h.readline()
    assert line == 'study_accession\tsample_accession\tsecondary_sample_accession\texperiment_accession\trun_accession\ttax_id\tscientific_name\tinstrument_model\tlibrary_layout\tfastq_ftp\tfastq_galaxy\tsubmitted_ftp\tsubmitted_galaxy\tsra_ftp\tsra_galaxy\tcram_index_ftp\tcram_index_galaxy\n', repr(line)
    for line in h:
        parts = line.rstrip("\n").split("\t")
        assert parts[0] == project
        urls = parts[9].split(";")
        for url in urls:
            if url.startswith("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR"):
                pass
            elif url.startswith("ftp.sra.ebi.ac.uk/vol1/fastq/ERR"):
                url = "ftp://" + url
            assert url.startswith("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR"), url
            filename = url[len("ftp://ftp.sra.ebi.ac.uk/"):]
            acc = parts[2]  # here using secondary_sample_accession
            if wanted and acc not in wanted:
                print("Not interested in %s from %s" % (filename, acc))
                continue
            if os.path.isfile(filename):
                print("Already have %s" % filename)
                continue
            if not filename.endswith(".fastq.gz"):
                print("Skipping %s" % filename)
                continue
            #Make directory...
            d = os.path.split(filename)[0]
            if not os.path.isdir(d):
                print("Making directory %s" % d)
                os.makedirs(d)
            #Download file...
            print("Downloading %s" % filename)
            rc = os.system("wget -nv -O %s %s" % (filename, url))
            assert not rc, rc
            #Now check the md5...
            print(filename)
    h.close()

process_fastq(project, fastq_file)
