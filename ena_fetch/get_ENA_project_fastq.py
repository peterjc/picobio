from __future__ import print_function

import os
try:
    from urllib.request import urlopen
except ImportError:
    from urllib import urlopen

project = "PRJEB2896"

# This was simple, but does not work anymore:
# fastq_url = "http://www.ebi.ac.uk/ena/data/view/reports/sra/fastq_files/internal/%s" % project
# This is the default column set:
# fastq_url = "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJEB2896&result=read_run&fields=study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp,fastq_galaxy,submitted_ftp,submitted_galaxy,sra_ftp,sra_galaxy,cram_index_ftp,cram_index_galaxy&download=txt"
# Or goto https://www.ebi.ac.uk/ena/data/view/PRJEB2896 and click on "TEXT" download link
#
# This includes important metadata including the fastq_md5 information,
fields = "study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_name,nominal_length,library_layout,read_count,experiment_title,study_title,study_alias,experiment_alias,run_alias,fastq_md5,fastq_ftp,submitted_md5,submitted_ftp,sra_md5,sra_ftp,cram_index_ftp".split(",")

fastq_url = 'https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=%s&result=read_run&fields=%s&download=txt' % (project, ",".join(fields))

fastq_file = "%s_metadata.tsv" % project

# The Gp life-stages:
wanted = """ERS091755
ERS092427
ERS001595
ERS092081
ERS092426
ERS092348
ERS092526
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

# A few populations of interest:
#wanted += ["ERR202431", "ERR202432", "ERR202433", "ERR202434"]
wanted += ["ERS092428", "ERS092429", "ERS092430", "ERS092431"]


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
    assert line == "\t".join(fields) + "\n", repr(line)
    for line in h:
        parts = line.rstrip("\n").split("\t")
        assert parts[0] == project
        urls = parts[fields.index("fastq_ftp")].split(";")
        md5s = parts[fields.index("fastq_md5")].split(";")
        for url, md5 in zip(urls, md5s):
            if url.startswith("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR"):
                pass
            elif url.startswith("ftp.sra.ebi.ac.uk/vol1/fastq/ERR"):
                url = "ftp://" + url
            assert url.startswith("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR"), url
            filename = url[len("ftp://ftp.sra.ebi.ac.uk/"):]
            pending = filename + ".tmp"
            acc = parts[2]  # here using secondary_sample_accession
            if wanted and acc not in wanted:
                print("Not interested in %s from %s" % (filename, acc))
                continue
            if not filename.endswith(".fastq.gz"):
                print("Skipping %s" % filename)
                continue
            # Make directory...
            d = os.path.split(filename)[0]
            if not os.path.isdir(d):
                print("Making directory %s" % d)
                os.makedirs(d)
            if os.path.isfile(filename):
                print("Already have %s" % filename)
                # Assume MD5 checked
                continue
            # Download file...
            print("Downloading %s" % filename)
            rc = os.system("wget -nv -O %s %s" % (pending, url))
            assert not rc, rc
            # Now check the md5...
            m = pending + ".md5"
            if not os.path.isfile(m):
                print("Creating %s with md5 %s" % (m, md5))
                with open(m, "w") as handle:
                    handle.write("%s  %s" % (md5, os.path.basename(pending)))
            print("Confirming %s has checksum %s" % (pending, md5))
            rc = os.system("cd %s && md5sum -c %s" % (d, os.path.split(m)[1]))
            assert not rc, rc
            # Rename files now that MD5 confirmed
            os.remove(m)
            m = filename + ".md5"
            with open(m, "w") as handle:
                handle.write("%s  %s" % (md5, os.path.basename(filename)))
            os.rename(pending, filename)
            print("Renamed %s to %s" % (pending, filename))
    h.close()

process_fastq(project, fastq_file)
