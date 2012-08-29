import os
import sys
import urllib

project = "ERP000297"

fastq_url = "http://www.ebi.ac.uk/ena/data/view/reports/sra/fastq_files/internal/%s" % project
fastq_file = "%s_fastq.tsv" % project


def download_in_one(url, filename):
    print "Fetching %s" % url
    n = urllib.urlopen(url)
    data = n.read()
    n.close()

    h = open(filename, "w")
    h.write(data)
    h.close()
    print "Saved as %s" % filename

if not os.path.isfile(fastq_file):
    download_in_one(fastq_url, fastq_file)

def process_fastq(project, fastq_filename):
    h = open(fastq_filename)
    line = h.readline()
    assert line == 'Study\tSample\tExperiment\tRun\tOrganism\tInstrument Platform\tInstrument Model\tLibrary Name\tLibrary Layout\tLibrary Source\tLibrary Selection\tRun Read Count\tRun Base Count\tFile Name\tFile Size\tmd5\tFtp\n', repr(line)
    for line in h:
        parts = line.rstrip("\n").split("\t")
        assert parts[0] == project
        url = parts[16]
        assert url.startswith("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR"), url
        filename = url[len("ftp://ftp.sra.ebi.ac.uk/"):]
        if os.path.isfile(filename):
            print "Already have %s" % filename
            continue
        if not filename.endswith(".fastq.gz"):
            print "Skipping %s" % filename
            continue
        #Make directory...
        d = os.path.split(filename)[0]
        if not os.path.isdir(d):
            print "Making directory %s" % d
            os.makedirs(d)
        #Download file...
        rc = os.system("wget -O %s %s" % (filename, url))
        assert not rc, rc
        #Now check the md5...
        print filename
    h.close()

process_fastq(project, fastq_file)
