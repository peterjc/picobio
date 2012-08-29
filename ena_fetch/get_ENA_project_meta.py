import os
import sys
import urllib

project = "ERP000297"
strain_file = "%s_strain.tsv" % project #output file

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

def get_strain(meta_xml_filename):
    h = open(meta_xml_filename)
    while True:
        line = h.readline()
        if not line:
            break
        if "<tag>strain</tag>" in line.lower():
            strain = h.readline().strip()
            assert strain.lower().startswith("<value>"), strain
            assert strain.lower().endswith("</value>"), strain
            h.close()
            return strain[7:-8]
    h.close()
    return None

def process_meta(project, fastq_filename, strain_file):
    h = open(fastq_filename)
    out = open(strain_file,  "w")
    line = h.readline()
    assert line == 'Study\tSample\tExperiment\tRun\tOrganism\tInstrument Platform\tInstrument Model\tLibrary Name\tLibrary Layout\tLibrary Source\tLibrary Selection\tRun Read Count\tRun Base Count\tFile Name\tFile Size\tmd5\tFtp\n', repr(line)
    out.write(line[:-1] + "\tStrain\n")
    for line in h:
        parts = line.rstrip("\n").split("\t")
        assert parts[0] == project
        url = parts[16]
        assert url.startswith("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR"), url

        sample = parts[1]
        assert sample.startswith("ERS")
        url = "http://www.ebi.ac.uk/ena/data/view/%s&display=xml" % sample
        url = "http://www.ebi.ac.uk/ena/data/view/%s&display=xml&download" % sample
        filename = "xml/%s.xml" % sample

        #Download file...
        if not os.path.isfile(filename):
            print url
            rc = os.system("wget -O %s '%s'" % (filename, url))
            assert not rc, rc

        strain = get_strain(filename)
        if not strain: strain = ""
        print filename, strain
        out.write(line[:-1] + "\t" + strain + "\n")
    h.close()
    out.close()

process_meta(project, fastq_file, strain_file)
