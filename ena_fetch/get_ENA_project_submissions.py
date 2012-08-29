import os
import sys
import urllib

project = "ERP000297"

submissions_url = "http://www.ebi.ac.uk/ena/data/view/reports/sra/submitted_files/internal/%s" % project
submissions_file = "%s_submissions.tsv" % project


def download_in_one(url, filename):
    print "Fetching %s" % url
    n = urllib.urlopen(url)
    data = n.read()
    n.close()

    h = open(filename, "w")
    h.write(data)
    h.close()
    print "Saved as %s" % filename

if not os.path.isfile(submissions_file):
    download_in_one(submissions_url, submissions_file)

def process_submissions(project, submissions_filename):
    h = open(submissions_filename)
    line = h.readline()
    assert line == 'Study\tSample\tExperiment\tRun\tOrganism\tInstrument Platform\tInstrument Model\tLibrary Name\tLibrary Layout\tLibrary Source\tLibrary Selection\tRun Read Count\tRun Base Count\tFile Name\tFile Size\tmd5\tFtp\n', repr(line)
    for line in h:
        parts = line.rstrip("\n").split("\t")
        assert parts[0] == project
        url = parts[16]
        assert url.startswith("ftp://ftp.sra.ebi.ac.uk/vol1/ERA")
        filename = url[len("ftp://ftp.sra.ebi.ac.uk/"):]
        if os.path.isfile(filename):
            print "Already have %s" % filename
            continue
        if filename.endswith(".srf"):
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

process_submissions(project, submissions_file)
