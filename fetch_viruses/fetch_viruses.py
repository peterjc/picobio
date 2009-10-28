import os
from Bio import Entrez
from Bio import SeqIO
from StringIO import StringIO
Entrez.email="peter.cock@scri.ac.uk"

name = "ssRnaViruses"
taxon_id = "439488"
search_text="txid%s[orgn] AND complete[prop]"%taxon_id
handle = Entrez.esearch("genome",term=search_text, usehistory=True)
search_results = Entrez.read(handle)
handle.close()
webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]
count = int(search_results["Count"])
print "%i hits" % count

data = Entrez.efetch("genome", rettype="brief", retstart=0, retmax=count,
                     webenv=webenv, query_key=query_key).read()
names = []
for line in data.split("\n") :
    line = line.rstrip()
    if not line : continue
    index,acc,rest = line.split(None,2)
    acc = acc.rstrip(".")
    if acc not in names :
        names.append(acc)
    else :
        print "Duplicated: %s" % line
assert len(names)==count

handle = open("GenBank/%s.txt" % name,"w")
handle.write("\n".join(names))
handle.close()

for index, name in enumerate(names):
    filename = "GenBank/%s.gbk" % name
    if not os.path.isfile(filename) :
        print "Fetching %i of %i" % (index+1, count)
        data = Entrez.efetch("genome", rettype="gb", retstart=index, retmax=1,
                             webenv=webenv, query_key=query_key).read()
        assert data.lstrip().startswith("LOCUS ")
        assert data.rstrip().endswith("//")
        #Test we can parse it:
        record = SeqIO.read(StringIO(data),"gb")
        assert names[index]==record.name
        #Save it:
        handle = open(filename, "w")
        handle.write(data)
        handle.close()
        print "%s saved" % record.name
    #Test we can parse it:
    record = SeqIO.read(open(filename),"gb")
    assert names[index]==record.name
    print "Got %s" % name
    
