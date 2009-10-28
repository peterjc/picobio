import os
from Bio import Entrez
from Bio import SeqIO
from StringIO import StringIO
Entrez.email="peter.cock@scri.ac.uk"

checked = set()
# dsDNA viruses, no RNA stage, Taxonomy ID: 35237
# dsRNA viruses, Taxonomy ID: 35325
# ssDNA viruses, Taxonomy ID: 29258
# ssRNA viruses, Taxonomy ID: 439488
# Viruses, Taxonomy ID: 10239
for name, taxon_id in [
    ("dsDnaViruses", "35237"),
    ("dsRnaViruses", "35325"),
    ("ssDnaViruses", "29258"),
    ("ssRnaViruses", "439488"),
    ("allViruses","10239")] :
    print "="*60
    print name
    print "="*60
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
    assert data, data
    names = []
    for line in data.split("\n") :
        line = line.rstrip()
        if not line : continue
        try :
            index,acc,rest = line.split(None,2)
        except ValueError, err :
            print line
            raise err
        acc = acc.rstrip(".")
        if acc not in names :
            names.append(acc)
        else :
            print "Duplicated: %s" % line
    assert len(names)==count, "%i vs %i" % (len(names), count)

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
            #print "%s saved" % record.name
        #Test we can parse it:
        if name not in checked :
            record = SeqIO.read(open(filename),"gb")
            assert names[index]==record.name
            print "Verified %s" % name
            checked.add(name)
    
