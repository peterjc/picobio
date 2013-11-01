import os
from Bio import Entrez
from Bio import TogoWS
from Bio import SeqIO
from StringIO import StringIO
Entrez.email="peter.cock@hutton.ac.uk"

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
    search_text="txid%s[orgn] AND complete[Properties]" % taxon_id
        
    handle = Entrez.esearch("nucleotide", term=search_text, usehistory=True)
    search_results = Entrez.read(handle)
    handle.close()
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    count = int(search_results["Count"])
    print "%i hits" % count
    names = []

    #Get the accessions...
    batch_size = 1000
    for start in range(0,count,batch_size):
        end = min(count, start+batch_size)
        print("Getting accessions for record %i to %i" % (start+1, end))
        fetch_handle = Entrez.efetch(db="nucleotide", rettype="acc", retmode="xml",
                                     retstart=start, retmax=batch_size,
                                     webenv=webenv, query_key=query_key)
        data = fetch_handle.read().strip().split()
        fetch_handle.close()
        assert len(data) == end - start
        names.extend(data)
    assert len(names) == count
    print("%s to %s" % (names[0], names[-1]))

    handle = open("GenBank/%s.txt" % name,"w")
    handle.write("\n".join(names))
    handle.close()

    #Get the sequences
    for acc in names:
        name, version = acc.split(".")
        filename = "GenBank/%s.gbk" % name
        if not os.path.isfile(filename) :
            print("Fetching %s" % (name))
            #Fails, seems efetch now requires using history :(
            #data = Entrez.efetch("genome", rettype="gb", id=acc).read()
            fetch_handle = TogoWS.entry("nuccore", acc)
            data = fetch_handle.read() # defaults to gb
            fetch_handle.close()
            assert data.lstrip().startswith("LOCUS "), data
            assert data.rstrip().endswith("//"), data
            #Test we can parse it:
            record = SeqIO.read(StringIO(data),"gb")
            assert name == record.name, "Got %r expected %r" % (record.name, name)
            #Save it:
            handle = open(filename, "w")
            handle.write(data)
            handle.close()
            #print "%s saved" % record.name
        #Test we can parse it:
        if name not in checked :
            record = SeqIO.read(open(filename),"gb")
            assert name == record.name, "Got %r from %r expected %r" % (record.name, filename, name)
            print("Verified %s" % name)
            checked.add(name)
