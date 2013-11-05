import os
try:
    from StringIO import StringIO # Python 2
except ImportError:
    from io import StringIO # Python 3

from Bio import SeqIO
from Bio import TogoWS
from Bio import Entrez

Entrez.email="peter.cock@hutton.ac.uk"

def download(acc, name, filename):
    fetch_handle = Entrez.efetch("nuccore", rettype="gbwithparts", id=acc)
    #fetch_handle = TogoWS.entry("nuccore", acc)
    data = fetch_handle.read() # defaults to gb                                                                                                                                      
    fetch_handle.close()
    assert data.lstrip().startswith("LOCUS "), data
    assert data.rstrip().endswith("//"), data
    #Test we can parse it:                                                                                                                                                           
    record = SeqIO.read(StringIO(data),"gb")
    assert name == record.name or record.id.startswith(name + "."), \
        "Got %r and %r expected %r" % (record.id, record.name, name)
    #Save it:                                                                                                                                                                        
    handle = open(filename, "w")
    handle.write(data)
    handle.close()

def download_batch(acc_list, check=False):
    missing = []
    for acc in acc_list:
        if "." in acc:
            name, version = acc.split(".")
        else:
            name = acc
        filename = "GenBank/%s.gbk" % name
        if not os.path.isfile(filename):
            missing.append(acc)
        elif check:
            check(acc, name, filename)
    count = len(missing)
    for index, acc in enumerate(missing):
        if "." in acc:
            name, version = acc.split(".")
        else:
            name = acc
        filename = "GenBank/%s.gbk" % name
        assert not os.path.isfile(filename)
        print("Fetching %s (%i of %i)" % (name, index+1, count))
        download(acc, name, filename)

def check(acc, name, filename):
    record = SeqIO.read(open(filename),"gb")
    assert name == record.name or record.id.startswith(name + "."), \
        "Got %r and %r expected %r" % (record.id, record.name, name)

# dsDNA viruses, no RNA stage, Taxonomy ID: 35237
# dsRNA viruses, Taxonomy ID: 35325
# ssDNA viruses, Taxonomy ID: 29258
# ssRNA viruses, Taxonomy ID: 439488
# Viruses, Taxonomy ID: 10239
for group, taxon_id in [
        ("dsDnaViruses", "35237"),
        ("dsRnaViruses", "35325"),
        ("ssDnaViruses", "29258"),
        ("ssRnaViruses", "439488"),
        ("allViruses","10239"),
    ] :
    print("="*60)
    print(group)
    print("="*60)

    if os.path.isfile("GenBank/%s.txt" % group):
        print("Pre-fetching any outstanding old search results...")
        handle = open("GenBank/%s.txt" % group)
        names =[line.strip() for line in handle]
        handle.close()
        download_batch(names)
    else:
        names = []

    print("Running NCBI search...")
    search_text = "txid%s[orgn] AND complete[Properties] AND genome" % taxon_id
    handle = Entrez.esearch("nucleotide", term=search_text, usehistory=True)
    search_results = Entrez.read(handle)
    handle.close()
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    count = int(search_results["Count"])
    print("%i hits" % count)

    if len(names) == count:
        print("Probably no new names...")
        continue

    #Get the accessions...
    names = []
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
    print("%i records, %s to %s" % (len(names), names[0], names[-1]))

    handle = open("GenBank/%s.txt" % group,"w")
    handle.write("\n".join(names))
    handle.close()

    #Get the sequences
    download_batch(names)
