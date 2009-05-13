#Prepares merged FASTA files to use for BLAST databases
#
# v000 - proteins only
# v001 - date based filename
#      - CDS nuc sequences too
# v002 - Use BLAST friendly names
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

tables = {"NC_008956":1, "NC_008954":1, "NC_008949":1, "NC_008948":1}
names = open("GenBank/dsDnaViruses.txt").read().split("\n")
protein_file = "BlastDB/dsDnaViruses_20090512_proteins.faa"
nuc_file = "BlastDB/dsDnaViruses_20090512_genes.ffn"

print "Looking at %i viruses" % len(names)

def get_nuc(seq, loc_string) :
    reverse = False
    if loc_string.startswith("complement(") :
        assert loc_string[-1]==")"
        loc_string = loc_string[11:-1]
        reverse = True
    start, end = map(int, loc_string.split(".."))
    nuc = seq[start-1:end]
    if reverse :
        return nuc.reverse_complement()
    else :
        return nuc

#TODO - Add this functionality to Biopython itself...
def get_feature_nuc(f, parent_seq) :
    if f.sub_features :
        if f.location_operator!="join":
            raise ValueError(f.location_operator)
        #TODO - This should recurse to cope with join(complement(...),...) properly
        #for mixed-strand features, BUT that is impossible with the current GenBank
        #parser due to how the strand is recorded on both the parent and subfeatures.
        f_subs = [parent_seq[f_sub.location.nofuzzy_start:f_sub.location.nofuzzy_end] \
                  for f_sub in f.sub_features]
        #f_subs = [get_feature_nuc(f_sub, parent_seq) for f_sub in f.sub_features]
        #TODO - Join support in Seq object?  But how to deal with alphabets...
        f_seq = Seq("".join(map(str,f_subs)),f_subs[0].alphabet)
    else :
        f_seq = parent_seq[f.location.nofuzzy_start:f.location.nofuzzy_end]
    if f.strand == -1 : f_seq = f_seq.reverse_complement()
    return f_seq

if os.path.isfile(protein_file) :
    print "Got %s" % protein_file
else :
    handle = open(protein_file,"w")
    bad = 0
    count = 0
    for index, name in enumerate(names):
        filename = "GenBank/%s.gbk" % name
        parent = None
        for record in SeqIO.parse(open(filename),"genbank-cds") :
            count+=1
            if record.seq is None :
                bad+=1
                print filename, record.annotations["raw_location"]
                if parent is None :
                    parent = SeqIO.read(open(filename),"gb")
                nuc = get_nuc(parent.seq, record.annotations["raw_location"])
                pro = nuc.translate(tables[name])
                assert pro.endswith("*") and pro.count("*")==1
                record.seq = pro
            if record.id.count("|")==0 :
                protein_id = record.annotations["protein_id"]
                gi = None
                for xref in record.dbxrefs :
                    if xref.lower().startswith("gi:") :
                        gi = xref[3:]
                        break
                assert gi and protein_id
                record.id = "gi|%s|ref|%s" % (gi, record.id)
            SeqIO.write([record], handle, "fasta")
        #print "%i: %i in %s" % (index+1, count, name)
    handle.close()
    print "Done"
    print "%i proteins" % count
    print "%i missing provided translation" % bad
    assert count == 73363

if os.path.isfile(nuc_file):
    print "Got %s" % nuc_file
else :
    handle = open(nuc_file,"w")
    count = 0
    for index, name in enumerate(names):
        filename = "GenBank/%s.gbk" % name
        print name
        parent = SeqIO.read(open(filename),"genbank")
        for f in parent.features :
            if f.type != "CDS" : continue
            nuc = get_feature_nuc(f, parent.seq)
            try :
                pro = nuc.translate(tables[name])
                assert pro.endswith("*") and pro.count("*")==1
            except KeyError, err :
                pass
            protein_id = f.qualifiers["protein_id"][0]
            gi = None
            for xref in f.qualifiers["db_xref"] :
                if xref.lower().startswith("gi:") :
                    gi = xref[3:]
                    break
            if not (gi and protein_id) :
                print f
                assert False
            #Bit of a hack, we are using the protein's ID here!
            record = SeqRecord(nuc, id="gi|%sref|%s" % (gi, protein_id),
                               description="; ".join(f.qualifiers.get("note",[])))
            SeqIO.write([record], handle, "fasta")
            count +=1
        #print "%i: %i in %s" % (index+1, count, name)
    handle.close()
    print "Done"
    print "%i genes" % count
    assert count == 73363
