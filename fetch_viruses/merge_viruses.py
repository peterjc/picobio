#Prepares merged FASTA files to use for BLAST databases
#
# v000 - proteins only
# v001 - date based filename
#      - CDS nuc sequences too
# v002 - Use BLAST friendly names
# v003 - multiple sets of viruses
# v004 - fixed missing | in fna names
# v005 - Handle feature extraction via Biopython
#      - Tested under Python 3
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

date_stamp = "20131105"

tables = {"NC_008956":1, "NC_008954":1, "NC_008949":1, "NC_008948":1,
          "NC_011452":1, "NC_008956":1}

def dedup(input_fasta, output_fasta):
    """Merge identical FASTA entries using NCBI NR Ctrl+A style"""
    #Try it in memory...
    print("Deduplicating %s..." % input_fasta)
    by_seq = dict()
    total = 0
    for record in SeqIO.parse(input_fasta, "fasta"):
        total += 1
        s = str(record.seq).upper()
        t = record.description
        try:
            by_seq[s].append(t)
        except KeyError:
            by_seq[s] = [t]
    print("Reduced %s from %i to %i records as %s (%0.1f%%)"
          % (input_fasta, total, len(by_seq), output_fasta,
             len(by_seq) * 100.0 / float(total)))
    handle = open(output_fasta, "w")
    for s in by_seq:
        titles = by_seq[s]
        #chr(1) = CTRL+A
        handle.write(">%s\n%s\n" % (chr(1).join(sorted(titles)), s))
    handle.close()


def get_nuc(seq, loc_string) :
    reverse = False
    if loc_string.startswith("complement(") :
        assert loc_string[-1]==")"
        loc_string = loc_string[11:-1]
        reverse = True
    start, end = [int(x.strip("<").strip(">")) for x in loc_string.split("..")]
    nuc = seq[start-1:end]
    if reverse :
        return nuc.reverse_complement()
    else :
        return nuc


for group in ["dsDnaViruses",
              "ssDnaViruses",
              "dsRnaViruses",
              "ssRnaViruses",
              "allViruses"] :
    print("="*len(group))
    print(group)
    print("="*len(group))
    names = open("GenBank/%s.txt" % group).read().split("\n")
    genomes_file = "%s_%s_genomes.fna" % (group, date_stamp)
    genomes_nr = "%s_%s_genomes_NR.fna" % (group, date_stamp)
    protein_file = "%s_%s_proteins.faa" % (group, date_stamp)
    protein_nr = "%s_%s_proteins_NR.faa" % (group, date_stamp)
    nuc_file = "%s_%s_genes.ffn" % (group, date_stamp)
    nuc_nr = "%s_%s_genes_NR.ffn" % (group, date_stamp)
    print("Looking at %i %s" % (len(names), group))


    if os.path.isfile(genomes_file):
        print("Got %s" % genomes_file)
    else:
        print("Writing %s..." % genomes_file)
        records = (SeqIO.read("GenBank/%s.gbk" % (acc+".").split(".")[0], "gb") for acc in names)
        count = SeqIO.write(records, genomes_file, "fasta")
        print("%i records in %s" % (count, genomes_file))

    if os.path.isfile(genomes_nr):
        print("Got %s" % genomes_nr)
    else:
        dedup(genomes_file, genomes_nr)


    if os.path.isfile(protein_file):
        print("Got %s" % protein_file)
    else :
        handle = open(protein_file,"w")
        bad = 0
        count = 0
        for index, name in enumerate(names):
            name = name.split(".", 1)[0]
            filename = "GenBank/%s.gbk" % name
            parent = None
            for record in SeqIO.parse(open(filename),"genbank-cds") :
                if "pseudo" in record.annotations : continue
                if "pseudogene" in record.annotations: continue
                count+=1
                try :
                    protein_id = record.annotations["protein_id"]
                except KeyError:
                    print(filename)
                    print(record)
                    assert False
                gi = None
                for xref in record.dbxrefs :
                    if xref.lower().startswith("gi:") :
                        gi = xref[3:]
                        break
                assert gi and protein_id, str(record)
                record.id = "gi|%s|ref|%s" % (gi, record.id)
                if record.description=="<unknown description>":
                    if "product" in record.annotations :
                        record.description = record.annotations["product"]
                    elif "note" in record.annotations :
                        record.description = record.annotations["note"]
                if record.seq is None :
                    bad+=1
                    print("%s %s" % (filename, record.annotations["raw_location"]))
                    if parent is None :
                        parent = SeqIO.read(open(filename),"gb")
                    nuc = get_nuc(parent.seq, record.annotations["raw_location"])
                    if "transl_table" in record.annotations :
                        table = int(record.annotations["transl_table"])
                    else :
                        table = tables[name]
                    pro = nuc.translate(table)
                    assert pro.endswith("*") and pro.count("*")==1
                    record.seq = pro[:-1] #remove stop
                SeqIO.write([record], handle, "fasta")
            #print("%i: %i in %s" % (index+1, count, name))
        handle.close()
        print("Done")
        print("%i proteins" % count)
        print("%i missing provided translation" % bad)

    if os.path.isfile(protein_nr):
        print("Got %s" % protein_nr)
    else:
        dedup(protein_file, protein_nr)
        
    if os.path.isfile(nuc_file):
        print("Got %s" % nuc_file)
    else :
        handle = open(nuc_file,"w")
        count = 0
        for index, name in enumerate(names):
            name = name.split(".", 1)[0]
            filename = "GenBank/%s.gbk" % name
            #print(name)
            parent = SeqIO.read(open(filename),"genbank")
            for f in parent.features :
                if f.type != "CDS" : continue
                if "pseudo" in f.qualifiers : continue
                if "pseudogene" in f.qualifiers : continue
                nuc = f.extract(parent.seq)
                protein_id = f.qualifiers["protein_id"][0]
                gi = None
                pro = nuc.translate(tables.get(name,1))
                if not (pro.endswith("*") and pro.count("*")==1) :
                    print("%s %s lacks stop codon" % (name, protein_id))
                for xref in f.qualifiers["db_xref"] :
                    if xref.lower().startswith("gi:") :
                        gi = xref[3:]
                        break
                if not (gi and protein_id) :
                    print(f)
                    assert False
                #Bit of a hack, we are using the protein's ID here!
                record = SeqRecord(nuc, id="gi|%s|ref|%s" % (gi, protein_id),
                                   description="; ".join(f.qualifiers.get("note",[])))
                SeqIO.write([record], handle, "fasta")
                count +=1
            #print "%i: %i in %s" % (index+1, count, name)
        handle.close()
        print("Done")
        print("%i genes" % count)

    if os.path.isfile(nuc_nr):
        print("Got %s" % nuc_nr)
    else:
        dedup(nuc_file, nuc_nr)
