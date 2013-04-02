import os
import time

to_profile = []

try:
    from Bio.Sequencing.SamBam import BamIterator
    def peter_iter(bam_filename, out_filename):
        """Peter's pure Python BAM iterator."""
        h = open(bam_filename, "rb")
        out_h = open(out_filename, "w")
        count = 0
        mapped = 0
        for read in BamIterator(h):
            count += 1
            if read.is_mapped:
                mapped += 1
                out_h.write("%s\t%s\n" % (read.rname, read.pos))
        h.close()
        out_h.close()
        return mapped, count
    to_profile.append(peter_iter)
except ImportError:
    pass

try:
    from pysam import Samfile
    def pysam_iter(bam_filename, out_filename):
        """PySam's Samfile as BAM iterator."""
        out_h = open(out_filename, "w")
        count = 0
        mapped = 0
        for read in Samfile(bam_filename, "rb"):
            count += 1
            if not read.is_unmapped:
                mapped += 1
                out_h.write("%s\t%s\n" % (read.rname, read.pos))
        out_h.close()
        return mapped, count
    to_profile.append(pysam_iter)
except ImportError:
    pass

print "Will profile %i functions:" % len(to_profile)
for p in to_profile:
    print p.__doc__
print
for f in os.listdir("."):
    if f.endswith(".bam"):
        print "Using %s" % f
        for p in to_profile:
            print "Profiling %s" % p.__doc__
            start = time.time()
            mapped, count = p(f, "/dev/null")
            taken = time.time() - start
            print "%s - %0.1fs giving %i/%i mapped" \
                % (p.__doc__, taken, mapped, count)
