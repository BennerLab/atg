import os
import pysam


def get_genome_stats(genome_fasta):
    reference_fasta_index = genome_fasta + '.fai'
    if not os.path.exists(reference_fasta_index):
        print("\nIndexing %s\n" % os.path.abspath(genome_fasta))
        pysam.faidx(genome_fasta)

    reference_genome = pysam.FastaFile(genome_fasta)
    total_bases = sum(reference_genome.lengths)

    return reference_genome.nreferences, total_bases

