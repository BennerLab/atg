import os
import subprocess
import multiprocessing
import tempfile
import shlex
import pysam
import pybedtools
import atg

BIGWIG_COMMAND = 'bedGraphToBigWig {bedgraph} {chrom_sizes} {bigwig}'
RPMM_SCALE_FACTOR = 10**6  # reads per million mapped reads

# TODO: add ChIP-Seq style
# TODO: parse ChIP-Seq/RNA-Seq from SAM/BAM header based on program used (PG?)
# TODO: parse genome from SAM/BAM header
# TODO: move get_chrom_sizes to some kind of data configuration module


def get_chrom_sizes(genome):
    """

    :param genome: a UCSC genome abbreviation, e.g. hg19
    :return: absolute path to chromosome sizes for specified genome
    """
    genome_dict = {'hg19': 'human', 'hg38': 'human', 'mm10': 'mouse', 'mm9': 'mouse', 'papAnu2': 'baboon'}
    data_root = os.path.expanduser(atg.config.settings['Data']['Root'])
    chrom_sizes_path = os.path.join(data_root, genome_dict[genome], 'Current', genome,
                                    'chrom.sizes')

    return os.path.abspath(chrom_sizes_path)


def get_genome_coverage(filename, scale, all_positions=False):
    if all_positions:
        return pybedtools.BedTool(filename).genome_coverage(bga=True, split=True, scale=scale).sort()
    else:
        return pybedtools.BedTool(filename).genome_coverage(bg=True, split=True, scale=scale).sort()


def prepend_chr(bedtool_feature):
    bedtool_feature.chrom = 'chr' + bedtool_feature.chrom
    return bedtool_feature


def translate_bedgraph(input_bedtool, output_filename, chrom_sizes_path):
    """
    Simple translation of a bedgraph from Ensembl chromosome numbering (e.g. 1) to UCSC (e.g. chr1). Any scaffolds or
    haplotype blocks will be ignored (not written to the output).

    :param input_bedtool: a BedTool object from pybedtools
    :param output_filename:
    :param chrom_sizes_path: path to chrom.sizes table
    :return:
    """
    chromosome_set = set()

    for line in open(chrom_sizes_path):
        chrom_name = line.strip().split()[0]
        chromosome_set.add(chrom_name)

    ucsc_bedgraph = input_bedtool.each(prepend_chr).filter(lambda x: x.chrom in chromosome_set)
    ucsc_bedgraph.saveas(output_filename)


class UnstrandedCoverageCalculator:
    def __init__(self, bam_filename):
        self.working_directory = tempfile.TemporaryDirectory()
        self.output_path = self.working_directory.name

        self.alignment_file = pysam.AlignmentFile(bam_filename)
        # index BAM file and re-load if there is no index already
        if not self.alignment_file.has_index():
            pysam.index(bam_filename)
            self.alignment_file = pysam.AlignmentFile(bam_filename)

        self.mapped_read_count = 0
        self.is_paired = False

        self.forward_strand_filename = os.path.join(self.output_path, 'fwd.bam')
        self.process_reads()

    def process_reads(self):
        """
        Create BAM file containing only uniquely mapped reads, necessary prior to bedgraph generation.
        :return:
        """
        forward_strand = pysam.AlignmentFile(self.forward_strand_filename, 'wb', template=self.alignment_file)
        for read in self.alignment_file.fetch():
            if read.mapq == 255:
                self.mapped_read_count += 1
                forward_strand.write(read)
                if read.is_read2:
                    self.is_paired = True

        # For paired-end runs, provide a rough count of fragments rather than reads
        if self.is_paired:
            self.mapped_read_count //= 2

    def write_bedgraph(self, forward_output=None, reverse_output=None, use_multiprocessing=False, scale=True,
                       all_positions=False):
        """

        :return:
        :param forward_output:
        :param reverse_output:
        :param use_multiprocessing:
        :param scale:
        :param all_positions:
        :return:
        """

        # scale coverage by million mapped reads/fragments
        if scale:
            scale_factor = RPMM_SCALE_FACTOR / self.mapped_read_count
        else:
            scale_factor = 1.0
        bedgraph_forward = get_genome_coverage(self.forward_strand_filename, scale_factor, all_positions)

        if forward_output:
            bedgraph_forward.moveto(forward_output)

        return bedgraph_forward

    def _check_chromosome_names(self, chrom_sizes_path):
        # check UCSC chromosome sizes file
        self.ucsc_chromosomes = set()
        for line in open(chrom_sizes_path):
            chrom_name = line.strip().split()[0]
            self.ucsc_chromosomes.add(chrom_name)

        # find chromosome names in BAM file
        alignment_reference_sequence_list = self.alignment_file.header['SQ']
        alignment_reference_seq_set = set()
        for reference_entry in alignment_reference_sequence_list:
            alignment_reference_seq_set.add(reference_entry['SN'])

        # any chromosomes in common?
        if len(alignment_reference_seq_set.intersection(self.ucsc_chromosomes)) == 0:
            return False

        return True

    def write_bigwig(self, forward_output, genome):
        chrom_sizes_path = get_chrom_sizes(genome)
        bedgraph_forward = self.write_bedgraph()

        if self._check_chromosome_names(chrom_sizes_path):
            bedgraph_forward_filename = bedgraph_forward.fn
        else:
            bedgraph_forward_filename = os.path.join(self.output_path, 'forward_ucsc.bedgraph')
            translate_bedgraph(bedgraph_forward, bedgraph_forward_filename, chrom_sizes_path)

        bigwig_command_forward = shlex.split(BIGWIG_COMMAND.format(bedgraph=bedgraph_forward_filename,
                                                                   chrom_sizes=chrom_sizes_path,
                                                                   bigwig=forward_output))
        subprocess.check_call(bigwig_command_forward)


class StrandedCoverageCalculator(UnstrandedCoverageCalculator):
    def __init__(self, bam_filename):
        super().__init__(bam_filename)
        self.reverse_strand_filename = os.path.join(self.output_path, 'rev.bam')
        self.process_stranded_reads()

    def process_reads(self):
        pass

    def process_stranded_reads(self):
        """
        Separate mapped reads into forward/reverse strand. For single-end runs, this only requires flipping the read
        strand. For paired-end runs, read1 is flipped, while read2 retains its strandedness.
        :return:
        """
        forward_strand = pysam.AlignmentFile(self.forward_strand_filename, 'wb', template=self.alignment_file)
        reverse_strand = pysam.AlignmentFile(self.reverse_strand_filename, 'wb', template=self.alignment_file)

        self.mapped_read_count = 0
        is_paired = False

        # sort reads into appropriate file
        # NOTE: only uniquely mapped reads (according to STAR) are retained
        for read in self.alignment_file.fetch():
            if read.mapq == 255:
                self.mapped_read_count += 1
                if read.is_read2:
                    is_paired = True
                    if read.is_reverse:
                        reverse_strand.write(read)
                    else:
                        forward_strand.write(read)
                else:
                    if read.is_reverse:
                        forward_strand.write(read)
                    else:
                        reverse_strand.write(read)

        # For paired-end runs, provide a rough count of fragments rather than reads
        if is_paired:
            self.mapped_read_count //= 2

        forward_strand.close()
        reverse_strand.close()

    def write_bedgraph(self, forward_output=None, reverse_output=None, use_multiprocessing=False, scale=True,
                       all_positions=False):
        """

        :return:
        :param forward_output:
        :param reverse_output:
        :param use_multiprocessing:
        :param scale:
        :param all_positions:
        :return:
        """

        # scale coverage by million mapped reads/fragments
        if scale:
            scale_factor = RPMM_SCALE_FACTOR/self.mapped_read_count
        else:
            scale_factor = 1.0
        task_list = [(self.forward_strand_filename, scale_factor, all_positions),
                     (self.reverse_strand_filename, -scale_factor, all_positions)]

        if use_multiprocessing:
            with multiprocessing.Pool() as pool:
                bedgraph_forward, bedgraph_reverse = pool.starmap(get_genome_coverage, task_list)
        else:
            bedgraph_forward, bedgraph_reverse = [get_genome_coverage(*i) for i in task_list]

        if forward_output:
            bedgraph_forward.moveto(forward_output)
        if reverse_output:
            bedgraph_reverse.moveto(reverse_output)

        return bedgraph_forward, bedgraph_reverse

    def write_bigwig(self, forward_output, reverse_output, genome):
        chrom_sizes_path = get_chrom_sizes(genome)
        bedgraph_forward, bedgraph_reverse = self.write_bedgraph()

        if self._check_chromosome_names(chrom_sizes_path):
            bedgraph_forward_filename = bedgraph_forward.fn
            bedgraph_reverse_filename = bedgraph_reverse.fn
        else:
            bedgraph_forward_filename = os.path.join(self.output_path, 'forward_ucsc.bedgraph')
            bedgraph_reverse_filename = os.path.join(self.output_path, 'reverse_ucsc.bedgraph')
            translate_bedgraph(bedgraph_forward, bedgraph_forward_filename, chrom_sizes_path)
            translate_bedgraph(bedgraph_reverse, bedgraph_reverse_filename, chrom_sizes_path)

        bigwig_command_forward = shlex.split(BIGWIG_COMMAND.format(bedgraph=bedgraph_forward_filename,
                                                                   chrom_sizes=chrom_sizes_path,
                                                                   bigwig=forward_output))
        bigwig_command_reverse = shlex.split(BIGWIG_COMMAND.format(bedgraph=bedgraph_reverse_filename,
                                                                   chrom_sizes=chrom_sizes_path,
                                                                   bigwig=reverse_output))
        subprocess.check_call(bigwig_command_forward)
        subprocess.check_call(bigwig_command_reverse)

    def __del__(self):
        pybedtools.cleanup()


def generate_bedgraph(namespace):
    bam_filename = namespace.bam_filename
    output_filename_list = namespace.output_filename
    scale = namespace.scale
    all_positions = namespace.all

    if len(output_filename_list) == 1:
        calculator = UnstrandedCoverageCalculator(bam_filename)
        calculator.write_bedgraph(output_filename_list[0], scale=scale, all_positions=all_positions)
    elif len(output_filename_list) == 2:
        calculator = StrandedCoverageCalculator(bam_filename)
        calculator.write_bedgraph(output_filename_list[0], output_filename_list[1], scale=scale,
                                  all_positions=all_positions)
    else:
        print("Please enter only 1 or 2 filenames.")
        return

def setup_subparsers(subparsers):
    coverage_parser = subparsers.add_parser('coverage', help='Generate a bedgraph from a BAM file.')

    coverage_parser.add_argument('bam_filename', help="STAR-mapped BAM file")
    coverage_parser.add_argument('output_filename', nargs='+',
                                 help="Name(s) for output bedgraph file(s). If one filename is given, output will be "
                                      "unstranded. If two filenames are given, the first will contain the positive"
                                      "strand and the second will contain the negative strand.")
    coverage_parser.add_argument('-s', '--scale', action='store_true', help="Scale coverage by million mapped reads")
    coverage_parser.add_argument('-a', '--all', action='store_true', help="Include all positions in bedgraph output"
                                                                          "(i.e., positions without coverage).")

    coverage_parser.set_defaults(func=generate_bedgraph)