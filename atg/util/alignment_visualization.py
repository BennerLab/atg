import os
import sys
import pandas
import argparse
import warnings
import atg.data
import atg.data.identifiers
from rpy2.robjects import r
from rpy2.robjects.packages import importr
from rpy2.rinterface import RRuntimeWarning


def visualize_locus(locus, alignment_list, flank_size, species, output, ucsc=False):
    if output.endswith('.png'):
        output_command = 'png("%s", height=600, width=800)' % output
    elif output.endswith('.pdf'):
        output_command = 'pdf("%s", height=6, width=8)' % output
    else:
        print('Output filename must end in .png or .pdf\n', file=sys.stderr)
        return

    # get relevant species information
    if not atg.data.config.has_section(species):
        print('%s is not an available species. Try one of the following:\n%s\n' %
              (species, ' '.join(atg.data.config.sections())), file=sys.stderr)
        return
    ensembl_shortname = atg.data.config.get(species, 'ensembl_shortname')
    ensembl_genome = atg.data.config.get(species, 'ensembl_genome')

    # get genome range from locus input
    if ':' in locus:
        chromosome, start, end = locus.replace(':', '-').split('-')
    else:
        translator = atg.data.identifiers.GeneIDTranslator(species)
        ensembl_id = translator.get_ensembl_id(locus)
        if not ensembl_id:
            print("Couldn't find a gene ID for %s in %s.\n" % (locus, species))
            return
        gene_annotation = pandas.read_csv(os.path.join(atg.data.genome_path[species], 'ensembl_gene.csv'), index_col=0)
        chromosome, start, end = tuple(gene_annotation.loc[ensembl_id].iloc[1:4])
        start -= flank_size
        end += flank_size

    # import libraries and set options
    # suppress warnings from Gviz and biomaRt
    warnings.filterwarnings('ignore', category=UserWarning)
    warnings.filterwarnings("ignore", category=RRuntimeWarning)
    biomart = importr('biomaRt')
    Gviz = importr('Gviz')
    warnings.resetwarnings()

    if not ucsc:
        r('options(ucscChromosomeNames=FALSE)')

    # get annotation for relevant region
    r('mart <- useMart(biomart="ensembl", dataset="%s_gene_ensembl")' % (ensembl_shortname))
    r('biomTrack <- BiomartGeneRegionTrack(genome = "%s", chromosome = "%s", start = %s, end = %s, '
      'name = "Genes", biomart = mart, collapseTranscripts="longest", transcriptAnnotation="symbol")' %
      (ensembl_genome, chromosome, start, end))

    # construct list of tracks
    r('track_list <- list(biomTrack)')

    for i, alignment_filename in enumerate(alignment_list):
        alignment_basename = os.path.splitext(os.path.basename(alignment_filename))[0]
        r('track_list[%d] <- AlignmentsTrack("%s", name="%s")' % (i+2, alignment_filename, alignment_basename))

    r(output_command)
    r('plotTracks(track_list, chromosome="%s", from=%s, to=%s)' % (chromosome, start, end))
    r('dev.off()')


def run_visualization(namespace):
    visualize_locus(namespace.locus, namespace.alignment_files, namespace.flank, namespace.organism, namespace.output,
                    namespace.ucsc)


def setup_subparsers(subparsers):
    retrieval_parser = subparsers.add_parser('gviz', help="Visualize alignments around a gene or genome range")
    retrieval_parser.add_argument('locus', help="a gene symbol (e.g. ACTB), Gene ID (ENSG00000075624), or\n"
                                                "genome range (e.g. 7:5527151-5563784)")
    retrieval_parser.add_argument('output', help="filename for output, ending with .png or .pdf")
    retrieval_parser.add_argument('alignment_files', nargs='*', help='Sorted BAM files')
    retrieval_parser.add_argument('-f', '--flank', type=int, default=2000,
                                  help="number of flanking nucleotides to show (default:2000)")
    retrieval_parser.add_argument('-o', '--organism', default="human", help="Common name of a species (default: human)")
    retrieval_parser.add_argument('-u', '--ucsc', action='store_true', help="Alignments use UCSC chromosome names")
    retrieval_parser.set_defaults(func=run_visualization)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    main_subparsers = parser.add_subparsers(dest="command", help='commands')
    main_subparsers.required = True
    setup_subparsers(main_subparsers)

    args = parser.parse_args()
    args.func(args)