import os
import pandas
import argparse
import atg.data.identifiers
from rpy2.robjects import r
from rpy2.robjects.packages import importr


def visualize_locus(locus, flank_size, species, output):
    # get genome range from locus input
    if ':' in locus:
        genome_location = locus
        chromosome, start, end = locus.replace(':', '-').split('-')
    else:
        translator = atg.data.identifiers.GeneIDTranslator(species)
        ensembl_id = translator.get_ensembl_id(locus)
        gene_annotation = pandas.read_csv(os.path.join(atg.data.genome_path[species], 'ensembl_gene.csv'), index_col=0)
        chromosome, start, end = tuple(gene_annotation.loc[ensembl_id].iloc[1:4])
        start -= flank_size
        end += flank_size

    print(chromosome, start, end)

    biomart = importr('biomaRt')
    Gviz = importr('Gviz')

    r('options(ucscChromosomeNames=FALSE)')
    r('pdf("/tmp/gviz_test.pdf", height=3, width=4)')

    r('mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")')
    r('biomTrack <- BiomartGeneRegionTrack(genome = "GRCh38", chromosome = "%s", start = %d, end = %d, '
      'name = "Genes", biomart = mart, collapseTranscripts="longest", transcriptAnnotation="symbol")' %
      (chromosome, start, end))

    r('plotTracks(list(biomTrack), chromosome="%s", from=%d, to=%d)' % (chromosome, start, end))
    r('dev.off()')


def run_visualization(namespace):
    visualize_locus(namespace.locus, namespace.flank, namespace.organism, namespace.output)


def setup_subparsers(subparsers):
    retrieval_parser = subparsers.add_parser('gviz', help="Visualize alignments around a gene or genome range")
    retrieval_parser.add_argument('locus', help="a gene symbol (e.g. ACTB), Gene ID (ENSG00000075624), or\n"
                                                "genome range (e.g. 7:5527151-5563784)")
    retrieval_parser.add_argument('output', help="filename for output, ending with .png or .pdf")
    retrieval_parser.add_argument('alignment_files', nargs='*', help='Sorted BAM files')
    retrieval_parser.add_argument('-f', '--flank', type=int, default=2000,
                                  help="number of flanking nucleotides to show (default:2000)")
    retrieval_parser.add_argument('-o', '--organism', default="human")
    retrieval_parser.set_defaults(func=run_visualization)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    main_subparsers = parser.add_subparsers(dest="command", help='commands')
    main_subparsers.required = True
    setup_subparsers(main_subparsers)

    args = parser.parse_args()
    args.func(args)