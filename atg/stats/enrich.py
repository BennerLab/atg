import pandas
import numpy
import os
import sys
import atg.data.identifiers
from scipy.stats import hypergeom
import seaborn as sns
import matplotlib.pyplot as plt

LOG_P_VALUE_MINIMUM = -80.0
LOG10_FACTOR = 1.0 / numpy.log(10)

GENE_COLUMN_LABEL_INDEX = 0
GO_COLUMN_LABEL_INDEX = 1
GO_DEFINITION_TERM_INDEX = 0
GO_DEFINITION_ACCESSION_INDEX = 1
GO_DEFINITION_DOMAIN_INDEX = 2
PLOTTED_CHARACTER_LIMIT = 40
DEFAULT_MAXIMUM_TERM_SIZE = 250
TERMINAL_OUTPUT_LINES = 20


def enrichment_significance(term_row):
    return LOG10_FACTOR * hypergeom.logsf(term_row['hit_count']-1, term_row['universe'],
                                          term_row['term_count'], term_row['list_size'])


def p_adjust_bh(p, n=None):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""

    p_array = numpy.asfarray(p)

    if n:
        p_array = numpy.append(p_array, numpy.ones(n - len(p)))

    by_descend = p_array.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p_array)) / numpy.arange(len(p_array), 0, -1)
    q = numpy.minimum(1, numpy.minimum.accumulate(steps * p_array[by_descend]))

    return q[by_orig][0:len(p)]


def check_significance_value(log_pvalue):
    if numpy.isneginf(log_pvalue) or log_pvalue < LOG_P_VALUE_MINIMUM:
        return LOG_P_VALUE_MINIMUM
    return log_pvalue


class EnrichmentCalculator:
    """
    calculate stats on statistical enrichment on gene categories, e.g. Gene Ontology
    """

    def __init__(self, gene_term_file, term_definition_file, maximum_size=DEFAULT_MAXIMUM_TERM_SIZE, minimum_size=3):
        self.term_definition = pandas.read_csv(term_definition_file).dropna()

        gene_term_df = pandas.read_csv(gene_term_file)
        self.gene_column_label = gene_term_df.columns[GENE_COLUMN_LABEL_INDEX]
        self.go_column_label = gene_term_df.columns[GO_COLUMN_LABEL_INDEX]

        term_count = gene_term_df[self.go_column_label].value_counts()
        relevant_terms = term_count[(term_count <= maximum_size) & (term_count >= minimum_size)]

        self.gene_term_df = gene_term_df[gene_term_df[self.go_column_label].isin(relevant_terms.index.values)]
        self.term_count = self.gene_term_df[self.go_column_label].value_counts()
        self.universe = self.gene_term_df[self.gene_column_label].nunique()

    def get_single_enrichment(self, gene_list, term, gene_universe=15000, gene_list_size=None):
        term_df = self.gene_term_df[self.gene_term_df[self.go_column_label] == term]
        n = len(term_df)
        x = sum(term_df[self.gene_column_label].isin(gene_list)) - 1  # matched genes
        N = gene_list_size if gene_list_size else len(gene_list)

        return hypergeom.logsf(x, gene_universe, n, N)

    def get_all_enrichment(self, gene_list, alternate_gene_list_size=None):
        """
        evaluate significance of the given gene list for all GO BP categories
        :param gene_list:
        :param alternate_gene_list_size: use a fixed number when calculating iterative enrichment
        :return: a DataFrame containing gene set counts and log p-values for all GO BP terms.
        """
        enrichment_df = pandas.DataFrame()
        enrichment_df['term_count'] = self.term_count
        enrichment_df['hit_count'] = (self.gene_term_df[self.gene_term_df[self.gene_column_label].isin(gene_list)]
                                                       [self.go_column_label].value_counts().astype(int))
        enrichment_df.dropna(inplace=True)

        enrichment_df['universe'] = self.universe
        enrichment_df['list_size'] = alternate_gene_list_size if alternate_gene_list_size else len(gene_list)
        enrichment_df['log_pvalue'] = enrichment_df.apply(enrichment_significance, axis=1, result_type='reduce')

        # replace negative infinite p-values with a small constant
        enrichment_df['log_pvalue'].replace(-numpy.inf, LOG_P_VALUE_MINIMUM, inplace=True)
        enrichment_df['log_adjusted_pvalue'] = numpy.log10(p_adjust_bh(numpy.power(10, enrichment_df['log_pvalue']),
                                                                       n=len(self.term_definition)))
        return enrichment_df.sort_values(['log_adjusted_pvalue', 'hit_count'], ascending=[True, False])

    def iterative_enrichment(self, gene_list, adj_log_p_value_threshold=-2.0):
        """

        :param gene_list:
        :param adj_log_p_value_threshold: stop performing enrichment calculations when the best adjusted p-value
        becomes larger than this value
        :return: a DataFrame containing raw and adjusted log p-values for GO BP terms.
        """
        current_gene_set = set(gene_list)
        full_gene_set_size = len(current_gene_set)
        enriched_term_results = []

        while len(current_gene_set) > 0:
            current_enrichment_df = self.get_all_enrichment(current_gene_set, full_gene_set_size)
            most_enriched_term_row = current_enrichment_df.iloc[0]
            if most_enriched_term_row['log_adjusted_pvalue'] > adj_log_p_value_threshold:
                break

            most_enriched_term = most_enriched_term_row.name
            enriched_term_results.append([most_enriched_term, most_enriched_term_row['log_pvalue'],
                                          most_enriched_term_row['log_adjusted_pvalue']])
            current_term_genes = self.gene_term_df.loc[self.gene_term_df[self.go_column_label] == most_enriched_term,
                                                       self.gene_column_label]
            # remove genes present in most enriched term
            current_gene_set.difference_update(current_term_genes.values)

        enrichment_df = pandas.DataFrame.from_records(enriched_term_results,
                                                      columns=['Go term name', 'log_pvalue', 'log_adjusted_pvalue'],
                                                      index='Go term name')
        return enrichment_df

    def plot_enrichment_single(self, gene_list, output_filename=None, log_pvalue_threshold=-2.0, num_terms=20,
                               iterative=False):

        if iterative:
            enrichment_df = self.iterative_enrichment(gene_list, log_pvalue_threshold)

        else:
            enrichment_df = self.get_all_enrichment(gene_list)

        annotated_enrichment_df = (enrichment_df.merge(self.term_definition, left_index=True,
                                                       right_on=self.term_definition.columns
                                                       [GO_DEFINITION_ACCESSION_INDEX])
                                                .query('log_adjusted_pvalue < @log_pvalue_threshold')
                                                .sort_values('log_adjusted_pvalue')
                                                .iloc[0:num_terms]
                                                .assign(display_log_pvalue=lambda x: x.log_adjusted_pvalue*-1))
        
        sns.set(style='whitegrid')
        plt.figure()
        plot = sns.stripplot(x='display_log_pvalue', y=self.term_definition.columns[GO_DEFINITION_TERM_INDEX],
                             data=annotated_enrichment_df, orient='h', edgecolor='gray', palette="Reds_r")
        plot.set_xlabel(r'$-log_{10}(p)$')
        plot.set_ylabel('')
        plot.xaxis.grid(False)
        plot.yaxis.grid(True)
        sns.despine(left=True, bottom=True)
        plt.tight_layout()

        if output_filename:
            plt.savefig(output_filename)
        else:
            plt.show()

    def iterative_enrichment_multilist(self, multi_gene_list, adj_log_p_value_threshold=-2.0):
        """
        perform iterative enrichment for multiple gene lists:
        1. 

        :param multi_gene_list:
        :param adj_log_p_value_threshold:
        :return:
        """

        # set up sets and track original list size
        gene_set_dict = {}
        gene_set_size = {}
        for set_name, gene_list in multi_gene_list.items():
            gene_set_dict[set_name] = set(gene_list)
            gene_set_size[set_name] = len(gene_list)

        # store results in a dictionary of lists, one list per input gene list and one list for enriched terms
        enrichment_results = {key: [] for key in gene_set_dict.keys()}
        enrichment_results['GO term accession'] = []

        while True:
            most_significant_log_pvalue = 0.0
            most_significant_term = None

            # find most significant term (MST) in any set
            for set_name, gene_set in gene_set_dict.items():
                current_enrichment_df = self.get_all_enrichment(gene_set, gene_set_size[set_name])
                if len(current_enrichment_df) == 0:
                    continue
                current_enriched_term_row = current_enrichment_df.iloc[0]
                current_enrichment_log_pvalue = current_enriched_term_row['log_adjusted_pvalue']

                # replace negative infinity with defined minimum value
                current_enrichment_log_pvalue = check_significance_value(current_enrichment_log_pvalue)

                if current_enrichment_log_pvalue < most_significant_log_pvalue:
                    most_significant_log_pvalue = current_enrichment_log_pvalue
                    most_significant_term = current_enriched_term_row.name

            if most_significant_log_pvalue >= adj_log_p_value_threshold:
                break

            enrichment_results['GO term accession'].append(most_significant_term)
            current_term_genes = self.gene_term_df.loc[self.gene_term_df[self.go_column_label] == most_significant_term,
                                                       self.gene_column_label]

            # evaluate MST in all sets
            for set_name, gene_set in gene_set_dict.items():
                current_term_significance = self.get_single_enrichment(gene_set, most_significant_term, self.universe,
                                                                       gene_set_size[set_name])
                current_term_significance = check_significance_value(current_term_significance)
                enrichment_results[set_name].append(current_term_significance)

            # remove MST genes from all sets
            for set_name, gene_set in gene_set_dict.items():
                gene_set_dict[set_name].difference_update(current_term_genes.values)

        enrichment_results_df = pandas.DataFrame.from_dict(enrichment_results).set_index('GO term accession')
        return enrichment_results_df

    def plot_enrichment_multiple(self, multi_gene_list, output_filename):
        enrichment_results_df = self.iterative_enrichment_multilist(multi_gene_list)

        # replace accession numbers with GO term descriptions
        formatted_enrichment_df = enrichment_results_df.merge(self.term_definition, left_index=True,
                                                              right_on='GO term accession')
        formatted_enrichment_df['GO term'] = (formatted_enrichment_df['GO term name'] + '\n[' +
                                              formatted_enrichment_df['GO term accession'] + ']')
        formatted_enrichment_df.drop(self.term_definition.columns, axis=1, inplace=True)
        formatted_enrichment_df.set_index('GO term', inplace=True)

        sns.set(style='whitegrid')
        plt.figure()
        plot = sns.clustermap(formatted_enrichment_df, cmap="Reds_r")
        plt.setp(plot.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
        # plot.set_xlabel(r'$-log_{10}(p)$')
        # plot.set_ylabel('')
        # plot.xaxis.grid(False)
        # plot.yaxis.grid(True)
        # sns.despine(left=True, bottom=True)

        if output_filename:
            plot.savefig(output_filename)
        else:
            plt.show()


def run_enrichment(namespace):
    species = namespace.species

    # get a translator ready even if it's not needed
    translator = atg.data.identifiers.GeneIDTranslator(species)
    go_term_path = os.path.join(atg.data.genome_path[species], 'go_biological_process.csv')
    go_definition_path = os.path.join(atg.data.genome_path[species], 'go_definition.csv')

    maximum_size = DEFAULT_MAXIMUM_TERM_SIZE
    if 'term_size_limit' in namespace:
        maximum_size = namespace.term_size_limit
    calculator = EnrichmentCalculator(go_term_path, go_definition_path, maximum_size)

    # for a single filename, read the gene list, assuming no header
    if len(namespace.filename) == 1:
        input_gene_series = pandas.read_csv(namespace.filename[0]).ix[:, 0]
        gene_list = translator.translate_identifiers(input_gene_series, input_type=None, output_type='ensembl')

        if namespace.plot:
            calculator.plot_enrichment_single(gene_list, output_filename=namespace.plot, iterative=not namespace.full)
        else:
            if namespace.full:
                enrichment_result = calculator.get_all_enrichment(gene_list)
            else:
                enrichment_result = calculator.iterative_enrichment(gene_list)

            if namespace.output:
                enrichment_result.reset_index().to_csv(namespace.output, index=False)
            else:
                enrichment_result.reset_index().ix[0:TERMINAL_OUTPUT_LINES, ].to_string(sys.stdout, index=False, )

    # read each file, using the filename as the gene set name
    else:
        gene_set_dict = {}
        for filename in namespace.filename:
            gene_set_name = os.path.splitext(os.path.basename(filename))[0]
            input_gene_series = pandas.read_csv(filename).ix[:, 0]
            gene_list = translator.translate_identifiers(input_gene_series, input_type=None, output_type='ensembl')
            gene_set_dict[gene_set_name] = gene_list

        if namespace.plot:
            calculator.plot_enrichment_multiple(gene_set_dict, namespace.plot)

        else:
            enrichment_result = calculator.iterative_enrichment_multilist(gene_set_dict)
            if namespace.output:
                enrichment_result.reset_index().to_csv(namespace.output, index=False)
            else:
                enrichment_result.reset_index().ix[0:TERMINAL_OUTPUT_LINES, ].to_string(sys.stdout, index=False)


def setup_subparsers(subparsers):
    enrichment_parser = subparsers.add_parser('enrich', help='Gene set enrichment')
    enrichment_parser.add_argument('filename', help="One or more single-column files containing gene identifiers",
                                   nargs='+')
    enrichment_parser.add_argument('-s', '--species', default="human", help="")
    enrichment_parser.add_argument('-n', dest='term_size_limit', default=DEFAULT_MAXIMUM_TERM_SIZE, type=int)
    enrichment_parser.add_argument('-p', '--plot', help="filename for graphical output")
    enrichment_parser.add_argument('-o', '--output', help="filename for tabular output")
    enrichment_parser.add_argument('-f', '--full', action="store_true", help="Report full enrichment analysis "
                                                                             "for single gene list input")

    enrichment_parser.set_defaults(func=run_enrichment)
