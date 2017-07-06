import pandas
import numpy
from scipy.stats import hypergeom

import seaborn as sns
import matplotlib.pyplot as plt

LOG_P_VALUE_MINIMUM = -40.0
LOG10_FACTOR = 1.0 / numpy.log(10)

GENE_COLUMN_LABEL_INDEX = 0
GO_COLUMN_LABEL_INDEX = 1
GO_DEFINITION_ACCESSION_INDEX = 1


def enrichment_significance(term_row):
    return hypergeom.logsf(term_row['hit_count'], term_row['universe'], term_row['term_count'], term_row['list_size'])


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


class EnrichmentCalculator:
    """
    calculate stats on statistical enrichment on gene categories, e.g. Gene Ontology
    """

    def __init__(self, gene_term_file, term_definition_file, maximum_size=1000, minimum_size=3):
        # TODO: filter by evidence code
        gene_term_df = pandas.read_csv(gene_term_file).dropna()

        self.gene_column_label = gene_term_df.columns[GENE_COLUMN_LABEL_INDEX]
        self.go_column_label = gene_term_df.columns[GO_COLUMN_LABEL_INDEX]

        gene_term_df.drop_duplicates([self.gene_column_label, self.go_column_label], inplace=True)

        term_count = gene_term_df[self.go_column_label].value_counts()
        relevant_terms = term_count[(term_count <= maximum_size) & (term_count >= minimum_size)]

        self.gene_term_df = gene_term_df[gene_term_df[self.go_column_label].isin(relevant_terms.index.values)]
        self.term_count = self.gene_term_df[self.go_column_label].value_counts()
        self.universe = self.gene_term_df[self.gene_column_label].nunique()

        self.term_description = pandas.read_csv(term_definition_file).dropna()

    def get_single_enrichment(self, gene_list, term, gene_universe=15000):
        term_df = self.gene_term_df[self.gene_term_df[self.go_column_label] == term]
        n = len(term_df)
        x = sum(term_df[self.gene_column_label].isin(gene_list))  # matched genes
        N = len(gene_list)

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
        enrichment_df['log_pvalue'] = enrichment_df.apply(enrichment_significance, axis=1, reduce=True) * LOG10_FACTOR

        # replace negative infinite p-values with a small constant
        enrichment_df['log_pvalue'].replace(-numpy.inf, LOG_P_VALUE_MINIMUM, inplace=True)
        enrichment_df['log_adjusted_pvalue'] = numpy.log10(p_adjust_bh(numpy.power(10, enrichment_df['log_pvalue']),
                                                                       n=len(self.term_description)))
        return enrichment_df.sort_values(['log_adjusted_pvalue', 'hit_count'], ascending=[True, False])

    def iterative_enrichment(self, gene_list, adj_log_p_value_threshold=-2.0, term_gene_limit=500):
        """

        :param gene_list:
        :param adj_log_p_value_threshold: stop performing enrichment calculations when the best adjusted p-value
        becomes larger than this value
        :param term_gene_limit: ignore any gene sets larger than the given size
        :return: a DataFrame containing raw and adjusted log p-values for GO BP terms.
        """
        current_gene_set = set(gene_list)
        full_gene_set_size = len(current_gene_set)
        enriched_term_results = []

        while len(current_gene_set) > 0:
            current_enrichment_df = self.get_all_enrichment(current_gene_set, full_gene_set_size)
            most_enriched_term_row = current_enrichment_df.query('term_count <= @term_gene_limit').iloc[0]
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

        annotated_enrichment_df = (enrichment_df.merge(self.term_description, left_index=True,
                                                       right_on=self.term_description.columns
                                                       [GO_DEFINITION_ACCESSION_INDEX])
                                                .query('log_adjusted_pvalue < @log_pvalue_threshold')
                                                .sort_values('log_adjusted_pvalue')
                                                .iloc[0:num_terms]
                                                .assign(display_log_pvalue=lambda x: x.log_adjusted_pvalue*-1))

        sns.set(style='whitegrid')
        plt.figure()
        plot = sns.stripplot(x='display_log_pvalue', y='GO term name', data=annotated_enrichment_df,
                             orient='h', edgecolor='gray', palette="Reds_r")
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
