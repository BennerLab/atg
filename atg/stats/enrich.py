import pandas
from scipy.stats import hypergeom
from numpy import log

import seaborn as sns
import matplotlib.pyplot as plt

LOG10_FACTOR = 1.0 / log(10)


def enrichment_significance(term_row):
    return hypergeom.logsf(term_row['hit_count'], term_row['universe'], term_row['term_count'], term_row['list_size'])


class EnrichmentCalculator:
    """
    calculate stats on statistical enrichment on gene categories, e.g. Gene Ontology
    """

    def __init__(self, gene_term_file, term_definition_file, maximum_size=1000, minimum_size=3):
        gene_term_df = pandas.read_csv(gene_term_file).dropna()
        term_count = gene_term_df['GO Term Accession'].value_counts()
        relevant_terms = term_count[(term_count <= maximum_size) & (term_count >= minimum_size)]

        self.gene_term_df = gene_term_df[gene_term_df['GO Term Accession'].isin(relevant_terms.index.values)]
        self.term_count = self.gene_term_df['GO Term Accession'].value_counts()
        self.universe = self.gene_term_df['Ensembl Gene ID'].nunique()

        self.term_description = pandas.read_csv(term_definition_file).dropna()

    def get_single_enrichment(self, gene_list, term, gene_universe=15000):
        term_df = self.gene_term_df[self.gene_term_df['GO Term Accession'] == term]
        n = len(term_df)
        x = sum(term_df['Ensembl Gene ID'].isin(gene_list))  # matched genes
        N = len(gene_list)

        return hypergeom.logsf(x, gene_universe, n, N)

    def get_all_enrichment(self, gene_list):
        enrichment_df = pandas.DataFrame()
        enrichment_df['term_count'] = self.term_count
        enrichment_df['hit_count'] = (self.gene_term_df[self.gene_term_df['Ensembl Gene ID'].isin(gene_list)]
                                                       ['GO Term Accession'].value_counts().astype(int))
        enrichment_df.dropna(inplace=True)

        enrichment_df['universe'] = self.universe
        enrichment_df['list_size'] = len(gene_list)
        enrichment_df['log_pvalue'] = enrichment_df.apply(enrichment_significance, axis=1, reduce=True) * LOG10_FACTOR

        return enrichment_df

    def iterative_enrichment(self, gene_list):
        pass

    def plot_enrichment_single(self, gene_list, output_filename=None, log_pvalue_threshold=-2.0, num_terms=20,
                               iterative=False):

        annotated_enrichment_df = (self.get_all_enrichment(gene_list)
                                       .merge(self.term_description, left_index=True, right_on="GO Term Accession")
                                       .query('log_pvalue < @log_pvalue_threshold')
                                       .sort_values('log_pvalue')
                                       .iloc[0:num_terms]
                                       .assign(display_log_pvalue=lambda x: x.log_pvalue*-1))

        sns.set(style='whitegrid')
        plt.figure()
        plot = sns.stripplot(x='display_log_pvalue', y='GO Term Name', data=annotated_enrichment_df,
                             orient='h', edgecolor='gray', palette="Reds_r")
        plot.set_xlabel(r'$-log_{10}(p)$')
        plot.set_ylabel('')
        plot.xaxis.grid(False)
        plot.yaxis.grid(True)
        sns.despine(left=True, bottom=True)
        plt.tight_layout()
        plt.savefig("/tmp/test.pdf")
