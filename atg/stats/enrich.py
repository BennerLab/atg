import pandas
from scipy.stats import hypergeom


class EnrichmentCalculator:
    """
    calculate stats on statistical enrichment on gene categories, e.g. Gene Ontology
    """

    def __init__(self, gene_term_file):
        self.gene_term_df = pandas.read_csv(gene_term_file).dropna()
        self.term_counts = self.gene_term_df['GO Term Accession'].value_counts()

    def get_single_enrichment(self, gene_list, term, gene_universe=15000):
        term_df = self.gene_term_df[self.gene_term_df['GO Term Accession'] == term]
        n = len(term_df)
        x = sum(term_df['Ensembl Gene ID'].isin(gene_list))  # matched genes
        N = len(gene_list)

        return hypergeom.logsf(x, gene_universe, n, N)

    def get_all_enrichment(self, gene_list):
        pass
