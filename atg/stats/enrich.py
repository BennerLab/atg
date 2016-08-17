import pandas
from scipy.stats import hypergeom


def enrichment_significance(term_row):
    return hypergeom.logsf(term_row['Hit count'], term_row['Universe'], term_row['Term count'], term_row['List size'])


class EnrichmentCalculator:
    """
    calculate stats on statistical enrichment on gene categories, e.g. Gene Ontology
    """

    def __init__(self, gene_term_file, maximum_size=1000, minimum_size=3):
        # self.gene_term_df = (pandas.read_csv(gene_term_file)
        #                      .dropna()
        #                      .groupby('GO Term Accession')
        #                      .filter(lambda x: maximum_size >= len(x) >= minimum_size))

        gene_term_df = pandas.read_csv(gene_term_file).dropna()
        term_count = gene_term_df['GO Term Accession'].value_counts()
        relevant_terms = term_count[(term_count <= maximum_size) & (term_count >= minimum_size)]

        self.gene_term_df = gene_term_df[gene_term_df['GO Term Accession'].isin(relevant_terms.index.values)]
        self.term_count = self.gene_term_df['GO Term Accession'].value_counts()
        self.universe = self.gene_term_df['Ensembl Gene ID'].nunique()

    def get_single_enrichment(self, gene_list, term, gene_universe=15000):
        term_df = self.gene_term_df[self.gene_term_df['GO Term Accession'] == term]
        n = len(term_df)
        x = sum(term_df['Ensembl Gene ID'].isin(gene_list))  # matched genes
        N = len(gene_list)

        return hypergeom.logsf(x, gene_universe, n, N)

    def get_all_enrichment(self, gene_list):
        enrichment_df = pandas.DataFrame()
        enrichment_df['Term count'] = self.term_count
        hit_count = self.gene_term_df[self.gene_term_df['Ensembl Gene ID'].isin(gene_list)]['GO Term Accession'].\
            value_counts()
        enrichment_df['Hit count'] = hit_count
        enrichment_df.dropna(inplace=True)

        enrichment_df['Universe'] = self.universe
        enrichment_df['List size'] = len(gene_list)
        enrichment_df['P-value'] = enrichment_df.apply(enrichment_significance, axis=1, reduce=True)
        enrichment_df.to_csv('/tmp/test.csv')

        return enrichment_df
