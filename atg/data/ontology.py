"""

http://purl.obolibrary.org/obo/go/go-basic.obo
"""

import pronto
import pandas


def process_ontology(ontology_file, gene_term_file):
    ont = pronto.Ontology(ontology_file)
    gene_term_df = pandas.read_csv(gene_term_file, usecols=[0, 1]).dropna().drop_duplicates()

    terms = gene_term_df.ix[:, 1].unique()
    df_list = []

    for term_id in terms:
        term = ont.get(term_id, False)

        if not term:
            continue

        if term.other.get('namespace')[0] == 'biological_process':
            complete_term_list = [term_id] + term.rchildren().id
            complete_gene_df = gene_term_df.ix[gene_term_df.ix[:, 1].isin(complete_term_list)].copy()
            complete_gene_df.iloc[:, 1] = term_id
            df_list.append(complete_gene_df.drop_duplicates())

    complete_gene_term_df = pandas.concat(df_list, ignore_index=True)
    complete_gene_term_df.to_csv('/tmp/go_biological_process.csv', index=False)


if __name__ == '__main__':
    process_ontology("/Users/mchang/Downloads/go-basic.obo", "/Users/mchang/ATGData/human/GRCh38/gene_go.csv")
    process_ontology("/Users/mchang/Downloads/go-basic.obo", "/Users/mchang/ATGData/mouse/GRCm38/gene_go.csv")
