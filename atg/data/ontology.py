"""

Process Gene Ontology
"""

import pronto
import pandas
import progress.bar


def process_ontology(gene_term_file, ontology_file='http://purl.obolibrary.org/obo/go/go-basic.obo'):
    ont = pronto.Ontology(ontology_file)
    gene_term_df = pandas.read_csv(gene_term_file, usecols=[0, 1]).dropna().drop_duplicates()

    terms = gene_term_df.ix[:, 1].unique()
    df_list = []

    progress_bar = progress.bar.Bar('Processing GO terms',
                                    suffix='%(index)d/%(max)d', max=len(terms))
    for term_id in terms:
        progress_bar.next()
        term = ont.get(term_id, False)

        if not term:
            continue

        if term.other.get('namespace')[0] == 'biological_process':
            complete_term_list = [term_id] + term.rchildren().id
            complete_gene_df = gene_term_df.ix[gene_term_df.ix[:, 1].isin(complete_term_list)].copy()
            complete_gene_df.iloc[:, 1] = term_id
            df_list.append(complete_gene_df.drop_duplicates())

    progress_bar.finish()

    complete_gene_term_df = pandas.concat(df_list, ignore_index=True)
    return complete_gene_term_df
