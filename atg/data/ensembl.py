"""
Find species data in Ensembl, recording genome and annotation URLs.
"""

import pandas
import ftplib
import string

ENSEMBL_GENOMES_DNA_BASE_LOCATION = string.Template('pub/current/$division/fasta/$species/dna/')
ENSEMBL_GENOMES_GTF_BASE_LOCATION = string.Template('pub/current/$division/gtf/$species/$assembly.$version.gtf.gz')

# ftp://ftp.ensemblgenomes.org//pub/release-35/plants/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz
# ftp://ftp.ensemblgenomes.org//pub/release-35/plants/gtf/zea_mays/Zea_mays.AGPv4.35.gtf.gz

def get_species_information(species):
    """

    :param species: genus and species (as named by Ensembl), e.g. zea_mays
    :return: dictionary containing URLs to genome fasta and gene annotation (GTF), if found
    """

    ensembl_genomes = pandas.read_table('/Users/mchang/Desktop/Ensembl species/species_Ensembl_Genomes.txt',
                                        index_col=False)
    ensembl_division = ensembl_genomes.ix[ensembl_genomes['species'] == species, 'division'].iat[0].\
        lstrip('Ensembl').lower()

    with ftplib.FTP('ftp.ensemblgenomes.org') as ftp:
        ftp.login()
        genome_listing = ftp.nlst(ENSEMBL_GENOMES_DNA_BASE_LOCATION.safe_substitute(division=ensembl_division,
                                                                                    species=species))
        genome_location = ''
        annotation_location = ''
        genome_assembly_version = ''

        # find toplevel unmasked genome
        for filename in genome_listing:
            if 'dna.toplevel' in filename:
                genome_location = filename
                break

        if genome_location != '':
            genome_filename = genome_location.split('/')[-1]
            genome_assembly = genome_filename.rstrip('.dna.toplevel.fa.gz')
            genome_assembly_version = genome_assembly.split('.', maxsplit=1)[1]

            annotation_listing = ftp.nlst(ENSEMBL_GENOMES_GTF_BASE_LOCATION.safe_substitute(division=ensembl_division,
                                                                                            species=species,
                                                                                            assembly=genome_assembly,
                                                                                            version=35))

            if len(annotation_listing) == 0:
                annotation_location = ''
            elif len(annotation_listing) == 1:
                annotation_location = annotation_listing[0]
            else:
                annotation_location = 'multiple'

        ftp.close()

    return {'genome': genome_location, 'annotation': annotation_location, 'version': genome_assembly_version}


if __name__ == '__main__':
    print(get_species_information('zea_mays'))
    print(get_species_information('schizosaccharomyces_pombe'))
    print(get_species_information('apis_mellifera'))
    print(get_species_information('daphnia_pulex'))