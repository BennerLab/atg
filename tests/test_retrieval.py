import pytest
import os
import tempfile
import pandas
import atg.data.retrieve
import atg.data.ensembl

UCSC_GZIPPED_FILE = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.agp.gz'
UCSC_CHROMOSOME_SIZE_FILE = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'
ASSEMBLY_FILE_ENTRIES_PER_LINE = 9
UCSC_CHROMOSOME_ENTRIES = 455
BIOMART_XML_FILE = os.path.join(os.path.dirname(__file__), 'data', 'GRCh38_chrY.xml')


def test_biomart_xml_fetch():
    with open(BIOMART_XML_FILE, 'r') as biomart_input:
        with tempfile.TemporaryDirectory() as working_directory:
            biomart_output = os.path.join(working_directory, "gene.csv")
            xml_string = ''.join(biomart_input.readlines())
            atg.data.retrieve.fetch_ensembl(xml_string, biomart_output)

            gene_df = pandas.read_csv(biomart_output)
            # at least 10 genes should be present
            assert gene_df.shape[0] > 10
            # the XML file restricts retrievals to chromosome Y
            assert gene_df['Chromosome/scaffold name'][0] == 'Y'


def test_ucsc_url_fetch():
    with tempfile.TemporaryDirectory() as working_directory:
        assembly_filename = os.path.join(working_directory, "test.agp")
        atg.data.retrieve.fetch_url(UCSC_GZIPPED_FILE, assembly_filename)

        with open(assembly_filename, 'r') as assembly:
            first_line = assembly.readline().split()
            assert first_line[0] == 'chr1'
            assert first_line[0] != 'chr2'
            assert len(first_line) == ASSEMBLY_FILE_ENTRIES_PER_LINE

        chromosome_filename = os.path.join(working_directory, "chrom.sizes")
        atg.data.retrieve.fetch_url(UCSC_CHROMOSOME_SIZE_FILE, chromosome_filename)
        chromosome_df = pandas.read_csv(chromosome_filename, names=['chrom', 'size'], sep='\t')
        assert chromosome_df['chrom'][0] == 'chr1'
        assert chromosome_df.shape[0] == UCSC_CHROMOSOME_ENTRIES


@pytest.fixture
def atg_data_tracker():
    with tempfile.TemporaryDirectory() as temp_dir:
        yield atg.data.retrieve.ATGDataTracker(temp_dir)


def test_atg_yeast_retrieval(atg_data_tracker):
    atg_data_tracker.retrieve_data('yeast')
    assert 'yeast' in os.listdir(atg_data_tracker.data_root)

    # check that files all exist and have non-zero size
    for genome_file in atg.data.retrieve.GENOME_FILES:
        current_genome_file = os.path.join(atg_data_tracker.data_root, "yeast", "R64-1-1", genome_file)
        assert os.path.exists(current_genome_file)
        assert os.path.getsize(current_genome_file) > 0

    # check that the GO biological process file is created
    go_bp_path = os.path.join(atg_data_tracker.data_root, "yeast", "R64-1-1", "go_biological_process.csv")
    assert not os.path.exists(go_bp_path)

    atg_data_tracker.derive_data('yeast')
    assert os.path.exists(go_bp_path)
    assert os.path.getsize(go_bp_path) > 0


@pytest.fixture
def ensembl_genomes():
    yield atg.data.ensembl.EnsemblSpecies()


organism_information = {
                        'zea_mays': {
                            'annotation': 'pub/release-35/plants/gtf/zea_mays/Zea_mays.AGPv4.35.gtf.gz',
                            'version': 'AGPv4',
                            'genome': 'pub/release-35/plants/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz',
                            'species': 'zea_mays'},
                        'apis_mellifera': {
                            'annotation': 'pub/release-35/metazoa/gtf/apis_mellifera/Apis_mellifera.'
                                          'GCA_000002195.1.35.gtf.gz',
                            'version': 'GCA_000002195.1',
                            'genome': 'pub/release-35/metazoa/fasta/apis_mellifera/dna/Apis_mellifera.'
                                       'GCA_000002195.1.dna.toplevel.fa.gz',
                            'species': 'apis_mellifera'},
                        'agaricus_bisporus_var_bisporus_h97': {
                            'genome': 'pub/release-35/fungi/fasta/fungi_basidiomycota1_collection/'
                                      'agaricus_bisporus_var_bisporus_h97/dna/'
                                      'Agaricus_bisporus_var_bisporus_h97.'
                                      'Agabi_varbisH97_2.dna.toplevel.fa.gz',
                            'version': 'Agabi_varbisH97_2',
                            'annotation': 'pub/release-35/fungi/gtf/fungi_basidiomycota1_collection/'
                                          'agaricus_bisporus_var_bisporus_h97/'
                                          'Agaricus_bisporus_var_bisporus_h97.'
                                          'Agabi_varbisH97_2.35.gtf.gz',
                            'species': 'agaricus_bisporus_var_bisporus_h97'}}


@pytest.mark.parametrize("organism,division", [('zea_mays', 'EnsemblPlants')])
def test_ensembl_genomes_table_read(ensembl_genomes, organism, division):
    one_record = ensembl_genomes.ensembl_species_df.loc[ensembl_genomes.ensembl_species_df['species'] == organism]
    # check that a single record is found
    assert one_record.shape[0] == 1
    assert one_record.iloc[0]['division'] == division


@pytest.mark.parametrize("organism,organism_values", organism_information.items())
def test_ensembl_genomes_species_information(ensembl_genomes, organism, organism_values):
    assert ensembl_genomes.get_species_information(organism) == organism_values


def test_ensembl_genomes_collect_species(ensembl_genomes):
    species_df = pandas.DataFrame.from_records([organism_information['apis_mellifera'],
                                                organism_information['zea_mays'],
                                                organism_information['agaricus_bisporus_var_bisporus_h97']])
    species_list = ['apis_mellifera', 'zea_mays', 'agaricus_bisporus_var_bisporus_h97']
    assert species_df.equals(ensembl_genomes.collect_species_information(species_list))


def test_ensembl_genomes_retrieval(ensembl_genomes):
    assert not ensembl_genomes.retrieve_species_data('homo_sapiens')
    assert ensembl_genomes.retrieve_species_data('agaricus_bisporus_var_bisporus_h97')
