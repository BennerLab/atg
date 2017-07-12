import os
import configparser
import atg.config

DATA_SOURCE_PATH = os.path.join(os.path.dirname(__file__), 'data_sources.ini')

data_root = os.path.expanduser(atg.config.settings['Data']['Root'])
config = configparser.ConfigParser(default_section="Common",
                                   interpolation=configparser.ExtendedInterpolation())
config.read(DATA_SOURCE_PATH)
genome_path = {}

for organism in config.sections():
    ensembl_genome = config.get(organism, 'ensembl_genome')
    current_genome_path = os.path.join(data_root, organism, ensembl_genome)
    genome_path[organism] = current_genome_path
