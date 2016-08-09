import os
import sys
import yaml
import multiprocessing
import atg.quantification.coverage

DEFAULT_TRACK_Y_LIMIT = 20
GENOMES_TXT = 'genome {genome}\ntrackDb {genome}/trackDb.txt'
HUB_TXT = ("hub {hub}\n"
           "shortLabel {hub}\n"
           "longLabel {hub}\n"
           "genomesFile genomes.txt\n"
           "email mchang@ucsd.edu")

HUB_TOPLEVEL = ("track {hub}\n"
                "superTrack on show\n"
                "shortLabel {hub}\n"
                "longLabel {hub}\n"
                "viewLimits {ymin}:{ymax}\n"
                "viewLimitsMax -100:100\n\n")

HUB_CONTAINER = ("\ttrack {name}\n"
                 "\tparent {parent}\n"
                 "\tcontainer multiWig\n"
                 "\tnoInherit on\n"
                 "\tshortLabel {name}\n"
                 "\tlongLabel {name}\n"
                 "\ttype bigWig\n"
                 "\tconfigurable on\n"
                 "\tvisibility full\n"
                 "\taggregate transparentOverlay\n"
                 "\tshowSubtrackColorOnUi on\n"
                 "\twindowingFunction maximum\n"
                 "\tpriority 1.4\n"
                 "\talwaysZero on\n"
                 "\tyLineMark 0\n"
                 "\tyLineOnOff on\n"
                 "\tmaxHeightPixels 125:125:11\n\n")

TRACK_ENTRY_STRANDED = ("\t\ttrack {track}-{strand}\n"
                        "\t\tbigDataUrl {path}\n"
                        "\t\tshortLabel {track}-{strand}\n"
                        "\t\tlongLabel {track}-{strand}\n"
                        "\t\ttype bigWig\n"
                        "\t\tparent {parent}\n"
                        "\t\tcolor {rgb}\n\n")

TRACK_ENTRY_UNSTRANDED = ("\t\ttrack {track}\n"
                          "\t\tbigDataUrl {path}\n"
                          "\t\tshortLabel {track}\n"
                          "\t\tlongLabel {track}\n"
                          "\t\ttype bigWig\n"
                          "\t\tparent {parent}\n"
                          "\t\tcolor {rgb}\n\n")

# UCSC Track documentation:
# https://genome.ucsc.edu/goldenpath/help/trackDb/trackDbDoc.html#wig_-_Signal_Graphing_Track_Settings

# Description of UCSC genome browser color transparency:
# https://groups.google.com/a/soe.ucsc.edu/d/msg/genome/zEUWA-UeyV8/eil11j1Zd8kJ

# Pastel1 from colorbrewer
COLOR_CYCLE = ('251,180,174', '179,205,227', '204,235,197', '222,203,228', '254,217,166', '255,255,204',
               '229,216,189', '253,218,236', '242,242,242')


def bam_to_bigwig(bam_path, genome, bigwig_forward_output, bigwig_reverse_output=None):
    if bigwig_reverse_output == 'None':
        coverage = atg.quantification.coverage.UnstrandedCoverageCalculator(bam_path)
        coverage.write_bigwig(bigwig_forward_output, genome)

    else:
        coverage = atg.quantification.coverage.StrandedCoverageCalculator(bam_path)
        coverage.write_bigwig(bigwig_forward_output, bigwig_reverse_output, genome)


class HubBuilder:
    """

    1. generate directory structure for output
    2. write metadata
        - genomes.txt
        - hub.txt
        - {genome}/trackDB.txt
    3. convert BAM files to stranded bigWig files

    """

    def __init__(self, hub_config_filename):
        self.hub_config = yaml.load(open(hub_config_filename))
        self.data_path = os.path.split(hub_config_filename)[0]
        self.hub = self.hub_config['hub']
        self.genome = self.hub_config['genome']
        self.base_output_path = self.hub_config['path']
        self.library = self.hub_config['library']

    def make_output_structure(self, overwrite=False):
        output_path = os.path.join(self.base_output_path, self.genome)
        try:
            os.makedirs(self.base_output_path, exist_ok=overwrite)
            os.makedirs(output_path, exist_ok=overwrite)
        except FileExistsError:
            print("Could not create output directory structure. [{}]".format(self.base_output_path), file=sys.stderr)
            return False

        return True

    def write_files(self):
        with open(os.path.join(self.base_output_path, 'genomes.txt'), 'w') as genomes_output:
            genomes_output.write(GENOMES_TXT.format(genome=self.genome))
        with open(os.path.join(self.base_output_path, 'hub.txt'), 'w') as hub_output:
            hub_output.write(HUB_TXT.format(hub=self.hub))

        track_queue = []

        with open(os.path.join(self.base_output_path, self.genome, 'trackDb.txt'), 'w') as trackDb_output:
            # determine y-axis scaling based on library type
            if self.library == 'unstranded':
                trackDb_output.write(HUB_TOPLEVEL.format(hub=self.hub, ymin=0, ymax=DEFAULT_TRACK_Y_LIMIT))
            elif self.library == 'rf':
                trackDb_output.write(HUB_TOPLEVEL.format(hub=self.hub, ymin=-1 * DEFAULT_TRACK_Y_LIMIT,
                                                         ymax=DEFAULT_TRACK_Y_LIMIT))
            else:
                print('Library type not recognized, should be either "unstranded" or "rf"')
                sys.exit(1)

            for multitrack in self.hub_config['multitracks']:
                for multitrack_name, track_list in multitrack.items():
                    trackDb_output.write(HUB_CONTAINER.format(parent=self.hub, name=multitrack_name))

                    for i, track_dict in enumerate(track_list):
                        bam_path = track_dict['path']
                        bam_filename = os.path.split(bam_path)[1]
                        bam_basename = os.path.splitext(bam_filename)[0]

                        output_filename_list = []

                        if self.library == 'unstranded':
                            current_path = bam_basename + os.path.extsep + 'bigwig'
                            output_filename_list.append(current_path)
                            trackDb_output.write(TRACK_ENTRY_UNSTRANDED.format(track=track_dict['track'],
                                                                               path=current_path,
                                                                               parent=multitrack_name,
                                                                               rgb=COLOR_CYCLE[i % len(COLOR_CYCLE)]))
                            track_queue.append([os.path.join(self.data_path, bam_path),
                                                self.genome,
                                                os.path.join(self.base_output_path,
                                                             self.genome, output_filename_list[0]),
                                                None])

                        elif self.library == 'rf':
                            # each sample will have a positive and negative track entry
                            for strand in ('pos', 'neg'):
                                current_path = bam_basename + '_' + strand + os.path.extsep + 'bigwig'
                                output_filename_list.append(current_path)
                                trackDb_output.write(TRACK_ENTRY_STRANDED.format(track=track_dict['track'],
                                                                                 path=current_path,
                                                                                 parent=multitrack_name,
                                                                                 rgb=COLOR_CYCLE[i % len(COLOR_CYCLE)],
                                                                                 strand=strand))

                            # queue up coverage calculation parameters
                            track_queue.append([os.path.join(self.data_path, bam_path),
                                                self.genome,
                                                os.path.join(self.base_output_path,
                                                             self.genome, output_filename_list[0]),
                                                os.path.join(self.base_output_path,
                                                             self.genome, output_filename_list[1])])

        # generate actual data tracks in parallel
        with multiprocessing.Pool(processes=8) as pool:
            pool.starmap(bam_to_bigwig, track_queue)


class TrackOrganizer:
    """
    Modify an existing Genome Browser Hub (generated with Homer), arranging tracks as specified by YAML file.
    """

    def __init__(self, hub_config_filename):
        self.hub_config = yaml.load(open(hub_config_filename))
        self.hub = self.hub_config['hub']
        self.genome = self.hub_config['genome']

    def write_track_db(self, output_filename='trackDb.txt'):
        """
        trackDb.txt
        :return:
        """

        with open(output_filename, 'w') as trackDb_output:
            trackDb_output.write(HUB_TOPLEVEL.format(hub=self.hub))

            for multitrack in self.hub_config['multitracks']:
                for multitrack_name, track_list in multitrack.items():
                    trackDb_output.write(HUB_CONTAINER.format(parent=self.hub, name=multitrack_name))

                    for i, track_dict in enumerate(track_list):
                        bigwig_path = track_dict['path']

                        trackDb_output.write(TRACK_ENTRY_UNSTRANDED.format(track=track_dict['track'], path=bigwig_path,
                                                                           parent=multitrack_name,
                                                                           rgb=COLOR_CYCLE[i % len(COLOR_CYCLE)]))

# Functions for main atg script


def make_hub(namespace):
    hub = HubBuilder(namespace.config)
    hub.make_output_structure(overwrite=True)
    hub.write_files()


def setup_hub(namespace):
    print('hub:', namespace.name)
    print('genome:', namespace.genome)
    print('path:', namespace.path)
    print('library:', namespace.library)
    print('multitracks:')
    print('  - default:')

    for filename in namespace.file_list:
        trackname = os.path.splitext(os.path.basename(filename))[0]
        print('    - track:', trackname)
        print('      path:', filename)


def reorganize_hub(namespace):
    hub = TrackOrganizer(namespace.config)
    hub.write_track_db(namespace.output)


def setup_subparsers(subparsers):
    make_hub_parser = subparsers.add_parser('make', help='Generate a Hub from YAML file, including conversion of BAM'
                                                         'to bigWig.')
    make_hub_parser.add_argument('config', help="YAML configuration file")
    make_hub_parser.set_defaults(func=make_hub)

    # Generate a basic YAML file for describing data and output
    setup_hub_parser = subparsers.add_parser('setup', help="Generate a simple YAML configuration file for "
                                                           "Hub grouping. Outputs to stdout.")
    setup_hub_parser.add_argument('file_list', nargs='+', help="BAM files")
    setup_hub_parser.add_argument('-n', '--name', help="Hub name", default="DEFAULT")
    setup_hub_parser.add_argument('-g', '--genome', help="UCSC genome, e.g. hg38", default="GENOME")
    setup_hub_parser.add_argument('-p', '--path', help="Path for hub output", default="PATH")
    setup_hub_parser.add_argument('-l', '--library', help="Strand specificity of sequencing library",
                                  choices=['unstranded', 'rf'], default='rf')
    setup_hub_parser.set_defaults(func=setup_hub)

    # New organization of existing bigwig files
    reorganize_hub_parser = subparsers.add_parser('reorganize', help="Generate a new trackDb.txt file given a YAML "
                                                                     "configuration. Will not create tracks.")
    reorganize_hub_parser.add_argument('config', help="YAML configuration file")
    reorganize_hub_parser.add_argument('-o', '--output', help="output file", default="trackDb.txt")
    reorganize_hub_parser.set_defaults(func=reorganize_hub)
