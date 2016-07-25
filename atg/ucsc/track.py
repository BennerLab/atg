import os
import sys
import yaml
import multiprocessing
from atg.quantification.coverage import CoverageCalculator

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
                "viewLimits -20:20\n"
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


def bam_to_bigwig(bam_path, bigwig_forward_output, bigwig_reverse_output, genome):
    coverage = CoverageCalculator(bam_path)
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
            trackDb_output.write(HUB_TOPLEVEL.format(hub=self.hub))

            for multitrack in self.hub_config['multitracks']:
                for multitrack_name, track_list in multitrack.items():
                    trackDb_output.write(HUB_CONTAINER.format(parent=self.hub, name=multitrack_name))

                    for i, track_dict in enumerate(track_list):
                        bam_path = track_dict['path']
                        bam_filename = os.path.split(bam_path)[1]
                        bam_basename = os.path.splitext(bam_filename)[0]

                        output_filename_list = []

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
                                            os.path.join(self.base_output_path, self.genome, output_filename_list[0]),
                                            os.path.join(self.base_output_path, self.genome, output_filename_list[1]),
                                            self.genome])

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
