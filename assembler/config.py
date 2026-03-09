import configparser
import os
import sys

class _Config:
    """A class to hold all configuration settings, loaded from config.ini."""
    def __init__(self):
        # This will be populated by load_config
        self.raw_config = None

    def _load_from_parser(self, config_parser):
        self.raw_config = config_parser

        # [debug]
        debug_section = self.raw_config['debug']
        self.DEBUG = debug_section.getboolean('DEBUG')
        self.PRINT_RUNTIME_LOGS = debug_section.getboolean('PRINT_RUNTIME_LOGS')
        self.OUTPUT_LOGGING_FILES = debug_section.getboolean('OUTPUT_LOGGING_FILES')

        # [anchor_params]
        anchor_params_section = self.raw_config['anchor_params']
        self.MIN_ANCHOR_LENGTH = anchor_params_section.getint('MIN_ANCHOR_LENGTH')
        self.MIN_NODES_IN_ANCHOR = anchor_params_section.getint('MIN_NODES_IN_ANCHOR')

        # [anchor_dict]
        anchor_dict_section = self.raw_config['anchor_dict']
        self.FORWARD_DICTIONARY = anchor_dict_section.getint('FORWARD_DICTIONARY')
        self.REVERSE_DICTIONARY = anchor_dict_section.getint('REVERSE_DICTIONARY')
        self.END_NODE_POS = anchor_dict_section.getint('END_NODE_POS')
        self.SNARL_ID_POS = anchor_dict_section.getint('SNARL_ID_POS')
        
        # [gaf_parsing]
        gaf_parsing_section = self.raw_config['gaf_parsing']
        self.EXPECTED_GAF_TAGS = gaf_parsing_section.getint('EXPECTED_GAF_TAGS')
        self.EXPECTED_MAP_Q = gaf_parsing_section.getint('EXPECTED_MAP_Q')
        self.MIN_CS_LEN = gaf_parsing_section.getint('MIN_CS_LEN')
        self.READ_NAME_ID = gaf_parsing_section.getint('READ_NAME_ID')
        self.READ_LEN = gaf_parsing_section.getint('READ_LEN')
        self.READ_START_ID = gaf_parsing_section.getint('READ_START_ID')
        self.RELATIVE_STRAND_ID = gaf_parsing_section.getint('RELATIVE_STRAND_ID')
        self.PATH_ID = gaf_parsing_section.getint('PATH_ID')
        self.PATH_START_ID = gaf_parsing_section.getint('PATH_START_ID')
        self.PATH_END_ID = gaf_parsing_section.getint('PATH_END_ID')
        self.DIV_ID = gaf_parsing_section.getint('DIV_ID')
        self.MAP_Q_ID = gaf_parsing_section.getint('MAP_Q_ID')
        self.CS_TAG_ID = gaf_parsing_section.getint('CS_TAG_ID')

        # [alignment_list]
        alignment_list_section = self.raw_config['alignment_list']
        self.READ_POSITION = alignment_list_section.getint('READ_POSITION')
        self.R_LEN_POSITION = alignment_list_section.getint('R_LEN_POSITION')
        self.READ_START_POS = alignment_list_section.getint('READ_START_POS')
        self.STRAND_POSITION = alignment_list_section.getint('STRAND_POSITION')
        self.START_POSITION = alignment_list_section.getint('START_POSITION')
        self.END_POSITION = alignment_list_section.getint('END_POSITION')
        self.NODE_POSITION = alignment_list_section.getint('NODE_POSITION')
        self.ORIENTATION_POSITION = alignment_list_section.getint('ORIENTATION_POSITION')
        self.CIGAR_POSITION = alignment_list_section.getint('CIGAR_POSITION')

        # [bp_matched_read_list]
        bp_matched_read_list_section = self.raw_config['bp_matched_read_list']
        self.READ_ID = bp_matched_read_list_section.getint('READ_ID')
        self.READ_STRAND = bp_matched_read_list_section.getint('READ_STRAND')
        self.ANCHOR_START = bp_matched_read_list_section.getint('ANCHOR_START')
        self.ANCHOR_END = bp_matched_read_list_section.getint('ANCHOR_END')
        self.MATCH_LIMIT = bp_matched_read_list_section.getint('MATCH_LIMIT')
        self.CS_LEFT_AVAIL = bp_matched_read_list_section.getint('CS_LEFT_AVAIL')
        self.CS_RIGHT_AVAIL = bp_matched_read_list_section.getint('CS_RIGHT_AVAIL')
        
        # [extension_merging]
        extension_merging_section = self.raw_config['extension_merging']
        self.MIN_ANCHOR_READS = extension_merging_section.getint('MIN_ANCHOR_READS')
        self.HET_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING = extension_merging_section.getfloat('HET_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING')
        self.HOMO_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING = extension_merging_section.getfloat('HOMO_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING')
        self.MIN_READS_REQUIRED_FOR_MERGING_R0 = extension_merging_section.getint('MIN_READS_REQUIRED_FOR_MERGING_R0')
        self.MIN_READS_REQUIRED_FOR_MERGING_R1 = extension_merging_section.getint('MIN_READS_REQUIRED_FOR_MERGING_R1')
        self.MIN_ANCHOR_READCOV = extension_merging_section.getint('MIN_ANCHOR_READCOV')
        self.MIN_ANCHOR_READCOV_FOR_INDEPENDENT_ANCHOR_EXTENSION = extension_merging_section.getint('MIN_ANCHOR_READCOV_FOR_INDEPENDENT_ANCHOR_EXTENSION')
        self.FRACTION_READS_FOR_SNARL_BOUNDARY_EXTENTION = extension_merging_section.getfloat('FRACTION_READS_FOR_SNARL_BOUNDARY_EXTENTION')
        self.MIN_READS_REQUIRED_FOR_BOUNDARY_EXTENSION = extension_merging_section.getint('MIN_READS_REQUIRED_FOR_BOUNDARY_EXTENSION')
        self.DROP_FRACTION = extension_merging_section.getfloat('DROP_FRACTION')

        # [reliability]
        reliability_section = self.raw_config['reliability']
        self.MIN_SNARL_LINKAGE_THRESHOLD = reliability_section.getint('MIN_SNARL_LINKAGE_THRESHOLD')
        self.ADD_BACK_HOMO_SNARLS = reliability_section.getboolean('ADD_BACK_HOMO_SNARLS')
        self.RELIABLE_SNARL_FRACTION_THRESHOLD = reliability_section.getfloat('RELIABLE_SNARL_FRACTION_THRESHOLD')
        self.ERROR_TOLERANCE_IN_COMPATIBILITY_CHECK = reliability_section.getint('ERROR_TOLERANCE_IN_COMPATIBILITY_CHECK')
        self.ENABLE_UNEQUAL_SET_COMPATIBILITY = reliability_section.getboolean('ENABLE_UNEQUAL_SET_COMPATIBILITY')
        self.MIN_READS_FOR_PARTITION_COMPATIBILITY = reliability_section.getint('MIN_READS_FOR_PARTITION_COMPATIBILITY')
        self.ENABLE_PROBABILISTIC_RELIABILITY_CHECKING = reliability_section.getboolean('ENABLE_PROBABILISTIC_RELIABILITY_CHECKING')
        self.ENABLE_REFINED_PROBABILISTIC_RELIABILITY_CHECKING = reliability_section.getboolean('ENABLE_REFINED_PROBABILISTIC_RELIABILITY_CHECKING')
        self.INVERSE_THRESHOLD = reliability_section.getint('INVERSE_THRESHOLD')
        self.ENABLE_BINOMIAL_RELIABILITY_CHECKING = reliability_section.getboolean('ENABLE_BINOMIAL_RELIABILITY_CHECKING')
        self.BINOMIAL_PVALUE_THRESHOLD = reliability_section.getfloat('BINOMIAL_PVALUE_THRESHOLD')
        self.MAX_POTENTIALLY_LINKED_SNARLS_TO_KEEP = reliability_section.getint('MAX_POTENTIALLY_LINKED_SNARLS_TO_KEEP')
        self.MAX_NEIGHBOURING_SNARLS_TO_PEEK_IN_READ = reliability_section.getint('MAX_NEIGHBOURING_SNARLS_TO_PEEK_IN_READ')

        # [gtest]
        gtest_section = self.raw_config['gtest']
        self.USE_GTEST_FOR_PARTITION_COMPATIBILITY = gtest_section.getboolean('USE_GTEST_FOR_PARTITION_COMPATIBILITY')
        self.DETANGLE_MAX_LOG_P = gtest_section.getint('DETANGLE_MAX_LOG_P')
        self.DETANGLE_MIN_LOG_P_DELTA = gtest_section.getint('DETANGLE_MIN_LOG_P_DELTA')
        self.DETANGLE_GTEST_EPSILON = gtest_section.getfloat('DETANGLE_GTEST_EPSILON')

settings = _Config()

def load_config(config_file=None):
    """
    Loads configuration from a given file or discovers it.
    This function populates the global 'settings' variable.
    """
    config = configparser.ConfigParser(inline_comment_prefixes="#")
    path_to_load = None
    
    if config_file:
        # If a specific file is provided, use it.
        if not os.path.exists(config_file):
            raise FileNotFoundError(f"Specified config file not found: {config_file}")
        path_to_load = config_file
    else:
        # Fallback logic if no file is provided via the command line.
        # 1. Check current working directory for a file named 'config.ini'
        if os.path.exists('config.ini'):
            path_to_load = 'config.ini'
        else:
            # 2. Check for the default config.ini bundled with the application
            if getattr(sys, 'frozen', False):
                # Path when running as a PyInstaller executable
                base_path = sys._MEIPASS
                default_path = os.path.join(base_path, 'config.ini')
            else:
                # Path when running from source
                base_path = os.path.dirname(os.path.abspath(__file__))
                default_path = os.path.join(base_path, '..', 'config.ini')

            if os.path.exists(default_path):
                path_to_load = default_path

    if path_to_load:
        config.read(path_to_load)
        settings._load_from_parser(config)
    else:
        raise FileNotFoundError(
            "config.ini not found. Please provide the path to one using the --config option, "
            "or place a 'config.ini' file in the current directory or the project's root directory."
        )
