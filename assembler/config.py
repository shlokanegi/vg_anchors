import configparser
import os
import sys

settings = None

def load_config(config_file=None):
    """
    Loads configuration from a given file or discovers it.
    This function populates the global 'settings' variable.
    """
    global settings
    
    config = configparser.ConfigParser()
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
        settings = config['DEFAULT']
    else:
        raise FileNotFoundError(
            "config.ini not found. Please provide the path to one using the --config option, "
            "or place a 'config.ini' file in the current directory or the project's root directory."
        )
