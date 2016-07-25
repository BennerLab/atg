import os
import sys
import configparser

USER_CONFIG_LOCATION = os.path.expanduser('~/.ATGConfig.txt')
INSTALL_CONFIG_LOCATION = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, "ATGConfig.txt"))

settings = configparser.ConfigParser()

if os.path.exists(USER_CONFIG_LOCATION):
    settings.read(USER_CONFIG_LOCATION)
elif os.path.exists(INSTALL_CONFIG_LOCATION):
    settings.read(INSTALL_CONFIG_LOCATION)
else:
    print("Could not load ATG configuration file.", file=sys.stderr)
    sys.exit(1)
