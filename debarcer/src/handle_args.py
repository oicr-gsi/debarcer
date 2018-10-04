import os
import sys
import argparse
import configparser


def handle_arg(var, alt, error):
    """Argument parsing and error handling"""
    
    if var is None:
 
        if alt is not None:
            var = alt
        
        if var is None:
            raise ValueError(error)
            sys.exit(1)
        
    return var




def arg_exists(arguments):
    for var in arguments:
        if ("/" in var and not os.path.exists(var)):
            print("ERR: Path or file, " + var + " does not exist.")




def config_validation(conf_paths):
    #conf_paths = dict(config.items('PATHS'))
    for key, val in conf_paths.items():
        if (not os.path.exists(val)):
            print("ERR: Path provided, for " + key + ", in the config file does not exist.")

