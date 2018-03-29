
import sys

def handle_arg(var, alt, error):
    """Argument parsing and error handling"""
    
    if var is None:
 
        if alt is not None:
            var = alt
        
        if var is None:
            raise ValueError(error)
            sys.exit(1)
        
    return var