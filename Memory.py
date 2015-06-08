#!/opt/local/bin/python-3.4.3
# -*- coding: utf-8 -*-

# Memory usage
import resource

def memory_usage_resource():

    rusage_denom = 1024.
    if sys.platform == 'darwin':
        # ... it seems that in OSX the output is different units ...
        rusage_denom = rusage_denom * rusage_denom
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
    return mem
