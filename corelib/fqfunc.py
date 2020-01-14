
import math
from collections import defaultdict
def fqFunc(cur_reads):
    read_info = {
        'reads' : 0,
        'length' : 0,
        'sites' : 0, 
        'read_lens' : defaultdict(int),
        'base_comp' : { "A" : 0, "T" : 0, "C" : 0, "G" : 0, "N" : 0 },
        'qual_pos' : defaultdict(int),
        'site_pos' : defaultdict(int)
        }
    # All the info possibily collected for a set of reads is compiled in this dictionary.

    for read in cur_reads:
        read_info['reads'] += 1;
        readlen = len(read['seq']);
        # Always get the number of reads and the read lengths.
    
        read_info['sites'] += readlen;
        cur_bin = math.floor(readlen/5)*5;
        read_info['read_lens'][cur_bin] += 1;
        # Get the read length bins if that option is specified.
    
    return read_info;
    