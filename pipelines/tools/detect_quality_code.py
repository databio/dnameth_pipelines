#!/usr/bin/env python
'''
	@author Brent Pedersen
    @author Nathan Sheffield
This script takes as input a fastq file and returns the quality encoding system
used by that fastq file.
The update makes it so by default, the script only tries to differentiate between
phred33 and phred64, instead of exact flavor, making it much faster for practical
application. 
'''
import sys
from argparse import ArgumentParser 


def get_qual_range(qual_str):
    """
    >>> get_qual_range("DLXYXXRXWYYTPMLUUQWTXTRSXSWMDMTRNDNSMJFJFFRMV")
    (68, 89)
    """
 
    vals = [ord(c) for c in qual_str]
    return min(vals), max(vals)
 
def get_encodings_in_range(rmin, rmax, ranges):
    valid_encodings = []
    for encoding, (emin, emax) in ranges.items():
        if rmin >= emin and rmax <= emax:
            valid_encodings.append(encoding)
    return valid_encodings
 

if __name__ == "__main__":
    parser = ArgumentParser(description='Detect quality codes.')
    parser.add_argument("-n", dest="n", help="number of qual lines to test default:-1"
                 " means test until end of file or until it it possible to "
                 " determine a single file-type", default=-1, type=int)

    parser.add_argument("-f", dest="filename", help="Filename to test", required=True)
    parser.add_argument("-c", dest="complex", help="Complex mode: Show the actual encoding system, not just the phred type", action="store_true", default=False)

    opts = parser.parse_args()
    if opts.complex:
        RANGES = {
    'phred33/Sanger': (33, 73),
    'phred33/Illumina-1.8': (33, 74),
    'phred33/Sanger2': (33, 75), # Another Illumina update
    'phred64/Solexa': (59, 104),
    'phred64/Illumina-1.3': (64, 104),
    'phred64/Illumina-1.5': (67, 104)
    }
    else:
        # Actually I don't care which exact one it is, I just want to know how to convert it. This will save you lots of time if you just need this much info.
        RANGES = {
	'phred33': (33,75),
	'phred64': (59,104)
}

    #print("# reading qualities from stdin")
    #print("RANGES:", RANGES)
    txt = open(opts.filename)
    gmin, gmax  = 99, 0
    valid = []
    for i, line in enumerate(txt):
        if ((i+1) % 4) != 0:
            continue
        #print i, i+1 % 4, line
        lmin, lmax = get_qual_range(line.rstrip())
        if lmin < gmin or lmax > gmax:
            gmin, gmax = min(lmin, gmin), max(lmax, gmax)
            valid = get_encodings_in_range(gmin, gmax, RANGES)
            if len(valid) == 0:
                print("no encodings for range: %s" % str((gmin, gmax)))
                raise SystemExit
            if len(valid) == 1 and opts.n == -1:
                print "\t".join(valid) + "\t" + str((gmin, gmax)) + "\t" + str(i)
                raise SystemExit
 
        if opts.n > 0 and i > opts.n:
            print("\t".join(valid) + "\t" + str((gmin, gmax)))
            raise SystemExit
 
    print("\t".join(valid) + "\t" + str((gmin, gmax)))


