#!/usr/bin/env python
"""Usage: remove_hdf5.py [-n] PATTERN HDF5FILE...
Remove datasets matching per Python's re.search(PATTERN, ...) within HDF5FILE.

Options:
    -n Don't actually change anything, just print what actions would be taken.
"""
import sys
import getopt
import h5py
import re

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def remove_hdf5(pattern, hdf5file, dryrun=False):
    print "Processing", hdf5file

    if dryrun:
        f = h5py.File(hdf5file, 'r')
    else:
        f = h5py.File(hdf5file, 'r+')

    # list(f.iterkeys()) as we mutate during iteration
    for name in list(f.iterkeys()):
        if re.search(pattern, name):
            print "Removing", name
            if not dryrun:
                del f[name]

    f.close()

def main(argv=None):

    # Permit interactive use
    if argv is None:
        argv = sys.argv

    # Parse and check incoming command line arguments
    dryrun = False
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hn", ["help"])
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o == "-n":
                dryrun = True
            elif o in ("-h", "--help"):
                print __doc__
                return 0
        if len(args) < 2:
            print >>sys.stderr, "Too few arguments.  See --help."
            return 1
    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

    # Process each file in turn
    pattern, hdf5files = args[0], args[1:]
    for hdf5file in hdf5files:
        remove_hdf5(pattern, hdf5file, dryrun)

if __name__ == "__main__":
    sys.exit(main())
