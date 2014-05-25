#!/usr/bin/env python
"""Usage: plot_status.py *.dat
Generate plots from suzerain status files.
Options:
  -f  --file_ext  Output file extension. Default is 'eps'.
  -h  --help      This help message.
      --quantity  Status quantity to plot.
  -s  --skip      Lines to skip in input files. Default is 20.
"""

import sys
import getopt
import h5py
import numpy as np
import gb
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

def getheader(infile):

    # Open file and read header
    f = open(infile,'r')

    # Read first line
    head = f.readline()
    head = head.split()

    # Return if file is empty
    if len(head) == 0:
      return 0

    # Delete unused columns
    del head[4]    # delete 'nt column'
    del head[0:3]  # delete first three columns

    return head


def getquantitykeys(quantity):

    keys = []
    if quantity == 'rms':
        keys.append('rms.fluct')
        keys.append('state.RMS')
    elif quantity == 'L2':
        keys.append('L2.mean')
        keys.append('state.L2')
    elif quantity == 'bulk':
        keys.append('bulk.state')
        keys.append('state.bulk')
    elif quantity == 'min':
        keys.append('state.min')
    elif quantity == 'xmin':
        keys.append('state.xmin')
    elif quantity == 'ymin':
        keys.append('state.ymin')
    elif quantity == 'zmin':
        keys.append('state.zmin')
    elif quantity == 'max':
        keys.append('state.max')
    elif quantity == 'xmax':
        keys.append('state.xmax')
    elif quantity == 'ymax':
        keys.append('state.ymax')
    elif quantity == 'zmax':
        keys.append('state.zmax')
    elif quantity == 'fneg':
        keys.append('state.fneg')
    elif quantity == 'lower.d0':
        keys.append('bc.lower.d0')
    elif quantity == 'lower.d1':
        keys.append('bc.lower.d1')
    elif quantity == 'lower.d2':
        keys.append('bc.lower.d2')
    elif quantity == 'upper.d0':
        keys.append('bc.upper.d0')
    elif quantity == 'upper.d1':
        keys.append('bc.upper.d1')
    elif quantity == 'upper.d2':
        keys.append('bc.upper.d2')

    return keys


def getstatus(infile,nskip,keys):

    data = []
    f = open(infile,'r')
    iskip = 0
    # walk all lines in file
    for line in f:
        # check if any of the keys for this quantity
        # is in the line
        for key in keys:
            if key in line:
	        row = line.split()
                # delete unused columns
	        del row[4]
	        row = row[3:]
	        try:
		    # try to converto to float
	            row = [float(i) for i in row]
                    # check if it should skip the line
	            if iskip % nskip == 0:
                        data.append(row)
	            iskip += 1
	        except ValueError:
		    # if the conversion falis,
                    # discard and continue
	            continue
		# if key was in line, continue with the next line
	        continue

    # convert to numpy array and return
    data = np.array(data)
    return data


def main(argv=None):

    # Permit interactive use
    if argv is None:
        argv = sys.argv

    # Default values
    nskip = 20
    fileext  = 'eps'
    quantity = 'rms'

    # Parse and check incoming command line arguments
    try:
        try:
	  opts, args = getopt.getopt(argv[1:], "hs:f:", ["help","skip=","file_type=","quantity="])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                return 0
            if o in ("-s", "--skip"):
                nskip = int(float(a))
            if o in ("-f", "--file_type"):
                fileext = a
            if o in ("--quantity"):
                quantity = a
        if len(args) < 1:
            print >>sys.stderr, "Incorrect number of arguments.  See --help."
            return 2
    except Usage as err:
        print >>sys.stderr, err.msg
        return 2

    # Get file names and number
    infiles = args
    nfiles  = len(infiles)

    # Get data from multiple files
    # and append the results to data_table
    data_table = np.empty([0])

    # Keep track of the number of rows of data_table
    rows = 0

    # First input file determines the file
    # kind (rms, L2, bulk), and the header

    # Get suzerain keys for the requested quantity
    keys = getquantitykeys(quantity)

    # Get header from first input file
    head_table = getheader(infiles[0])

    # Process each file in turn
    for infile in infiles:
        print "Processing", infile
        head  = getheader(infile)
        # read data if the header is not empty
        if (head != 0 and head==head_table):
            data  = getstatus(infile,nskip,keys)
            rows += data.shape[0]
            cols  = data.shape[1]
            data_table = np.append(data_table, data)
            head_table = head

    # Reshape table and sort by time
    data_table = np.array(data_table).reshape(rows, cols)
    data_table = data_table[data_table[:,0].argsort()]

    # Plot all files
    for figid in xrange(1,cols):
        pyplot.figure(figid)
        key = quantity + '_' + head_table[figid]
        pyplot.plot(data_table[:,0], data_table[:,figid], '-o',linewidth=3, label=key)
        pyplot.legend(loc=0)
        pyplot.savefig(key + '.' + fileext, bbox_inches='tight')

if __name__ == "__main__":
    sys.exit(main())
