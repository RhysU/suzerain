#!/usr/bin/env python
"""Usage: plot_status.py *.dat
Generate plots from L2.dat, rms.fluct.dat, or bulk.dat files.
Options:
  -f  --file_ext  Output file extension. Default is 'eps'.
  -h  --help      This help message.
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


def getfilekind(infile):

    if infile.find('rms.fluct.dat')!=-1:
       kind = 'rms_'
    elif infile.find('L2.dat')!=-1:
       kind = 'L2_'
    elif infile.find('bulk.dat')!=-1:
       kind = 'bulk_'
    else:
       kind = ''

    return kind


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


def getstatus(infile,nskip):

    # Load data from file, except the header
    data = np.genfromtxt(infile, skip_header=1)

    # Resize if there is only one line
    if data.shape[0] == data.size:
      data = np.array(data).reshape(1, data.size)

    # Delete unused columns
    data = np.delete(data,[0,1,2,4],axis=1)

    # Keep every (nskip)th row
    r     = data.shape[0]
    data  = data[np.s_[1:r:nskip]]

    return data
    

def main(argv=None):

    # Permit interactive use
    if argv is None:
        argv = sys.argv

    # Default values
    nskip = 20
    fileext = 'eps'

    # Parse and check incoming command line arguments
    try:
        try:
	  opts, args = getopt.getopt(argv[1:], "hs:f:", ["help","skip=","file_type="])
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                return 0
            if o in ("-s", "--skip="):
                nskip = int(float(a))
            if o in ("-f", "--file_type="):
                fileext = a
        if len(args) < 1:
            print >>sys.stderr, "Incorrect number of arguments.  See --help."
            return 2
    except Usage, err:
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

    # TODO: Add sanity check for first input file
    # Get file "kind" from first input file
    kind_table = getfilekind(infiles[0])

    # Get file "kind" from first input file
    head_table = getheader(infiles[0])

    # Process each file in turn  
    for infile in infiles:
        kind = getfilekind(infile)
	if kind != kind_table:
	   print "Skipping ", infile, " (different kind)"
	else:
           print "Processing", infile
           head  = getheader(infile)
           # read data if the header is not empty
	   if (head != 0 and head==head_table):
              data  = getstatus(infile,nskip)
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
        key = kind_table + head_table[figid]
        pyplot.plot(data_table[:,0], data_table[:,figid], '-o',linewidth=3, label=key)
        pyplot.legend(loc=0)
        pyplot.savefig(key + '.' + fileext, bbox_inches='tight')
  
if __name__ == "__main__":
    sys.exit(main())
