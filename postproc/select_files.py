#!/bin/py
"""Usage: python get_files.py OPTIONS HDF5FILE
Given the input files (in general not sorted),
this script selects the input files such that the sorted time
series (t0, t1, ..., tn) contains the ref_file (if given), and
the time difference between files is dt_target, with a relative
dt_tolerance smaller than dt_target*tolerance. To set the
processing parameters refer to the Options description below.
Also:
* The files with "time" out of tolerance are discarded.
* The files selected for processing are those with smaller
  time discrepancy with respect to the exact time series.
* If there are gaps in the final time series the
  output file is not written, unless the option
  --ignore_gaps is set.
* The results are stored in the "files.txt"
  file. If the file exists the output is not written
  unless the --clobber option is set (in which
  case the existing output file is deleted).
Options:
      -h  --help        This help message.
      -r  --ref_file=   File to use for time reference.
                        The time of the first input file
                        is used by default.
          --clobber     Overwrite old output file.
          --dt_target=  Target time difference between files.
                        The time difference of the first
                        two files is used by default.
          --ignore_gaps Process files even if there
                        are gaps in the time series of
                        files to process. The computation
                        doesn't proceed by default.
          --tolerance=  Maximum fraction of delta t
                        abs(target time-file time)/dt_target
                        for a file to be processed. The
                        tolerance value should be between
                        0 and 0.5. Default is 0.1.
          --min_t=      Minimum value of t to include
          --max_t=      Maximum value of t to include

"""

import getopt
import glob
import h5py
import math
import numpy as np
import os
import sys


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def get_files_to_process(files,
                         t_ref_file = "",
                         dt_target  = 0.0,
                         tolerance  = 0.1,
                         min_t      = float("-inf"),
                         max_t      = float("inf")):

    # Set target delta t and tolerance
    dt_tgt    = float(dt_target)
    tolerance = float(tolerance)

    # Pick the first file in input to set the time
    # reference if a reference file is not provided explicitly
    if t_ref_file == "":
        t_ref_file = files[0]

    # Grab reference time
    ff  = h5py.File(t_ref_file,'r')
    t_ref=ff['t'].value[0]
    ff.close()
    #print "Reference time file: ", t_ref_file, t_ref

    # Reset dt_tgt if zero to the time difference
    # of the first two input files
    if dt_tgt == 0:
        ff  = h5py.File(files[0],'r+')
        t0  =ff['t'].value[0]
        ff.close()
        ff  = h5py.File(files[1],'r+')
        t1  =ff['t'].value[0]
        ff.close()
        dt_tgt = abs(t0-t1)
        print dt_tgt

    # Pre-select files to process
    tt    = []
    files_to_process = []
    print "Selecting files to process"
    for fi in files:
        if not os.path.isfile(fi):
            print "File ", fi, " not found"
            continue

        f  = h5py.File(fi,'r')

        # Grab time
        time=f['t'].value[0]

        # Decide if it should process it or skip it
        # If time is within tolerance and a file
        # has been selected for that slot, keep the
        # one with less difference wrt the ideal
        # time for that slot
        skip_this_file = False
        is_time_accounted_for = False
        replace_for_this_file = False
        time_mod = ((time-t_ref) % dt_tgt) / dt_tgt
        if time_mod > 0.5: time_mod = 1 - time_mod
        # Check if time is within cuttoffs:
        if time<min_t or time>max_t:
            skip_this_file = True
        # Check if time is within tolerance:
        elif (time_mod < tolerance):
            # This is a candidate file,
            # now check if this time has been accounted for
            for i in xrange(len(tt)):
                time_accounted_for = tt[i]
                if ((np.abs(time-time_accounted_for)/dt_tgt) < 0.5):
                    # This time slot has been accounted for, skip
                    # adding it to the list
                    is_time_accounted_for = True
                    skip_this_file = True
                    # ... but check if it should replace the existing one
                    # because of being closer to the corresponding
                    # target time
                    time_mod0 = ((time_accounted_for-t_ref) % dt_tgt) / dt_tgt
                    if time_mod0 > 0.5: time_mod0 = 1 - time_mod0
                    # print time_mod, time_mod0, time_mod < time_mod0
                    if time_mod < time_mod0:
                        # .. replace the existing one by the current candidate
                        del files_to_process[i]
                        del tt[i]
                        files_to_process.insert(i, os.path.realpath(fi).rstrip())
                        tt.insert(i, time)
                        replace_for_this_file = True
                        break
        else:
            # Skip this file because is not within tolerance
            skip_this_file = True

        # print "time                  = ", time
        # print "time_mod              = ", time_mod
        # print "tolerance             = ", tolerance
        # print "Is time mod within tolerance? ", (time_mod < tolerance)
        # print "skip_this_file        = ", skip_this_file
        # print "is_time_accounted_for = ", is_time_accounted_for
        # print "replace_for_this_file = ", replace_for_this_file

        # Now ignore the file or add it to the list
        if skip_this_file:
            # Continue with the next file without adding this one
            # to the list
            continue
        else:
            # Add this file to the list of files to process
            files_to_process.append(os.path.realpath(fi).rstrip())
            tt.append(time)
    #print files_to_process

    # Check for gaps in the time series
    gaps = False
    dts = np.diff(sorted(tt))
    for dt in dts:
        if dt > dt_tgt * (1+2*tolerance):
            gaps = True
            #print 'Found gap: ', dt, dt_tgt
            break

    # Sort files in increasing simulation time order
    all_files = np.array(tt).reshape(len(tt),1)
    all_files = np.concatenate((all_files, np.array(files_to_process).reshape(len(files_to_process),1)), axis=1)

    all_files_sort = all_files[all_files[:,0].argsort()]
    all_files = all_files.transpose()
    all_files_sort = all_files_sort.transpose()
    #print all_files_sort

    # Return values
    retval = {}
    retval['files'] = all_files_sort
    retval['gaps']  = gaps

    return retval


def main(argv=None):

    # Set defaults
    clobber     = False
    dt_target   = 0.0
    ignore_gaps = False
    outfile     = "files.txt"
    ref_file    = ""
    tolerance   = 0.1
    min_t       = float("-inf")
    max_t       = float("inf")

    # Permit interactive use
    if argv is None:
        argv = sys.argv

    # Parse and check incoming command line arguments
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hr:",
        ["help", "ref_file=", "dt_target=", "tolerance=",
         "ignore_gaps","clobber",
         "min_t=", "max_t="
        ])
        except getopt.error, msg:
            raise Usage(msg)

        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                return 0
            if o in ("-r", "--ref_file"):
                ref_file = a
            if o in ("--dt_target"):
                dt_target=float(a)
            if o in ("--tolerance"):
                tolerance=float(a)
            if o in ("--ignore_gaps"):
                ignore_gaps=True
            if o in ("--clobber"):
                clobber = True
            if o in ("--min_t"):
                min_t = a
            if o in ("--max_t"):
                max_t = a

        # Check for proper number of arguments
        if len(args) < 1:
            print args, len(args)
            print >> sys.stderr, "Incorrect number of arguments. See --help."
            return 2

        # Check for reference time file
        if ref_file != "" and not os.path.isfile(ref_file):
            print >> sys.stderr, "Reference file not found"
            return 2

        # Check for output file
        if os.path.isfile(outfile):
            if clobber:
                os.system("rm " + outfile)
            else:
                print >> sys.stderr, ("Output file already exists.\n"
                "Use option --clobber to overwrite.")
                return 2

        # Check for tolerance range
        if tolerance < 0.0 or tolerance > 0.5:
            print >> sys.stderr, "Option 'tolerance' out of range. See --help."
            return 2

        # Check time cuttofs
        if not is_number(min_t):
            print >> sys.stderr, "min_t is not a number. See --help."
            return 2
        else:
            min_t = float(min_t)

        if not is_number(max_t):
            print >> sys.stderr, "max_t is not a number. See --help."
            return 2
        else:
            max_t = float(max_t)

        if min_t > max_t:
            print >> sys.stderr, "min_t is larger than max_t. See --help."
            return 2

    except Usage, err:
        print >> sys.stderr, err.msg
        return 2

    # Input files
    hdf5files = args

    # Get sorted list of files
    files_dict = get_files_to_process(hdf5files, ref_file, dt_target, tolerance,
                                      min_t, max_t)
    files = files_dict['files']
    gaps  = files_dict['gaps']

    # Print sorted files and dt
    print files[1,0], files[0,0]
    for i in xrange(1,len(files[0,:])):
        print files[1,i], files[0,i], float(files[0,i]) - float(files[0,i-1])

    # Output to text file
    if not gaps:
        print "There are no gaps in the time series"
        np.savetxt(outfile, files[1,:], fmt='%s')
    else: # there are gaps in the time series
        if ignore_gaps:
            print ("WARNING: There are gaps in the time series."
                   " Generating ", outfile, " due to --ignore_gaps.")
            np.savetxt(outfile, files[1,:], fmt='%s')
        else: # ignore gaps and compute anyway
            print >>sys.stderr, "There are gaps in the time series."
            return 2

    sys.exit("Done")


if __name__ == "__main__":
    sys.exit(main())
