//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Suzerain is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------

/** @file
 * Implementation of \ref driver_base::summary_run.
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/support/driver_base.hpp>

namespace suzerain {

namespace support {

const char * const driver_base::summary_description =
"Invocable in four distinct ways:\n"
"\n"
"  1) perfect_summary                INFILE.h5 ...\n"
"\n"
"     This first way processes each INFILE.h5 in turn outputting a\n"
"     corresponding INFILE.mean containing a whitespace-separated table\n"
"     of means from the first samples in the file.  Useful primarily for\n"
"     quick plotting of a single snapshot.\n"
"\n"
"  2) perfect_summary -s             INFILE.h5 ...\n"
"\n"
"     This second way (-s) sends the data from all samples to standard\n"
"     output sorted according to the simulation time with a blank line\n"
"     separating adjacent times.  Useful primarily for quick plotting of\n"
"     multiple snapshots.\n"
"\n"
"  3) perfect_summary -f OUTFILE.dat INFILE.h5 ...\n"
"\n"
"     This third way (-f) is identical to the second except the output is\n"
"     automatically sent to the file named OUTFILE.dat.\n"
"\n"
"  4) perfect_summary -o OUTFILE.h5  INFILE.h5 ...\n"
"\n"
"     This fourth way (-o) outputs a single HDF5 file called OUTFILE.h5\n"
"     combining all samples.  Additionally, automatic autocorrelation\n"
"     analysis using autoregresive modeling techniques is run on the\n"
"     combined samples and output as HDF5 attributes.\n"
"\n"
"Options -s, -f, and -o may be specified simultaneously.\n";

const char * const driver_base::summary_argument_synopsis
    = "RESTART-OR-SAMPLE-HDF5-FILE...";

void
driver_base::summary_run(int argc, char **argv)
{
    // FIXME Implement
    return;
}

} // end namespace support

} // end namespace suzerain
