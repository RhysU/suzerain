/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Suzerain is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with <APP/LIBRARY>.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * driver_sine.cc: A P3DFFT test driver based on work by Dmitry Pekurovsky
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include "config.h"

#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp> 
#include <boost/program_options.hpp>
#include <boost/shared_array.hpp>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <log4cxx/logger.h>
#include <mpi.h>
#include <numeric>
#include "p3dfft_d.h"
#include <sstream>
#include <string>
#include <vector>

// TODO: Place into pecos::suzerain namespace

namespace ublas = boost::numeric::ublas; // Shorthand, TODO remove in some way

#define ONLYPROC0(expr) { if (!procid) { expr ; } };

double FORTNAME(t1),FORTNAME(t2),FORTNAME(t3),FORTNAME(t4),FORTNAME(tp1);

void print_all(double *A,long int nar,int procid,long int Nglob)
{
  int x,y,z,Fstart[3],Fsize[3],Fend[3];

  get_dims(Fstart,Fend,Fsize,2);
  Fsize[0] *= 2;
  Fstart[0] = 1 + (Fstart[0]-1)*2;
  for (long int i=0; i < nar; i++)
    if (abs(A[i]) > Nglob *1.25e-4)
      {
        z = i/(Fsize[0]*Fsize[1]);
        y = i/Fsize[0] - z*Fsize[1];
        x = i-z*Fsize[0]*Fsize[1] - y*Fsize[0];
        printf("(%d,%d,%d) %f\n",x+Fstart[0],y+Fstart[1],z+Fstart[2],A[i]);
      }
}

// Print usage information
template<typename charT, typename traits>
void print_help(std::basic_ostream<charT, traits>& out,
                const std::string application_name,
                const boost::program_options::options_description options)
{
  using namespace std;

  out << endl
      << "Usage: " << application_name << " [OPTION] [FILE]..." << endl
// TODO: Provide a legitimate description for the help message
//      << endl
//      << "Description: " << endl
//      << endl
      << options
      << endl;
}

// Print version information
template<typename charT, typename traits>
void print_version(std::basic_ostream<charT, traits>& out)
{
  out << PACKAGE_STRING
      << " (built " __DATE__ " " __TIME__ ")" << std::endl;
}

int main(int argc,char **argv)
{
  int nproc;        // Number of processors in MPI environment
  int procid;       // This processor's global processor ID
  int nx, ny, nz;   // Domain dimensions in x, y, and z directions
  int dims[2];      // Processor grid dimensions in 1st, 2nd directions
  int nrep;         // Number of times to repeat the test

  MPI_Init(&argc, &argv);                   // Initialize MPI on startup
  atexit((void (*) ()) MPI_Finalize);       // Finalize down MPI at exit
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);    // TODO Const-ness
  MPI_Comm_rank(MPI_COMM_WORLD, &procid);   // TODO Const-ness

  // TODO Compute width from magnitude of nproc
  std::ostringstream procname;
  procname << "proc" << std::setfill('0') << std::setw(3) << procid;
  log4cxx::LoggerPtr logger = log4cxx::Logger::getLogger(procname.str());

  // Find default processor grid size based on nproc
  dims[0] = dims[1] = 0;                // Zeroed for MPI_Dims_create
  MPI_Dims_create(nproc, 2, dims);      // Find a Cartesian grid
  if (dims[0] > dims[1])                // Ensure first dimension smaller
    {
      std::swap(dims[0], dims[1]);
    }

  namespace po = boost::program_options;

  // Options accepted on command line and in configuration file
  po::options_description desc_config("Configuration options");
  desc_config.add_options()
    ("nx",  po::value<int>(&nx)->default_value(16),
        "Domain grid size in X direction")
    ("ny",  po::value<int>(&ny)->default_value(16),
        "Domain grid size in Y direction")
    ("nz",  po::value<int>(&nz)->default_value(16),
        "Domain grid size in Z direction")
    ("rep", po::value<int>(&nrep)->default_value(1),
        "Number of repetitions to perform for timing purposes")
    ("pg1", po::value<int>(&dims[0])->default_value(dims[0]),
        "Processor grid size in first direction.")
    ("pg2", po::value<int>(&dims[1])->default_value(dims[1]),
        "Processor grid size in second direction.")
    ;

  // Options allowed only on command line
  po::options_description desc_clionly("Program information");
  desc_clionly.add_options()
    ("help,h",    "show usage information")
    ("version,v", "print version string")
    ;

  // Options allowed on command line and in configuration file
  // Not shown to the user
  po::options_description desc_hidden("Hidden options");
  desc_hidden.add_options()
    ("input-file", po::value< std::vector<std::string> >(), "input file")
    ;

  // Build the options acceptable on the CLI, in a file, and in help message
  po::options_description opts_cli;
  opts_cli.add(desc_config).add(desc_hidden).add(desc_clionly);
  po::options_description opts_file;
  opts_file.add(desc_config).add(desc_hidden);
  po::options_description opts_visible;
  opts_visible.add(desc_config).add(desc_clionly);

  // Have positional parameters act like input-file
  po::positional_options_description opts_positional;
  opts_positional.add("input-file", -1);

  // Parse all the command line options
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
    options(opts_cli).positional(opts_positional).run(), vm);

  // Process command-line only parameters
  if (vm.count("help"))
    {
      ONLYPROC0(print_help(std::cout, argv[0], opts_visible));
      exit(0);
    }
  if (vm.count("version"))
    {
      ONLYPROC0(print_version(std::cout));
      exit(0);
    }

  // Parse any input files provided on the command line
  if (vm.count("input-file"))
    {
      BOOST_FOREACH(std::string filename,
                    vm["input-file"].as< std::vector<std::string> >())
        {
          ONLYPROC0(LOG4CXX_DEBUG(logger, "Reading input file " << filename));
          std::ifstream ifs( (vm["input-file"].as< std::string >()).c_str() );
          po::store(po::parse_config_file(ifs, opts_file), vm);
        }
    }

  // Perform options callbacks now that we're done parsing options
  po::notify(vm);

  ONLYPROC0(LOG4CXX_DEBUG(logger, "Number of processors: " << nproc));
  ONLYPROC0(LOG4CXX_DEBUG(logger, "Physical grid dimensions: "
      << boost::format("(%d, %d, %d)") % nx % ny %nz));
  ONLYPROC0(LOG4CXX_DEBUG(logger, "Processor grid dimensions: "
      << boost::format("(%d, %d)") % dims[0] % dims[1]));
  if (dims[0]*dims[1] != nproc)
    {
      ONLYPROC0(LOG4CXX_WARN(logger,
          "Processor grid dimensions incompatible with number of processors"));
    }

  /* Initialize P3DFFT */
  p3dfft_setup(dims, nx, ny, nz, 1 /* safe to overwrite btrans */);

  /* Get dimensions for input and output arrays */
  ublas::c_vector<int, 3> istart(3), isize(3), iend(3);
  ublas::c_vector<int, 3> fstart(3), fsize(3), fend(3);

  get_dims(istart.data(), iend.data(), isize.data(), 1 /* physical pencil */);
  get_dims(fstart.data(), fend.data(), fsize.data(), 2 /* wave pencil */);

  ublas::vector<double> sinx(nx);
  ublas::vector<double> siny(ny);
  ublas::vector<double> sinz(nz);

  for (ublas::vector<double>::size_type i=0; i < isize[0]; ++i) {
    sinx[i] = sin((i+istart[0]-1)*2.0*M_PI/nx);
  }
  for (ublas::vector<double>::size_type i=0; i < isize[1]; ++i) {
    siny[i] = sin((i+istart[1]-1)*2.0*M_PI/ny);
  }
  for (ublas::vector<double>::size_type i=0; i < isize[2]; ++i) {
    sinz[i] = sin((i+istart[2]-1)*2.0*M_PI/nz);
  }

  /* Allocate and initialize state space */
  ublas::shallow_array_adaptor<double> 
    dataA(std::max( isize[0]*isize[1]*isize[2], fsize[0]*fsize[1]*fsize[2]*2));
  ublas::vector<double, ublas::shallow_array_adaptor<double> > 
    vecA(dataA.size(), dataA);
  double * const A = &dataA[0];

  {
    double *p = A;
    for (int z = 0; z < isize[2]; z++)
      for (int y=0; y < isize[1]; y++)
        {
          const double sinyz = siny[y]*sinz[z];
          for (int x=0; x < isize[0]; x++)
            *p++ = sinx[x]*sinyz;
        }
  }

  const long int Ntot  = fsize[0]*fsize[1]*fsize[2]*2;
  const long int Nglob = nx*ny*nz;
  const double factor = 1.0/Nglob;

  double rtime1 = 0.0;
  for (int m=0; m < nrep; m++)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      rtime1 = rtime1 - MPI_Wtime();
      ONLYPROC0(LOG4CXX_DEBUG(logger, "Iteration " << m));

      p3dfft_ftran_r2c(A,A);          // Physical to wave transform
      rtime1 = rtime1 + MPI_Wtime();

      ONLYPROC0(LOG4CXX_DEBUG(logger, "Forward transform results "));
      print_all(A,Ntot,procid,Nglob);
      vecA *= factor;                 // normalize for grid size

      MPI_Barrier(MPI_COMM_WORLD);
      rtime1 = rtime1 - MPI_Wtime();
      p3dfft_btran_c2r(A,A);          // Wave to physical transfrom
      rtime1 = rtime1 + MPI_Wtime();

    }
  p3dfft_clean();   // Free work space

  /* Check results */
  double cdiff = 0.0;
  {
    double *p = A;
    for (int z=0; z < isize[2]; z++)
      for (int y=0; y < isize[1]; y++)
        {
          const double sinyz = siny[y]*sinz[z];
          for (int x=0; x < isize[0]; x++)
            {
              const double ans = sinx[x]*sinyz;
              if (cdiff < abs(*p - ans))
                cdiff = abs(*p - ans);
              p++;
            }
        }
  }

  // Gather error indicator
  double ccdiff = 0.0;
  MPI_Reduce(&cdiff,&ccdiff,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  ONLYPROC0(LOG4CXX_DEBUG(logger, 
      "Maximum difference: " << std::scientific << ccdiff));

  // Gather timing statistics
  double rtime2 = 0.0;
  MPI_Reduce(&rtime1,&rtime2,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  ONLYPROC0(LOG4CXX_DEBUG(logger, "Time per loop: " << rtime2/((double)nrep)));
}
