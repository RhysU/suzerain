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
 * along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * driver_sine.cc: A P3DFFT test driver based on work by Dmitry Pekurovsky
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <suzerain/config.h>
#include <suzerain/common.hpp>
#pragma hdrstop

#include <log4cxx/logger.h>
#include <mpi.h>
#include <p3dfft_d.h>

#include <suzerain/pencil_grid.hpp>
#include <suzerain/pencil.hpp>

#define ONLYPROC0(expr) { if (!procid) { expr ; } };

double real_data(const double x, const double y, const double z) {
    return   1.0* 1.0
           + 2.0* 3.0*sin(x) // + 2.0* 6.0*sin(2*x) + 2.0* 9.0*sin(3*x)
           + 2.0* 5.0*sin(y) // + 2.0*10.0*sin(2*y) + 2.0*15.0*sin(3*y)
           + 2.0* 7.0*sin(z) // + 2.0*14.0*sin(2*z) + 2.0*21.0*sin(3*z)
           + 4.0*11.0*sin(x)*sin(y)
           + 4.0*13.0*sin(x)*sin(z)
           + 4.0*17.0*sin(y)*sin(z)
           + 8.0*19.0*sin(x)*sin(y)*sin(z)
           ;
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

int main(int argc, char **argv)
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

    if (dims[0] > dims[1]) {              // Ensure first dimension smaller
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
    if (vm.count("help")) {
        ONLYPROC0(print_help(std::cout, argv[0], opts_visible));
        exit(0);
    }

    if (vm.count("version")) {
        ONLYPROC0(print_version(std::cout));
        exit(0);
    }

    // Parse any input files provided on the command line
    if (vm.count("input-file")) {
        BOOST_FOREACH(std::string filename,
                      vm["input-file"].as< std::vector<std::string> >()) {
            ONLYPROC0(LOG4CXX_DEBUG(logger, "Reading input file " << filename));
            std::ifstream ifs( (vm["input-file"].as< std::string >()).c_str() );
            po::store(po::parse_config_file(ifs, opts_file), vm);
        }
    }

    // Perform options callbacks now that we're done parsing options
    po::notify(vm);

    ONLYPROC0(LOG4CXX_INFO(logger, "Number of processors: " << nproc));

    ONLYPROC0(LOG4CXX_INFO(logger, "Physical grid dimensions: "
                           << boost::format("(% 4d, % 4d, % 4d)") % nx % ny % nz));

    ONLYPROC0(LOG4CXX_INFO(logger, "Processor grid dimensions: "
                           << boost::format("(%d, %d)") % dims[0] % dims[1]));

    if (dims[0]*dims[1] != nproc) {
        ONLYPROC0(LOG4CXX_WARN(logger,
                "Processor grid dimensions incompatible with number of processors"));
    }

    /* Initialize P3DFFT */
    p3dfft_setup(dims, nx, ny, nz, 1 /* safe to overwrite btrans */);

    /* Get dimensions for input and output arrays */
    boost::array<int, 3> istart, isize, iend, fstart, fsize, fend;
    get_dims(istart.data(), iend.data(), isize.data(), 1/* physical pencil */);
    get_dims(fstart.data(), fend.data(), fsize.data(), 2/* wave pencil */);

    /* P3DFFT STRIDE1 wave space sizes return in (Y, X, Z) ordering */
    /* Suzerain's pencils require (X, Y, Z) ordering */
    std::swap(fstart[0], fstart[1]);
    std::swap(fend[0], fend[1]);
    std::swap(fsize[0], fsize[1]);

    LOG4CXX_DEBUG(logger, "Physical space pencil sizes from P3DFFT:           "
                          << boost::format("(% 4d, % 4d, % 4d)")
                          % isize[0] % isize[1] % isize[2]);

    LOG4CXX_DEBUG(logger, "Wave space pencil sizes from P3DFFT:               "
                          << boost::format("(% 4d, % 4d, % 4d)")
                          % fsize[0] % fsize[1] % fsize[2]);

    /* Transform indices for C conventions, ranges like [istart, iend) */
    std::transform(istart.begin(), istart.end(), istart.begin(),
            std::bind2nd(std::minus<int>(),1));
    std::transform(fstart.begin(), fstart.end(), fstart.begin(),
            std::bind2nd(std::minus<int>(),1));


    LOG4CXX_DEBUG(logger, "Physical space pencil sizes modified after P3DFFT: "
                          << boost::format("(% 4d, % 4d, % 4d)")
                          % isize[0] % isize[1] % isize[2]);

    LOG4CXX_DEBUG(logger, "Wave space pencil sizes modified after P3DFFT:     "
                          << boost::format("(% 4d, % 4d, % 4d)")
                          % fsize[0] % fsize[1] % fsize[2]);

    /* Create a uniform tensor product grid */
    std::valarray<double> gridx(isize[0]), gridy(isize[1]), gridz(isize[2]);
    for (size_t i = 0; i < isize[0]; ++i) {
        gridx[i] = (i+istart[0]) * 2*M_PI/nx;
        LOG4CXX_TRACE(logger, boost::format("gridx[%3d] = % 6g") % i % gridx[i]);
    }
    for (size_t j = 0; j < isize[1]; ++j) {
        gridy[j] = (j+istart[1]) * 2*M_PI/ny;
        LOG4CXX_TRACE(logger, boost::format("gridy[%3d] = % 6g") % j % gridy[j]);
    }
    for (size_t k = 0; k < isize[2]; ++k) {
        gridz[k] = (k+istart[2]) * 2*M_PI/nz;
        LOG4CXX_TRACE(logger, boost::format("gridz[%3d] = % 6g") % k % gridz[k]);
    }

    /* Allocate and initialize state space */
    using pecos::suzerain::pencil;
    pencil<> A(istart.data(), isize.data(), fstart.data(), fsize.data());

    LOG4CXX_INFO(logger,
                 "Physical space pencil start and end: "
                 << boost::format("[(%3d, %3d, %3d) ... (%3d, %3d, %3d))")
                 % A.physical.start_x
                 % A.physical.start_y
                 % A.physical.start_z
                 % A.physical.end_x
                 % A.physical.end_y
                 % A.physical.end_z);

    LOG4CXX_INFO(logger,
                 "Wave space pencil start and end:     "
                 << boost::format("[(%3d, %3d, %3d) ... (%3d, %3d, %3d))")
                 % A.wave.start_x
                 % A.wave.start_y
                 % A.wave.start_z
                 % A.wave.end_x
                 % A.wave.end_y
                 % A.wave.end_z);

    MPI_Barrier(MPI_COMM_WORLD);

    for (pencil<>::size_type j = 0; j < A.physical.size_y; ++j) {
        for (pencil<>::size_type k = 0; k < A.physical.size_z; ++k) {
            for (pencil<>::size_type i = 0; i < A.physical.size_x; ++i) {
                const double value = real_data(gridx[i], gridy[j], gridz[k]);
                LOG4CXX_TRACE(logger, boost::format(
                              "Physical space (% 6.4f, % 6.4f, % 6.4f) = % 6g for index (%3d, %3d, %3d)")
                              % gridx[i] % gridy[j] % gridz[k] % value
                              % i % j % k);
                A.physical(i,j,k) = value;
            }
        }
    }

    const long int Ntot  = fsize[0] * fsize[1] * fsize[2] * 2;

    const long int Nglob = nx * ny * nz;

    const double factor = 1.0 / Nglob;

    double rtime1 = 0.0;

    for (int m = 0; m < nrep; m++) {
        MPI_Barrier(MPI_COMM_WORLD);
        rtime1 = rtime1 - MPI_Wtime();
        ONLYPROC0(LOG4CXX_INFO(logger, "Iteration " << m));

        p3dfft_ftran_r2c(A.data(), A.data()); // Physical to wave
        rtime1 = rtime1 + MPI_Wtime();

        std::transform(A.begin(), A.end(), A.begin(),
                       std::bind1st(std::multiplies<pencil<>::real_type>(),factor));

        ONLYPROC0(LOG4CXX_INFO(logger, "Forward transform results "));

        for (pencil<>::complex_iterator it = A.wave.begin();
             it != A.wave.end();
             ++it) {
            if (abs(*it) > 1e-8) {
                pencil<>::size_type i, j, k;
                A.wave.inverse_global_offset(it - A.wave.begin(), &i, &j, &k);
                LOG4CXX_INFO(logger,
                        boost::format("(%3d, %3d, %3d) = (%12g, %12g) at index %3d")
                        % i % j % k
                        % it->real() % it->imag()
                        % (it-A.wave.begin()) );
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        rtime1 = rtime1 - MPI_Wtime();
        p3dfft_btran_c2r(A.data(), A.data()); // Wave to physical
        rtime1 = rtime1 + MPI_Wtime();

    }

    MPI_Barrier(MPI_COMM_WORLD);
    p3dfft_clean();   // Free work space

    /* Check results */
    // FIXME Error results are fishy, compare P3DFFT's test_sine_inplace.x
    double cdiff = 0.0;
    for (pencil<>::size_type j = 0; j < A.physical.size_y; ++j) {
        for (pencil<>::size_type k = 0; k < A.physical.size_z; ++k) {
            for (pencil<>::size_type i = 0; i < A.physical.size_x; ++i) {
                const double answer = real_data(gridx[i], gridy[j], gridz[k]);
                const double abserr = fabs(A.physical(i,j,k) - answer);
                cdiff               = std::max(cdiff, abserr);
            }
        }
    }

    // Gather error indicator
    double ccdiff = 0.0;
    MPI_Reduce(&cdiff, &ccdiff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    ONLYPROC0(LOG4CXX_INFO(logger,
                           "Maximum difference: " << std::scientific << ccdiff));

    // Gather timing statistics
    double rtime2 = 0.0;
    MPI_Reduce(&rtime1, &rtime2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    ONLYPROC0(LOG4CXX_INFO(logger, "Time per loop: " << rtime2 / ((double)nrep)));
}
