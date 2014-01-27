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
 * Application executing \ref suzerain::perfect::driver_summary::run.
 */

// FIXME Support loading multiple sample collections per file
// FIXME Allow excluding particular time ranges in the output

#include <esio/esio.h>

#include <suzerain/common.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/math.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/pre_gsl.h>
#include <suzerain/rholut.hpp>
#include <suzerain/support/definition_grid.hpp>
#include <suzerain/support/definition_isothermal.hpp>
#include <suzerain/support/definition_time.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/program_options.hpp>
#include <suzerain/support/support.hpp>
#include <suzerain/validation.hpp>

#include "driver.hpp"
#include "perfect.hpp"
#include "quantities.hpp"
#include "definition_scenario.hpp"

#pragma warning(disable:1419)

namespace suzerain {

namespace perfect {

/** Summarize mean statistics for one or more channel restart files. */
struct driver_summary : public driver
{
    driver_summary(const std::string& revstr)
        : driver("Compressible, perfect gas simulation summarization",
                 "RESTART-OR-SAMPLE-HDF5-FILE...",
"Invocable in four distinct ways:\n"
"\t1) perfect_summary                INFILE.h5 ...\n"
"\t2) perfect_summary -s             INFILE.h5 ...\n"
"\t3) perfect_summary -f OUTFILE.dat INFILE.h5 ...\n"
"\t4) perfect_summary -o OUTFILE.h5  INFILE.h5 ...\n"
"\n"
"The first way processes each INFILE.h5 in turn outputting a corresponding\n"
"INFILE.mean containing a whitespace-separated table of means from the first\n"
"sample collection in the file.  The second way (-s) sends the data from all\n "
"sample collections to standard output sorted according to the simulation\n "
"time with a blank line separating adjacent times.  The third way (-f)\n"
"is identical to the second except the output is automatically sent to the\n"
"file named OUTFILE.dat.  The fourth way (-o) outputs a single HDF5 file\n"
"called OUTFILE.h5 containing all sample collections. Options -s, -f,  and\n"
"-o may be specified simultaneously.\n",
                 revstr)
        , who("summary")
    {}

    /** Implementation below in this file */
    int run(int argc, char **argv);

private:

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;
};

/**
 * Details on the sampled and computed quantities
 */
namespace sample {

/** A Boost.Preprocessor sequence of tuples of grid-related details */
#define SEQ_GRID                                                                                      \
    ((t,            "Simulation time"))                                                               \
    ((y,            "Wall-normal collocation point locations"))                                       \
    ((bulk_weights, "Take dot product of these weights against any quantity to find the bulk value"))

/** A Boost.Preprocessor sequence of tuples of directly sampled quantities.  */
#define SEQ_SAMPLED                                                                                                       \
    ((bar_rho,              "Reynolds-averaged density"))                                                                 \
    ((bar_rho_u,            "Reynolds-averaged streamwise momentum"))                                                     \
    ((bar_rho_v,            "Reynolds-averaged wall-normal momentum"))                                                    \
    ((bar_rho_w,            "Reynolds-averaged spanwise momentum"))                                                       \
    ((bar_rho_E,            "Reynolds-averaged total (intrinsic plus kinetic) energy"))                                   \
    ((bar_E,                "Reynolds-averaged total (intrinsic plus kinetic) energy per unit mass"))                     \
    ((bar_T,                "Reynolds-averaged temperature"))                                                             \
    ((bar_a,                "Reynolds-averaged speed of sound"))                                                          \
    ((bar_h0,               "Reynolds-averaged stagnation enthalpy"))                                                     \
    ((bar_H0,               "Reynolds-averaged stagnation enthalpy per unit mass"))                                       \
    ((bar_mu,               "Reynolds-averaged dynamic viscosity"))                                                       \
    ((bar_nu,               "Reynolds-averaged kinematic viscosity"))                                                     \
    ((bar_u,                "Reynolds-averaged streamwise velocity"))                                                     \
    ((bar_v,                "Reynolds-averaged wall-normal velocity"))                                                    \
    ((bar_w,                "Reynolds-averaged spanwise velocity"))                                                       \
    ((bar_symxx_grad_u,     "Symmetric part (x,x)-component of Reynolds-averaged velocity gradient"))                     \
    ((bar_symxy_grad_u,     "Symmetric part (x,y)-component of Reynolds-averaged velocity gradient"))                     \
    ((bar_symxz_grad_u,     "Symmetric part (x,z)-component of Reynolds-averaged velocity gradient"))                     \
    ((bar_symyy_grad_u,     "Symmetric part (y,y)-component of Reynolds-averaged velocity gradient"))                     \
    ((bar_symyz_grad_u,     "Symmetric part (y,z)-component of Reynolds-averaged velocity gradient"))                     \
    ((bar_symzz_grad_u,     "Symmetric part (z,z)-component of Reynolds-averaged velocity gradient"))                     \
    ((bar_symxx_rho_grad_u, "Symmetric part (x,x)-component of Reynolds-averaged density times velocity gradient"))       \
    ((bar_symxy_rho_grad_u, "Symmetric part (x,y)-component of Reynolds-averaged density times velocity gradient"))       \
    ((bar_symxz_rho_grad_u, "Symmetric part (x,z)-component of Reynolds-averaged density times velocity gradient"))       \
    ((bar_symyy_rho_grad_u, "Symmetric part (y,y)-component of Reynolds-averaged density times velocity gradient"))       \
    ((bar_symyz_rho_grad_u, "Symmetric part (y,z)-component of Reynolds-averaged density times velocity gradient"))       \
    ((bar_symzz_rho_grad_u, "Symmetric part (z,z)-component of Reynolds-averaged density times velocity gradient"))       \
    ((bar_gradx_T,          "Reynolds-averaged x-component of temperature gradient"))                                     \
    ((bar_grady_T,          "Reynolds-averaged y-component of temperature gradient"))                                     \
    ((bar_gradz_T,          "Reynolds-averaged z-component of temperature gradient"))                                     \
    ((bar_rho_gradx_T,      "Reynolds-averaged x-component of density times temperature gradient"))                       \
    ((bar_rho_grady_T,      "Reynolds-averaged y-component of density times temperature gradient"))                       \
    ((bar_rho_gradz_T,      "Reynolds-averaged z-component of density times temperature gradient"))                       \
    ((bar_tau_colon_grad_u, "Reynolds-averaged contraction of the viscous stress tensor against the velocity gradient"))  \
    ((bar_tauxx,            "Reynolds-averaged (x,x)-component of the viscous stress tensor"))                            \
    ((bar_tauxy,            "Reynolds-averaged (x,y)-component of the viscous stress tensor"))                            \
    ((bar_tauxz,            "Reynolds-averaged (x,z)-component of the viscous stress tensor"))                            \
    ((bar_tauyy,            "Reynolds-averaged (y,y)-component of the viscous stress tensor"))                            \
    ((bar_tauyz,            "Reynolds-averaged (y,z)-component of the viscous stress tensor"))                            \
    ((bar_tauzz,            "Reynolds-averaged (z,z)-component of the viscous stress tensor"))                            \
    ((bar_tauux,            "Reynolds-averaged x-component of the viscous stress tensor applied to the velocity"))        \
    ((bar_tauuy,            "Reynolds-averaged y-component of the viscous stress tensor applied to the velocity"))        \
    ((bar_tauuz,            "Reynolds-averaged z-component of the viscous stress tensor applied to the velocity"))        \
    ((bar_p_div_u,          "Reynolds-averaged pressure times divergence of the velocity"))                               \
    ((bar_rho_u_u,          "Reynolds-averaged (x,x)-component of the momentum times the velocity"))                      \
    ((bar_rho_u_v,          "Reynolds-averaged (x,y)-component of the momentum times the velocity"))                      \
    ((bar_rho_u_w,          "Reynolds-averaged (x,z)-component of the momentum times the velocity"))                      \
    ((bar_rho_v_v,          "Reynolds-averaged (y,y)-component of the momentum times the velocity"))                      \
    ((bar_rho_v_w,          "Reynolds-averaged (y,z)-component of the momentum times the velocity"))                      \
    ((bar_rho_w_w,          "Reynolds-averaged (z,z)-component of the momentum times the velocity"))                      \
    ((bar_rho_u_u_u,        "Reynolds-averaged (x,x,x)-component of the momentum times the velocity times the velocity")) \
    ((bar_rho_u_u_v,        "Reynolds-averaged (x,x,y)-component of the momentum times the velocity times the velocity")) \
    ((bar_rho_u_u_w,        "Reynolds-averaged (x,x,z)-component of the momentum times the velocity times the velocity")) \
    ((bar_rho_u_v_v,        "Reynolds-averaged (x,y,y)-component of the momentum times the velocity times the velocity")) \
    ((bar_rho_u_v_w,        "Reynolds-averaged (x,y,z)-component of the momentum times the velocity times the velocity")) \
    ((bar_rho_u_w_w,        "Reynolds-averaged (x,z,z)-component of the momentum times the velocity times the velocity")) \
    ((bar_rho_v_v_v,        "Reynolds-averaged (y,y,y)-component of the momentum times the velocity times the velocity")) \
    ((bar_rho_v_v_w,        "Reynolds-averaged (y,y,z)-component of the momentum times the velocity times the velocity")) \
    ((bar_rho_v_w_w,        "Reynolds-averaged (y,z,z)-component of the momentum times the velocity times the velocity")) \
    ((bar_rho_w_w_w,        "Reynolds-averaged (z,z,z)-component of the momentum times the velocity times the velocity")) \
    ((bar_rho_T_u,          "Reynolds-averaged x-component of the temperature times the velocity"))                       \
    ((bar_rho_T_v,          "Reynolds-averaged y-component of the temperature times the velocity"))                       \
    ((bar_rho_T_w,          "Reynolds-averaged z-component of the temperature times the velocity"))                       \
    ((bar_rho_mu,           "Reynolds-averaged dynamic viscosity times the density"))                                     \
    ((bar_mu_Sxx,           "Reynolds-averaged (x,x)-component of the deviatoric portion of the strain rate"))            \
    ((bar_mu_Sxy,           "Reynolds-averaged (x,y)-component of the deviatoric portion of the strain rate"))            \
    ((bar_mu_Sxz,           "Reynolds-averaged (x,z)-component of the deviatoric portion of the strain rate"))            \
    ((bar_mu_Syy,           "Reynolds-averaged (y,y)-component of the deviatoric portion of the strain rate"))            \
    ((bar_mu_Syz,           "Reynolds-averaged (y,z)-component of the deviatoric portion of the strain rate"))            \
    ((bar_mu_Szz,           "Reynolds-averaged (z,z)-component of the deviatoric portion of the strain rate"))            \
    ((bar_mu_div_u,         "Reynolds-averaged dynamic viscosity times divergence of the velocity"))                      \
    ((bar_mu_gradx_T,       "Reynolds-averaged x-component of dynamic viscosity times the temperature gradient"))         \
    ((bar_mu_grady_T,       "Reynolds-averaged y-component of dynamic viscosity times the temperature gradient"))         \
    ((bar_mu_gradz_T,       "Reynolds-averaged z-component of dynamic viscosity times the temperature gradient"))         \
    ((bar_fx,               "Reynolds-averaged x-component of the momentum forcing"))                                     \
    ((bar_fy,               "Reynolds-averaged y-component of the momentum forcing"))                                     \
    ((bar_fz,               "Reynolds-averaged z-component of the momentum forcing"))                                     \
    ((bar_qb,               "Reynolds-averaged volumetric energy forcing"))                                               \
    ((bar_f_dot_u,          "Reynolds-averaged energy contribution due to momentum forcing work"))                        \
    ((bar_Srho,             "Reynolds-averaged mass contributions due to slow growth forcing"))                           \
    ((bar_Srhou,            "Reynolds-averaged streamwise momentum contributions due to slow growth forcing"))            \
    ((bar_Srhov,            "Reynolds-averaged wall-normal momentum contributions due to slow growth forcing"))           \
    ((bar_Srhow,            "Reynolds-averaged spanwise momentum contributions due to slow growth forcing"))              \
    ((bar_SrhoE,            "Reynolds-averaged total energy contributions due to slow growth forcing"))                   \
    ((bar_Srhou_dot_u,      "Reynolds-averaged energy contribution due to slow growth forcing work"))                     \
    ((bar_Crho,             "Reynolds-averaged mass contributions due to various integral constraints"))                  \
    ((bar_Crhou,            "Reynolds-averaged streamwise momentum contributions due to various integral constraints"))   \
    ((bar_Crhov,            "Reynolds-averaged wall-normal momentum contributions due to various integral constraints"))  \
    ((bar_Crhow,            "Reynolds-averaged spanwise momentum contributions due to various integral constraints"))     \
    ((bar_CrhoE,            "Reynolds-averaged total energy contributions due to various integral constraints"))          \
    ((bar_Crhou_dot_u,      "Reynolds-averaged energy contribution due to work by various integral constraints"))

// Helpers for working with sequences of tuples
#define NAME(tuple)  BOOST_PP_TUPLE_ELEM(2,0,tuple)
#define SNAME(tuple) BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(2,0,tuple))
#define DESC(tuple)  BOOST_PP_TUPLE_ELEM(2,1,tuple)

// Prepare sequences of first wall-normal derivatives
#define TRANSFORM_Y(r, data, tuple)                                                 \
    (BOOST_PP_CAT(NAME(tuple),__y),  "Wall-normal first derivative of "DESC(tuple))
#define SEQ_SAMPLED_Y BOOST_PP_SEQ_TRANSFORM(TRANSFORM_Y,,SEQ_SAMPLED)

// Prepare sequences of second wall-normal derivatives
#define TRANSFORM_YY(r, data, tuple)                                                 \
    (BOOST_PP_CAT(NAME(tuple),__yy), "Wall-normal second derivative of "DESC(tuple))
#define SEQ_SAMPLED_YY BOOST_PP_SEQ_TRANSFORM(TRANSFORM_YY,,SEQ_SAMPLED)

// Building a Boost.Preprocessor sequence of all data of interest
//   #define SEQ_ALL
//       SEQ_GRID SEQ_SAMPLED    SEQ_DERIVED    SEQ_LOCALS SEQ_RESIDUAL
//                SEQ_SAMPLED_Y  SEQ_DERIVED_Y
//                SEQ_SAMPLED_YY SEQ_DERIVED_YY
// appears to be possible as the sequence has more than 256 elements.
// Instead, we have to invoke on each component sequence in turn.

    /** Number of scalar quantities processed as a wall-normal function */
    static const std::size_t count = BOOST_PP_SEQ_SIZE(SEQ_GRID)
                                   + BOOST_PP_SEQ_SIZE(SEQ_SAMPLED)
                                   + BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y)
                                   + BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_YY);

#define OP(s, data, tuple) NAME(tuple)
    /** Provides named index constants for each quantity */
    enum index {
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_GRID))      ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_SAMPLED))   ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_RESIDUAL))  ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_SAMPLED_Y)) ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_SAMPLED_YY)),
    };
#undef OP

#define OP(s, data, tuple) SNAME(tuple)
    /** Provides names indexed on \ref index */
    static const char * name[count] = {
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_GRID))      ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_SAMPLED))   ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_SAMPLED_Y)) ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_SAMPLED_YY))
    };
#undef OP

#define OP(s, data, tuple) DESC(tuple)
    /** Provides human-readable descriptions indexed on \ref index */
    static const char * desc[count] = {
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_GRID))      ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_SAMPLED))   ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_SAMPLED_Y)) ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_SAMPLED_YY))
    };
#undef OP

    /** Type used to store all quantities in a single contiguous region */
    typedef Eigen::Array<real_t, Eigen::Dynamic, count> storage_type;

    /**
     * Map type used to manage and sort samples across a time series.
     */
    typedef boost::ptr_map<real_t, storage_type> storage_map_type;

    /** Output names in a manner suitable for columns output by \ref iofmt */
    static void write_names(std::ostream &out)
    {
        for (size_t i = 0; i < sample::count; ++i) {  // Headings
            out << std::setw(std::numeric_limits<real_t>::digits10 + 11)
                << sample::name[i];
            if (i < sample::count - 1) out << " ";
        }
        out << std::endl;
    }

    /** Used for formatting output data to match \ref sample::write_names. */
    static const Eigen::IOFormat iofmt(
            Eigen::FullPrecision, 0, "     ", "\n", "    ");

} // namespace sample

// Provided by main_summary_svnrev.{c,h} to speed recompilation
#pragma warning(push,disable:1419)
extern "C" const char revstr[];
#pragma warning(pop)

#pragma warning(disable:383 1572)

/**
 * Compute the integration weights necessary to compute a bulk quantity from
 * the quantity's value at collocation points using a dot product.
 */
static VectorXr compute_bulk_weights(
        real_t Ly,
        bspline& b,
        bsplineop_lu& boplu)
{
    // Obtain coefficient -> bulk quantity weights
    VectorXr bulkcoeff(b.n());
    b.integration_coefficients(0, bulkcoeff.data());
    bulkcoeff /= Ly;

    // Form M^-1 to map from collocation point values to coefficients
    MatrixXXr mat = MatrixXXr::Identity(b.n(),b.n());
    boplu.solve(b.n(), mat.data(), 1, b.n());

    // Dot the coefficients with each column of M^-1
    VectorXr retval(b.n());
    for (int i = 0; i < b.n(); ++i) {
        retval[i] = bulkcoeff.dot(mat.col(i));
    }

    return retval;
}

} // namespace perfect

} // namespace suzerain

/**
 * Compute all quantities from namespace \ref sample using the sample
 * collections present in \c filename using the wall-normal discretization from
 * \c filename.
 *
 * @param filename   To be loaded.
 * @param i_scenario If <tt>!i_scenario</tt>,
 *                   populated with the definition_scenario from the file.
 * @param i_grid     Handled identically to <tt>i_scenario</tt>.
 * @param i_timedef  Handled identically to <tt>i_scenario</tt>.
 * @param i_b        If <tt>!!i_b</tt> on entry, after computation interpolate
 *                   the results onto the collocation points given by \c i_b.
 *                   Otherwise, perform no additional interpolation and
 *                   update \c i_b with the basis in \c filename.
 * @param i_bop      Handled identically to \c i_b.
 * @param i_boplu    Handled identically to \c i_b.
 *
 * @return A map of quantities keyed on the nondimensional simulation time.
 */
static sample::storage_map_type process(
        const std::string& filename,
        shared_ptr<definition_scenario>& i_scenario,
        shared_ptr<definition_grid    >& i_grid,
        shared_ptr<definition_time    >& i_timedef,
        shared_ptr<bspline            >& i_b,
        shared_ptr<bsplineop          >& i_bop,
        shared_ptr<bsplineop_lu       >& i_boplu);

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);                         // Initialize MPI
    atexit((void (*) ()) MPI_Finalize);             // Finalize MPI at exit
    logging::initialize(MPI_COMM_WORLD,             // Initialize logging
                        support::log4cxx_config_console);

    DEBUG0("Establishing floating point environment from GSL_IEEE_MODE");
    mpi_gsl_ieee_env_setup(suzerain::mpi::comm_rank(MPI_COMM_WORLD));

    // Process incoming arguments
    std::vector<std::string> restart_files;
    std::string datfile;
    std::string hdffile;
    bool use_stdout = false;
    bool use_dat    = false;
    bool use_hdf5   = false;
    bool describe   = false;
    {
        suzerain::support::program_options options(
                "Suzerain-based channel mean quantity computations",
                "RESTART-OR-SAMPLE-HDF5-FILE...",
"Invocable in four distinct ways:\n"
"\t1) perfect_summary                INFILE.h5 ...\n"
"\t2) perfect_summary -s             INFILE.h5 ...\n"
"\t3) perfect_summary -f OUTFILE.dat INFILE.h5 ...\n"
"\t4) perfect_summary -o OUTFILE.h5  INFILE.h5 ...\n"
"\n"
"The first way processes each INFILE.h5 in turn outputting a corresponding\n"
"INFILE.mean containing a comma-separated table of means from the first\n"
"sample collection in the file.  The second way (-s) sends the data from all\n "
"sample collections to standard output sorted according to the simulation\n "
"time with a blank line separating adjacent times.  The third way (-f)\n"
"is identical to the second except the output is automatically sent to the\n"
"file named OUTFILE.dat.  The fourth way (-o) outputs a single HDF5 file\n"
"called OUTFILE.h5 containing all sample collections. Options -s, -f,  and\n"
"-o may be specified simultaneously.\n"
                , revstr);
        options.add_options()
            ("stdout,s",   "Write results to standard output?")
            ("datfile,f",   boost::program_options::value(&datfile),
                           "Write results to a textual output file")
            ("hdffile,o",   boost::program_options::value(&hdffile),
                           "Write results to an HDF5 output file")
            ("describe,d", "Dump all sample details to standard output")
            ;
        restart_files = options.process(argc, argv);
        switch (options.verbose()) {
            case 0:                   break;
            case 1:  DEBUG0_ENABLE(); break;
            default: TRACE0_ENABLE(); break;
        }
        switch (options.verbose_all()) {
            case 0:                   break;
            case 1:  DEBUG_ENABLE();  break;
            default: TRACE_ENABLE();  break;
        }
        use_stdout = options.variables().count("stdout");
        use_dat    = options.variables().count("datfile");
        use_hdf5   = options.variables().count("hdffile");
        describe   = options.variables().count("describe");
    }

    // Ensure that we're running in a single processor environment
    if (suzerain::mpi::comm_size(MPI_COMM_WORLD) > 1) {
        FATAL(argv[0] << " only intended to run on single rank");
        return EXIT_FAILURE;
    }

    // Dump a banner containing one-indexed columns, names, and descriptions
    if (describe) {
        boost::io::ios_all_saver ias(std::cout);

        const std::size_t ndxwidth = 1 + static_cast<std::size_t>(
                std::floor(std::log10(static_cast<real_t>(sample::count))));

        std::size_t namewidth = 0;
        for (std::size_t i = 0; i < sample::count; ++i) {
            namewidth = std::max(namewidth, strlen(sample::name[i]));
        }

        for (size_t i = 0; i < sample::count; ++i) {
            std::cout << "# "
                      << std::setw(ndxwidth) << std::right << i
                      << "  "
                      << std::setw(namewidth) << std::left << sample::name[i]
                      << "  "
                      << std::left << sample::desc[i]
                      << '\n';
        }
        std::cout << std::flush;
    }

    // Scenario and grid details provided to process(...)
    shared_ptr<definition_scenario> scenario;
    shared_ptr<definition_grid    > grid;
    shared_ptr<definition_time    > timedef;
    shared_ptr<bspline            > b;
    shared_ptr<bsplineop          > cop;
    shared_ptr<bsplineop_lu       > boplu;

    // Processing differs slightly when done file-by-file versus
    // aggregated across multiple files...
    if (!use_stdout && !use_dat && !use_hdf5) {

        BOOST_FOREACH(const std::string& filename, restart_files) {

            // Load data from filename
            sample::storage_map_type data = process(
                    filename, scenario, grid, timedef, b, cop, boplu);

            // Save quantities to `basename filename .h5`.mean
            static const char suffix[] = ".h5";
            const std::size_t suffix_len = sizeof(suffix) - 1;
            std::string outname;
            if (filename.rfind(suffix) == filename.length() - suffix_len) {
                outname = filename.substr(
                        0, filename.length() - suffix_len) + ".mean";
            } else {
                outname = filename + ".mean";
            }
            DEBUG("Saving nondimensional quantities to " << outname);

            // Write header followed by data values separated by blanks
            std::ofstream ofs(outname.c_str());
            sample::write_names(ofs);
            BOOST_FOREACH(sample::storage_map_type::value_type i, data) {
                ofs << i->second->format(sample::iofmt) << std::endl
                    << std::endl;
            }
            ofs.close();

            // Numerics details reset to avoid carrying grid across files.
            scenario.reset();
            grid.reset();
            timedef.reset();
            b.reset();
            cop.reset();
            boplu.reset();
        }

    } else {

        // A single map of data is stored across all files.  Because the map
        // key is the simulation time, we automatically get a well-ordered,
        // unique set of data across all files.
        sample::storage_map_type pool;

        // Scenario and grid details preserved across multiple files!
        // The last file on the command line determines the projection target
        // grid because of the BOOST_REVERSE_FOREACH below.  That causes
        // date-sorting input files to behave sensibly.
        BOOST_REVERSE_FOREACH(const std::string& filename, restart_files) {

            if (!scenario) INFO0 ("Output file has scenario per " << filename);
            if (!grid)     DEBUG0("Output file has grid per "     << filename);
            if (!timedef)  DEBUG0("Output file has timedef per "  << filename);

            // Load data from filename
            sample::storage_map_type data = process(
                    filename, scenario, grid, timedef, b, cop, boplu);

            // Output status to the user so they don't thing we're hung.
            BOOST_FOREACH(sample::storage_map_type::value_type i, data) {
                INFO0("Read sample for t = " << i->first
                       << " from " << filename);
            }

            // Transfer data into larger pool (which erases it from data)
            pool.transfer(data);

            // Warn on any duplicate values which were not transfered
            BOOST_FOREACH(sample::storage_map_type::value_type i, data) {
                WARN0("Duplicate sample time "
                      << i->first << " from " << filename << " ignored");
            }

        }

        if (use_stdout) {
            // Write header followed by data values separated by blanks
            sample::write_names(std::cout);
            BOOST_FOREACH(sample::storage_map_type::value_type i, pool) {
                std::cout << i->second->format(sample::iofmt) << std::endl
                          << std::endl;
            }
        }

        if (use_dat) {
            INFO0("Writing file " << datfile);
            std::ofstream outf(datfile.c_str());
            sample::write_names(outf);
            BOOST_FOREACH(sample::storage_map_type::value_type i, pool) {
                outf << i->second->format(sample::iofmt) << std::endl
                     << std::endl;
            }
        }

        if (use_hdf5) {
            // Create a file-specific ESIO handle using RAII
            shared_ptr<boost::remove_pointer<esio_handle>::type> h(
                    esio_handle_initialize(MPI_COMM_WORLD),
                    esio_handle_finalize);

            // Create output file
            DEBUG("Creating file " << hdffile);
            esio_file_create(h.get(), hdffile.c_str(), 1 /* overwrite */);

            // Store the scenario and numerics metadata
            scenario->save(h.get());
            grid->save(h.get());
            shared_ptr<suzerain::bsplineop> gop(new suzerain::bsplineop(
                        *b, 0, SUZERAIN_BSPLINEOP_GALERKIN_L2));
            support::save(h.get(), b, cop, gop);
            gop.reset();
            timedef->save(h.get());

            // Determine how many time indices and collocation points we have.
            // We'll build a vector of time values to write after iteration.
            const int Nt = pool.size();
            const int Ny = grid->N.y();
            std::vector<real_t> t;
            t.reserve(Nt);

            // Loop over each entry in pool...
            BOOST_FOREACH(sample::storage_map_type::value_type i, pool) {

                // ...writing every wall-normal pencil of data to file...
                esio_plane_establish(h.get(), Nt, t.size(), 1, Ny, 0, Ny);
                for (std::size_t j = 0; j < sample::count; ++j) {
                    // ...skipping those which do not vary in time...
                    if (    j == sample::t
                         || j == sample::y
                         || j == sample::bulk_weights) {
                        continue;
                    }
                    esio_plane_write(h.get(), sample::name[j],
                                     i->second->col(j).data(), 0, 0,
                                     sample::desc[j]);
                }

                // ...and adding the time value to the running vector of times.
                t.push_back(i->first);

                // Output status to the user so they don't thing we're hung.
                INFO0("Wrote sample " << t.size() << " of " << Nt
                      << " for t = " << t.back());

            }

            // Set "/t" to be the one-dimensional vector containing all times.
            esio_line_establish(h.get(), Nt, 0, t.size());
            esio_line_write(h.get(), sample::name[sample::t],
                            t.size() ? &t.front() : NULL,
                            0, sample::desc[sample::t]);

            // Set "/y" to be the one-dimensional vector of collocation points.
            // Strictly speaking unnecessary, but useful shorthand for scripts.
            t.resize(Ny);
            esio_line_establish(h.get(), t.size(), 0, t.size());
            for (int i = 0; i < Ny; ++i) t[i] = b->collocation_point(i);
            esio_line_write(h.get(), sample::name[sample::y], &t.front(),
                            0, sample::desc[sample::y]);

            // (Re-) compute the bulk weights and then output those as well.
            const suzerain::VectorXr bulk_weights
                    = compute_bulk_weights(grid->L.y(), *b, *boplu);
            esio_line_establish(h.get(), bulk_weights.size(),
                                0, bulk_weights.size());
            esio_line_write(h.get(), sample::name[sample::bulk_weights],
                            bulk_weights.data(), 0,
                            sample::desc[sample::bulk_weights]);
        }
    }

    return EXIT_SUCCESS;
}

static sample::storage_map_type process(
        const std::string& filename,
        shared_ptr<definition_scenario>& i_scenario,
        shared_ptr<definition_grid    >& i_grid,
        shared_ptr<definition_time    >& i_timedef,
        shared_ptr<bspline            >& i_b,
        shared_ptr<bsplineop          >& i_bop,
        shared_ptr<bsplineop_lu       >& i_boplu)
{
    using sample::storage_type;
    using sample::storage_map_type;

    sample::storage_map_type retval;

    // Create a file-specific ESIO handle using RAII
    shared_ptr<boost::remove_pointer<esio_handle>::type> h(
            esio_handle_initialize(MPI_COMM_WORLD), esio_handle_finalize);

    DEBUG("Loading file " << filename);
    esio_file_open(h.get(), filename.c_str(), 0 /* read-only */);

    // Load time, scenario, grid, timedef, and B-spline details from file.
    // The definition_time defaults are ignored but required as that
    // class lacks a default constructor (by design).
    real_t time;
    definition_scenario scenario;
    definition_grid grid;
    definition_time timedef(/* advance_dt */ 0,
                 /* advance_nt */ 0,
                 /* advance_wt */ 0,
                 /* status_dt  */ 0,
                 /* status_nt  */ 0,
                 /* min_dt     */ 0,
                 /* max_dt     */ 0);
    shared_ptr<bspline> b;
    shared_ptr<bsplineop> cop;
    support::load_time(h.get(), time);
    scenario.load(h.get());
    grid.load(h.get());
    timedef.load(h.get());
    support::load(h.get(), b, cop);
    assert(b->n() == grid.N.y());

    // Return the scenario, grid, and timedef to the caller if not already set
    if (!i_scenario) i_scenario.reset(new definition_scenario(scenario));
    if (!i_grid)     i_grid    .reset(new definition_grid    (grid    ));
    if (!i_timedef)  i_timedef .reset(new definition_time    (timedef ));

    // Compute factorized mass matrix for current operators in use
    shared_ptr<suzerain::bsplineop_lu> boplu
        = suzerain::make_shared<suzerain::bsplineop_lu>(*cop.get());
    boplu->factor_mass(*cop.get());

    // Likewise, use b and friends if i_b was not supplied by the caller
    if (!i_b) {
        i_b     = b;
        i_bop   = cop;
        i_boplu = boplu;
    }

    // Load samples as coefficients
    std::auto_ptr<perfect::quantities> q(new perfect::quantities(time, b->n()));
    q->load(h.get());
    if (q->t >= 0) {
        DEBUG0("Successfully loaded sample collection from " << filename);
    } else {
        WARN0("No valid sample collection found in " << filename);
        return retval;
    }

    // Convert samples into collocation point values in s
    std::auto_ptr<storage_type> s(new storage_type(b->n(),
                (storage_type::Index) storage_type::ColsAtCompileTime));
    s->fill(std::numeric_limits<real_t>::quiet_NaN());  // ++paranoia

#define ACCUMULATE(coeff_name, coeff_col, point_name)                 \
    cop->accumulate(0, 1.0, q->coeff_name().col(coeff_col).data(), 1, \
                    0.0, s->col(sample::point_name).data(),   1)
    ACCUMULATE(rho,              0, bar_rho               );
    ACCUMULATE(rho_u,            0, bar_rho_u             );
    ACCUMULATE(rho_u,            1, bar_rho_v             );
    ACCUMULATE(rho_u,            2, bar_rho_w             );
    ACCUMULATE(rho_E,            0, bar_rho_E             );
    ACCUMULATE(E,                0, bar_E                 );
    ACCUMULATE(T,                0, bar_T                 );
    ACCUMULATE(a,                0, bar_a                 );
    ACCUMULATE(h0,               0, bar_h0                );
    ACCUMULATE(H0,               0, bar_H0                );
    ACCUMULATE(mu,               0, bar_mu                );
    ACCUMULATE(nu,               0, bar_nu                );
    ACCUMULATE(u,                0, bar_u                 );
    ACCUMULATE(u,                1, bar_v                 );
    ACCUMULATE(u,                2, bar_w                 );
    ACCUMULATE(sym_grad_u,       0, bar_symxx_grad_u      );
    ACCUMULATE(sym_grad_u,       1, bar_symxy_grad_u      );
    ACCUMULATE(sym_grad_u,       2, bar_symxz_grad_u      );
    ACCUMULATE(sym_grad_u,       3, bar_symyy_grad_u      );
    ACCUMULATE(sym_grad_u,       4, bar_symyz_grad_u      );
    ACCUMULATE(sym_grad_u,       5, bar_symzz_grad_u      );
    ACCUMULATE(sym_rho_grad_u,   0, bar_symxx_rho_grad_u  );
    ACCUMULATE(sym_rho_grad_u,   1, bar_symxy_rho_grad_u  );
    ACCUMULATE(sym_rho_grad_u,   2, bar_symxz_rho_grad_u  );
    ACCUMULATE(sym_rho_grad_u,   3, bar_symyy_rho_grad_u  );
    ACCUMULATE(sym_rho_grad_u,   4, bar_symyz_rho_grad_u  );
    ACCUMULATE(sym_rho_grad_u,   5, bar_symzz_rho_grad_u  );
    ACCUMULATE(grad_T,           0, bar_gradx_T           );
    ACCUMULATE(grad_T,           1, bar_grady_T           );
    ACCUMULATE(grad_T,           2, bar_gradz_T           );
    ACCUMULATE(rho_grad_T,       0, bar_rho_gradx_T       );
    ACCUMULATE(rho_grad_T,       1, bar_rho_grady_T       );
    ACCUMULATE(rho_grad_T,       2, bar_rho_gradz_T       );
    ACCUMULATE(tau_colon_grad_u, 0, bar_tau_colon_grad_u  );
    ACCUMULATE(tau,              0, bar_tauxx             );
    ACCUMULATE(tau,              1, bar_tauxy             );
    ACCUMULATE(tau,              2, bar_tauxz             );
    ACCUMULATE(tau,              3, bar_tauyy             );
    ACCUMULATE(tau,              4, bar_tauyz             );
    ACCUMULATE(tau,              5, bar_tauzz             );
    ACCUMULATE(tau_u,            0, bar_tauux             );
    ACCUMULATE(tau_u,            1, bar_tauuy             );
    ACCUMULATE(tau_u,            2, bar_tauuz             );
    ACCUMULATE(p_div_u,          0, bar_p_div_u           );
    ACCUMULATE(rho_u_u,          0, bar_rho_u_u           );
    ACCUMULATE(rho_u_u,          1, bar_rho_u_v           );
    ACCUMULATE(rho_u_u,          2, bar_rho_u_w           );
    ACCUMULATE(rho_u_u,          3, bar_rho_v_v           );
    ACCUMULATE(rho_u_u,          4, bar_rho_v_w           );
    ACCUMULATE(rho_u_u,          5, bar_rho_w_w           );
    ACCUMULATE(rho_u_u_u,        0, bar_rho_u_u_u         );
    ACCUMULATE(rho_u_u_u,        1, bar_rho_u_u_v         );
    ACCUMULATE(rho_u_u_u,        2, bar_rho_u_u_w         );
    ACCUMULATE(rho_u_u_u,        3, bar_rho_u_v_v         );
    ACCUMULATE(rho_u_u_u,        4, bar_rho_u_v_w         );
    ACCUMULATE(rho_u_u_u,        5, bar_rho_u_w_w         );
    ACCUMULATE(rho_u_u_u,        6, bar_rho_v_v_v         );
    ACCUMULATE(rho_u_u_u,        7, bar_rho_v_v_w         );
    ACCUMULATE(rho_u_u_u,        8, bar_rho_v_w_w         );
    ACCUMULATE(rho_u_u_u,        9, bar_rho_w_w_w         );
    ACCUMULATE(rho_T_u,          0, bar_rho_T_u           );
    ACCUMULATE(rho_T_u,          1, bar_rho_T_v           );
    ACCUMULATE(rho_T_u,          2, bar_rho_T_w           );
    ACCUMULATE(rho_mu,           0, bar_rho_mu            );
    ACCUMULATE(mu_S,             0, bar_mu_Sxx            );
    ACCUMULATE(mu_S,             1, bar_mu_Sxy            );
    ACCUMULATE(mu_S,             2, bar_mu_Sxz            );
    ACCUMULATE(mu_S,             3, bar_mu_Syy            );
    ACCUMULATE(mu_S,             4, bar_mu_Syz            );
    ACCUMULATE(mu_S,             5, bar_mu_Szz            );
    ACCUMULATE(mu_div_u,         0, bar_mu_div_u          );
    ACCUMULATE(mu_grad_T,        0, bar_mu_gradx_T        );
    ACCUMULATE(mu_grad_T,        1, bar_mu_grady_T        );
    ACCUMULATE(mu_grad_T,        2, bar_mu_gradz_T        );
    ACCUMULATE(f,                0, bar_fx                );
    ACCUMULATE(f,                1, bar_fy                );
    ACCUMULATE(f,                2, bar_fz                );
    ACCUMULATE(qb,               0, bar_qb                );
    ACCUMULATE(f_dot_u,          0, bar_f_dot_u           );
    ACCUMULATE(Srho,             0, bar_Srho              );
    ACCUMULATE(Srhou,            0, bar_Srhou             );
    ACCUMULATE(Srhou,            1, bar_Srhov             );
    ACCUMULATE(Srhou,            2, bar_Srhow             );
    ACCUMULATE(SrhoE,            0, bar_SrhoE             );
    ACCUMULATE(Srhou_dot_u,      0, bar_Srhou_dot_u       );
    ACCUMULATE(Crho,             0, bar_Crho              );
    ACCUMULATE(Crhou,            0, bar_Crhou             );
    ACCUMULATE(Crhou,            1, bar_Crhov             );
    ACCUMULATE(Crhou,            2, bar_Crhow             );
    ACCUMULATE(CrhoE,            0, bar_CrhoE             );
    ACCUMULATE(Crhou_dot_u,      0, bar_Crhou_dot_u       );
#undef ACCUMULATE

    // Store time and collocation points into s.
    // Not strictly necessary, but very useful for textual output
    // and as a sanity check of any later grid projection.
    s->col(sample::t).fill(q->t);
    for (int i = 0; i < b->n(); ++i)
        s->col(sample::y)[i] = b->collocation_point(i);

    // Free coefficient-related resources
    q.reset();

    // Introduce shorthand for constants
    const real_t Ma    = scenario.Ma;
    const real_t Re    = scenario.Re;
    const real_t Pr    = scenario.Pr;
    const real_t gamma = scenario.gamma;

    // Shorthand for referring to a particular column
#define C(name) s->col(sample::name)

    // Shorthand for computing derivatives within a particular __y, __yy
#define D(name)                                      \
    C(name##__y) = C(name);                          \
    boplu->solve(1, C(name##__y).data(), 1, b->n()); \
    C(name##__yy) = C(name##__y);                    \
    cop->apply(1, 1.0, C(name##__y).data(),  1);     \
    cop->apply(2, 1.0, C(name##__yy).data(), 1)

    // Computations following "Sampling logistics" in writeups
    C(tilde_u) = C(bar_rho_u)/C(bar_rho);
    C(tilde_v) = C(bar_rho_v)/C(bar_rho);
    C(tilde_w) = C(bar_rho_w)/C(bar_rho);
    C(tilde_E) = C(bar_rho_E)/C(bar_rho);
    C(tilde_u_u) = C(bar_rho_u_u)/C(bar_rho);
    C(tilde_u_v) = C(bar_rho_u_v)/C(bar_rho);
    C(tilde_u_w) = C(bar_rho_u_w)/C(bar_rho);
    C(tilde_v_v) = C(bar_rho_v_v)/C(bar_rho);
    C(tilde_v_w) = C(bar_rho_v_w)/C(bar_rho);
    C(tilde_w_w) = C(bar_rho_w_w)/C(bar_rho);
    C(tilde_upp_upp) = C(tilde_u_u) - C(tilde_u).square();
    C(tilde_upp_vpp) = C(tilde_u_v) - C(tilde_u)*C(tilde_v);
    C(tilde_upp_wpp) = C(tilde_u_w) - C(tilde_u)*C(tilde_w);
    C(tilde_vpp_vpp) = C(tilde_v_v) - C(tilde_v).square();
    C(tilde_vpp_wpp) = C(tilde_v_w) - C(tilde_v)*C(tilde_w);
    C(tilde_wpp_wpp) = C(tilde_w_w) - C(tilde_w).square();
    C(tilde_k) = (C(tilde_upp_upp) + C(tilde_vpp_vpp) + C(tilde_wpp_wpp)) / 2;
    C(tilde_T) = (gamma*(gamma-1))*(C(tilde_E) - Ma*Ma*(
                  (  C(tilde_u).square()
                   + C(tilde_v).square()
                   + C(tilde_w).square() ) / 2
               +  C(tilde_k)
               ));
    C(tilde_H) = C(tilde_E) + C(tilde_T)/gamma;
    C(bar_p)   = C(bar_rho) * C(tilde_T)/gamma;
    D(tilde_u);  // Form derivatives
    D(tilde_v);  // Form derivatives
    D(tilde_w);  // Form derivatives
    C(bar_tau_colon_grad_upp) = C(bar_tau_colon_grad_u)
                              - C(bar_tauxy)*C(tilde_u__y)
                              - C(bar_tauyy)*C(tilde_v__y)
                              - C(bar_tauyz)*C(tilde_w__y);
    C(tilde_epsilon) = C(bar_tau_colon_grad_upp)/C(bar_rho);
    C(bar_rhop_up) = C(bar_rho_u) - C(bar_rho)*C(bar_u);
    C(bar_rhop_vp) = C(bar_rho_v) - C(bar_rho)*C(bar_v);
    C(bar_rhop_wp) = C(bar_rho_w) - C(bar_rho)*C(bar_w);
    C(bar_upp) = C(bar_u) - C(tilde_u);
    C(bar_vpp) = C(bar_v) - C(tilde_v);
    C(bar_wpp) = C(bar_w) - C(tilde_w);
    C(bar_f_dot_upp) = C(bar_f_dot_u)
                     - C(bar_fx) * C(tilde_u)
                     - C(bar_fy) * C(tilde_v)
                     - C(bar_fz) * C(tilde_w);
    C(bar_Srhou_dot_upp) = C(bar_Srhou_dot_u)
                         - C(bar_Srhou) * C(tilde_u)
                         - C(bar_Srhov) * C(tilde_v)
                         - C(bar_Srhow) * C(tilde_w);
    C(bar_Crhou_dot_upp) = C(bar_Crhou_dot_u)
                         - C(bar_Crhou) * C(tilde_u)
                         - C(bar_Crhov) * C(tilde_v)
                         - C(bar_Crhow) * C(tilde_w);
    C(bar_tauuppx) = C(bar_tauux)
                   - C(bar_tauxx)*C(tilde_u)
                   - C(bar_tauxy)*C(tilde_v)
                   - C(bar_tauxz)*C(tilde_w);
    C(bar_tauuppy) = C(bar_tauuy)
                   - C(bar_tauxy)*C(tilde_u)
                   - C(bar_tauyy)*C(tilde_v)
                   - C(bar_tauyz)*C(tilde_w);
    C(bar_tauuppz) = C(bar_tauuz)
                   - C(bar_tauxz)*C(tilde_u)
                   - C(bar_tauyz)*C(tilde_v)
                   - C(bar_tauzz)*C(tilde_w);
    D(bar_v);  // Form derivatives
    C(bar_pp_div_upp) = C(bar_p_div_u) - C(bar_p)*C(bar_v__y);
    C(tilde_u_u_u) = C(bar_rho_u_u_u)/C(bar_rho);
    C(tilde_u_u_v) = C(bar_rho_u_u_v)/C(bar_rho);
    C(tilde_u_u_w) = C(bar_rho_u_u_w)/C(bar_rho);
    C(tilde_u_v_v) = C(bar_rho_u_v_v)/C(bar_rho);
    C(tilde_u_v_w) = C(bar_rho_u_v_w)/C(bar_rho);
    C(tilde_u_w_w) = C(bar_rho_u_w_w)/C(bar_rho);
    C(tilde_v_v_v) = C(bar_rho_v_v_v)/C(bar_rho);
    C(tilde_v_v_w) = C(bar_rho_v_v_w)/C(bar_rho);
    C(tilde_v_w_w) = C(bar_rho_v_w_w)/C(bar_rho);
    C(tilde_w_w_w) = C(bar_rho_w_w_w)/C(bar_rho);
    C(tilde_u2u) = C(tilde_u_u_u) + C(tilde_u_v_v) + C(tilde_u_w_w);
    C(tilde_u2v) = C(tilde_u_u_v) + C(tilde_v_v_v) + C(tilde_v_w_w);
    C(tilde_u2w) = C(tilde_u_u_w) + C(tilde_v_v_w) + C(tilde_w_w_w);
    C(tilde_upp2upp) = C(tilde_u2u)
                     - (C(tilde_u_u) + C(tilde_v_v) + C(tilde_w_w))*C(tilde_u)
                     - 2*(C(tilde_u_u)*C(tilde_u) + C(tilde_u_v)*C(tilde_v) + C(tilde_u_w)*C(tilde_w))
                     + 2*(C(tilde_u).square() + C(tilde_v).square() + C(tilde_w).square())*C(tilde_u);
    C(tilde_upp2vpp) = C(tilde_u2v)
                     - (C(tilde_u_u) + C(tilde_v_v) + C(tilde_w_w))*C(tilde_v)
                     - 2*(C(tilde_u_v)*C(tilde_u) + C(tilde_v_v)*C(tilde_v) + C(tilde_v_w)*C(tilde_w))
                     + 2*(C(tilde_u).square() + C(tilde_v).square() + C(tilde_w).square())*C(tilde_v);
    C(tilde_upp2wpp) = C(tilde_u2w)
                     - (C(tilde_u_u) + C(tilde_v_v) + C(tilde_w_w))*C(tilde_w)
                     - 2*(C(tilde_u_w)*C(tilde_u) + C(tilde_v_w)*C(tilde_v) + C(tilde_w_w)*C(tilde_w))
                     + 2*(C(tilde_u).square() + C(tilde_v).square() + C(tilde_w).square())*C(tilde_w);
    C(tilde_T_u) = C(bar_rho_T_u)/C(bar_rho);
    C(tilde_T_v) = C(bar_rho_T_v)/C(bar_rho);
    C(tilde_T_w) = C(bar_rho_T_w)/C(bar_rho);
    C(tilde_Tpp_upp) = C(tilde_T_u) - C(tilde_T)*C(tilde_u);
    C(tilde_Tpp_vpp) = C(tilde_T_v) - C(tilde_T)*C(tilde_v);
    C(tilde_Tpp_wpp) = C(tilde_T_w) - C(tilde_T)*C(tilde_w);
    C(tilde_mu) = C(bar_rho_mu)/C(bar_rho);
    C(bar_mupp) = C(bar_mu) - C(tilde_mu);
    C(tilde_nu) = C(bar_mu)/C(bar_rho);
    C(bar_nupp) = C(bar_nu) - C(tilde_nu);
    C(tilde_symxx_grad_u) = C(bar_symxx_rho_grad_u)/C(bar_rho);
    C(tilde_symxy_grad_u) = C(bar_symxy_rho_grad_u)/C(bar_rho);
    C(tilde_symxz_grad_u) = C(bar_symxz_rho_grad_u)/C(bar_rho);
    C(tilde_symyy_grad_u) = C(bar_symyy_rho_grad_u)/C(bar_rho);
    C(tilde_symyz_grad_u) = C(bar_symyz_rho_grad_u)/C(bar_rho);
    C(tilde_symzz_grad_u) = C(bar_symzz_rho_grad_u)/C(bar_rho);
    C(tilde_Sxx) = C(tilde_symxx_grad_u)
                 - (C(tilde_symxx_grad_u) + C(tilde_symyy_grad_u) + C(tilde_symzz_grad_u)) / 3;
    C(tilde_Sxy) = C(tilde_symxy_grad_u);
    C(tilde_Sxz) = C(tilde_symxz_grad_u);
    C(tilde_Syy) = C(tilde_symyy_grad_u)
                 - (C(tilde_symxx_grad_u) + C(tilde_symyy_grad_u) + C(tilde_symzz_grad_u)) / 3;
    C(tilde_Syz) = C(tilde_symyz_grad_u);
    C(tilde_Szz) = C(tilde_symzz_grad_u)
                 - (C(tilde_symxx_grad_u) + C(tilde_symyy_grad_u) + C(tilde_symzz_grad_u)) / 3;
    C(tilde_nupp_Sppxx) = C(bar_mu_Sxx)/C(bar_rho) - C(tilde_nu)*C(tilde_Sxx);
    C(tilde_nupp_Sppxy) = C(bar_mu_Sxy)/C(bar_rho) - C(tilde_nu)*C(tilde_Sxy);
    C(tilde_nupp_Sppxz) = C(bar_mu_Sxz)/C(bar_rho) - C(tilde_nu)*C(tilde_Sxz);
    C(tilde_nupp_Sppyy) = C(bar_mu_Syy)/C(bar_rho) - C(tilde_nu)*C(tilde_Syy);
    C(tilde_nupp_Sppyz) = C(bar_mu_Syz)/C(bar_rho) - C(tilde_nu)*C(tilde_Syz);
    C(tilde_nupp_Sppzz) = C(bar_mu_Szz)/C(bar_rho) - C(tilde_nu)*C(tilde_Szz);
    C(tilde_nupp_div_upp) = C(bar_mu_div_u)/C(bar_rho)
                          - C(tilde_nu)*(C(tilde_symxx_grad_u) + C(tilde_symyy_grad_u) + C(tilde_symzz_grad_u));
    C(tilde_nupp_gradxTpp) = (C(bar_mu_gradx_T) - C(tilde_nu)*C(bar_rho_gradx_T))/C(bar_rho);
    C(tilde_nupp_gradyTpp) = (C(bar_mu_grady_T) - C(tilde_nu)*C(bar_rho_grady_T))/C(bar_rho);
    C(tilde_nupp_gradzTpp) = (C(bar_mu_gradz_T) - C(tilde_nu)*C(bar_rho_gradz_T))/C(bar_rho);

    // Differentiate SAMPLED
    // Uses that bar_rho{,__y,__yy} is the first entry in SAMPLED{,_Y,_YY}
    s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y)>(sample::bar_rho__y)
        = s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED)>(sample::bar_rho);
    boplu->solve(BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y),
            s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y)>(sample::bar_rho__y).data(),
            1, b->n());
    s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_YY)>(sample::bar_rho__yy)
        = s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y)>(sample::bar_rho__y);
    cop->apply(1, BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y), 1.0,
            s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y)>(sample::bar_rho__y).data(),
            1, b->n());
    cop->apply(2, BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y), 1.0,
            s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_YY)>(sample::bar_rho__yy).data(),
            1, b->n());

    // Differentiate DERIVED
    // Uses that tilde_u{,__y,__yy} is the first entry in DERIVED{,_Y,_YY}
    s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_DERIVED_Y)>(sample::tilde_u__y)
        = s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_DERIVED)>(sample::tilde_u);
    boplu->solve(BOOST_PP_SEQ_SIZE(SEQ_DERIVED_Y),
            s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_DERIVED_Y)>(sample::tilde_u__y).data(),
            1, b->n());
    s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_DERIVED_YY)>(sample::tilde_u__yy)
        = s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_DERIVED_Y)>(sample::tilde_u__y);
    cop->apply(1, BOOST_PP_SEQ_SIZE(SEQ_DERIVED_Y), 1.0,
            s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_DERIVED_Y)>(sample::tilde_u__y).data(),
            1, b->n());
    cop->apply(2, BOOST_PP_SEQ_SIZE(SEQ_DERIVED_Y), 1.0,
            s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_DERIVED_YY)>(sample::tilde_u__yy).data(),
            1, b->n());

    // Computations of local quantities (see descriptions for definitions).
    // This must occur after differentiation of SAMPLED and DERIVED.
    // Note the following:
    //
    // 1)  In local_Mat computation, ".abs()" is present to avoid taking the
    // square root of very small, negative tilde_k arising from negative
    // tilde_{upp_upp,vpp_vpp_vpp} in laminar situations due to round off
    // errors (i.e. tilde_u_u - tilde_u**2 ~= -eps).
    //
    // 2)  In local_Re computation, the coefficient scenario.Re arises because
    // (bar_rho_u * L / bar_mu) are already nondimensional.  Multiplying by Re
    // re-incorporates the reference quantities rho_0, u_0, L_0, and mu_0 to
    // cause the nondimensional local_Re to be correctly formed from
    // dimensional quantities.  Ditto for Re in the eddy viscosity. Ditto for
    // Ma in local_Ma and local_Mat.
    C(local_Ma)  = Ma * C(bar_u) / C(bar_a);
    C(local_Mat) = Ma * (std::sqrt(real_t(2))*C(tilde_k).abs().sqrt()) / C(bar_a);
    C(local_Prt) = (C(tilde_upp_vpp) * C(tilde_T__y)) / (C(tilde_Tpp_vpp) * C(tilde_u__y));
    C(local_nut) = - Re * C(tilde_upp_vpp) / C(tilde_u__y);
    C(local_Re)  = Re * C(bar_rho_u) /* L = 1 */ / C(bar_mu);

    // Computation of Favre-averaged density equation residual following writeup
    C(bar_rho__t) =
        // - \nabla\cdot\bar{\rho}\tilde{u}
           - C(bar_rho_v__y)
        // + \overline{\mathscr{S}_{\rho}}
           + C(bar_Srho)
        // + \overline{\mathscr{C}_{\rho}}
           + C(bar_Crho)
        ;

    // Computation of Favre-averaged streamwise momentum residual following writeup
    C(bar_rho_u__t) =
        // - \nabla\cdot(\tilde{u}\otimes\bar{\rho}\tilde{u})
           - C(tilde_v)*C(bar_rho_u__y) - C(bar_rho_u)*C(tilde_v__y)
        // - \frac{1}{\Mach^2}\nabla{}\bar{p}
           - 0
        // + \nabla\cdot\left( \frac{\bar{\tau}}{\Reynolds} \right)
           + C(bar_tauxy__y) / Re
        // - \nabla\cdot\left( \bar{\rho} \widetilde{u''\otimes{}u''} \right)
           - C(bar_rho)*C(tilde_upp_vpp__y) - C(tilde_upp_vpp)*C(bar_rho__y)
        // + \bar{f}
           + C(bar_fx)
        // + \overline{\mathscr{S}_{\rho{}u}}
           + C(bar_Srhou)
        // + \overline{\mathscr{C}_{\rho{}u}}
           + C(bar_Crhou)
        ;

    // Computation of Favre-averaged wall-normal momentum residual following writeup
    C(bar_rho_v__t) =
        // - \nabla\cdot(\tilde{u}\otimes\bar{\rho}\tilde{u})
           - C(tilde_v)*C(bar_rho_v__y) - C(bar_rho_v)*C(tilde_v__y)
        // - \frac{1}{\Mach^2}\nabla{}\bar{p}
           - C(bar_p__y)
        // + \nabla\cdot\left( \frac{\bar{\tau}}{\Reynolds} \right)
           + C(bar_tauyy__y) / Re
        // - \nabla\cdot\left( \bar{\rho} \widetilde{u''\otimes{}u''} \right)
           - C(bar_rho)*C(tilde_vpp_vpp__y) - C(tilde_vpp_vpp)*C(bar_rho__y)
        // + \bar{f}
           + C(bar_fy)
        // + \overline{\mathscr{S}_{\rho{}u}}
           + C(bar_Srhov)
        // + \overline{\mathscr{C}_{\rho{}u}}
           + C(bar_Crhov)
        ;

    // Computation of Favre-averaged spanwise momentum residual following writeup
    C(bar_rho_w__t) =
        // - \nabla\cdot(\tilde{u}\otimes\bar{\rho}\tilde{u})
           - C(tilde_v)*C(bar_rho_w__y) - C(bar_rho_w)*C(tilde_v__y)
        // - \frac{1}{\Mach^2}\nabla{}\bar{p}
           - 0
        // + \nabla\cdot\left( \frac{\bar{\tau}}{\Reynolds} \right)
           + C(bar_tauyz__y) / Re
        // - \nabla\cdot\left( \bar{\rho} \widetilde{u''\otimes{}u''} \right)
           - C(bar_rho)*C(tilde_vpp_wpp__y) - C(tilde_vpp_wpp)*C(bar_rho__y)
        // + \bar{f}
           + C(bar_fz)
        // + \overline{\mathscr{S}_{\rho{}u}}
           + C(bar_Srhow)
        // + \overline{\mathscr{C}_{\rho{}u}}
           + C(bar_Crhow)
        ;

    // Computation of Favre-averaged total energy residual following writeup
    C(bar_rho_E__t) =
    // - \nabla\cdot\bar{\rho}\tilde{H}\tilde{u}
       - C(bar_rho_v)*C(tilde_H__y) - C(tilde_H)*C(bar_rho_v__y)
    // + \Mach^{2} \nabla\cdot\left( \frac{\bar{\tau}}{\Reynolds} \right) \tilde{u}
       + (Ma*Ma/Re)*( C(tilde_u)*C(bar_tauxy__y) + C(bar_tauxy)*C(tilde_u__y)
                    + C(tilde_v)*C(bar_tauyy__y) + C(bar_tauyy)*C(tilde_v__y)
                    + C(tilde_w)*C(bar_tauyz__y) + C(bar_tauyz)*C(tilde_w__y) )
    // - \Mach^{2} \nabla\cdot\left( \bar{\rho} \widetilde{u''\otimes{}u''} \right) \tilde{u}
       - (Ma*Ma)*( C(bar_rho_u)*C(tilde_upp_vpp__y) + C(tilde_upp_vpp)*C(bar_rho_u__y)
                 + C(bar_rho_v)*C(tilde_vpp_vpp__y) + C(tilde_vpp_vpp)*C(bar_rho_v__y)
                 + C(bar_rho_w)*C(tilde_vpp_wpp__y) + C(tilde_vpp_wpp)*C(bar_rho_w__y) )
    // - \frac{1}{2}\Mach^{2} \nabla\cdot\left( \bar{\rho}\widetilde{{u''}^{2}u''} \right)
       - (Ma*Ma/2)*( C(bar_rho)*C(tilde_upp2vpp__y) + C(tilde_upp2vpp)*C(bar_rho__y) )
    // + \Mach^{2} \nabla\cdot\left( \frac{\overline{\tau{}u''}}{\Reynolds} \right)
       + (Ma*Ma/Re)*( C(bar_tauuppy__y) )

    // + \frac{1}{\gamma-1} \nabla\cdot\left(
    //     \frac{ \bar{\mu} \widetilde{\nabla{}T} }{\Reynolds\Prandtl}
    //   \right)
       + (
            C(tilde_nu)*C(bar_rho_grady_T__y) + C(bar_rho_grady_T)*C(tilde_nu__y)
         ) / ((gamma-1)*Re*Pr)
    // + \frac{1}{\gamma-1} \nabla\cdot\left(
    //     \frac{ \bar{\rho} \widetilde{\nu'' \left(\nabla{}T\right)''} }{\Reynolds\Prandtl}
    //   \right)
       + (
            C(bar_rho)*C(tilde_nupp_gradyTpp) + C(tilde_nupp_gradyTpp)*C(bar_rho__y)
         ) / ((gamma-1)*Re*Pr)
    // - \frac{1}{\gamma-1} \nabla\cdot\left( \bar{\rho} \widetilde{T''u''} \right)
       - (
            C(bar_rho)*C(tilde_Tpp_vpp__y) + C(tilde_Tpp_vpp)*C(bar_rho__y)
         ) / (gamma-1)
    // + \Mach^{2} \bar{f}\cdot\tilde{u}
       + Ma*Ma*(C(bar_fx)*C(tilde_u) + C(bar_fy)*C(tilde_v) + C(bar_fz)*C(tilde_w))
    // + \Mach^{2} \overline{f\cdot{}u''}
       + Ma*Ma*C(bar_f_dot_upp)
    // + \bar{q}_b
       + C(bar_qb)
    // + \overline{\mathscr{S}_{\rho{}E}}
       + C(bar_SrhoE)
    // + \overline{\mathscr{C}_{\rho{}E}}
       + C(bar_CrhoE)
        ;

    // Computation of Favre-averaged turbulent kinetic energy residual following writeup
    C(bar_rho_k__t) =
    // - \nabla\cdot\bar{\rho}k\tilde{u}
       - C(bar_rho_v)*C(tilde_k__y) - C(tilde_k)*C(bar_rho_v__y)
    // - \bar{\rho} \widetilde{u''\otimes{}u''} : \nabla\tilde{u}
       - C(bar_rho)*( C(tilde_upp_vpp)*C(tilde_u__y)
                    + C(tilde_vpp_vpp)*C(tilde_v__y)
                    + C(tilde_vpp_wpp)*C(tilde_w__y) )
    // - \frac{\bar{\rho} \epsilon}{\Reynolds}
       - C(bar_rho)*C(tilde_epsilon) / Re
    // - \frac{1}{2}\nabla\cdot\left( \bar{\rho} \widetilde{{u''}^{2}u''} \right)
       - C(bar_rho)*C(tilde_upp2vpp__y)/2 - C(tilde_upp2vpp)*C(bar_rho__y)/2
    // + \nabla\cdot\left( \frac{\overline{\tau{}u''}}{\Reynolds} \right)
       + C(bar_tauuppy__y)/Re
    // + \frac{1}{\Mach^2} \left(
    //       \bar{p}\nabla\cdot\overline{u''}
    //     + \overline{p' \nabla\cdot{}u''}
    //     - \frac{1}{\gamma} \nabla\cdot\bar{\rho} \widetilde{T''u''}
    //   \right)
       + 1/(Ma*Ma)*( C(bar_p)*C(bar_vpp__y)
                   + C(bar_pp_div_upp)
                   - C(bar_rho)*C(tilde_Tpp_vpp__y)/gamma
                   - C(tilde_Tpp_vpp)*C(bar_rho__y)/gamma )
    // + \overline{f\cdot{}u''}
       + C(bar_f_dot_upp)
    // + \overline{\mathscr{S}_{\rho{}u}\cot{}upp}
       + C(bar_Srhou_dot_upp)
    // + \overline{\mathscr{C}_{\rho{}u}\cot{}upp}
       + C(bar_Crhou_dot_upp)
        ;

#undef C

    const real_t bsplines_dist = b->distance_to(*i_b);
    if (bsplines_dist <= suzerain_bspline_distance_distinct) {

        // Compute bulk integration weights
        s->col(sample::bulk_weights)
                = compute_bulk_weights(grid.L.y(), *b, *boplu);

        // Results match target numerics to within acceptable tolerance.
        retval.insert(time, s);

    } else {

        // Results do not match target numerics.
        // Must project onto target collocation points.

        // Convert all results in s to coefficients
        boplu->solve(sample::count, s->data(), 1, b->n());

        // Obtain target collocation points
        suzerain::ArrayXr buf(i_b->n());
        for (int i = 0; i < i_b->n(); ++i) buf[i] = i_b->collocation_point(i);

        // Evaluate coefficients onto the target collocation points
        std::auto_ptr<storage_type> r(new storage_type(i_b->n(),
                    (storage_type::Index) storage_type::ColsAtCompileTime));
        for (std::size_t i = 0; i < sample::count; ++i) {
            b->linear_combination(0, s->col(i).data(),
                                  buf.size(), buf.data(), r->col(i).data());
        }

        // Notice that sample::t, being a constant, and sample::y, being a
        // linear, should have been converted to the target collocation points
        // without more than epsilon-like floating point loss.

        // Compute bulk integration weights (which will not translate directly)
        r->col(sample::bulk_weights)
                = compute_bulk_weights(grid.L.y(), *b, *boplu);

        retval.insert(time, r);

    }

    return retval;
}
