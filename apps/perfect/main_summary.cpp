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

// FIXME Support loading multiple samples per file
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
#include <suzerain/samples.hpp>
#include <suzerain/support/definition_grid.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/program_options.hpp>
#include <suzerain/support/support.hpp>
#include <suzerain/validation.hpp>

#include "ar.hpp"
#include "driver.hpp"
#include "perfect.hpp"
#include "definition_scenario.hpp"

#pragma warning(disable:1419)

namespace suzerain {

namespace perfect {

/** Summarize mean statistics for one or more restart files. */
struct driver_summary : public driver
{
    driver_summary(const std::string& revstr)
        : driver("Compressible, perfect gas simulation summarization",
                 "RESTART-OR-SAMPLE-HDF5-FILE...",
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
"Options -s, -f, and -o may be specified simultaneously.\n",
                 revstr)
        , who("summary")
        , boplu()
    {
        // Almost none of the common application/driver infrastructure is used:
        fftwdef.reset();     // No FFTs
        grid.reset();        // Grid taken from input files only
        isothermal.reset();  // No monkeying with boundary conditions...
        rad.reset();         // ...or inviscid radial flow parameters
        restartdef.reset();  // ...or writing restart files
        scenario.reset();    // Scenario taken from input files only
        sg.reset();          // No monkeying with slow growth...
        statsdef.reset();    // ...or writing statistics files
        timedef.reset();     // ...or advancing time
    }

    /** Logging requirements are simpler than what superclass provides. */
    virtual std::string log4cxx_config() { return support::log4cxx_config; }

    /** Implementation below in this file */
    int run(int argc, char **argv);

private:

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;

    /** Discrete operators controlling the grid used for output. */
    shared_ptr<bsplineop_lu> boplu;
};

/**
 * Details on the sampled and computed quantities
 */
namespace summary {

/** A Boost.Preprocessor sequence of tuples of grid-related details */
#define SEQ_GRID                                                                                      \
    ((t,            "Simulation time"))                                                               \
    ((y,            "Wall-normal collocation point locations"))                                       \
    ((bulk_weights, "Take dot product of these weights against any quantity to find the bulk value"))

/**
 * A Boost.Preprocessor sequence of tuples of directly sampled quantities.
 * Automatically synchronized with \ref SUZERAIN_SAMPLES in \ref samples.h.
 */
#define SEQ_SAMPLED \
    SUZERAIN_SAMPLES_COMPONENTS(bar_)

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
//       SEQ_GRID SEQ_SAMPLED    SEQ_DERIVED
//                SEQ_SAMPLED_Y  SEQ_DERIVED_Y
//                SEQ_SAMPLED_YY SEQ_DERIVED_YY
// appears to be impossible as the sequence has more than 256 elements.
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
        for (size_t i = 0; i < summary::count; ++i) {  // Headings
            out << std::setw(std::numeric_limits<real_t>::digits10 + 11)
                << summary::name[i];
            if (i < summary::count - 1) out << " ";
        }
        out << std::endl;
    }

    /** Used for formatting output data to match \ref summary::write_names. */
    static const Eigen::IOFormat iofmt(
            Eigen::FullPrecision, 0, "     ", "\n", "    ");

    /**
     * Compute all quantities from namespace \ref summary using the sample
     * collections present in \c filename using the wall-normal discretization
     * from \c filename.
     *
     * @param filename   To be loaded.
     * @param i_scenario If <tt>!i_scenario</tt>,
     *                   populated with the definition_scenario from the file.
     * @param i_grid     Handled identically to <tt>i_scenario</tt>.
     * @param i_b        If <tt>!!i_b</tt> on entry, after computation interpolate
     *                   the results onto the collocation points given by \c i_b.
     *                   Otherwise, perform no additional interpolation and
     *                   update \c i_b with the basis in \c filename.
     * @param i_bop      Handled identically to \c i_b.
     * @param i_boplu    Handled identically to \c i_b.
     *
     * @return A map of quantities keyed on the nondimensional simulation time.
     */
    storage_map_type process(
            const std::string& filename,
            shared_ptr<definition_scenario     >& i_scenario,
            shared_ptr<support::definition_grid>& i_grid,
            shared_ptr<bspline                 >& i_b,
            shared_ptr<bsplineop               >& i_bop,
            shared_ptr<bsplineop_lu            >& i_boplu);

} // namespace summary

#pragma warning(disable:383 1572)

} // namespace perfect

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

} // namespace suzerain

// Provided by main_summary_svnrev.{c,h} to speed recompilation
extern "C" const char revstr[];

/** Instantiate and invoke the application. */
int main(int argc, char **argv)
{
    suzerain::perfect::driver_summary app(revstr);
    return app.run(argc, argv);
}

/** The logic executed by main(). */
int
suzerain::perfect::driver_summary::run(int argc, char **argv)
{
    namespace po = boost::program_options;

    // Establish general, binary-specific options
    bool clobber;
    std::string datfile;
    std::string hdffile;
    options.add_options()
        ("clobber",    po::bool_switch(&clobber)->default_value(false),
         "Overwrite any existing HDF5 output files?")
        ("stdout,s",
         "Write results to standard output?")
        ("datfile,f",  po::value(&datfile),
         "Write results to a textual output file")
        ("hdffile,o",  po::value(&hdffile),
         "Write results to an HDF5 output file")
        ("describe,d",
         "Dump all sample descriptions to standard output")
        ;

    // Establish options related to autoregressive model processing
    po::options_description ar_options(
        "Automatic autocorrelation analysis by AR(p) models");
    std::string ar_criterion = "CIC";
    std::size_t ar_minorder  = 0;
    std::size_t ar_maxorder  = 512;
    bool        ar_absrho    = false;
    double      ar_wlenT0    = 7;
    ar_options.add_options()
        ("criterion",
         po::value(&ar_criterion)->default_value(ar_criterion),
         "Use the specified model selection criterion")
        ("minorder",
         po::value(&ar_minorder)->default_value(ar_minorder)->value_name("MIN"),
         "Consider only models of at least order AR(p=MIN)")
        ("maxorder",
         po::value(&ar_maxorder)->default_value(ar_maxorder)->value_name("MAX"),
         "Consider only models of at most order AR(p=MAX)")
        ("absolute_rho",
         po::bool_switch(&ar_absrho)->default_value(ar_absrho),
         "Use absolute autocorrelation when integrating for T0")
        ("window_T0",
         po::value(&ar_wlenT0)->default_value(ar_wlenT0)->value_name("WLEN"),
         "Integrate for T0 until WLEN times the input length")
        ;
    options.options().add(ar_options);  // Yeeeech

    // Initialize application and then process binary-specific options
    // (henceforth suzerain::support::logging macros become usable)
    const std::vector<std::string> restart_files = initialize(argc, argv);
    const bool use_stdout = options.variables().count("stdout");
    const bool use_dat    = options.variables().count("datfile");
    const bool use_hdf5   = options.variables().count("hdffile");
    const bool describe   = options.variables().count("describe");

    // Ensure that we're running in a single processor environment
    if (mpi::comm_size(MPI_COMM_WORLD) > 1) {
        FATAL(argv[0] << " only intended to run on single rank");
        return EXIT_FAILURE;
    }

    // Die if the non-HDF5 datfile argument ended with '.h5' like an HDF5 file
    // Defends (somewhat) against clobbering potentially valuable results
    if (use_dat && boost::algorithm::ends_with(datfile, ".h5")) {
        FATAL("Cowardly refusing to output a 'datfile' ending in '.h5'");
        return EXIT_FAILURE;
    }

    // Dump a banner containing one-indexed columns, names, and descriptions
    if (describe) {
        boost::io::ios_all_saver ias(std::cout);

        const std::size_t ndxwidth = 1 + static_cast<std::size_t>(
                std::floor(std::log10(static_cast<real_t>(summary::count))));

        std::size_t namewidth = 0;
        for (std::size_t i = 0; i < summary::count; ++i) {
            namewidth = std::max(namewidth, strlen(summary::name[i]));
        }

        for (size_t i = 0; i < summary::count; ++i) {
            std::cout << "# "
                      << std::setw(ndxwidth) << std::right << i
                      << "  "
                      << std::setw(namewidth) << std::left << summary::name[i]
                      << "  "
                      << std::left << summary::desc[i]
                      << '\n';
        }
        std::cout << std::flush;
    }

    // Processing differs slightly when done file-by-file versus
    // aggregated across multiple files...
    if (!use_stdout && !use_dat && !use_hdf5) {

        BOOST_FOREACH(const std::string& filename, restart_files) {

            // Load data from filename
            summary::storage_map_type data = summary::process(
                    filename, scenario, grid, b, cop, boplu);

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
            summary::write_names(ofs);
            BOOST_FOREACH(summary::storage_map_type::value_type i, data) {
                ofs << i->second->format(summary::iofmt) << std::endl
                    << std::endl;
            }
            ofs.close();

            // Numerics details reset to avoid carrying grid across files.
            scenario.reset();
            grid.reset();
            b.reset();
            cop.reset();
            boplu.reset();
        }

    } else {

        // A single map of data is stored across all files.  Because the map
        // key is the simulation time, we automatically get a well-ordered,
        // unique set of data across all files.
        summary::storage_map_type pool;

        // Scenario and grid details preserved across multiple files!
        // The last file on the command line determines the projection target
        // grid because of the BOOST_REVERSE_FOREACH below.  That causes
        // date-sorting input files to behave sensibly.
        BOOST_REVERSE_FOREACH(const std::string& filename, restart_files) {

            if (!scenario) INFO0 ("Output file has scenario per " << filename);
            if (!grid)     DEBUG0("Output file has grid per "     << filename);

            // Load data from filename
            summary::storage_map_type data = summary::process(
                    filename, scenario, grid, b, cop, boplu);

            // Output status to the user so they don't think we're hung.
            BOOST_FOREACH(summary::storage_map_type::value_type i, data) {
                INFO0("Read sample for t = " << i->first
                       << " from " << filename);
            }

            // Transfer data into larger pool (which erases it from data)
            pool.transfer(data);

            // Warn on any duplicate values which were not transfered
            BOOST_FOREACH(summary::storage_map_type::value_type i, data) {
                WARN0("Duplicate sample time "
                      << i->first << " from " << filename << " ignored");
            }

        }

        if (use_stdout) {
            // Write header followed by data values separated by blanks
            summary::write_names(std::cout);
            BOOST_FOREACH(summary::storage_map_type::value_type i, pool) {
                std::cout << i->second->format(summary::iofmt) << std::endl
                          << std::endl;
            }
        }

        if (use_dat) {
            INFO0("Writing file " << datfile);
            std::ofstream outf(datfile.c_str());
            summary::write_names(outf);
            BOOST_FOREACH(summary::storage_map_type::value_type i, pool) {
                outf << i->second->format(summary::iofmt) << std::endl
                     << std::endl;
            }
        }

        if (use_hdf5) {
            // Create a file-specific ESIO handle using RAII
            shared_ptr<boost::remove_pointer<esio_handle>::type> h(
                    esio_handle_initialize(MPI_COMM_WORLD),
                    esio_handle_finalize);

            // Create output file and store metadata
            DEBUG("Creating file " << hdffile);
            esio_file_create(h.get(), hdffile.c_str(), clobber);
            save_metadata(h.get());

            // Determine how many time indices and collocation points we have.
            // We'll build a vector of time values to write after iteration.
            const int Nt = pool.size();
            const int Ny = grid->N.y();
            std::vector<real_t> t;
            t.reserve(Nt);

            // Loop over each entry in pool...
            BOOST_FOREACH(summary::storage_map_type::value_type i, pool) {

                // ...writing every wall-normal pencil of data to file...
                esio_plane_establish(h.get(), Nt, t.size(), 1, Ny, 0, Ny);
                for (std::size_t j = 0; j < summary::count; ++j) {
                    // ...skipping those which do not vary in time...
                    if (    j == summary::t
                         || j == summary::y
                         || j == summary::bulk_weights) {
                        continue;
                    }
                    esio_plane_write(h.get(), summary::name[j],
                                     i->second->col(j).data(), 0, 0,
                                     summary::desc[j]);
                }

                // ...and adding the time value to the running vector of times.
                t.push_back(i->first);

                // Output status to the user so they don't think we're hung.
                INFO0("Wrote sample " << t.size() << " of " << Nt
                      << " for t = " << t.back());

            }

            // Set "/t" to be the one-dimensional vector containing all times.
            esio_line_establish(h.get(), Nt, 0, t.size());
            esio_line_write(h.get(), summary::name[summary::t],
                            t.size() ? &t.front() : NULL,
                            0, summary::desc[summary::t]);

            // Set "/y" to be the one-dimensional vector of collocation points.
            // Strictly speaking unnecessary, but useful shorthand for scripts.
            t.resize(Ny);
            esio_line_establish(h.get(), t.size(), 0, t.size());
            for (int i = 0; i < Ny; ++i) t[i] = b->collocation_point(i);
            esio_line_write(h.get(), summary::name[summary::y], &t.front(),
                            0, summary::desc[summary::y]);

            // (Re-) compute the bulk weights and then output those as well.
            const VectorXr bulk_weights
                    = compute_bulk_weights(grid->L.y(), *b, *boplu);
            esio_line_establish(h.get(), bulk_weights.size(),
                                0, bulk_weights.size());
            esio_line_write(h.get(), summary::name[summary::bulk_weights],
                            bulk_weights.data(), 0,
                            summary::desc[summary::bulk_weights]);
        }
    }

    return EXIT_SUCCESS;
}

suzerain::perfect::summary::storage_map_type
suzerain::perfect::summary::process(
        const std::string& filename,
        shared_ptr<definition_scenario     >& i_scenario,
        shared_ptr<support::definition_grid>& i_grid,
        shared_ptr<bspline                 >& i_b,
        shared_ptr<bsplineop               >& i_bop,
        shared_ptr<bsplineop_lu            >& i_boplu)
{
    storage_map_type retval;

    // Create a file-specific ESIO handle using RAII
    shared_ptr<boost::remove_pointer<esio_handle>::type> h(
            esio_handle_initialize(MPI_COMM_WORLD), esio_handle_finalize);

    DEBUG("Loading file " << filename);
    esio_file_open(h.get(), filename.c_str(), 0 /* read-only */);

    // Load time, scenario, grid, time, and B-spline details from file
    // (cannot use driver_base::load_metadata as method is freestanding)
    real_t time;
    definition_scenario scenario;
    support::definition_grid grid;
    shared_ptr<bspline> b;
    shared_ptr<bsplineop> cop;
    support::load_time(h.get(), time);
    scenario.load(h.get());
    grid.load(h.get());
    support::load_bsplines(h.get(), b, cop);
    assert(b->n() == grid.N.y());

    // Return the scenario and grid to the caller if not already set
    if (!i_scenario) i_scenario.reset(new definition_scenario     (scenario));
    if (!i_grid)     i_grid    .reset(new support::definition_grid(grid    ));

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
    std::auto_ptr<samples> q(new samples(time, b->n()));
    support::load_samples(h.get(), *q);
    if (q->t >= 0) {
        DEBUG0("Successfully loaded samples from " << filename);
    } else {
        WARN0("No valid samples found in " << filename);
        return retval;
    }

    // Convert samples into collocation point values in s
    std::auto_ptr<storage_type> s(new storage_type(b->n(),
                (storage_type::Index) storage_type::ColsAtCompileTime));
    s->fill(std::numeric_limits<real_t>::quiet_NaN());  // ++paranoia

#define ACCUMULATE(quantity, component, offset, description)                          \
    cop->accumulate(0, 1.0, q->quantity().col(offset).data(),                     1,  \
                       0.0, s->col(summary::BOOST_PP_CAT(bar_,component)).data(), 1);

    SUZERAIN_SAMPLES_COMPONENTS_FOR_EACH(ACCUMULATE, SUZERAIN_SAMPLES)

#undef ACCUMULATE

    // Store time and collocation points into s.
    // Not strictly necessary, but very useful for textual output
    // and as a sanity check of any later grid projection.
    s->col(summary::t).fill(q->t);
    for (int i = 0; i < b->n(); ++i)
        s->col(summary::y)[i] = b->collocation_point(i);

    // Free coefficient-related resources
    q.reset();

    // Shorthand for referring to a particular column
#define C(name) s->col(summary::name)

    // Differentiate SAMPLED
    // Uses that bar_rho{,__y,__yy} is the first entry in SAMPLED{,_Y,_YY}
    s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y)>(summary::bar_rho__y)
        = s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED)>(summary::bar_rho);
    boplu->solve(BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y),
            s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y)>(summary::bar_rho__y).data(),
            1, b->n());
    s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_YY)>(summary::bar_rho__yy)
        = s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y)>(summary::bar_rho__y);
    cop->apply(1, BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y), 1.0,
            s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y)>(summary::bar_rho__y).data(),
            1, b->n());
    cop->apply(2, BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y), 1.0,
            s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_YY)>(summary::bar_rho__yy).data(),
            1, b->n());

#undef C

    const real_t bsplines_dist = b->distance_to(*i_b);
    if (bsplines_dist <= suzerain_bspline_distance_distinct) {

        // Compute bulk integration weights
        s->col(summary::bulk_weights)
                = compute_bulk_weights(grid.L.y(), *b, *boplu);

        // Results match target numerics to within acceptable tolerance.
        retval.insert(time, s);

    } else {

        // Results do not match target numerics.
        // Must project onto target collocation points.

        // Convert all results in s to coefficients
        boplu->solve(summary::count, s->data(), 1, b->n());

        // Obtain target collocation points
        suzerain::ArrayXr buf(i_b->n());
        for (int i = 0; i < i_b->n(); ++i) buf[i] = i_b->collocation_point(i);

        // Evaluate coefficients onto the target collocation points
        std::auto_ptr<storage_type> r(new storage_type(i_b->n(),
                    (storage_type::Index) storage_type::ColsAtCompileTime));
        for (std::size_t i = 0; i < summary::count; ++i) {
            b->linear_combination(0, s->col(i).data(),
                                  buf.size(), buf.data(), r->col(i).data());
        }

        // Notice that summary::t, being a constant, and summary::y, being a
        // linear, should have been converted to the target collocation points
        // without more than epsilon-like floating point loss.

        // Compute bulk integration weights (which will not translate directly)
        r->col(summary::bulk_weights)
                = compute_bulk_weights(grid.L.y(), *b, *boplu);

        retval.insert(time, r);

    }

    return retval;
}
