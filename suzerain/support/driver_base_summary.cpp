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

#include <esio/esio.h>

#include <suzerain/arsel.hpp>
#include <suzerain/common.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/running_statistics.hpp>
#include <suzerain/summary.hpp>
#include <suzerain/support/definition_arsel.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/program_options.hpp>
#include <suzerain/support/support.hpp>

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

int
driver_base::summary_run(int argc, char **argv)
{
    static const char who[] = "summary";

    // Almost none of the common application/driver infrastructure is used:
    fftwdef.reset();     // No FFTs
    grid.reset();        // Grid taken from input files only
    restartdef.reset();  // No writing restart files
    statsdef.reset();    // ...or writing statistics files
    timedef.reset();     // ...or advancing time

    namespace po = boost::program_options;
    using namespace std;

    // Establish binary-specific options
    support::definition_arsel arspec;
    options.add_definition(arspec);
    bool clobber;
    string datfile;
    string hdffile;
    string tgtfile;
    options.add_options()
        ("clobber",    po::bool_switch(&clobber)->default_value(false),
         "Overwrite any existing HDF5 output files?")
        ("describe,d",
         "Dump all sample descriptions to standard output")
        ("stdout,s",
         "Write results to standard output?")
        ("datfile,f",  po::value(&datfile),
         "Write results to a textual output file")
        ("hdffile,o",  po::value(&hdffile),
         "Write results to an HDF5 output file")
        ("target,t",   po::value(&tgtfile),
         "Project output (per --datfile, --hdffile) onto target grid"
         " and scenario from this file."
         " If omitted, the final positional argument is used.")
        ;

    // Initialize application and then process binary-specific options
    // (henceforth suzerain::support::logging macros become usable)
    const vector<string> restart_files = initialize(argc, argv);
    const bool use_stdout = options.variables().count("stdout");
    const bool use_dat    = options.variables().count("datfile");
    const bool use_hdf5   = options.variables().count("hdffile");
    const bool describe   = options.variables().count("describe");

    // Ensure that we're running in a single processor environment
    if (mpi::comm_size(MPI_COMM_WORLD) > 1) {
        FATAL0(argv[0] << " only intended to run on single rank");
        return EXIT_FAILURE;
    }

    // Die if the non-HDF5 datfile argument ended with '.h5' like an HDF5 file
    // Defends (somewhat) against clobbering potentially valuable results
    if (use_dat && boost::algorithm::ends_with(datfile, ".h5")) {
        FATAL0("Cowardly refusing to output a 'datfile' ending in '.h5'");
        return EXIT_FAILURE;
    }

    // Dump a banner containing one-indexed columns, names, and descriptions
    if (describe) {
        boost::io::ios_all_saver ias(cout);

        const size_t ndxwidth = 1 + static_cast<size_t>(
                floor(log10(static_cast<real_t>(summary::nscalars::total))));
        size_t namewidth = 0;
        for (size_t i = 0; i < summary::nscalars::total; ++i) {
            namewidth = max(namewidth, strlen(summary::name[i]));
        }

        for (size_t i = 0; i < summary::nscalars::total; ++i) {
            cout << "# "
                 << setw(ndxwidth) << right << i
                 << "  "
                 << setw(namewidth) << left << summary::name[i]
                 << "  "
                 << left << summary::description[i]
                 << '\n';
        }
        cout << flush;
    }

    // Processing differs slightly when done file-by-file versus
    // aggregated across multiple files...
    typedef boost::ptr_map<real_t, summary> pool_type;
    if (!use_stdout && !use_dat && !use_hdf5) {

        for (vector<string>::const_iterator it = restart_files.begin();
             it != restart_files.end();
             ++it) {

            // Load data from filename using a clean ESIO handle
            const string& filename = *it;
            shared_ptr<boost::remove_pointer<esio_handle>::type> h(
                    esio_handle_initialize(MPI_COMM_WORLD),
                    esio_handle_finalize);
            esio_file_open(h.get(), filename.c_str(), /*read-only*/0);
            pool_type data(support::load_summary(h.get()));

            // Save quantities to `basename filename .h5`.mean
            static const char suffix[] = ".h5";
            const size_t suffix_len = sizeof(suffix) - 1;
            string outname;
            if (filename.rfind(suffix) == filename.length() - suffix_len) {
                outname = filename.substr(
                        0, filename.length() - suffix_len) + ".mean";
            } else {
                outname = filename + ".mean";
            }
            DEBUG0("Saving nondimensional quantities to " << outname);

            // Write header followed by data values separated by blanks
            ofstream ofs(outname.c_str());
            for (pool_type::const_iterator i = data.begin();
                 i != data.end(); ++i) {
                i->second->write(ofs, data.begin() == i) << '\n' << endl;
            }
            ofs.close();

            // Numerics details reset to avoid carrying grid across files.
            grid.reset();
            b.reset();
            cop.reset();
        }

    } else {

        // Load the target grid and scenario per --target or positionals
        // Notice load_metadata gets everything the base class desires
        {
            shared_ptr<boost::remove_pointer<esio_handle>::type> h(
                    esio_handle_initialize(MPI_COMM_WORLD),
                    esio_handle_finalize);

            if (options.variables().count("target")) {
                esio_file_open(h.get(), tgtfile.c_str(), /*read-only*/0);
            } else if (restart_files.size() > 0) {
                esio_file_open(h.get(),
                               restart_files.back().c_str(), /*read-only*/0);
            } else {
                FATAL0("One or more positional arguments required when "
                       " --target is not supplied");
                return EXIT_FAILURE;
            }

            load_metadata(h.get());
        }

        // A single map of data is stored across all files.  Because the map
        // key is the simulation time, we automatically get a well-ordered,
        // unique set of data across all files.
        pool_type pool;

        for (size_t i = 0; i < restart_files.size(); ++i) {
            const string& filename = restart_files[i];

            // Load data from filename using a clean ESIO handle
            shared_ptr<boost::remove_pointer<esio_handle>::type> h(
                    esio_handle_initialize(MPI_COMM_WORLD),
                    esio_handle_finalize);
            esio_file_open(h.get(), filename.c_str(), /*read-only*/0);
            pool_type data(support::load_summary(h.get(), b));

            // Output status to the user so they don't think we're hung.
            for (pool_type::const_iterator i = data.begin();
                 i != data.end(); ++i) {
                INFO0("Read sample for t = " << i->first
                      << " from " << filename);
            }

            // Transfer data into pool (which erases it from data)
            pool.transfer(data);

            // Warn on any duplicate values which were not transfered
            for (pool_type::const_iterator i = data.begin();
                 i != data.end(); ++i) {
                WARN0("Duplicate sample time "
                      << i->first << " from " << filename << " ignored");
            }

        }

        if (use_stdout) {
            // Write header followed by data values separated by blanks
            for (pool_type::const_iterator i = pool.begin();
                 i != pool.end(); ++i) {
                i->second->write(cout, pool.begin() == i) << '\n' << endl;
            }
        }

        if (use_dat) {
            INFO0("Writing file " << datfile);
            ofstream outf(datfile.c_str());
            for (pool_type::const_iterator i = pool.begin();
                 i != pool.end(); ++i) {
                i->second->write(outf, pool.begin() == i) << '\n' << endl;
            }
        }

        if (use_hdf5) {

            // Create output file and store metadata per base class
            DEBUG0("Creating file " << hdffile);
            shared_ptr<boost::remove_pointer<esio_handle>::type> h(
                    esio_handle_initialize(MPI_COMM_WORLD),
                    esio_handle_finalize);
            esio_file_create(h.get(), hdffile.c_str(), clobber);
            save_metadata(h.get());

            // Determine vector of unique times in the pool and save as "/t"
            std::vector<real_t> t;
            t.reserve(pool.size());
            for (pool_type::const_iterator i = pool.begin();
                 i != pool.end(); ++i) {
                t.push_back(i->first);
            }

            esio_line_establish(h.get(), t.size(), 0, t.size());
            esio_line_write(h.get(), summary::name[summary::offset::t],
                    t.data(), 0, summary::description[summary::offset::t]);

            // Set "/y" to be the one-dimensional vector of collocation points.
            // Strictly speaking unnecessary, but useful shorthand for scripts.
            std::vector<real_t> y;
            y.resize(b->n());
            for (size_t i = 0; i < y.size(); ++i) {
                y[i] = b->collocation_point(i);
            }
            esio_line_establish(h.get(), y.size(), 0, y.size());
            esio_line_write(h.get(), summary::name[summary::offset::y],
                    y.data(), 0, summary::description[summary::offset::y]);

            // Compute the bulk weights and then output those as well.
            suzerain::bsplineop_lu masslu(*cop.get());
            masslu.factor_mass(*cop.get());
            const VectorXr bulk_weights
                    = support::compute_bulk_weights(*b, masslu);
            esio_line_establish(h.get(), bulk_weights.size(),
                                0, bulk_weights.size());
            esio_line_write(h.get(), "bulk_weights", bulk_weights.data(), 0,
                            "Take dot product of these weights against any"
                            " quantity to find the bulk value");

            // Process each component's spatiotemporal trace...
            running_statistics<real_t,1> dtstats;
            esio_plane_establish(h.get(), t.size(), 0, t.size(),
                                          y.size(), 0, y.size());
#           pragma omp parallel default(shared)
            {
                ArrayXXr data;
                std::vector<real_t> eff_N, eff_var, mu, mu_sigma, p, T;
#               pragma omp for schedule(dynamic) lastprivate(dtstats)
                for (size_t c = summary::offset::nongrid;
                    c < summary::nscalars::total;
                    ++c) {

                    INFO0("Processing component " << summary::name[c]);

                    // ...assembling from the pool into contiguous memory
                    data.setConstant(y.size(), t.size(),
                                    numeric_limits<real_t>::quiet_NaN());
                    size_t off = 0;
                    for (pool_type::const_iterator i = pool.begin();
                        i != pool.end(); ++i) {
                        data.col(off++) = i->second->storage.col(c);
                    }

                    // ...running the automatic autocorrelation analysis
                    dtstats = arsel(t, data, arspec,
                                    eff_N, eff_var, mu, mu_sigma, p, T);

                    // ...saving the contiguous spatiotemporal plane to disk
                    // ...and writing arsel results as additional attributes
                    // TODO Parallelize IO after isolating/fixing segfault
#                   pragma omp critical
                    {
                        esio_handle esioh = h.get();
                        const char * const name = summary::name[c];
                        esio_plane_write(esioh, name, data.data(),
                                         data.outerStride(),
                                         data.innerStride());
                        esio_attribute_writev(esioh, name, "eff_N",
                                eff_N.data(), eff_N.size());
                        esio_attribute_writev(esioh, name, "eff_var",
                                eff_var.data(), eff_var.size());
                        esio_attribute_writev(esioh, name, "mu",
                                mu.data(), mu.size());
                        esio_attribute_writev(esioh, name, "mu_sigma",
                                mu_sigma.data(), mu_sigma.size());
                        esio_attribute_writev(esioh, name, "p",
                                p.data(), p.size());
                        esio_attribute_writev(esioh, name, "T",
                                T.data(),        T.size());
                    }
                }
            }

            INFO0("Finished processing aggregated data");
            arspec.save(h.get());

            // Summarize the temporal content of the data iff nontrivial
            if (t.size() == 1) {
                INFO0("Collection contains a single sample at time " << t[0]);
            } else if (t.size() > 1) {
                ostringstream msg;
                msg.precision(static_cast<int>(
                        numeric_limits<real_t>::digits10*0.75));
                msg << "Collection contains " << t.size()
                    << " samples spanning times ["
                    << t[0] << ", " << t[t.size() - 1] << "]";
                INFO0(who, msg.str());
                msg.str("");
                msg.precision(static_cast<int>(
                        numeric_limits<real_t>::digits10*0.50));
                msg << "Min/avg/max/std of time between samples: "
                    << dtstats.min(0) << ", "
                    << dtstats.avg(0) << ", "
                    << dtstats.max(0) << ", "
                    << dtstats.std(0);
                INFO0(who, msg.str());
            }
        }

    }

    return EXIT_SUCCESS;
}

} // end namespace support

} // end namespace suzerain
