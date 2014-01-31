//--------------------------------------------------------------------------
//
// Copyright (C) 2008-2014 Rhys Ulerich
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

#ifndef SUZERAIN_SUPPORT_APPLICATION_BASE_HPP
#define SUZERAIN_SUPPORT_APPLICATION_BASE_HPP

/** @file
 * Building blocks for basic Suzerain-based applications.
 */

#include <suzerain/common.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/state_fwd.hpp>
#include <suzerain/support/esio_fwd.hpp>
#include <suzerain/support/definition_fftw.hpp>
#include <suzerain/support/definition_grid.hpp>
#include <suzerain/support/program_options.hpp>

namespace suzerain {

namespace support {

/**
 * An concrete base class for managing a Suzerain application.
 * Instantiate from within \c main().
 *
 * Signal handling capabilities may misbehave if multiple instances
 * are executing within the same operating system process.
 */
class application_base
{
public:

    /**
     * Constructor providing details necessary for <tt>--help</tt> output.
     *
     * @param application_synopsis Application synopsis for <tt>--help</tt>.
     * @param argument_synopsis    Argument synopsis to display for
     *                             <tt>--help</tt> option.  For example,
     *                             "[FILE]..." or "SOURCE... DIRECTORY".
     * @param description          Extent description of the application to
     *                             be displayed at the bottom of
     *                             <tt>--help</tt>.
     * @param revstr               Version information to be displayed
     *                             when <tt>--version</tt> is used.
     */
    application_base(const std::string &application_synopsis,
                     const std::string &argument_synopsis = "",
                     const std::string &description = "",
                     const std::string &revstr = "");

    virtual ~application_base();

    /**
     * Initialize everything, including MPI, necessary for the application.
     * Changing default values, e.g. \ref fftwdef, or adding additional options
     * to \ref options must be completed \e prior to invoking this method.
     *
     * @param argc Incoming arguments per <code>main(argc, ...)</code>
     * @param argv Incoming arguments per <code>main(..., argv)</code>
     *
     * @return A vector of the positional arguments.
     */
    virtual std::vector<std::string> initialize(
            int argc,
            char **argv);

    /**
     * The log4cxx configuration to use, which is built upon but may differ
     * from support::log4cxx_config.
     *
     * <tt>${FOO}</tt> syntax may be used to pick up environment variables in
     * addition to properties See
     * http://logging.apache.org/log4j/1.2/apidocs/org/apache/log4j/PatternLayout.html
     * and
     * http://logging.apache.org/log4j/1.2/apidocs/org/apache/log4j/PropertyConfigurator.html
     */
    virtual std::string log4cxx_config();

    /** The application's revision string which is reported by \c --version. */
    std::string revstr;

    /** Command line options received by \c initialize. */
    program_options options;

    /** A global grid extent definition. */
    shared_ptr<definition_grid> grid;

    /** Details of how FFTs are performed within the application. */
    shared_ptr<definition_fftw> fftwdef;

    /** A thread-unsafe, mutable B-spline workspace for building operators. */
    shared_ptr<bspline> b;

    /** Collocation-based B-spline operators used by the application. */
    shared_ptr<bsplineop> cop;

    /** Galerkin-based B-spline operators used by the application. */
    shared_ptr<bsplineop> gop;

    /** Decomposition details used for MPI-readiness and parallel FFTs. */
    shared_ptr<pencil_grid> dgrid;

    /**
     * Storage type always kept in Fourier wave space.  Name arises because
     * the state is most closely associated with linear operator #L used by
     * subclasses.
     */
    typedef interleaved_state<4, complex_t> state_linear_type;

    /**
     * Storage type transformable to/from Fourier physical space.  Name arises
     * because the state is most closely associated with nonlinear operator #N
     * within driver subclasses.
     */
    typedef contiguous_state<4, complex_t> state_nonlinear_type;

    /**
     * An ancestor common to \ref state_linear_type and \ref
     * state_nonlinear_type.  A common ancestor is necessary for
     * interoperability of #L and #N.
     */
    typedef multi_array::ref<complex_t, 4> state_common_type;

    /** Linear state storage.  See \ref #state_linear_type. */
    shared_ptr<state_linear_type> state_linear;

    /** Nonlinear state storage.  See \ref #state_nonlinear_type. */
    shared_ptr<state_nonlinear_type> state_nonlinear;

    /**
     * Load a grid and discrete operators from a file or initialize
     * them per #grid.  The following are modified:
     * \li #grid
     * \li #b
     * \li #cop
     * \li #gop
     *
     * @param esioh An ESIO handle pointing to an open restart file.
     *              If \c NULL, no details are loaded from disk
     *              and the discrete operators are formed using
     *              the contents of #grid.
     */
    virtual void load_grid_and_operators(
            const esio_handle esioh);

    /**
     * Save a grid and discrete operators to a file.
     * The following are modified:
     * \li #grid
     * \li #b
     * \li #cop
     * \li #gop
     *
     * @param esioh An ESIO handle pointing to an open restart file.
     */
    virtual void save_grid_and_operators(
            const esio_handle esioh);

    /**
     * Establish the parallel decomposition per #grid.  The following are
     * modified:
     * \li #dgrid
     *
     * @param output_size Display metrics on the global degrees of freedom?
     * @param output_plan Display metrics on the FFTW planning process?
     * @param output_load Display metrics on the load balancing?
     */
    virtual void establish_decomposition(
            const bool output_size = true,
            const bool output_plan = true,
            const bool output_load = true);

    /**
     * Establish the state storage per #dgrid and a supplied number of scalar
     * fields.  If storage with the supplied number of fields is already
     * present, invocation is a NOP.  If storage with a different number of
     * fields is present, it is released and new storage allocated.  Requesting
     * zero scalar fields releases any associated storage.
     *
     * The following may be modified:
     * \li #state_linear
     * \li #state_nonlinear
     *
     * @param linear_nfields    Number of scalar fields of wave-only
     *                          storage to be allocated.
     * @param nonlinear_nfields Number of scalar fields of transformable
     *                          storage to be allocated.
     */
    virtual void establish_state_storage(
            const std::size_t linear_nfields,
            const std::size_t nonlinear_nfields);

    /**
     * Modify IEEE floating point environment per <a
     * href="http://www.gnu.org/software/gsl/manual/html_node/Setting-up-your-IEEE-environment.html">GSL_IEEE_MODE</a>.
     *
     * Should usually be called after "startup" completes as startup processing
     * relies on NaNs.
     */
    virtual void establish_ieee_mode();

    /**
     * Compute and display some discretization quality metrics.  These include
     * the B-spline discrete conservation error and the B-spline mass matrix
     * condition number.
     */
    virtual void log_discretization_quality();

protected:

    /** Wall time at which MPI_Init completed */
    double wtime_mpi_init;

    /** Wall time elapsed during FFTW planning */
    double wtime_fftw_planning;

private:

    /** Use P3DFFT for parallel FFT operations */
    bool use_p3dfft;

    /** Use underling for parallel FFT operations */
    bool use_underling;

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;

};

} // end namespace support

} // end namespace suzerain

#endif // SUZERAIN_SUPPORT_APPLICATION_BASE_HPP
