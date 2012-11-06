//--------------------------------------------------------------------------
//
// Copyright (C) 2008, 2009, 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
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
// application_base.hpp: building blocks for basic application logic
// $Id$

#ifndef SUZERAIN_SUPPORT_APPLICATION_BASE_HPP
#define SUZERAIN_SUPPORT_APPLICATION_BASE_HPP

#include <esio/esio.h>

#include <suzerain/common.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/fftw_definition.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/program_options.hpp>
#include <suzerain/restart_definition.hpp>
#include <suzerain/signal_definition.hpp>
#include <suzerain/state.hpp>
#include <suzerain/statistics_definition.hpp>
#include <suzerain/time_definition.hpp>
#include <suzerain/timestepper.hpp>

namespace suzerain {

namespace support {

class field;

/**
 * An abstract driver base class for managing a Suzerain application.
 * Intended for time-varying, three-dimensional problems.
 * Instantiate from within \c main().
 *
 * Signal handling capabilities may misbehave if multiple instances
 * are executing within the same process.
 */
class application_base
{
public:

    typedef interleaved_state<4,complex_t> linear_state_type;

    typedef contiguous_state<4,complex_t> nonlinear_state_type;

    application_base(const std::string &application_synopsis,
                     const std::string &description = "",
                     const std::string &revstr = "");

    /**
     * Initialize everything, including MPI, necessary for the application.
     * Changes to default values, e.g. \ref statsdef, or adding of additional
     * options to \ref options must be completed prior to invoking this
     * method.
     *
     * @param argc Incoming arguments per <code>main(argc, ...)</code>
     * @param argv Incoming arguments per <code>main(..., argv)</code>
     *
     * @return A vector of the positional arguments.
     */
    virtual std::vector<std::string> initialize(int argc, char **argv);

    /**
    * The log4cxx configuration to use, which is built upon but may differ from
    * support::log4cxx_config.
    *
    * <tt>${FOO}</tt> syntax may be used to pick up environment variables in
    * addition to properties See
    * http://logging.apache.org/log4j/1.2/apidocs/org/apache/log4j/PatternLayout.html
    * and
    * http://logging.apache.org/log4j/1.2/apidocs/org/apache/log4j/PropertyConfigurator.html
    */
    virtual std::string log4cxx_config();

    virtual ~application_base();

    std::string revstr;

    std::vector<support::field> fields;

    grid_definition grid;

    fftw_definition fftwdef;

    program_options options;

    boost::shared_ptr<bspline> b;

    boost::shared_ptr<bsplineop> bop; // Collocation operators

    boost::shared_ptr<bsplineop> gop; // Galerkin L2 operators

    boost::shared_ptr<pencil_grid> dgrid;

    boost::shared_ptr<linear_state_type> state_linear;

    boost::shared_ptr<nonlinear_state_type> state_nonlinear;

    /**
     * Routine to save a restart file, generally called via a timecontroller.
     *
     * The restart saves the data in \ref state_linear.
     * The data in \ref state_nonlinear is destroyed by this call.
     */
    virtual bool save_restart(
            real_t t,
            size_t nt);

protected:

private:

    /** Wall time at which MPI_Init completed */
    double wtime_mpi_init;

    /** Wall time elapsed during FFTW planning */
    double wtime_fftw_planning;

#if defined(SUZERAIN_HAVE_P3DFFT) && defined(SUZERAIN_HAVE_UNDERLING)
    /** Use P3DFFT for parallel FFT operations */
    bool use_p3dfft;

    /** Use underling for parallel FFT operations */
    bool use_underling;
#endif

};

} // end namespace support

} // end namespace suzerain

#endif // SUZERAIN_SUPPORT_APPLICATION_BASE_HPP
