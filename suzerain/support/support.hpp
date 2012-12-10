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
// support.hpp: Support logic spanning potentially many applications
// $Id$

#ifndef SUZERAIN_SUPPORT_SUPPORT_HPP
#define SUZERAIN_SUPPORT_SUPPORT_HPP

#include <esio/esio.h>
#ifdef HAVE_UNDERLING
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <underling/error.h>
#include <underling/underling.hpp>
#include <underling/underling_fftw.hpp>
#endif

#include <suzerain/common.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/diffwave.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/grid_specification.hpp>
#include <suzerain/inorder.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/state_fwd.hpp>
#include <suzerain/support/grid_definition.hpp>
#include <suzerain/support/time_definition.hpp>
#include <suzerain/timestepper.hpp>

namespace suzerain {

/**
 * Contains cross-cutting functionality used within various Suzerain
 * applications.
 */
namespace support {

/**
 * Default log4cxx configuration to use when none found in environment.
 * Appends output to the console and to a file.
 */
extern const char log4cxx_config[];

/**
 * A log4cxx configuration similar to channel::log4cxx_config for use by
 * console-only applications where file logging is unnecessary.
 */
extern const char log4cxx_config_console[];

/** Log-and-abort handler for errors originating in the GSL */
void mpi_abort_on_error_handler_gsl(const char * reason,
                                    const char * file,
                                    int line,
                                    int error_code);

/** Log-and-abort handler for errors originating in Suzerain */
void mpi_abort_on_error_handler_suzerain(const char * reason,
                                         const char * file,
                                         int line,
                                         int error_code);

/** Log-and-abort handler for errors originating in ESIO */
void mpi_abort_on_error_handler_esio(const char * reason,
                                     const char * file,
                                     int line,
                                     int error_code);

/** Log-and-abort handler for errors originating in underling */
void mpi_abort_on_error_handler_underling(const char * reason,
                                          const char * file,
                                          int line,
                                          int error_code);

/** Common logic for all error handlers */
void mpi_abort_on_error_handler(const char * reason,
                                const char * file,
                                int line,
                                int error_code,
                                const char * origin,
                                const char * strerror);

/** If wisdom_file is not empty, read wisdom on rank zero and broadcast */
void wisdom_broadcast(const std::string& wisdom_file);

/** If wisdom_file is not empty, gather wisdom to rank zero and write it */
void wisdom_gather(const std::string& wisdom_file);

// FIXME Remove and place logic into application_base
/**
 * Create a B-spline workspace on [left,right] per ndof, k, and htdelta.
 * @return the absolute error in reproducing prescribed abscissae.
 */
real_t create(const int ndof,
              const int k,
              const double left,
              const double right,
              const double htdelta,
              shared_ptr<bspline>& b,
              shared_ptr<bsplineop>& cop);

// FIXME Remove and place logic into application_base
/** Save a \ref bspline workspace in a file. */
void save(const esio_handle h,
          const shared_ptr<bspline>& b,
          const shared_ptr<bsplineop>& cop,
          const shared_ptr<bsplineop>& gop);

// FIXME Remove and place logic into application_base
/**
 * Load a \ref bspline workspace from a file.
 * @return the absolute error in reproducing prescribed abscissae.
 */
real_t load(const esio_handle h,
            shared_ptr<bspline>& b,
            shared_ptr<bsplineop>& cop);

/** Save the current simulation time information */
void save_time(const esio_handle h,
               real_t time);

/** Load the current simulation time information */
void load_time(const esio_handle h,
               real_t &time);

/**
 * Forward declaration to allocate state padded for transformation to/from
 * physical space.  Accounts for parallel decomposition details.  Emphatically
 * \em NOT thread safe.  The caller is responsible for <tt>delete</tt>-ing the
 * returned pointer.  No guarantees are made about the memory contents.
 */
template<class StateType>
StateType* allocate_padded_state(
           const std::size_t howmany_fields,
           const pencil_grid& dgrid);

/**
 * Specialization of allocate_padded_state for contiguous_state.  Emphatically
 * \em NOT thread safe.  The caller is responsible for <tt>delete</tt>-ing the
 * returned pointer.  No guarantees are made about the memory contents.
 */
template<>
contiguous_state<4,complex_t>* allocate_padded_state(
           const std::size_t howmany_fields,
           const pencil_grid& dgrid);

/**
 * Parses "min:max", "min:[defaultmax]", or "[defaultmin]:max" into valmin, \c
 * valmax where \c absmin <= \c valmin <= \c valmax <= \c absmax is enforced
 * with the outer two inequalities being considered a validation failure.
 */
template<typename T>
void parse_range(const std::string& s,
                 T *valmin, T *valmax,
                 const T defaultmin, const T defaultmax,
                 const T absmin, const T absmax,
                 const char *name);

// Prior forward declaration suppresses Intel warnings
template<typename T>
void parse_range(const std::string& s,
                 T *valmin, T *valmax,
                 const T defaultmin, const T defaultmax,
                 const T absmin, const T absmax,
                 const char *name)
{
    assert(absmin <= defaultmin);
    assert(defaultmin <= defaultmax);
    assert(defaultmax <= absmax);

    // Split s on a mandatory colon into whitespace-trimmed s_{min,max}
    const size_t colonpos = s.find_first_of(':');
    if (colonpos == std::string::npos) {
        throw std::invalid_argument(std::string(name)
            + " not in format \"low:high\", \"[low]:high\", or low:[high].");
    }
    std::string s_min(s, 0, colonpos);
    std::string s_max(s, colonpos + 1);
    boost::algorithm::trim(s_min);
    boost::algorithm::trim(s_max);

    // Parse recognized formats into valmin and valmax
    if (s_min.length() == 0 && s_max.length() == 0) {
        throw std::invalid_argument(std::string(name)
            + " not in format \"low:high\", \"[low]:high\", or low:[high].");
    } else if (s_min.length() == 0) {
        *valmin = defaultmin;
        *valmax = exprparse<T>(s_max, name);
    } else if (s_max.length() == 0) {
        *valmin = exprparse<T>(s_min, name);
        *valmax = defaultmax;
    } else {
        *valmin = exprparse<T>(s_min, name);
        *valmax = exprparse<T>(s_max, name);
    }

    // Ensure valmin <= valmax
    if (*valmin > *valmax) std::swap(*valmin, *valmax);

    // Validate range is within [absmin, absmax]
    if (*valmin < absmin || absmax < *valmax) {
        std::ostringstream oss;
        oss << name << " value [" << *valmin << ":" << *valmax
            << "] is outside valid range [" << absmin << ":" << absmax <<  "]";
        throw std::invalid_argument(oss.str());
    }
}

} // end namespace support

} // end namespace suzerain

#endif // SUZERAIN_SUPPORT_SUPPORT_HPP
