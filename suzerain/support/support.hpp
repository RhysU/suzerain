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

#ifndef SUZERAIN_SUPPORT_SUPPORT_HPP
#define SUZERAIN_SUPPORT_SUPPORT_HPP

/** @file
 * Support logic with complex prerequisites spanning multiple applications
 */

#include <suzerain/common.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/extrema.hpp>
#include <suzerain/state_fwd.hpp>
#include <suzerain/support/esio_fwd.hpp>

namespace suzerain {

// Forward declarations
class bspline;
class bsplineop;
class bsplineop_lu;
class pencil_grid;
class samples;
class specification_grid;
class summary;

/**
 * Contains cross-cutting functionality used within various Suzerain
 * applications.
 */
namespace support {

// Forward declarations
class field;

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
void
mpi_abort_on_error_handler_gsl(const char* reason,
                               const char* file,
                               int line,
                               int error_code);

/** Log-and-abort handler for errors originating in Suzerain */
void
mpi_abort_on_error_handler_suzerain(const char* reason,
                                    const char* file,
                                    int line,
                                    int error_code);

/** Log-and-abort handler for errors originating in ESIO */
void
mpi_abort_on_error_handler_esio(const char* reason,
                                const char* file,
                                int line,
                                int error_code);

/** Log-and-abort handler for errors originating in underling */
void
mpi_abort_on_error_handler_underling(const char* reason,
                                     const char* file,
                                     int line,
                                     int error_code);

/** Common logic for all error handlers */
void
mpi_abort_on_error_handler(const char* reason,
                           const char* file,
                           int line,
                           int error_code,
                           const char* origin,
                           const char* strerror);

/**
 * If wisdom_file is not empty, read wisdom on rank zero and broadcast.
 *
 * @return True if wisdom was successfully broadcast.  False otherwise.
 */
bool
wisdom_broadcast(const std::string& wisdom_file);

/**
 * If wisdom_file is not empty, gather wisdom to rank zero and write it.
 *
 * @return True if wisdom was successfully gathered.  False otherwise.
 */
bool
wisdom_gather(const std::string& wisdom_file);

/**
 * Create a B-spline workspace on [left,right] per ndof, k, and htdelta.
 * @return the absolute error in reproducing prescribed abscissae.
 */
real_t
create_bsplines(const int ndof,
                const int k,
                const double left,
                const double right,
                const double htdelta,
                shared_ptr<bspline>& b,
                shared_ptr<bsplineop>& cop,
                const bool verbose = false);

/** Save a \ref bspline workspace and associated operators in a file. */
void
save_bsplines(const esio_handle h,
              bspline&    b,
              const bsplineop&  cop);

/** Save a \ref bspline workspace and associated operators in a file. */
void
save_bsplines(const esio_handle h,
              bspline&    b,
              const bsplineop&  cop,
              const bsplineop&  gop);

/**
 * Load a \ref bspline workspace from a file.
 * @return the absolute error in reproducing prescribed abscissae.
 */
real_t
load_bsplines(const esio_handle      h,
              shared_ptr<bspline>&   b);

/**
 * Load a \ref bspline and \ref bsplineop workspace from a file.
 * @return the absolute error in reproducing prescribed abscissae.
 */
real_t
load_bsplines(const esio_handle      h,
              shared_ptr<bspline>&   b,
              shared_ptr<bsplineop>& cop);

/**
 * Compute the integration weights necessary to compute a bulk quantity from the
 * quantity's value at collocation points using only a dot product.
 */
VectorXr
compute_bulk_weights(bspline& b,
                     const bsplineop_lu& masslu);

/** Save the current simulation time information */
void
save_time(const esio_handle h,
          real_t time);

/** Load the current simulation time information */
void
load_time(const esio_handle h,
          real_t& time);

/**
 * Save sampled quantities in a file using name \c prefix.
 *
 * @return True if all sampled quantities could be saved.  False otherwise.
 */
bool
save_samples(const esio_handle h,
             const samples& s,
             const char* const prefix = "bar_");

/**
 * Load sampled quantities from a file using name \c prefix. Statistics not
 * present in the file are considered to be all NaNs.  Member #t, which is not
 * modified by this routine, is presumably set in some other fashion.
 *
 * @return True if all sampled quantities could be loaded.  False otherwise.
 */
bool
load_samples(const esio_handle h,
             const shared_ptr<bspline> target,
             samples& s,
             const char* const prefix  = "bar_",
             const char* const sizeper = "bar_rho");

/**
 * Save min and max for each field in a file using name \c min_prefix
 * and \c max_prefix, as well as their x and z coordinates with \c minx_prefix,
 * \c maxx_prefix, \c minz_prefix, and \c maxc_prefix. Save as well the
 * global fraction of negative values using \c fneg_prefix.
 *
 * @return True if all sampled quantities could be saved.  False otherwise.
 */
bool
save_extrema(const esio_handle h,
             const std::vector<field>& fields,
             const std::vector<field_extrema_xz>& extrema,
             const specification_grid& grid,
             const char* const min_prefix  = "min_",
             const char* const minx_prefix = "minx_",
             const char* const minz_prefix = "minz_",
             const char* const max_prefix  = "max_",
             const char* const maxx_prefix = "maxx_",
             const char* const maxz_prefix = "maxz_",
             const char* const fneg_prefix = "fneg_");

/**
 * Load min and max for each field in a file using name \c min_prefix
 * and \c max_prefix, as well as their x and z coordinates with \c minx_prefix,
 * \c maxx_prefix, \c minz_prefix, and \c maxc_prefix. Load as well the
 * global fraction of negative values using \c fneg_prefix. Statistics not
 * present in the file are considered to be all NaNs. 
 *
 * @return True if all sampled quantities could be loaded.  False otherwise.
 */
bool
load_extrema(const esio_handle h,
             const std::vector<field>& fields,
             std::vector<field_extrema_xz>& extrema,
             const char* const min_prefix  = "min_",
             const char* const minx_prefix = "minx_",
             const char* const minz_prefix = "minz_",
             const char* const max_prefix  = "max_",
             const char* const maxx_prefix = "maxx_",
             const char* const maxz_prefix = "maxz_",
             const char* const fneg_prefix = "fneg_");

// TODO Separate into prepare_summary(samples&) and load_summary(esio_handle)
// TODO Recognize and load multiple summaries from the file
// TODO At the very least, this default argument stuff is ugly
/**
 * Convert all \ref samples in the file into \ref summary, indexed by time,
 * detailing the data and its derivatives at the collocation points of the \c
 * target.
 *
 * @return A resource-owning container mapping unique simulation times to
 * summarized statistics. See <a
 * href="http://www.boost.org/doc/libs/release/libs/ptr_container/"> Boost
 * Pointer Container</a> for the ownership semantics.
 */
std::auto_ptr<boost::ptr_map<real_t, summary> >
load_summary(const esio_handle h,
             shared_ptr<bspline>& target);

/**
 * Forward declaration to allocate state padded for transformation to/from
 * physical space.  Accounts for parallel decomposition details.  Emphatically
 * \em NOT thread safe.  The caller is responsible for <tt>delete</tt>-ing the
 * returned pointer.  No guarantees are made about the memory contents.
 */
template<class StateType>
StateType*
allocate_padded_state(
    const std::size_t howmany_fields,
    const pencil_grid& dgrid);

/**
 * Specialization of allocate_padded_state for contiguous_state.  Emphatically
 * \em NOT thread safe.  The caller is responsible for <tt>delete</tt>-ing the
 * returned pointer.  No guarantees are made about the memory contents.
 */
template<>
contiguous_state<4, complex_t>*
allocate_padded_state(
    const std::size_t howmany_fields,
    const pencil_grid& dgrid);

} // end namespace support

} // end namespace suzerain

#endif // SUZERAIN_SUPPORT_SUPPORT_HPP
