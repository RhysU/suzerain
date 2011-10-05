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
 * channel.hpp: Channel-related functionality spanning binaries
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef CHANNEL_HPP
#define CHANNEL_HPP

#include <Eigen/Core>
#include <esio/esio.h>
#include <suzerain/bspline.hpp>
#include <suzerain/diffwave.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/inorder.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/scenario_definition.hpp>
#include <suzerain/state.hpp>
#include <suzerain/time_definition.hpp>
#include <suzerain/timestepper.hpp>

#include "precision.hpp"
#include "nsctpl_rholut.hpp"

/**
 * Contains cross-cutting functionality used within the channel binaries.
 */
namespace channel {

/** Default log4cxx configuration to use when none found in environment */
extern const char log4cxx_config[];

/** Contains basic details about the scalar state fields employed */
namespace field {

/** Contains state variable indices within state storage */
namespace ndx {

// Anonymous enum to declare our state variable storage indices.
// Update count just below if you modify this enum!
enum {
    rho,  /**< Nondimensional density */
    rhou, /**< Nondimensional streamwise momentum */
    rhov, /**< Nondimensional wall-normal momentum */
    rhow, /**< Nondimensional spanwise momentum */
    rhoe  /**< Nondimensional total energy */
};

} // end namespace ndx;

/** Contains the number of distinct state variables we track */
const std::size_t count = static_cast<std::size_t>(ndx::rhoe) + 1;

/** Field names as stored in restart files */
extern const boost::array<const char *, count> name;

/** Field descriptions as stored in restart files */
extern const boost::array<const char *, count> description;

} // end namespace field

/** Manufactured solution employed throughout the channel code */
typedef nsctpl_rholut::manufactured_solution<real_t> manufactured_solution;

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

/** Store a ScenarioDefinition in a restart file */
void store(const esio_handle h,
           const suzerain::problem::ScenarioDefinition<real_t>& scenario);

/** Load a ScenarioDefinition from a restart file */
void load(const esio_handle h,
          suzerain::problem::ScenarioDefinition<real_t>& scenario);

/** Store a GridDefinition in a restart file */
void store(const esio_handle h,
           const suzerain::problem::GridDefinition& grid,
           const real_t Lx,
           const real_t Lz);

/** Load a GridDefinition from a restart file */
void load(const esio_handle h,
          suzerain::problem::GridDefinition& grid);

/** Store a TimeDefinition in a restart file */
void store(const esio_handle h,
           const suzerain::problem::TimeDefinition<real_t>& timedef);

/** Load a TimeDefinition from a restart file */
void load(const esio_handle h,
          suzerain::problem::TimeDefinition<real_t>& timedef);

/**
 * Store manufactured solution parameters in a restart file.
 * Parameters are only stored when \c msoln evaluates as true.
 */
void store(const esio_handle h,
           const suzerain::problem::ScenarioDefinition<real_t>& scenario,
           const boost::shared_ptr<manufactured_solution> & msoln);

/**
 * Load manufactured solution parameters from a restart file.
 * If the restart file contains active manufactured solution parameters, \c
 * msoln will be modified to contain an appropriate instance.  If it does not,
 * \c msoln will be reset.
 */
void load(const esio_handle h,
          const suzerain::problem::ScenarioDefinition<real_t>& scenario,
          boost::shared_ptr<manufactured_solution>& msoln);

/**
 * Create a B-spline workspace on [left,right] per ndof, k, and htdelta.
 * @return the absolute error in reproducing prescribed abscissae.
 */
real_t create(const int ndof,
              const int k,
              const double left,
              const double right,
              const double htdelta,
              boost::shared_ptr<suzerain::bspline>& b,
              boost::shared_ptr<suzerain::bsplineop>& bop);

/** Store a suzerain::bspline workspace in a restart file */
void store(const esio_handle h,
           const boost::shared_ptr<suzerain::bspline>& b,
           const boost::shared_ptr<suzerain::bsplineop>& bop,
           const boost::shared_ptr<suzerain::bsplineop>& gop);

/**
 * Load a suzerain::bspline workspace from a restart file.
 * @return the absolute error in reproducing prescribed abscissae.
 */
real_t load(const esio_handle h,
            boost::shared_ptr<suzerain::bspline>& b,
            boost::shared_ptr<suzerain::bsplineop>& bop);

/** Store the current simulation time information */
void store_time(const esio_handle h,
                real_t time);

/** Load the current simulation time information */
void load_time(const esio_handle h,
               real_t &time);

/**
 * Forward declaration to allocate state padded for transformation to/from
 * physical space.  Accounts for parallel decomposition details.
 */
template<class StateType>
StateType* allocate_padded_state(
           const std::size_t howmany_fields,
           const suzerain::pencil_grid& dgrid);

/** Specialization of allocate_padded_state for ContiguousState. */
template<>
suzerain::ContiguousState<4,complex_t>* allocate_padded_state(
           const std::size_t howmany_fields,
           const suzerain::pencil_grid& dgrid);

/**
 * Store the current simulation conserved state as expansion coefficients into
 * an open restart file.   Only non-dealiased, conserved state is saved as
 * "wave space" coefficients.  This is the most efficient and flexible way to
 * save state to disk.
 */
void store_coefficients(
        const esio_handle h,
        const suzerain::ContiguousState<4,complex_t> &state,
        suzerain::ContiguousState<4,complex_t> &scratch,
        const suzerain::problem::ScenarioDefinition<real_t>& scenario,
        const suzerain::problem::GridDefinition& grid,
        const suzerain::pencil_grid& dgrid);

/**
 * Store the current simulation primitive state as collocation point values
 * into an open restart file.  Note that <tt>state</tt>'s contents are
 * destroyed.  Collocation point values required only for dealiasing purposes
 * <i>are</i> stored.  This method is less efficient and the restart data less
 * flexible than that produced by store_coefficients().
 */
void store_collocation_values(
        const esio_handle h,
        const suzerain::ContiguousState<4,complex_t> &state,
        suzerain::ContiguousState<4,complex_t>& scratch,
        const suzerain::problem::ScenarioDefinition<real_t>& scenario,
        const suzerain::problem::GridDefinition& grid,
        const suzerain::pencil_grid& dgrid,
        suzerain::bspline& b,
        const suzerain::bsplineop& bop);

/**
 * Load the current simulation state from an open coefficient-based restart
 * file.  Handles the very non-trivial task of adjusting the restart to match
 * the provided \c grid, \c dgrid, \c b, and \c bop.
 */
void load_coefficients(const esio_handle h,
                       suzerain::ContiguousState<4,complex_t> &state,
                       const suzerain::problem::GridDefinition& grid,
                       const suzerain::pencil_grid& dgrid,
                       const suzerain::bspline& b,
                       const suzerain::bsplineop& bop);


/**
 * Load the current simulation state from an open collocation point value
 * restart file.  Cannot handle interpolating onto a different grid.
 */
void load_collocation_values(
        const esio_handle h,
        suzerain::ContiguousState<4,complex_t>& state,
        const suzerain::problem::ScenarioDefinition<real_t>& scenario,
        const suzerain::problem::GridDefinition& grid,
        const suzerain::pencil_grid& dgrid,
        suzerain::bspline& b,
        const suzerain::bsplineop& bop);

/**
 * Interrogate an open restart file and invoke either load_coefficients()
 * or load_collocation_values() as necessary.
 */
void load(const esio_handle h,
          suzerain::ContiguousState<4,complex_t>& state,
          const suzerain::problem::ScenarioDefinition<real_t>& scenario,
          const suzerain::problem::GridDefinition& grid,
          const suzerain::pencil_grid& dgrid,
          suzerain::bspline& b,
          const suzerain::bsplineop& bop);

/** Options definitions for adding random noise to momentum fields */
class NoiseDefinition : public suzerain::problem::IDefinition {

public:

    /** Construct an instance with the given default values */
    explicit NoiseDefinition(real_t fluct_percent = 0,
                             unsigned long fluct_seed = 12345);

    /**
     * Maximum fluctuation magnitude to add as a percentage
     * of centerline streamwise momentum.
     */
    real_t percent;

    /**
     * Fraction of the X direction wavenumbers in [0,1] below
     * which fluctuations will not be added.
     */
    real_t kxfrac_min;

    /**
     * Fraction of the X direction wavenumbers in [0,1] above
     * which fluctuations will not be added.
     */
    real_t kxfrac_max;

    /**
     * Fraction of the Z direction wavenumbers in [0,1] below
     * which fluctuations will not be added.
     */
    real_t kzfrac_min;

    /**
     * Fraction of the Z direction wavenumbers in [0,1] above
     * which fluctuations will not be added.
     */
    real_t kzfrac_max;

    /** RngStream generator seed (see L'Ecuyer et al. 2002) */
    unsigned long seed;

};

/**
 * Add random momentum field perturbations ("noise") according to
 * the provided NoiseDefinition.
 */
void
add_noise(suzerain::ContiguousState<4,complex_t> &state,
          const NoiseDefinition& noisedef,
          const suzerain::problem::ScenarioDefinition<real_t>& scenario,
          const suzerain::problem::GridDefinition& grid,
          const suzerain::pencil_grid& dgrid,
          suzerain::bspline &b,
          const suzerain::bsplineop& bop);

/** Read a complex-valued field via ESIO */
template< typename I >
inline
void complex_field_read(esio_handle h, const char *name, complex_t *field,
                        I cstride = 0, I bstride = 0, I astride = 0)
{
    using boost::numeric_cast;

    esio_field_readv(h, name, reinterpret_cast<real_t *>(field),
                     2*numeric_cast<int>(cstride),
                     2*numeric_cast<int>(bstride),
                     2*numeric_cast<int>(astride),
                     2);
}

/** Read a complex-valued field via ESIO */
inline
void complex_field_read(esio_handle h, const char *name, complex_t *field)
{
    // When no strides are provided, we must specify the stride type.
    return complex_field_read<int>(h, name, field);
}

/** Write a complex-valued field via ESIO */
template< typename I >
inline
void complex_field_write(esio_handle h,
                         const char *name, const complex_t *field,
                         I cstride = 0, I bstride = 0, I astride = 0,
                         const char * comment = 0)
{
    using boost::numeric_cast;

    esio_field_writev(h, name, reinterpret_cast<const real_t *>(field),
                      2*numeric_cast<int>(cstride),
                      2*numeric_cast<int>(bstride),
                      2*numeric_cast<int>(astride),
                      2, comment);
}

/** Write a complex-valued field via ESIO */
inline
void complex_field_write(esio_handle h,
                         const char *name, const complex_t *field)
{
    // When no strides are provided, we must specify the stride type.
    return complex_field_write<int>(h, name, field);
}

/**
 * A template typedef for how to view multiple state fields in physical space,
 * Including a convenient method for constructing such an instance.
 */
template <int NFields>
struct physical_view {

    /**
     * In physical space, we'll employ a view to reshape the 4D row-major (F,
     * Y, Z, X) with contiguous (Y, Z, X) into a 2D (F, Y*Z*X) layout where we
     * know F a priori.  Reducing the dimensionality encourages linear access
     * and eases indexing overhead.
     */
    typedef Eigen::Map<
                    Eigen::Array<real_t, NFields,
                                 Eigen::Dynamic, Eigen::RowMajor>,
                    Eigen::Unaligned, // FIXME Defensive but likely unnecessary
                    Eigen::OuterStride<Eigen::Dynamic>
                > type;

    /**
     * Create a view instance given state storage and sufficient information
     * about the parallel decomposition.
     * */
    static inline type create(
            const suzerain::pencil_grid &dgrid,
            suzerain::ContiguousState<4,complex_t> &state)
    {
        type retval(reinterpret_cast<real_t *>(state.origin()),
                    NFields,                            // F
                    dgrid.local_physical_extent.prod(), // Y*Z*X
                    Eigen::OuterStride<>(  state.strides()[0]
                                         * sizeof(complex_t)/sizeof(real_t)));

        return retval;
    }

};


/** Holds information on the \f$L^2\f$ norm of a scalar field */
struct L2 {
    real_t mean2;
    real_t fluctuating2;
    real_t total2()      const { return mean2 + fluctuating2;    };
    real_t total()       const { return std::sqrt(total2());     };
    real_t mean()        const { return std::sqrt(mean2);        };
    real_t fluctuating() const { return std::sqrt(fluctuating2); };
};

/**
 * Compute information about the \f$L^2\f$ norm of all scalar fields.
 * See writeup/L2.tex for full details.
 */
boost::array<L2,field::count>
field_L2(const suzerain::ContiguousState<4,complex_t> &state,
         const suzerain::problem::ScenarioDefinition<real_t>& scenario,
         const suzerain::problem::GridDefinition& grid,
         const suzerain::pencil_grid& dgrid,
         const suzerain::bsplineop& gop);

/**
 * Accumulate the result of adding \c alpha times the manufactured solution \c
 * msoln times \c beta times the given wave-space state \c swave.  Setting
 * <tt>alpha=1</tt> and <tt>beta=0</tt> may be used to initialize a
 * manufactured solution field.  Setting <tt>alpha=-1</tt> and <tt>beta=1</tt>
 * may be used to compute error against the manufactured solution.
 */
void accumulate_manufactured_solution(
        const real_t alpha,
        const manufactured_solution &msoln,
        const real_t beta,
        suzerain::ContiguousState<4,complex_t> &swave,
        const suzerain::problem::ScenarioDefinition<real_t> &scenario,
        const suzerain::problem::GridDefinition &grid,
        const suzerain::pencil_grid &dgrid,
        suzerain::bspline &b,
        const suzerain::bsplineop &bop,
        const real_t simulation_time);

} // end namespace channel

#endif // CHANNEL_HPP
