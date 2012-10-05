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

#ifndef SUPPORT_HPP
#define SUPPORT_HPP

#ifdef HAVE_UNDERLING
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <underling/error.h>
#include <underling/underling.hpp>
#include <underling/underling_fftw.hpp>
#endif

#include <esio/esio.h>
#include <suzerain/common.hpp>
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
           const suzerain::problem::ScenarioDefinition& scenario);

/** Load a ScenarioDefinition from a restart file */
void load(const esio_handle h,
          suzerain::problem::ScenarioDefinition& scenario);

/** Store a GridDefinition in a restart file */
void store(const esio_handle h,
           const suzerain::problem::GridDefinition& grid);

/** Load a GridDefinition from a restart file */
void load(const esio_handle h,
          suzerain::problem::GridDefinition& grid);

/** Store a TimeDefinition in a restart file */
void store(const esio_handle h,
           const suzerain::problem::TimeDefinition& timedef);

/** Load a TimeDefinition from a restart file */
void load(const esio_handle h,
          suzerain::problem::TimeDefinition& timedef);

/**
 * Store manufactured solution parameters in a restart file.
 * Parameters are only stored when \c msoln evaluates as true.
 */
void store(const esio_handle h,
           const suzerain::problem::ScenarioDefinition& scenario,
           const suzerain::problem::GridDefinition& grid,
           const boost::shared_ptr<manufactured_solution> & msoln);

/**
 * Load manufactured solution parameters from a restart file.
 * If the restart file contains active manufactured solution parameters, \c
 * msoln will be modified to contain an appropriate instance.  If it does not,
 * \c msoln will be reset.
 */
void load(const esio_handle h,
          const suzerain::problem::ScenarioDefinition& scenario,
          const suzerain::problem::GridDefinition& grid,
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

/**
 * Compute the "distance" between two B-spline bases.  Distance is "huge" if
 * any of the order, number of degrees of freedom, or number of knots differ.
 * When all those criteria match the distance becomes the maximum absolute
 * difference between the knot vectors.
 */
real_t distance(const suzerain::bspline& a,
                const suzerain::bspline& b);

/**
 * Common constant used to define distinct B-spline bases per
 * bspline_bases_distance() in the presence of floating point error.
 */
extern const real_t bsplines_distinct_distance;

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
 * physical space.  Accounts for parallel decomposition details.  Emphatically
 * \em NOT thread safe.  The caller is responsible for <tt>delete</tt>-ing the
 * returned pointer.  No guarantees are made about the memory contents.
 */
template<class StateType>
StateType* allocate_padded_state(
           const std::size_t howmany_fields,
           const suzerain::pencil_grid& dgrid);

/**
 * Specialization of allocate_padded_state for ContiguousState.  Emphatically
 * \em NOT thread safe.  The caller is responsible for <tt>delete</tt>-ing the
 * returned pointer.  No guarantees are made about the memory contents.
 */
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
        const suzerain::ContiguousState<4,complex_t> &swave,
        const suzerain::problem::ScenarioDefinition& scenario,
        const suzerain::problem::GridDefinition& grid,
        const suzerain::pencil_grid& dgrid);

/**
 * Store the current simulation primitive state as collocation point values
 * into an open restart file.  Note that <tt>state</tt>'s contents are
 * destroyed.  Collocation point values required only for dealiasing purposes
 * <i>are</i> stored but only those points informed by non-dealiased state.
 * This method is less efficient and the restart data less flexible than that
 * produced by store_coefficients().
 */
void store_collocation_values(
        const esio_handle h,
        suzerain::ContiguousState<4,complex_t>& swave,
        const suzerain::problem::ScenarioDefinition& scenario,
        const suzerain::problem::GridDefinition& grid,
        const suzerain::pencil_grid& dgrid,
        suzerain::bspline& b,
        const suzerain::bsplineop& bop);

/**
 * Load the current simulation state from an open coefficient-based restart
 * file.  Handles the non-trivial task of adjusting the restart to match the
 * provided \c grid, \c dgrid, \c b, and \c bop.
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
        const suzerain::problem::ScenarioDefinition& scenario,
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
          const suzerain::problem::ScenarioDefinition& scenario,
          const suzerain::problem::GridDefinition& grid,
          const suzerain::pencil_grid& dgrid,
          suzerain::bspline& b,
          const suzerain::bsplineop& bop);

/**
 * Hold temperature and density constant while changing the Mach number and
 * ratio of specific heats.  On input, \c state should contain total energy
 * fields using \c old_Ma and \c old_gamma.  On output \c state will contain
 * total energy fields using <tt>scenario.Ma</tt> and <tt>scenario.gamma</tt>.
 */
void
adjust_scenario(suzerain::ContiguousState<4,complex_t> &swave,
                const suzerain::problem::ScenarioDefinition& scenario,
                const suzerain::problem::GridDefinition& grid,
                const suzerain::pencil_grid& dgrid,
                suzerain::bspline &b,
                const suzerain::bsplineop& bop,
                const real_t old_Ma,
                const real_t old_gamma);

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
          const suzerain::problem::ScenarioDefinition& scenario,
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
    esio_field_readv(h, name, reinterpret_cast<real_t *>(field),
                     2*boost::numeric_cast<int>(cstride),
                     2*boost::numeric_cast<int>(bstride),
                     2*boost::numeric_cast<int>(astride),
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
    esio_field_writev(h, name, reinterpret_cast<const real_t *>(field),
                      2*boost::numeric_cast<int>(cstride),
                      2*boost::numeric_cast<int>(bstride),
                      2*boost::numeric_cast<int>(astride),
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
 * Including a convenient method for constructing such an instance.  The
 * optional first template parameter may be specified to provide a number of
 * fields known at compile time.
 */
template <int NFields = Eigen::Dynamic>
struct physical_view {

    BOOST_STATIC_ASSERT(NFields == Eigen::Dynamic || NFields >= 0);

    /**
     * In physical space, we'll employ a view to reshape the 4D row-major (F,
     * Y, Z, X) with contiguous (Y, Z, X) into a 2D (F, Y*Z*X) layout where we
     * know F a priori.  Reducing the dimensionality encourages linear access
     * and eases indexing overhead.
     */
    typedef Eigen::Map<
                Eigen::Array<real_t, NFields, Eigen::Dynamic, Eigen::RowMajor>,
                Eigen::Unaligned, // FIXME Defensive but likely unnecessary
                Eigen::OuterStride<Eigen::Dynamic>
            > type;

    /**
     * Create a view instance given state storage and sufficient information
     * about the parallel decomposition.  The default value of \c nfields may
     * only be used when the template parameter \c NFields was not
     * Eigen::Dynamic.
     */
    static inline type create(
            const suzerain::pencil_grid &dgrid,
            suzerain::ContiguousState<4,complex_t> &state,
            const int nfields = NFields)
    {
        if (NFields == Eigen::Dynamic || NFields == nfields) {
            using Eigen::OuterStride;
            type retval(reinterpret_cast<real_t *>(state.origin()),
                        nfields,                            // F
                        dgrid.local_physical_extent.prod(), // Y*Z*X
                        OuterStride<>(  state.strides()[0]
                                     * sizeof(complex_t)/sizeof(real_t)));

            return retval;
        }

        throw std::invalid_argument(
                "NFields, nfields mismatch in physical_view::create");
    }

};

/**
 * Accumulate the result of adding \c alpha times the manufactured solution \c
 * msoln times \c beta times the given wave-space state \c swave.  Setting
 * <tt>alpha=1</tt> and <tt>beta=0</tt> may be used to initialize a
 * manufactured solution field.  Setting <tt>alpha=-1</tt> and <tt>beta=1</tt>
 * may be used to compute error against the manufactured solution.  The
 * manufactured solution lives on \e only the non-dealiased, non-Nyquist modes.
 */
void accumulate_manufactured_solution(
        const real_t alpha,
        const manufactured_solution &msoln,
        const real_t beta,
        suzerain::ContiguousState<4,complex_t> &swave,
        const suzerain::problem::ScenarioDefinition &scenario,
        const suzerain::problem::GridDefinition &grid,
        const suzerain::pencil_grid &dgrid,
        suzerain::bspline &b,
        const suzerain::bsplineop &bop,
        const real_t simulation_time);

/**
 * Encapsulate the mean quantities detailed in the "Sampling logistics" section
 * of <tt>writeups/derivation.tex</tt>.
 *
 * Samples of each quantity are made available through a two-dimensional,
 * column-major arrays.  The row index iterates over wall-normal collocation
 * point locations and the column index iterates over tensor indices.  Scalars
 * have only a single tensor index.  Vector quantities have three indices
 * corresponding to the streamwise x, wall-normal y, and spanwise z directions.
 * Symmetric tensors (for example, \f$\overline{\mu{}S}\f$}) have six entries
 * corresponding to the <tt>xx</tt>, <tt>xy</tt>, <tt>xz</tt>, <tt>yy</tt>,
 * <tt>yz</tt>, and <tt>zz</tt> indices.  Rank one triple products (for example
 * \f$\overline{\rho{}u\otimes{}u\otimes{}u}\f$) have ten entries corresponding
 * to the <tt>xxx</tt>, <tt>xxy</tt>, <tt>xxz</tt>, <tt>xyy</tt>, <tt>xyz</tt>,
 * <tt>xzz</tt>, <tt>yyy</tt>, <tt>yyz</tt>, <tt>yzz</tt>, and <tt>zzz</tt>
 * indices.
 *
 * \internal Many memory overhead and implementation consistency issues have
 * been traded for the headache of reading Boost.Preprocessor-based logic.  So
 * it goes.
 */
class mean
{
public:

/**
 * A Boost.Preprocessor sequence of tuples of quantities computed in wave
 * space.
 */
#define CHANNEL_MEAN_WAVE                                  \
    ((rho,                      1)) /* scalar           */ \
    ((rhou,                     3)) /* vector           */ \
    ((rhoe,                     1)) /* scalar           */

/**
 * A Boost.Preprocessor sequence of tuples of quantities computed in physical
 * space.
 */
#define CHANNEL_MEAN_PHYSICAL                        \
    ((mu,                1))  /* scalar           */ \
    ((nu,                1))  /* scalar           */ \
    ((u,                 3))  /* vector           */ \
    ((sym_rho_grad_u,    6))  /* symmetric tensor */ \
    ((rho_grad_T,        3))  /* vector           */ \
    ((tau_colon_grad_u,  1))  /* scalar           */ \
    ((tau,               6))  /* symmetric tensor */ \
    ((tau_u,             3))  /* vector           */ \
    ((p_div_u,           1))  /* scalar           */ \
    ((rho_u_u,           6))  /* symmetric tensor */ \
    ((rho_u_u_u,        10))  /* symmetric tensor */ \
    ((rho_T_u,           3))  /* vector           */ \
    ((rho_mu,            1))  /* scalar           */ \
    ((mu_S,              6))  /* symmetric tensor */ \
    ((mu_div_u,          1))  /* scalar           */ \
    ((mu_grad_T,         3))  /* vector           */

/**
 * A Boost.Preprocessor sequence of tuples of quantities computed
 * through implicit forcing.
 */
#define CHANNEL_MEAN_IMPLICIT                               \
    ((f,                        3))  /* vector           */ \
    ((qb,                       1))  /* scalar           */ \
    ((f_dot_u,                  1))  /* scalar           */

/** A Boost.Preprocessor sequence of tuples of all sampled quantities. */
#define CHANNEL_MEAN \
    CHANNEL_MEAN_WAVE CHANNEL_MEAN_PHYSICAL CHANNEL_MEAN_IMPLICIT

    /* Compile-time totals of the number of scalars sampled at each point */
    struct nscalars { enum {
#define EXTRACT(r, data, tuple) BOOST_PP_TUPLE_ELEM(2, 1, tuple)
#define SUM(s, state, x) BOOST_PP_ADD(state, x)

        wave = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0,
                BOOST_PP_SEQ_TRANSFORM(EXTRACT,,CHANNEL_MEAN_WAVE)),

        physical = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0,
                BOOST_PP_SEQ_TRANSFORM(EXTRACT,,CHANNEL_MEAN_PHYSICAL)),

        implicit = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0,
                BOOST_PP_SEQ_TRANSFORM(EXTRACT,,CHANNEL_MEAN_IMPLICIT)),

        total = BOOST_PP_SEQ_FOLD_LEFT(SUM, 0,
                BOOST_PP_SEQ_TRANSFORM(EXTRACT,,CHANNEL_MEAN))

#undef EXTRACT
#undef SUM
    }; };

    /** Simulation time when mean quantities were obtained */
    real_t t;

    /** Type of the contiguous storage used to house all scalars */
    typedef Eigen::Array<real_t, Eigen::Dynamic, nscalars::total> storage_type;

    /** Contiguous storage used to house all means */
    storage_type storage;

    /**
     * Constructor setting <tt>this->t = NaN</tt>.
     * Caller will need to resize <tt>this->storage</tt> prior to use.
     */
    mean()
        : t(std::numeric_limits<real_t>::quiet_NaN())
    {}

    /**
     * Constructor setting <tt>this->t = t</tt>.
     * Caller will need to resize <tt>this->storage</tt> prior to use.
     */
    explicit mean(real_t t)
        : t(t)
    {}

    /**
     * Constructor setting <tt>this->t = t</tt> and preparing a zero-filled \c
     * storage containing \c Ny rows.
     */
    mean(real_t t, storage_type::Index Ny)
        : t(t),
          storage(storage_type::Zero(Ny, storage_type::ColsAtCompileTime))
    {}

#define OP(r, data, tuple)                                              \
    BOOST_PP_TUPLE_ELEM(2, 0, tuple) = BOOST_PP_TUPLE_ELEM(2, 1, tuple)

    /** Compile-time offsets for each quantity within \c storage */
    struct start { enum {
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(
                OP,,SUZERAIN_SHIFTED_SUM(CHANNEL_MEAN)))
    }; };

    /** Compile-time sizes for each quantity within \c storage */
    struct size { enum {
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP,,CHANNEL_MEAN))
    }; };

#undef OP

    // Declare a named, mutable "view" into storage for each quantity
#define DECLARE(r, data, tuple)                                               \
    storage_type::NColsBlockXpr<size::BOOST_PP_TUPLE_ELEM(2, 0, tuple)>::Type \
    BOOST_PP_TUPLE_ELEM(2, 0, tuple)()                                        \
    {                                                                         \
        return storage.middleCols<size::BOOST_PP_TUPLE_ELEM(2, 0, tuple)>(    \
                start::BOOST_PP_TUPLE_ELEM(2, 0, tuple));                     \
    }
    BOOST_PP_SEQ_FOR_EACH(DECLARE,,CHANNEL_MEAN)
#undef DECLARE

    // Declare a named, immutable "view" into storage for each quantity
#define DECLARE(r, data, tuple)                                                    \
    storage_type::ConstNColsBlockXpr<size::BOOST_PP_TUPLE_ELEM(2, 0, tuple)>::Type \
    BOOST_PP_TUPLE_ELEM(2, 0, tuple)() const                                       \
    {                                                                              \
        return storage.middleCols<size::BOOST_PP_TUPLE_ELEM(2, 0, tuple)>(         \
                start::BOOST_PP_TUPLE_ELEM(2, 0, tuple));                          \
    }
    BOOST_PP_SEQ_FOR_EACH(DECLARE,,CHANNEL_MEAN)
#undef DECLARE

    /**
     * A foreach operation iterating over all mutable quantities in \c storage.
     * The functor is invoked as <tt>f(std::string("foo",
     * storage_type::NColsBlockXpr<size::foo>::Type))</tt> for a quantity named
     * "foo".  See Eigen's "Writing Functions Taking Eigen Types as Parameters"
     * for suggestions on how to write a functor, especially the \c const_cast
     * hack details therein.  See <tt>boost::ref</tt> for how to use a stateful
     * functor.
     */
    template <typename BinaryFunction>
    void foreach(BinaryFunction f) {
#define INVOKE(r, data, tuple) \
        f(::std::string(BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(2, 0, tuple))), \
          this->BOOST_PP_TUPLE_ELEM(2, 0, tuple)());
        BOOST_PP_SEQ_FOR_EACH(INVOKE,,CHANNEL_MEAN)
    }
#undef INVOKE

    /**
     * A foreach operation iterating over all immutable quantities in \c
     * storage.  The functor is invoked as <tt>f(std::string("foo",
     * storage_type::NColsBlockXpr<size::foo>::Type))</tt> for a quantity named
     * "foo".  See Eigen's "Writing Functions Taking Eigen Types as Parameters"
     * for suggestions on how to write a functor.  See <tt>boost::ref</tt> for
     * how to use a stateful functor.
     */
    template <typename BinaryFunction>
    void foreach(BinaryFunction f) const {
#define INVOKE(r, data, tuple) \
        f(::std::string(BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(2, 0, tuple))), \
          this->BOOST_PP_TUPLE_ELEM(2, 0, tuple)());
        BOOST_PP_SEQ_FOR_EACH(INVOKE,,CHANNEL_MEAN)
    }
#undef INVOKE
};

/**
 * Using the provided state, sample the mean quantities declared in \ref mean
 * with the notable exceptions of \f$\bar{f}\f$, \f$\overline{\rho{}q_b}\f$,
 * and \f$\overline{f\cdot{}u}\f$.  This is an expensive, collective method
 * producing valid results <em>only on rank zero</em>.
 *
 * @param[in]     scenario Scenario parameters.
 * @param[in]     grid     Grid parameters.
 * @param[in]     dgrid    Pencil decomposition parameters.
 * @param[in,out] b        B-spline basis workspace.
 * @param[in]     bop      B-spline operator workspace.
 * @param[in,out] swave    Destroyed in the computation
 * @param[in]     t        Current simulation time
 *
 * @return Mean quantities as B-spline coefficients.
 */
mean sample_mean_quantities(
        const suzerain::problem::ScenarioDefinition &scenario,
        const suzerain::problem::GridDefinition &grid,
        const suzerain::pencil_grid &dgrid,
        suzerain::bspline &b,
        const suzerain::bsplineop &bop,
        suzerain::ContiguousState<4,complex_t> &swave,
        const real_t t);

/** Store a \ref mean instance in a restart file */
void store(const esio_handle h, const mean& m);

/**
 * Load a \ref mean instance from a restart file.  Statistics not present in
 * the restart file are considered to be all NaNs.  On utter failure,
 * <tt>m.t</tt> will be NaN as well.
 */
void load(const esio_handle h, mean& m);

} // end namespace channel

#endif // SUPPORT_HPP
