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
// perfect.hpp: Support logic for the Suzerain perfect gas application
// $Id$

#ifndef SUZERAIN_PERFECT_HPP
#define SUZERAIN_PERFECT_HPP

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
#include <suzerain/state.hpp>
#include <suzerain/time_definition.hpp>
#include <suzerain/timestepper.hpp>

#include "nsctpl_rholut.hpp"
#include "scenario_definition.hpp"

namespace suzerain {

/**
 * Functionality used throughout the Suzerain perfect gas application.
 */
namespace perfect {

/** Store a ScenarioDefinition in a restart file */
void store(const esio_handle h,
           const ScenarioDefinition& scenario);

/** Load a ScenarioDefinition from a restart file */
void load(const esio_handle h,
          ScenarioDefinition& scenario);

/** Manufactured solution employed throughout the channel code */
typedef nsctpl_rholut::manufactured_solution<real_t> manufactured_solution;

/**
 * Store manufactured solution parameters in a restart file.
 * Parameters are only stored when \c msoln evaluates as true.
 */
void store(const esio_handle h,
           const ScenarioDefinition& scenario,
           const problem::GridDefinition& grid,
           const boost::shared_ptr<manufactured_solution> & msoln);

/**
 * Load manufactured solution parameters from a restart file.
 * If the restart file contains active manufactured solution parameters, \c
 * msoln will be modified to contain an appropriate instance.  If it does not,
 * \c msoln will be reset.
 */
void load(const esio_handle h,
          const ScenarioDefinition& scenario,
          const problem::GridDefinition& grid,
          boost::shared_ptr<manufactured_solution>& msoln);

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
        ContiguousState<4,complex_t>& swave,
        const ScenarioDefinition& scenario,
        const problem::GridDefinition& grid,
        const pencil_grid& dgrid,
        bspline& b,
        const bsplineop& bop);

/**
 * Load the current simulation state from an open collocation point value
 * restart file.  Cannot handle interpolating onto a different grid.
 */
void load_collocation_values(
        const esio_handle h,
        ContiguousState<4,complex_t>& state,
        const ScenarioDefinition& scenario,
        const problem::GridDefinition& grid,
        const pencil_grid& dgrid,
        bspline& b,
        const bsplineop& bop);

/**
 * Interrogate an open restart file and invoke either load_coefficients()
 * or load_collocation_values() as necessary.
 */
void load(const esio_handle h,
          ContiguousState<4,complex_t>& state,
          const ScenarioDefinition& scenario,
          const problem::GridDefinition& grid,
          const pencil_grid& dgrid,
          bspline& b,
          const bsplineop& bop);

/**
 * Hold temperature and density constant while changing the Mach number and
 * ratio of specific heats.  On input, \c state should contain total energy
 * fields using \c old_Ma and \c old_gamma.  On output \c state will contain
 * total energy fields using <tt>scenario.Ma</tt> and <tt>scenario.gamma</tt>.
 */
void
adjust_scenario(ContiguousState<4,complex_t> &swave,
                const ScenarioDefinition& scenario,
                const problem::GridDefinition& grid,
                const pencil_grid& dgrid,
                bspline &b,
                const bsplineop& bop,
                const real_t old_Ma,
                const real_t old_gamma);

/** Options definitions for adding random noise to momentum fields */
class NoiseDefinition : public problem::IDefinition {

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
add_noise(ContiguousState<4,complex_t> &state,
          const NoiseDefinition& noisedef,
          const ScenarioDefinition& scenario,
          const problem::GridDefinition& grid,
          const pencil_grid& dgrid,
          bspline &b,
          const bsplineop& bop);

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
        ContiguousState<4,complex_t> &swave,
        const problem::GridDefinition &grid,
        const pencil_grid &dgrid,
        bspline &b,
        const bsplineop &bop,
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
    ((rho_u,                    3)) /* vector           */ \
    ((rho_E,                    1)) /* scalar           */

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
    ((mu_grad_T,         3))  /* vector           */ \
    ((Srho,              1))  /* scalar           */ \
    ((Srhou,             3))  /* vector           */ \
    ((SrhoE,             1))  /* scalar           */ \
    ((Srhou_dot_u,       1))  /* scalar           */

/**
 * A Boost.Preprocessor sequence of tuples of quantities computed
 * through implicit forcing.
 */
#define CHANNEL_MEAN_IMPLICIT                               \
    ((f,                 3))  /* vector           */ \
    ((qb,                1))  /* scalar           */ \
    ((f_dot_u,           1))  /* scalar           */ \
    ((Crho,              1))  /* scalar           */ \
    ((Crhou,             3))  /* vector           */ \
    ((CrhoE,             1))  /* scalar           */ \
    ((Crhou_dot_u,       1))  /* scalar           */

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
        const ScenarioDefinition &scenario,
        const problem::GridDefinition &grid,
        const pencil_grid &dgrid,
        bspline &b,
        const bsplineop &bop,
        ContiguousState<4,complex_t> &swave,
        const real_t t);

/** Store a \ref mean instance in a restart file */
void store(const esio_handle h, const mean& m);

/**
 * Load a \ref mean instance from a restart file.  Statistics not present in
 * the restart file are considered to be all NaNs.  On utter failure,
 * <tt>m.t</tt> will be NaN as well.
 */
void load(const esio_handle h, mean& m);

} // end namespace perfect

} // end namespace suzerain

#endif // SUZERAIN_PERFECT_HPP
