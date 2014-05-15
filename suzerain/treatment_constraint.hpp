//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
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

#ifndef SUZERAIN_TREATMENT_CONSTRAINT_HPP
#define SUZERAIN_TREATMENT_CONSTRAINT_HPP

/** @file
 * Provides \ref constraint::treatment.
 */

#include <suzerain/common.hpp>
#include <suzerain/lowstorage.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/state_fwd.hpp>

namespace suzerain {

// Forward declarations
class bspline;
class pencil_grid;

namespace constraint {

// Forward declarations
class base;
class disabled;

/**
 * A wrapper applying integral constraint treatment atop any linear operator.
 * During \ref invert_mass_plus_scaled_operator implicit momentum forcing is
 * applied following the section of <tt>writeups/treatment_channel.tex</tt>
 * titled "Enforcing a target bulk momentum via the linear operator".
 */
class treatment
    : public lowstorage::linear_operator<
          multi_array::ref<complex_t,4>,
          contiguous_state<4,complex_t>
      >
    , public boost::noncopyable
{
public:

    /**
    * Abstract interface for tracking instantaneous collocation point
    * profiles taken to be inputs to constraint treatment logic.
    */
    class inputs
    {
    public:
        virtual ~inputs();                            ///< Virtual for interface
        virtual ArrayXXr::ConstColXpr u()  const = 0; ///< Streamwise velocity
        virtual ArrayXXr::ConstColXpr v()  const = 0; ///< Wall-normal velocity
        virtual ArrayXXr::ConstColXpr w()  const = 0; ///< Spanwise velocity
        virtual ArrayXXr::ConstColXpr uu() const = 0; ///< Reynolds stress component uu
        virtual ArrayXXr::ConstColXpr uv() const = 0; ///< Reynolds stress component uv
        virtual ArrayXXr::ConstColXpr uw() const = 0; ///< Reynolds stress component uw
        virtual ArrayXXr::ConstColXpr vv() const = 0; ///< Reynolds stress component vv
        virtual ArrayXXr::ConstColXpr vw() const = 0; ///< Reynolds stress component vw
        virtual ArrayXXr::ConstColXpr ww() const = 0; ///< Reynolds stress component ww
    };

    /**
    * Abstract interface for tracking forcing collocation point profile
    * outputs maintained by constraint treatment logic.
    */
    class outputs
    {
    public:
        virtual ~outputs();                          ///< Virtual for interface
        virtual ArrayXXr::ColXpr fx()           = 0; ///< Streamwise momentum
        virtual ArrayXXr::ColXpr fy()           = 0; ///< Wall-normal momentum
        virtual ArrayXXr::ColXpr fz()           = 0; ///< Spanwise momentum
        virtual ArrayXXr::ColXpr f_dot_u()      = 0; ///< Forcing work
        virtual ArrayXXr::ColXpr qb()           = 0; ///< Heating
        virtual ArrayXXr::ColXpr CrhoE()        = 0; ///< Constrained total energy
        virtual ArrayXXr::ColXpr Crhou()        = 0; ///< Constrained streamwise momentum
        virtual ArrayXXr::ColXpr Crhov()        = 0; ///< Constrained wall-normal momentum
        virtual ArrayXXr::ColXpr Crhow()        = 0; ///< Constrained spanwise momentum
        virtual ArrayXXr::ColXpr Crhou_dot_u()  = 0; ///< Constrained forcing work
        virtual ArrayXXr::ColXpr Crho()         = 0; ///< Constrained density
        virtual ArrayXXr::ColXpr C2rhoE()       = 0; ///< Squared \ref CrhoE
        virtual ArrayXXr::ColXpr C2rhou()       = 0; ///< Squared \ref Crhou
        virtual ArrayXXr::ColXpr C2rhov()       = 0; ///< Squared \ref Crhov
        virtual ArrayXXr::ColXpr C2rhow()       = 0; ///< Squared \ref Crhow
        virtual ArrayXXr::ColXpr C2rhou_dot_u() = 0; ///< Squared \ref Crhou_dot_u
        virtual ArrayXXr::ColXpr C2rho()        = 0; ///< Squared \ref Crho
    };

    /**
     * Constructor.  After construction, #L must be provided.  Generally, one
     * or more constraints will be supplied via #physical or #numerical before
     * usage.
     *
     * @param Ma     Nondimensional Mach number for kinetic energy computations.
     *               In a dimensional setting, this should be one.  Notice that
     *               it is a \e reference permitting the value to track another
     *               setting.  For example, one in a \ref definition_scenario.
     * @param dgrid  Parallel decomposition details used to determine which
     *               rank houses the "zero-zero" Fourier modes.
     * @param inp    Profiles to be used during constraint application.
     * @param out    Profiles to be tracked due to constraint application.
     */
    treatment(const real_t& Ma,
              const pencil_grid& dgrid,
              bspline& b,
              inputs&  inp,
              outputs& out);

    /** Delegates invocation to #L */
    virtual void apply_mass_plus_scaled_operator(
             const complex_t &phi,
             multi_array::ref<complex_t,4> &state,
             const std::size_t substep_index) const;

    /** Delegates invocation to #L */
    virtual void accumulate_mass_plus_scaled_operator(
             const complex_t &phi,
             const multi_array::ref<complex_t,4> &input,
             const complex_t &beta,
             contiguous_state<4,complex_t> &output,
             const std::size_t substep_index) const;

    /**
     * Force the channel problem delegating to #L when appropriate.
     * #L is responsible for enforcing all boundary conditions.
     */
    virtual void invert_mass_plus_scaled_operator(
            const complex_t& phi,
            multi_array::ref<complex_t,4>& state,
            const lowstorage::method_interface<complex_t>& method,
            const real_t delta_t,
            const std::size_t substep_index,
            multi_array::ref<complex_t,4> *ic0 = NULL) const;

    /** The operator whose behavior is modified by this instance. */
    shared_ptr<lowstorage::linear_operator<
                multi_array::ref<complex_t,4>,
                contiguous_state<4,complex_t>
            > > L;

    /**
     * Catalog of physically-oriented constraints to be indexed by equation
     * number ndx::e, ndx::mx, ndx::my, ndx::mz, or ndx::rho.  These are
     * constraints which interpretable as physics-related source terms which
     * may include contributions to multiple equations.  For example, an
     * ndx::mx momentum constraint will also contribute forcing work terms to
     * the ndx::e total energy.
     *
     * \warning The ndx::rho entry is not treated "physically" in the sense
     * that density added by the constraint updates neither the momentum nor
     * total energy equations.  Therefore, when the physical[ndx::rho]
     * constraint modifies the density it also adjusts the velocity and
     * specific total energy.
     *
     * Means of the implicit ndx::mx, ndx::my, ndx::mz, and ndx::e forcing are
     * maintained across each individual time step for sampling the statistics
     * \c /bar_f, \c /bar_f_dot_u, and \c /bar_qb using #com.  The mean
     * implicit numerical[ndx::rho] forcing is combined with
     * physical[ndx::rho]'s effect and stored within \c /bar_Crho.
     */
    array<shared_ptr<constraint::base>, 5> physical;

    /**
     * Catalog of numerically-oriented constraints to be indexed by equation
     * number ndx::e, ndx::mx, ndx::my, ndx::mz, or ndx::rho.  These are
     * numerical implementation artifacts to control the behavior of single
     * equations rather than physics-related sources.  For example, an
     * ndx::mx momentum constraint will not not contribute forcing work terms
     * to the ndx::e total energy.
     *
     * Means of the implicit ndx::mx, ndx::my, ndx::mz, and ndx::e forcing are
     * maintained across each individual time step for sampling the statistics
     * \c /bar_Crhou, \c /bar_Crhou_dot_u, and \c /bar_CrhoE using #com.
     * The mean implicit numerical[ndx::rho] forcing is combined with
     * physical[ndx::rho]'s effect and stored within \c /bar_Crho.
     */
    array<shared_ptr<constraint::base>, 5> numerical;

    /** An appropriately-size, do-nothing constraint usable by callers. */
    const shared_ptr<constraint::disabled> none;

protected:

    /**
     * What Mach number should be used to scale kinetic
     * energy contributions to the total energy equation?
     */
    const real_t& Ma;

    /** Does this rank possess the "zero-zero" Fourier modes? */
    const bool rank_has_zero_zero_modes;

    /** Profiles to be used during constraint application. */
    inputs&  inp;

    /** Profiles to be tracked due to constraint application. */
    outputs& out;

private:

    /**
     * Constraint data passed to #L indexed by ndx::type.  Entries <tt>[ndx::e,
     * ndx::rho]</tt> pertain to #physical.  Entries <tt>[ndx::rho_(1) +
     * ndx::e, ndx::rho_(1) + ndx::rho]</tt> pertain to #numerical.  Mutable
     * member avoids repeated allocation/deallocation.
     */
    mutable MatrixXAc cdata;

    /**
     * Least squares constraint solver.
     * Mutable member avoids repeated allocation/deallocation.
     */
    mutable Eigen::JacobiSVD<
            MatrixXXr, Eigen::FullPivHouseholderQRPreconditioner
        > jacobiSvd;

};

} // namespace constraint

} // namespace suzerain

#endif /* SUZERAIN_TREATMENT_CONSTRAINT_HPP */
