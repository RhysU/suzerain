/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
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
 * bspline_definition.hpp: classes handling B-spline parameters
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_BSPLINE_DEFINITION_HPP
#define __SUZERAIN_BSPLINE_DEFINITION_HPP

#include <suzerain/common.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/types.hpp>
#include <suzerain/validation.hpp>

/** @file
 * Provides classes handling B-spline basis definitions, which are runtime
 * arguments or scenario parameters used to perform calculations.
 */

namespace suzerain {

namespace problem {

// TODO Add parameter for different collocation point choices when appropriate.

/**
 * Holds basic B-spline basis details, including the order and parameter(s)
 * describing breakpoints.  Basis order is the piecewise polynomial order
 * plus one.  For example, piecewise cubics are order four.
 *
 */
template< typename FPT = double >
class BsplineDefinition : public IDefinition, public integral_types
{
public:
    /**
     * Construct an instance with the given default order k and breakpoint
     * stretching parameter alpha.
     *
     * @param default_k Default B-spline basis order.
     * @param default_alpha Default grid stretching parameter.
     */
    BsplineDefinition(size_type default_k = 16, FPT default_alpha = 5.0);

    /**
     * Retrieve the B-spline basis order.
     * order plus one.  For example, piecewise cubics
     *
     * @return the B-spline basis order.
     */
    size_type k() const { return k_; }

    /**
     * Retrieve the breakpoint stretching parameter.  The parameter is used to
     * cluster breakpoints near solid boundaries, e.g. at a flat plate or at
     * the top and bottom wall in a channel flow.
     *
     * @return the breakpoint stretching parameter.
     * @see For information on how alpha is used, see suzerain::math::stretchspace.
     *
     */
    FPT alpha() const { return alpha_; }

    /*! @copydoc IDefinition::options */
    const boost::program_options::options_description& options() {
        return options_;
    }

private:

    /** Stores the program options processing information */
    boost::program_options::options_description options_;

    size_type k_;  /**< Stores the B-spline basis order */
    FPT alpha_;    /**< Stores the breakpoint stretching parameter */
};

template< typename FPT >
BsplineDefinition<FPT>::BsplineDefinition(size_type default_k,
                                          FPT default_alpha)
    : options_("B-spline basis definition"),
      k_(default_k),
      alpha_(default_alpha)
{
    namespace po = ::boost::program_options;

    using ::std::bind2nd;
    using ::std::ptr_fun;
    using ::suzerain::validation::ensure_nonnegative;
    using ::suzerain::validation::ensure_positive;

    // Created to solve ambiguous type issues below
    ::std::pointer_to_binary_function<FPT,const char*,void>
        ptr_fun_ensure_positive_FPT(ensure_positive<FPT>);

    options_.add_options()
        ("k", po::value<size_type>(&global_extents_[0])
            ->notifier(bind2nd(ptr_fun(ensure_positive<size_type>),"k"))
            ->default_value(default_size),
        "B-spline basis order where k = 4 indicates piecewise cubics")
        ("alpha", po::value<FPT>(&alpha_)
            ->notifier(bind2nd(ptr_fun_ensure_positive_FPT,"alpha"))
            ->default_value(default_length),
        "B-spline breakpoint stretching parameter")
    ;
}

} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_BSPLINE_DEFINITION_HPP
