/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * functional.hpp: functional tools
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_FUNCTIONAL_HPP
#define __SUZERAIN_FUNCTIONAL_HPP

#include <suzerain/common.hpp>

namespace suzerain {

/**
 * Provides general functional programming tools.
 */
namespace functional {

/**
 * A functor that performs assignment to type \c Target from type \c Source.
 * \c Target must be assignable from \c Source.
 *
 * The \c Enable template parameter allows using <tt>boost::enable_if</tt> to
 * specialize or partially specialize the functor per <tt>enable_if</tt>'s
 * documentation section <a
 * href="http://www.boost.org/doc/libs/1_41_0/libs/utility/enable_if.html">
 * 3.1: Enabling template class specializations</a>.
 */
template<class Target, class Source, class Enable = void>
struct assign {

    /**
     * Create an instance which assigns \c s when applied.
     *
     * @param s source of assignment operation occurring via
     *          <tt>operator()</tt>.
     */
    assign(const Source &s) : s_(s) {}

    /**
     * Assign the value provided at construction to \c t.
     *
     * @param t to be assigned.
     */
    void operator()(Target& t) const { t = s_; }

private:
    const Source &s_; /**< Source for assignment operations */
};

} // namespace functional

} // namespace suzerain

#endif // __SUZERAIN_FUNCTIONAL_HPP
