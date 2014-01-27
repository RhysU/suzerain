//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
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

#ifndef SUZERAIN_SUPPORT_DEFINITION_ISOTHERMAL_HPP
#define SUZERAIN_SUPPORT_DEFINITION_ISOTHERMAL_HPP

/** @file
 * Provides \ref definition_isothermal.
 */

#include <suzerain/common.hpp>
#include <suzerain/specification_isothermal.hpp>
#include <suzerain/support/definition_base.hpp>
#include <suzerain/support/esio_fwd.hpp>
#include <suzerain/support/loadable.hpp>
#include <suzerain/support/overridable.hpp>
#include <suzerain/support/populatable.hpp>
#include <suzerain/support/savable.hpp>

namespace suzerain {

namespace support {

// FIXME: Declare other constructors per specification_isothermal.

/**
 * Holds parameters defining isothermal boundary cases.
 */
class definition_isothermal
    : public virtual definition_base
    , public virtual loadable
    , public virtual overridable<specification_isothermal>
    , public virtual populatable<specification_isothermal>
    , public virtual savable
    , public specification_isothermal
{
public:

    /**
     * Construct an instance with all parameters set to NaN.
     * Clients can use NaN as a not-yet-specified or use-the-default value.
     */
    definition_isothermal();

    /**
     * @copydoc specification_isothermal(real_t)
     */
    definition_isothermal(real_t wall_T);

    /**
     * @copydoc specification_isothermal(real_t,const std::vector<real_t>&)
     */
    definition_isothermal(real_t wall_T,
                          const std::vector<real_t>& wall_cs);

    /**
     * @copydoc specification_isothermal(real_t,real_t)
     */
    definition_isothermal(real_t wall_T,
                          real_t inflow_velocity);

    /**
     * @copydoc specification_isothermal(real_t,real_t,const std::vector<real_t>&)
     */
    definition_isothermal(real_t wall_T,
                          real_t inflow_velocity,
                          const std::vector<real_t>& wall_cs);

    /**
     * @copydoc specification_isothermal(real_t,real_t,real_t,real_t)
     */
    definition_isothermal(real_t lower_T,
                          real_t lower_v,
                          real_t upper_T,
                          real_t upper_v);

    /**
     * @copydoc specification_isothermal(real_t,real_t,const std::vector<real_t>&,real_t,real_t,const std::vector<real_t>&)
     */
    definition_isothermal(real_t lower_T,
                          real_t lower_v,
                          const std::vector<real_t>& lower_cs,
                          real_t upper_T,
                          real_t upper_v,
                          const std::vector<real_t>& upper_cs);

    /**
     * @copydoc specification_isothermal(real_t,real_t,real_t,const std::vector<real_t>&,real_t,real_t,real_t,const std::vector<real_t>&)
     */
    definition_isothermal(real_t lower_T,
                          real_t lower_v,
                          real_t lower_rho,
                          const std::vector<real_t>& lower_cs,
                          real_t upper_T,
                          real_t upper_v,
                          real_t upper_rho,
                          const std::vector<real_t>& upper_cs);

    /** @copydoc populatable::populate */
    virtual void populate(
            const specification_isothermal& that,
            const bool verbose = false);

    /** @copydoc overridable::override */
    virtual void override(
            const specification_isothermal& that,
            const bool verbose = false);

    /** @copydoc savable::save */
    virtual void save(
            const esio_handle h) const;

    /** @copydoc loadable::load */
    virtual void load(
            const esio_handle h,
            const bool verbose = true);

    /** @copydoc definition_base::options_description() */
    virtual boost::program_options::options_description options_description();

};

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_DEFINITION_ISOTHERMAL_HPP
