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

#ifndef SUZERAIN_SLOWGROWTH_HPP
#define SUZERAIN_SLOWGROWTH_HPP

/** @file
 * Provides \ref slowgrowth to facilitate slow growth forcing usage.
 */

#include <suzerain/l2.hpp>
#include <suzerain/largo_state.hpp>
#include <suzerain/state_fwd.hpp>

namespace suzerain {

// Forward declarations
class operator_tools;
class specification_largo;

// TODO Doxygen for class slowgrowth

/**
 * Provides relatively clean API for applying slow growth forcing.
 * Largo's APIs are good, but gathering the required information to
 * use them is cumbersome, error-prone, and best not duplicated.
 */
class slowgrowth
{

public:

    slowgrowth(const specification_largo& sg,
               const real_t code_Ma);

    void
    initialize(const std::size_t substep_index = 0);

    void
    gather_wavexz(const operator_tools &otool,
                  const contiguous_state<4,complex_t> &swave);

    /**
     * Abstract data source for \ref gather_physical_cons.  Each column
     * represents a mean quantity indexed by collocation point.
     */
    class physical_cons
    {
    public:

        /** Virtual for interface */
        virtual ~physical_cons();

        /** Total energy */
        virtual ArrayXXr::ConstColXpr rhoE()  const = 0;

        /** Streamwise momentum */
        virtual ArrayXXr::ConstColXpr rhou()  const = 0;

        /** Wall-normal momentum */
        virtual ArrayXXr::ConstColXpr rhov()  const = 0;

        /** Spanwise momentum */
        virtual ArrayXXr::ConstColXpr rhow()  const = 0;

        /** Density */
        virtual ArrayXXr::ConstColXpr rho ()  const = 0;

        /** Pressure */
        virtual ArrayXXr::ConstColXpr p   ()  const = 0;

        /** Pressure squared */
        virtual ArrayXXr::ConstColXpr p2  ()  const = 0;
    };

    void
    gather_physical_cons(const operator_tools &otool,
                         const physical_cons &data);

    /**
     * Abstract data source for \ref gather_physical_rqq.  Each column
     * represents a mean quantity indexed by collocation point.
     */
    class physical_rqq
    {
    public:
        /** Virtual for interface */
        virtual ~physical_rqq();

        /** Density times specific energy squared */
        virtual ArrayXXr::ConstColXpr rhoEE() const = 0;

        /** Density times streamwise velocity squared */
        virtual ArrayXXr::ConstColXpr rhouu() const = 0;

        /** Density times wall-normal velocity squared */
        virtual ArrayXXr::ConstColXpr rhovv() const = 0;

        /** Density times spanwise velocity squared */
        virtual ArrayXXr::ConstColXpr rhoww() const = 0;

        /** Density */
        virtual ArrayXXr::ConstColXpr rho  () const = 0;
    };

    void
    gather_physical_rqq(const operator_tools &otool,
                        const physical_rqq &data);

    void
    inner_y(const int j,
            const real_t y_j);

    void
    inner_xz(largo_state &local_state,
             largo_state &slowgrowth_forcing);

protected:

    void
    calculate_baseflow(const real_t y,
                       largo_state &base,
                       largo_state &dy,
                       largo_state &dt,
                       largo_state &dx,
                       largo_state &src) const;

    const specification_largo &sg;

    const real_t inv_codeMa2;

private:

    largo_state basewall;
    std::vector<field_L2xz> meanrms;
    std::vector<field_L2xz> meanrms_y;
    ArrayX5r rqq;
    ArrayX5r rqq_y;

};

} // namespace suzerain

#endif  /* SUZERAIN_SLOWGROWTH_HPP */
