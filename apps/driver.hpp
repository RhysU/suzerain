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
// driver.hpp: Application driver logic spanning multiple applications
// $Id$

#ifndef SUZERAIN_SUPPORT_DRIVER_HPP
#define SUZERAIN_SUPPORT_DRIVER_HPP

#include <suzerain/common.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/fftw.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/restart_definition.hpp>
#include <suzerain/signal_definition.hpp>
#include <suzerain/state.hpp>
#include <suzerain/statistics_definition.hpp>
#include <suzerain/time_definition.hpp>
#include <suzerain/timestepper.hpp>

namespace suzerain {

namespace support {

class Driver
{
public:

    Driver();

    virtual ~Driver();

protected:

    typedef InterleavedState<4,complex_t> linear_state_type;

    typedef ContiguousState<4,complex_t>  nonlinear_state_type;

    const problem::GridDefinition grid;

    const fftw::FFTWDefinition fftwdef;

    const problem::RestartDefinition restart;

    const problem::StatisticsDefinition statsdef;

    const problem::TimeDefinition timedef;

    const problem::SignalDefinition sigdef;

private:
};

} // end namespace support

} // end namespace suzerain

#endif // SUZERAIN_SUPPORT_DRIVER_HPP
