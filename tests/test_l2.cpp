//--------------------------------------------------------------------------
//
// Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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

/** @file
 * A suzerain::compute_field_L2xz test using suzerain::support infrastructure.
 */

#include <suzerain/l2.hpp>

#include <suzerain/common.hpp>
#include <suzerain/error.h>
#include <suzerain/operator_base.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/support/application_base.hpp>
#include <suzerain/support/logging.hpp>

// One dimensional test function suggested by S. Johnson.  Magically has mean
// zero, total RMS 1/Sqrt(2), but no closed Fourier representation.  See
// http://sci.tech-archive.net/Archive/sci.math/2008-05/msg00401.html
template< typename Scalar >
static inline
Scalar sample_data(const Scalar& x, const Scalar& L)
{
    using std::cos;
    const Scalar twopi_L = 2*boost::math::constants::pi<Scalar>()/L;
    return cos(twopi_L*(x + 2*cos(twopi_L*3*x)));
}

// Tensor product of 2 univariate test functions.
// Again, zero mean but now total RMS 1/2.
template< typename Scalar >
static inline
Scalar sample_data(const Scalar& x, const Scalar& Lx,
                   const Scalar& z, const Scalar& Lz)
{
    return sample_data(x, Lx) * sample_data(z, Lz);
}

using namespace suzerain;

int main(int argc, char **argv)
{
    // Instantiate application_base instance and set desired defaults.
    // Defaults designed to increase test coverage by tickling edge cases.
    support::application_base app("Tests suzerain::compute_field_L2xz");
    app.grid->htdelta = 0.5;
    app.grid->k       = 4;
    app.grid->L.x()   = 8 * boost::math::constants::pi<real_t>();
    app.grid->L.y()   = 2;
    app.grid->L.z()   = 4 * boost::math::constants::pi<real_t>();
    app.grid->Nx(32);
    app.grid->Ny(16);
    app.grid->Nz(64);
    app.grid->DAFx(1.5);
    app.grid->DAFz(2.0);

    // Add additional command line options
    std::size_t nfields = 3;
    bool check          = true;
    app.options.add_options()
        ("nfields,n", boost::program_options::value(&nfields)
         ->default_value(nfields),
         "Number of independent scalar fields")
        ("check,c", boost::program_options::value(&check)
         ->default_value(check)->zero_tokens(),
         "Check results against expected values")
    ;

    // Initialize the parallel decomposition
    app.initialize(argc, argv);
    app.establish_ieee_mode();
    app.load_grid_and_operators(NULL);
    app.establish_decomposition();

    // Allocate necessary storage and prepare physical space view
    app.establish_state_storage(/* linear state is wave-only        */ 0,
                                /* nonlinear state is transformable */ nfields);
    physical_view<> p(*app.dgrid, *app.state_nonlinear);

    // Initialize physical-space manufactured field using operator_base.
    // p is a 2D, real-valued view of (NFIELDS, X*Z*Y) leftmost fastest
    // We need x_i, y_j, z_k positions making structure below more complicated
    suzerain::operator_base o(*app.grid, *app.dgrid, *app.cop, *app.b);
    for (int offset = 0, j = app.dgrid->local_physical_start.y();  // Y
         j < app.dgrid->local_physical_end.y();
         ++j) {
        const real_t y = o.y(j);

        for (int k = app.dgrid->local_physical_start.z();          // Z
             k < app.dgrid->local_physical_end.z();
             ++k) {
            const real_t z = o.z(k);

            for (int i = app.dgrid->local_physical_start.x();      // X
                 i < app.dgrid->local_physical_end.x();
                 ++i, /* NB */ ++offset) {
                const real_t x = o.x(i);

                const real_t data = y * sample_data(x, app.grid->L.x(),
                                                    z, app.grid->L.z());
                for (std::size_t f = 0; f < nfields; ++f) {        // Fields
                    p(f, offset) = (f + 1)*data;
                }

            }

        }

    }

    // Convert the physical space values to wave space
    for (std::size_t f = 0; f < nfields; ++f) {
        app.dgrid->transform_physical_to_wave(&p.coeffRef(0, 0));  // X, Z
        o.zero_dealiasing_modes(*app.state_nonlinear, 0);
        o.bop_solve(*o.massluz(), *app.state_nonlinear, 0);        // Y
    }

    // Compute the L^2 norm over the X and Z directions at collocation points
    std::vector<field_L2xz> L2xz = compute_field_L2xz(
            *app.state_nonlinear, *app.grid, *app.dgrid, *app.cop);

    // TODO Output results

    // Stop processing if correctness checking not requested
    if (!check) return EXIT_SUCCESS;

    // TODO Error if results not good to some tolerance

    return EXIT_SUCCESS;
}
