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
#include <suzerain/format.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/support/application_base.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/support.hpp>

using namespace suzerain;

/** Test application using suzerain::support framework. */
struct test : public support::application_base
{
    test()
        : support::application_base("Tests suzerain::compute_field_L2xz")
        , who("test")
    {}

    std::string log4cxx_config() { return support::log4cxx_config_console; }

    // Invoked from main(...) and implemented below.
    int run(int argc, char **argv);

private:

    std::string who;
};

// One dimensional test function suggested by S. Johnson.  Magically has mean
// zero, total RMS 1/Sqrt(2), but no closed Fourier representation.  See
// http://sci.tech-archive.net/Archive/sci.math/2008-05/msg00401.html
template< typename Scalar >
static inline
Scalar sample_data(const Scalar& x, const Scalar& L)
{
    using boost::math::constants::pi;
    using std::cos;
    const Scalar twopi_L = 2*pi<Scalar>()/L;
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

int main(int argc, char **argv)
{
    test t;
    return t.run(argc, argv);
}

int test::run(int argc, char **argv)
{
    // Establish default grid and domain extents
    grid.reset(new support::grid_definition( 5              // Lx
                                           , 32             // Nx
                                           , 2              // DAFx
                                           , 3              // Ly
                                           , 6              // Ny
                                           , 4              // k
                                           , 0.5            // htdelta
                                           , 7              // Lz
                                           , 64             // Nz
                                           , real_t(3) / 2  // DAFz
            ));

    // Add additional command line options
    std::size_t nfields = 3;
    bool check          = true;
    bool constant       = false;
    options.add_options()
        ("nfields,n", boost::program_options::value(&nfields)
         ->default_value(nfields),
         "Number of independent scalar fields")
        ("check,c", boost::program_options::value(&check)
         ->default_value(check)->zero_tokens(),
         "Check results against expected values")
        ("constant,C", boost::program_options::value(&constant)
         ->default_value(constant)->zero_tokens(),
         "Employ simpler constant-valued scalar manufactured fields")
    ;

    // Initialize the parallel decomposition, state storage, and operators
    std::vector<std::string> positional = initialize(argc, argv);
    if (positional.size() != 0) {
        FATAL0(who, "No positional arguments accepted");
        return EXIT_FAILURE;
    }
    establish_ieee_mode();
    load_grid_and_operators(NULL);
    establish_decomposition();
    suzerain::operator_base o(*grid, *dgrid, *cop, *b);
    establish_state_storage(/* linear state is wave-only        */ 0,
                            /* nonlinear state is transformable */ nfields);
    physical_view<> p(*dgrid, *state_nonlinear);
    p.setConstant(std::numeric_limits<real_t>::quiet_NaN());  // ++paranoia

    // Either constant-valued or spatially-varying fields may be initialized.
    // The former is good for checking grid-, dealiasing-, and normalization.
    // The latter is good for checking fluctuating magnitude computations.
    if (constant) {
        INFO0(who, "Initializing constant-valued scalar fields");
        for (std::size_t f = 0; f < nfields; ++f)
            p.row(f).setConstant(f + 1);
    } else {
        INFO0(who, "Initializing spatially-varying scalar fields");
        for (int offset = 0, j = dgrid->local_physical_start.y();     // Y
             j < dgrid->local_physical_end.y();
             ++j) {
            const real_t y = o.y(j);

            for (int k = dgrid->local_physical_start.z();             // Z
                 k < dgrid->local_physical_end.z();
                 ++k) {
                const real_t z = o.z(k);

                for (int i = dgrid->local_physical_start.x();         // X
                     i < dgrid->local_physical_end.x();
                     ++i, /* NB */ ++offset) {
                    const real_t x = o.x(i);

                    // Notice scaling by wall-normal coordinate and field index
                    const real_t data = y * sample_data(x, grid->L.x(),
                                                        z, grid->L.z());
                    for (std::size_t f = 0; f < nfields; ++f) {       // Fields
                        p(f, offset) = (f + 1)*data;
                    }
                }
            }
        }
    }
    for (std::size_t f = 0; f < nfields; ++f) {
        dgrid->transform_physical_to_wave(&p.coeffRef(f, 0));  // X, Z
        o.zero_dealiasing_modes(*state_nonlinear, f);
        o.bop_solve(*o.massluz(), *state_nonlinear, f);        // Y
    }
    p *= dgrid->chi();

    INFO0(who, "Compute L^2_{xz} at each collocation point");
    std::vector<field_L2xz> L2xz = compute_field_L2xz(
            *state_nonlinear, *grid, *dgrid, *cop);

    // Coefficient to convert L^2_xz results into RMS results per l2.hpp
    const real_t rms_adjust = 1 / std::sqrt(grid->L.x() * grid->L.z());

    INFO0(who, "Mean RMS for each field as a function of collocation point:");
    for (int j = 0; j < grid->N.y(); ++j) {
        std::ostringstream msg;
        msg << fullprec<>(o.y(j));
        for (std::size_t f = 0; f < nfields; ++f) {
            msg << ' ' << fullprec<>(rms_adjust * L2xz[f].mean(j));
        }
        INFO0(who, msg.str());
    }

    INFO0(who,
          "Fluctuating RMS for each field as a function of collocation point:");
    for (int j = 0; j < grid->N.y(); ++j) {
        std::ostringstream msg;
        msg << fullprec<>(o.y(j));
        for (std::size_t f = 0; f < nfields; ++f) {
            msg << ' ' << fullprec<>(rms_adjust * L2xz[f].fluctuating(j));
        }
        INFO0(who, msg.str());
    }

    // Stop processing if correctness checking not requested
    if (!check) return EXIT_SUCCESS;

    // TODO Error if results not good to some tolerance

    return EXIT_SUCCESS;
}
