//--------------------------------------------------------------------------
//
// Copyright (C) 2008-2014 Rhys Ulerich
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

/** @file
 * A suzerain::pencil_grid test based on Dmitry Pekurovsky's P3DFFT work.
 */

#include <suzerain/common.hpp>
#include <suzerain/error.h>
#include <suzerain/operator_base.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/support/application_base.hpp>
#include <suzerain/support/logging.hpp>

using namespace suzerain;

static real_t sample_field(const real_t x, const real_t y, const real_t z)
{
    using std::sin;
    return   1* 1
           +    2*y
           + 2* 3*sin(x)
           + 2* 5*sin(z)
           + 2* 7*sin(x)*y
           + 2*11*y*sin(z)
           + 4*13*sin(x)*sin(z)
           + 4*17*sin(x)*y*sin(z)
           ;
}

int main(int argc, char **argv)
{
    // Instantiate application_base instance and set desired defaults
    support::application_base app(
            "Benchmark parallel FFTs implemented using suzerain::pencil_grid",
            /* No arguments */ "", /* No further description*/ "", REVISIONSTR);
    app.grid->htdelta = 0;
    app.grid->k       = 4;
    app.grid->L.x()   = 8 * boost::math::constants::pi<real_t>();
    app.grid->L.y()   = 2;
    app.grid->L.z()   = 4 * boost::math::constants::pi<real_t>();
    app.grid->Nx(16);
    app.grid->Ny(16);
    app.grid->Nz(16);
    app.grid->DAFx(1);
    app.grid->DAFz(1);

    // Add additional command line options
    std::size_t repeat  = 1;
    std::size_t nfields = 1;
    bool check          = false;
    app.options.add_options()
        ("repeat,r", boost::program_options::value(&repeat)
         ->default_value(repeat),
         "Number of repetitions to perform")
        ("nfields,n", boost::program_options::value(&nfields)
         ->default_value(nfields),
         "Number of independent fields")
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

    // Initialize manufactured field using operator_base utilities
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

                const real_t data = sample_field(x, y, z);
                for (std::size_t f = 0; f < nfields; ++f) {        // Fields
                    p(f, offset) = f + data;
                }

            }

        }

    }

    SUZERAIN_MPICHKQ(MPI_Barrier(MPI_COMM_WORLD));

    // Perform transformations accumulating parallel FFT timing in time
    // Barriers are used to ensure each timing starts from synchronized state
    double time = 0;
    for (std::size_t m = 0; m < repeat; m++) {
        INFO0("Iteration " << (m + 1));

        // Physical to wave
        SUZERAIN_MPICHKQ(MPI_Barrier(MPI_COMM_WORLD));
        time -= MPI_Wtime();
        for (std::size_t l = 0; l < nfields; ++l) {
            app.dgrid->transform_physical_to_wave(p.data());
        }
        time += MPI_Wtime();

        // Wave space normalization uses only dealiased Fourier X and Z extents
        app.state_nonlinear->scale(app.dgrid->chi());

        // Wave to physical
        SUZERAIN_MPICHKQ(MPI_Barrier(MPI_COMM_WORLD));
        time -= MPI_Wtime();
        for (std::size_t l = 0; l < nfields; ++l) {
            app.dgrid->transform_wave_to_physical(p.data());
        }
        time += MPI_Wtime();

    }

    // Report the worst time observed on any rank
    const int rank = mpi::comm_rank(MPI_COMM_WORLD);
    struct { double time; int rank; } worst = { time, rank };
    SUZERAIN_MPICHKQ(MPI_Allreduce(MPI_IN_PLACE, &worst, 1, MPI_DOUBLE_INT,
                                   MPI_MAXLOC, MPI_COMM_WORLD));
    if (worst.rank == rank) {
        INFO("Seconds per iteration " << time / (double) repeat);
    }

    // Stop processing if correctness checking not requested
    if (!check) return EXIT_SUCCESS;

    // Quantify discrepancies between source data and computed results
    namespace acc = boost::accumulators;
    acc::accumulator_set<
            real_t, acc::stats< acc::tag::min, acc::tag::max >
        > track;
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

                const real_t data = sample_field(x, y, z);
                for (std::size_t f = 0; f < nfields; ++f) {        // Fields
                    track((data + f) - p(f,offset));
                }

            }

        }

    }

    // Report the minimum discrepancy observed on any rank
    struct { double discrepancy; int rank; } min = { acc::min(track), rank };
    SUZERAIN_MPICHKQ(MPI_Allreduce(MPI_IN_PLACE, &min, 1, MPI_DOUBLE_INT,
                                   MPI_MINLOC, MPI_COMM_WORLD));
    if (min.rank == rank) {
        INFO("Minimum discrepancy " << std::scientific << min.discrepancy);
    }

    // Report the maximum discrepancy observed on any rank
    struct { double discrepancy; int rank; } max = { acc::max(track), rank };
    SUZERAIN_MPICHKQ(MPI_Allreduce(MPI_IN_PLACE, &max, 1, MPI_DOUBLE_INT,
                                   MPI_MAXLOC, MPI_COMM_WORLD));
    if (max.rank == rank) {
        INFO("Maximum discrepancy " << std::scientific << max.discrepancy);
    }

    return EXIT_SUCCESS;
}
