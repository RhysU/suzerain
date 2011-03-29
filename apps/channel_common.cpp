/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
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
 * channel_common.cpp: Channel-related functionality spanning binaries
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <suzerain/htstretch.h>
#include "logger.hpp"
#include "channel_common.hpp"

using boost::numeric_cast;

const boost::array<const char *,5> field_names = {{
    "rho", "rhou", "rhov", "rhow", "rhoe"
}};

const boost::array<const char *,5> field_descriptions = {{
    "Nondimensional density coefficients stored row-major ZXY",
    "Nondimensional X momentum coefficients stored row-major ZXY",
    "Nondimensional Y momentum coefficients stored row-major ZXY",
    "Nondimensional Z momentum coefficients stored row-major ZXY",
    "Nondimensional total energy coefficients stored row-major ZXY"
}};

void store(const esio_handle h,
           const suzerain::problem::ScenarioDefinition<real_t>& scenario)
{
    DEBUG0("Storing ScenarioDefinition parameters");

    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);

    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));

    esio_line_write(h, "Re", &scenario.Re, 0,
            scenario.options().find("Re",false).description().c_str());

    esio_line_write(h, "Pr", &scenario.Pr, 0,
            scenario.options().find("Pr",false).description().c_str());

    esio_line_write(h, "gamma", &scenario.gamma, 0,
            scenario.options().find("gamma",false).description().c_str());

    esio_line_write(h, "beta", &scenario.beta, 0,
            scenario.options().find("beta",false).description().c_str());

    esio_line_write(h, "Lx", &scenario.Lx, 0,
            scenario.options().find("Lx",false).description().c_str());

    esio_line_write(h, "Ly", &scenario.Ly, 0,
            scenario.options().find("Ly",false).description().c_str());

    esio_line_write(h, "Lz", &scenario.Lz, 0,
            scenario.options().find("Lz",false).description().c_str());
}

void load(const esio_handle h,
          suzerain::problem::ScenarioDefinition<real_t>& scenario)
{
    DEBUG0("Loading ScenarioDefinition parameters");

    esio_line_establish(h, 1, 0, 1); // All ranks load

    if (scenario.Re) {
        INFO0("Overriding scenario using Re = " << scenario.Re);
    } else {
        esio_line_read(h, "Re", &scenario.Re, 0);
    }

    if (scenario.Pr) {
        INFO0("Overriding scenario using Pr = " << scenario.Pr);
    } else {
        esio_line_read(h, "Pr", &scenario.Pr, 0);
    }

    if (scenario.gamma) {
        INFO0("Overriding scenario using gamma = " << scenario.gamma);
    } else {
        esio_line_read(h, "gamma", &scenario.gamma, 0);
    }

    if (scenario.beta) {
        INFO0("Overriding scenario using beta = " << scenario.beta);
    } else {
        esio_line_read(h, "beta", &scenario.beta, 0);
    }

    if (scenario.Lx) {
        INFO0("Overriding scenario using Lx = " << scenario.Lx);
    } else {
        esio_line_read(h, "Lx", &scenario.Lx, 0);
    }

    if (scenario.Ly) {
        INFO0("Overriding scenario using Ly = " << scenario.Ly);
    } else {
        esio_line_read(h, "Ly", &scenario.Ly, 0);
    }

    if (scenario.Lz) {
        INFO0("Overriding scenario using Lz = " << scenario.Lz);
    } else {
        esio_line_read(h, "Lz", &scenario.Lz, 0);
    }
}

void store(const esio_handle h,
           const suzerain::problem::GridDefinition& grid,
           const real_t Lx,
           const real_t Lz)
{
    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);

    DEBUG0("Storing GridDefinition parameters");

    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));

    esio_line_write(h, "Nx", &grid.N.x(), 0,
            grid.options().find("Nx",false).description().c_str());

    esio_line_write(h, "DAFx", &grid.DAF.x(), 0,
            grid.options().find("DAFx",false).description().c_str());

    esio_line_write(h, "Ny", &grid.N.y(), 0,
            grid.options().find("Ny",false).description().c_str());

    esio_line_write(h, "k", &grid.k, 0,
            grid.options().find("k",false).description().c_str());

    esio_line_write(h, "htdelta", &grid.htdelta, 0,
            grid.options().find("htdelta",false).description().c_str());

    esio_line_write(h, "Nz", &grid.N.z(), 0,
            grid.options().find("Nz",false).description().c_str());

    esio_line_write(h, "DAFz", &grid.DAF.z(), 0,
            grid.options().find("DAFz",false).description().c_str());

    DEBUG0("Storing wavenumber vectors for Fourier bases");

    Eigen::ArrayXc buf(std::max(grid.N.x(), grid.N.z()));

    // Obtain wavenumbers via computing 1*(i*kx)/i
    buf.fill(complex_t(1,0));
    suzerain::diffwave::apply(1, 0, complex_t(0,-1), buf.data(),
            Lx, Lz, 1, grid.N.x(), grid.N.x(), 0, grid.N.x(), 1, 1, 0, 1);
    esio_line_establish(h, grid.N.x(), 0, (procid == 0 ? grid.N.x() : 0));
    esio_line_write(h, "kx", reinterpret_cast<real_t *>(buf.data()),
            2, "Wavenumbers in streamwise X direction"); // Re(buf)

    // Obtain wavenumbers via computing 1*(i*kz)/i
    buf.fill(complex_t(1,0));
    suzerain::diffwave::apply(1, 0, complex_t(0,-1), buf.data(),
            Lx, Lz, 1, 1, 1, 0, 1, grid.N.z(), grid.N.z(), 0, grid.N.z());
    esio_line_establish(h, grid.N.z(), 0, (procid == 0 ? grid.N.z() : 0));
    esio_line_write(h, "kz", reinterpret_cast<real_t *>(buf.data()),
            2, "Wavenumbers in spanwise Z direction"); // Re(buf)
}

void load(const esio_handle h,
          suzerain::problem::GridDefinition& grid)
{
    DEBUG0("Loading GridDefinition parameters");

    esio_line_establish(h, 1, 0, 1); // All ranks load

    if (grid.N.x()) {
        INFO0("Overriding grid using Nx = " << grid.N.x());
    } else {
        int value;
        esio_line_read(h, "Nx", &value, 0);
        grid.Nx(value);
    }

    if (grid.DAF.x()) {
        INFO0("Overriding grid using DAFx = " << grid.DAF.x());
    } else {
        double factor;
        esio_line_read(h, "DAFx", &factor, 0);
        grid.DAFx(factor);
    }

    if (grid.N.y()) {
        INFO0("Overriding grid using Ny = " << grid.N.y());
    } else {
        int value;
        esio_line_read(h, "Ny", &value, 0);
        grid.Ny(value);
    }

    if (grid.k) {
        INFO0("Overriding grid using k = " << grid.k);
    } else {
        esio_line_read(h, "k", &grid.k, 0);
    }

    if (boost::math::signbit(grid.htdelta)) { // (-0.0 "<=" htdelta)
        esio_line_read(h, "htdelta", &grid.htdelta, 0);
    } else {
        INFO0("Overriding grid using htdelta = " << grid.htdelta);
    }

    if (grid.N.z()) {
        INFO0("Overriding grid using Nz = " << grid.N.z());
    } else {
        int value;
        esio_line_read(h, "Nz", &value, 0);
        grid.Nz(value);
    }

    if (grid.DAF.z()) {
        INFO0("Overriding grid using DAFz = " << grid.DAF.z());
    } else {
        double factor;
        esio_line_read(h, "DAFz", &factor, 0);
        grid.DAFz(factor);
    }
}

void create(const int ndof,
            const int k,
            const double a,
            const double b,
            const double htdelta,
            boost::shared_ptr<const suzerain::bspline>& bspw)
{
    // Compute breakpoint locations
    Eigen::ArrayXd breakpoints(ndof + 2 - k);
    suzerain::math::linspace(0.0, 1.0, breakpoints.size(), breakpoints.data());
    for (int i = 0; i < breakpoints.size(); ++i) {
        breakpoints[i] = suzerain_htstretch2(htdelta, 1.0, breakpoints[i]);
    }
    breakpoints = (b - a) * breakpoints + a;

    // Generate the B-spline workspace based on order and breakpoints
    // Maximum non-trivial derivative operators included
    bspw = boost::make_shared<const suzerain::bspline>(
            k, k - 2, breakpoints.size(), breakpoints.data());
    assert(bspw->ndof() == ndof);
}

void store(const esio_handle h,
           const boost::shared_ptr<const suzerain::bspline>& bspw)
{
    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);

    DEBUG0("Storing B-spline knot details");

    Eigen::ArrayXr buf(bspw->nknots());

    bspw->knots(buf.data(), 1);
    esio_line_establish(h, bspw->nknots(),
            0, (procid == 0 ? bspw->nknots() : 0));
    esio_line_write(h, "knots", buf.data(), 0,
            "Knots used to build B-spline basis");

    bspw->breakpoints(buf.data(), 1);
    esio_line_establish(h, bspw->nbreakpoints(),
            0, (procid == 0 ? bspw->nbreakpoints() : 0));
    esio_line_write(h, "breakpoints", buf.data(), 0,
            "Breakpoint locations used to build B-spline basis");

    bspw->collocation_points(buf.data(), 1);
    esio_line_establish(h, bspw->ndof(),
            0, (procid == 0 ? bspw->ndof() : 0));
    esio_line_write(h, "collocation_points", buf.data(), 0,
            "Collocation points used to build discrete operators");

    DEBUG0("Storing B-spline derivative operators");

    char name[8];
    char comment[127];
    for (int k = 0; k <= bspw->nderivatives(); ++k) {
        snprintf(name, sizeof(name)/sizeof(name[0]), "Dy%d", k);
        snprintf(comment, sizeof(comment)/sizeof(comment[0]),
                "Wall-normal derivative Dy%d(i,j) = D%d[j,ku+i-j] for"
                " 0 <= j < n, max(0,j-ku-1) <= i < min(m,j+kl)", k, k);
        const int lda = bspw->ku(k) + 1 + bspw->kl(k);
        esio_plane_establish(h,
                bspw->ndof(), 0, (procid == 0 ? bspw->ndof() : 0),
                lda,          0, (procid == 0 ? lda          : 0));
        esio_plane_write(h, name, bspw->D(k), 0, 0, comment);
        esio_attribute_write(h, name, "kl", bspw->kl(k));
        esio_attribute_write(h, name, "ku", bspw->ku(k));
        esio_attribute_write(h, name, "m",  bspw->ndof());
        esio_attribute_write(h, name, "n",  bspw->ndof());
    }
}

void load(const esio_handle h,
          boost::shared_ptr<const suzerain::bspline>& bspw)
{
    DEBUG0("Loading B-spline workspace based on order and breakpoints");

    // All ranks load B-spline order
    int k;
    esio_line_establish(h, 1, 0, 1);
    esio_line_read(h, "k", &k, 0);

    // htdelta is ignored

    // All ranks load B-spline breakpoints
    int nbreak;
    esio_line_size(h, "breakpoints", &nbreak);
    esio_line_establish(h, nbreak, 0, nbreak);
    Eigen::ArrayXr buf(nbreak);
    esio_line_read(h, "breakpoints", buf.data(), 0);

    // Collocation points are ignored

    // Construct B-spline workspace
    bspw.reset(new suzerain::bspline(k, k - 2, nbreak, buf.data()));
}

void store_time(const esio_handle h,
                real_t time)
{
    // Root writes details
    int rank;
    esio_handle_comm_rank(h, &rank);
    esio_line_establish(h, 1, 0, (rank == 0) ? 1 : 0);

    esio_line_write(h, "t", &time, 0, "Simulation physical time");

    DEBUG0("Stored simulation time " << time);
}

void load_time(const esio_handle h,
               real_t &time)
{
    // All ranks read details
    esio_line_establish(h, 1, 0, 1);

    esio_line_read(h, "t", &time, 0);

    DEBUG0("Loaded simulation time " << time);
}

void store(const esio_handle h,
           const suzerain::NoninterleavedState<4,complex_t> &state,
           const suzerain::problem::GridDefinition& grid,
           const suzerain::pencil_grid& dgrid)
{
    typedef suzerain::NoninterleavedState<4,complex_t> store_type;

    // Ensure state storage meets this routine's assumptions
    assert(                  state.shape()[0]  == field_names.size());
    assert(numeric_cast<int>(state.shape()[1]) == dgrid.local_wave_extent.y());
    assert(numeric_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
    assert(numeric_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());

    // Compute wavenumber translation logistics for X direction
    int fxb[2], fxe[2], mxb[2], mxe[2];
    suzerain::inorder::wavenumber_translate(grid.N.x(),
                                            grid.dN.x(),
                                            dgrid.local_wave_start.x(),
                                            dgrid.local_wave_end.x(),
                                            fxb[0], fxe[0], fxb[1], fxe[1],
                                            mxb[0], mxe[0], mxb[1], mxe[1]);
    // X contains only positive wavenumbers => second range must be empty
    assert(fxb[1] == fxe[1]);
    assert(mxb[1] == mxe[1]);

    // Compute wavenumber translation logistics for Z direction
    // One or both ranges may be empty
    int fzb[2], fze[2], mzb[2], mze[2];
    suzerain::inorder::wavenumber_translate(grid.N.z(),
                                            grid.dN.z(),
                                            dgrid.local_wave_start.z(),
                                            dgrid.local_wave_end.z(),
                                            fzb[0], fze[0], fzb[1], fze[1],
                                            mzb[0], mze[0], mzb[1], mze[1]);

    // Save each scalar field in turn...
    for (size_t i = 0; i < field_names.static_size; ++i) {

        // Create a view of the state for just the i-th scalar
        // Not strictly necessary, but aids greatly in debugging
        boost::multi_array_types::index_range all;
        boost::const_array_view_gen<store_type,3>::type field
            = state[boost::indices[i][all][all][all]];

        // ...which requires two writes per field, once per Z range
        for (int j = 0; j < 2; ++j) {

            // Source of write is NULL for empty WRITE operations
            // Required since MultiArray triggers asserts on invalid indices
            const complex_t * src = NULL;
            if (mxb[0] != mxe[0] && mzb[j] != mze[j]) {
                src = &field[0]
                            [mxb[0] - dgrid.local_wave_start.x()]
                            [mzb[j] - dgrid.local_wave_start.z()];
            }

            // Collectively establish size of read across all ranks
            esio_field_establish(h,
                                 grid.N.z(),     fzb[j], (fze[j] - fzb[j]),
                                 grid.N.x()/2+1, fxb[0], (fxe[0] - fxb[0]),
                                 grid.N.y(),     0,      (     grid.N.y()));

            // Perform collective write operation from state_linear
            complex_field_write(h, field_names[i], src,
                                field.strides()[2],
                                field.strides()[1],
                                field.strides()[0],
                                field_descriptions[i]);
        }
    }
}

void load(const esio_handle h,
          suzerain::NoninterleavedState<4,complex_t> &state,
          const suzerain::problem::GridDefinition& grid,
          const suzerain::pencil_grid& dgrid,
          const suzerain::bspline& bspw)
{
    typedef suzerain::NoninterleavedState<4,complex_t> load_type;

    // Ensure local state storage meets this routine's assumptions
    assert(                  state.shape()[0]  == field_names.size());
    assert(numeric_cast<int>(state.shape()[1]) == dgrid.global_wave_extent.y());
    assert(numeric_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
    assert(numeric_cast<int>(state.shape()[3]) == dgrid.global_wave_extent.z());

    // Obtain details on the restart field's global sizes
    int Fz, Fx, Fy, ncomponents;
    esio_field_sizev(h, field_names[0], &Fz, &Fx, &Fy, &ncomponents);
    assert(ncomponents == 2);

    // Prepare a file-specific B-spline basis
    boost::shared_ptr<const suzerain::bspline> Fbspw;
    load(h, Fbspw);
    assert(Fy == Fbspw->ndof());

    // Check if the B-spline basis in the file differs from ours.  Use strict
    // equality as minor knot differences magnify once collocation points
    // and operators are computed.
    bool bsplines_same =    bspw.order()  == Fbspw->order()
                         && bspw.ndof()   == Fbspw->ndof()
                         && bspw.nknots() == Fbspw->nknots();
    for (int j = 0; bsplines_same && j < bspw.nknots(); ++j) {
        double x_j, y_j;
        bspw.knot(j, &x_j);
        Fbspw->knot(j, &y_j);
        bsplines_same = x_j == y_j;
    }

    // Compute wavenumber translation logistics for X direction.
    // Requires turning a C2R FFT complex-valued coefficient count into a
    // real-valued coefficient count.  Further, need to preserve even- or
    // odd-ness of the coefficient count to handle, for example, Fx == 1.
    int fxb[2], fxe[2], mxb[2], mxe[2];
    suzerain::inorder::wavenumber_translate(2 * (Fx - 1) + (Fx & 1),
                                            grid.dN.x(),
                                            dgrid.local_wave_start.x(),
                                            dgrid.local_wave_end.x(),
                                            fxb[0], fxe[0], fxb[1], fxe[1],
                                            mxb[0], mxe[0], mxb[1], mxe[1]);
    // X contains only positive wavenumbers => second range must be empty
    assert(fxb[1] == fxe[1]);
    assert(mxb[1] == mxe[1]);

    // Compute wavenumber translation logistics for Y direction
    // One or both ranges may be empty
    int fzb[2], fze[2], mzb[2], mze[2];
    suzerain::inorder::wavenumber_translate(Fz,
                                            grid.dN.z(),
                                            dgrid.local_wave_start.z(),
                                            dgrid.local_wave_end.z(),
                                            fzb[0], fze[0], fzb[1], fze[1],
                                            mzb[0], mze[0], mzb[1], mze[1]);

    // Possibly prepare a temporary buffer into which to read each scalar field
    // and a factorization of bspw's mass matrix.  Used only when
    // !bsplines_same.
    typedef boost::multi_array<
        complex_t, 3, suzerain::blas::allocator<complex_t>::type
    > tmp_type;
    boost::scoped_ptr<tmp_type> tmp;
    boost::scoped_ptr<suzerain::bspline_luz> mass;
    if (!bsplines_same) {
        DEBUG0("Differences in B-spline basis require restart projection");
        const boost::array<tmp_type::index,3> extent = {{
            Fy, state.shape()[2], state.shape()[3]
        }};
        tmp.reset(new tmp_type(extent, boost::fortran_storage_order()));
        mass.reset(new suzerain::bspline_luz(bspw));
        mass->form_mass(bspw);
    }

    DEBUG0("Started loading simulation fields");

    // Load each scalar field in turn
    for (size_t i = 0; i < field_names.static_size; ++i) {

        // Create a view of the state for just the i-th scalar
        boost::multi_array_types::index_range all;
        boost::array_view_gen<load_type, 3>::type field
                = state[boost::indices[i][all][all][all]];

        // Clear storage prior to load to zero not-loaded coefficents
        if (bsplines_same) {
            suzerain::multi_array::fill(field, 0);
        } else {
            suzerain::multi_array::fill(*tmp, 0);
        }

        // Two ESIO read operations per field (once per Z range)
        for (int j = 0; j < 2; ++j) {

            // Collectively establish size of read across all ranks
            esio_field_establish(h, Fz, fzb[j], (fze[j] - fzb[j]),
                                        Fx, fxb[0], (fxe[0] - fxb[0]),
                                        Fy,      0, (            Fy));


            // Destination of read is NULL for empty READ operations.
            // Otherwise find starting point for coefficient loading.
            complex_t * dst = NULL;
            boost::array<load_type::index,3> dst_strides = {{ 0, 0, 0 }};
            if (mxb[0] != mxe[0] && mzb[j] != mze[j]) {
                const boost::array<load_type::index,3> index_list = {{
                        0,
                        mxb[0] - dgrid.local_wave_start.x(),
                        mzb[j] - dgrid.local_wave_start.z()
                }};
                if (bsplines_same) {
                    dst = &field(index_list);
                    std::copy(field.strides(), field.strides() + 3,
                              dst_strides.begin());
                } else {
                    dst = &(*tmp)(index_list);
                    std::copy(tmp->strides(), tmp->strides() + 3,
                              dst_strides.begin());
                }
            }

            // Perform collective read operation into dst
            complex_field_read(h, field_names[i], dst,
                               dst_strides[2], dst_strides[1], dst_strides[0]);
        }


        // If necessary, interpolate between B-spline bases.
        // Relies heavily on both bases being collocation-based.
        // This will change dramatically if we ever go the L_2 route.
        if (!bsplines_same) {

            // Step 0: Obtain collocation points for new basis
            Eigen::ArrayXd points(bspw.ndof());
            bspw.collocation_points(points.data(), 1);

            // Step 1: Affine transformation of points into old basis domain
            // Allows stretching of old range onto new range.
            double oldmin, oldmax, newmin, newmax;
            bspw.collocation_point(                  0, &newmin);
            bspw.collocation_point(    bspw.ndof() - 1, &newmax);
            Fbspw->collocation_point(                0, &oldmin);
            Fbspw->collocation_point(Fbspw->ndof() - 1, &oldmax);
            // Only pay for floating point loss if strictly necessary
            if (oldmin != newmin || oldmax != newmax) {
                points = (points - newmin)
                    * ((oldmax - oldmin) / (newmax - newmin)) + oldmin;
            }

            // Step 2: Evaluate old basis + coefficients at points into new
            const int jmax = numeric_cast<int>(field.shape()[1]);
            const int kmax = numeric_cast<int>(field.shape()[2]);
            for (int k = 0; k < kmax; ++k) {
                for (int j = 0; j < jmax; ++j) {
                    Fbspw->zevaluate(0, &(*tmp)[0][j][k], bspw.ndof(),
                                     points.data(), &field[0][j][k], 0);
                }
            }

            // Step 3: Invert for new coefficients using factored mass matrix
            for (int k = 0; k < kmax; ++k) {
                mass->solve(jmax, &field[0][0][k],
                            numeric_cast<int>(field.strides()[0]),
                            numeric_cast<int>(field.strides()[1]));
            }
        }
    }

    DEBUG0("Finished loading simulation fields");
}
