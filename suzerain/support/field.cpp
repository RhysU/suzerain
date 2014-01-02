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
 * @copydoc field.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/support/field.hpp>

#include <esio/esio.h>
#include <esio/error.h>

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/error.h>
#include <suzerain/inorder.hpp>
#include <suzerain/operator_tools.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/state.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/support.hpp>

using boost::numeric_cast;
using std::size_t;

namespace suzerain {

namespace support {

/** Read a complex-valued field via ESIO */
static inline
void complex_field_read(esio_handle h, const char *name, complex_t *field,
                        int cstride = 0, int bstride = 0, int astride = 0)
{
    esio_field_readv(h, name, reinterpret_cast<real_t *>(field),
                     2*boost::numeric_cast<int>(cstride),
                     2*boost::numeric_cast<int>(bstride),
                     2*boost::numeric_cast<int>(astride),
                     2);
}

/** Write a complex-valued field via ESIO */
static inline
void complex_field_write(esio_handle h,
                         const char *name, const complex_t *field,
                         int cstride = 0, int bstride = 0, int astride = 0,
                         const char * comment = 0)
{
    esio_field_writev(h, name, reinterpret_cast<const real_t *>(field),
                      2*boost::numeric_cast<int>(cstride),
                      2*boost::numeric_cast<int>(bstride),
                      2*boost::numeric_cast<int>(astride),
                      2, comment);
}

void save_coefficients(
        const esio_handle h,
        const std::vector<field> &fields,
        const contiguous_state<4,complex_t> &swave,
        const grid_specification& grid,
        const pencil_grid& dgrid)
{
    // Ensure swave meets this routine's assumptions
    SUZERAIN_ENSURE(                  swave.shape()[0]  == fields.size()              );
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[1]) == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[3]) == dgrid.local_wave_extent.z());

    // Compute wavenumber translation logistics for X direction
    int fxb[2], fxe[2], mxb[2], mxe[2];
    inorder::wavenumber_translate(grid.N.x(),
                                  grid.dN.x(),
                                  dgrid.local_wave_start.x(),
                                  dgrid.local_wave_end.x(),
                                  fxb[0], fxe[0], fxb[1], fxe[1],
                                  mxb[0], mxe[0], mxb[1], mxe[1]);
    // X contains only positive wavenumbers => second range must be empty
    SUZERAIN_ENSURE(fxb[1] == fxe[1]);
    SUZERAIN_ENSURE(mxb[1] == mxe[1]);

    // Compute wavenumber translation logistics for Z direction
    // One or both ranges may be empty
    int fzb[2], fze[2], mzb[2], mze[2];
    inorder::wavenumber_translate(grid.N.z(),
                                  grid.dN.z(),
                                  dgrid.local_wave_start.z(),
                                  dgrid.local_wave_end.z(),
                                  fzb[0], fze[0], fzb[1], fze[1],
                                  mzb[0], mze[0], mzb[1], mze[1]);

    // Save each scalar field in turn...
    for (size_t i = 0; i < fields.size(); ++i) {

        // ...first generate a metadata comment...
        std::string comment = fields[i].description;
        comment += " stored row-major ZXY using a Fourier basis in Z stored"
                   " in-order per /kz; a Fourier basis in X stored using"
                   " Hermitian symmetry per /kx; and a B-spline basis in Y"
                   " defined by /k, /breakpoints_y, and /knots";

        // ...followed by two collective writes per field (once per Z range)
        for (int j = 0; j < 2; ++j) {

            // Collectively establish size of write across all ranks
            esio_field_establish(h,
                                 grid.N.z(),     fzb[j], (fze[j] - fzb[j]),
                                 grid.N.x()/2+1, fxb[0], (fxe[0] - fxb[0]),
                                 grid.N.y(),     0,      (     grid.N.y()));

            // Source of write is NULL for empty WRITE operations
            // Required since MultiArray triggers asserts on invalid indices
            const complex_t * src = NULL;
            if (mxb[0] != mxe[0] && mzb[j] != mze[j]) {
                src = &swave[i][0][mxb[0] - dgrid.local_wave_start.x()]
                                  [mzb[j] - dgrid.local_wave_start.z()];
            }

            // Perform collective write operation
            complex_field_write(h, fields[i].location.c_str(), src,
                                swave.strides()[3],
                                swave.strides()[2],
                                swave.strides()[1],
                                comment.c_str());
        }
    }
}

void load_coefficients(const esio_handle h,
                       const std::vector<field> &fields,
                       contiguous_state<4,complex_t> &state,
                       const grid_specification& grid,
                       const pencil_grid& dgrid,
                       const bsplineop& cop,
                       const bspline& b)
{
    typedef contiguous_state<4,complex_t> load_type;

    // Ensure local state storage meets this routine's assumptions
    SUZERAIN_ENSURE(                  state.shape()[0]  == fields.size());
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[1]) == dgrid.global_wave_extent.y());
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());

    // Obtain details on the restart field's global sizes
    int Fz, Fx, Fy, ncomponents;
    esio_field_sizev(h, fields[0].location.c_str(), &Fz, &Fx, &Fy, &ncomponents);
    SUZERAIN_ENSURE(ncomponents == 2);

    // Prepare a file-specific B-spline basis
    shared_ptr<bspline> Fb;
    shared_ptr<bsplineop> Fbop;
    load(h, Fb, Fbop);
    SUZERAIN_ENSURE(Fy == Fb->n());

    // Check if the B-spline basis in the file differs from ours.
    const double bsp_dist = b.distance_to(*Fb);
    const bool   bsp_same = bsp_dist < suzerain_bspline_distance_distinct;

    // Compute wavenumber translation logistics for X direction.
    // Requires turning a C2R FFT complex-valued coefficient count into a
    // real-valued coefficient count.  Further, need to preserve even- or
    // odd-ness of the coefficient count to handle, for example, Fx == 1.
    int fxb[2], fxe[2], mxb[2], mxe[2];
    inorder::wavenumber_translate(2 * (Fx - 1) + (Fx & 1),
                                  grid.dN.x(),
                                  dgrid.local_wave_start.x(),
                                  dgrid.local_wave_end.x(),
                                  fxb[0], fxe[0], fxb[1], fxe[1],
                                  mxb[0], mxe[0], mxb[1], mxe[1]);
    // X contains only positive wavenumbers => second range must be empty
    SUZERAIN_ENSURE(fxb[1] == fxe[1]);
    SUZERAIN_ENSURE(mxb[1] == mxe[1]);

    // Compute wavenumber translation logistics for Z direction
    // One or both ranges may be empty
    int fzb[2], fze[2], mzb[2], mze[2];
    inorder::wavenumber_translate(Fz,
                                  grid.dN.z(),
                                  dgrid.local_wave_start.z(),
                                  dgrid.local_wave_end.z(),
                                  fzb[0], fze[0], fzb[1], fze[1],
                                  mzb[0], mze[0], mzb[1], mze[1]);

    // Possibly prepare a tmp buffer into which to read each scalar field and a
    // factorization of b's mass matrix.  Used only when !bsp_same.
    typedef boost::multi_array<
        complex_t, 3, blas::allocator<complex_t>::type
    > tmp_type;
    scoped_ptr<tmp_type> tmp;
    scoped_ptr<bsplineop_luz> mass;
    if (!bsp_same) {
        INFO0("Differences in B-spline basis require restart projection ("
              << bsp_dist << " >= "
              << suzerain_bspline_distance_distinct << ")");
        const array<tmp_type::size_type,3> extent = {{
            static_cast<tmp_type::size_type>(Fy),
            state.shape()[2],
            state.shape()[3]
        }};
        tmp.reset(new tmp_type(extent, boost::fortran_storage_order()));
        mass.reset(new bsplineop_luz(cop));
        mass->factor_mass(cop);
    }

    // Load each scalar field in turn
    for (size_t i = 0; i < fields.size(); ++i) {

        // Create a view of the state for just the i-th scalar
        boost::multi_array_types::index_range all;
        boost::array_view_gen<load_type, 3>::type field
                = state[boost::indices[i][all][all][all]];

        // Clear storage prior to load to zero not-loaded coefficients
        if (bsp_same) {
            multi_array::fill(field, 0);
        } else {
            multi_array::fill(*tmp, 0);
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
            array<load_type::index,3> dst_strides = {{ 0, 0, 0 }};
            if (mxb[0] != mxe[0] && mzb[j] != mze[j]) {
                const array<load_type::index,3> index_list = {{
                        0,
                        mxb[0] - dgrid.local_wave_start.x(),
                        mzb[j] - dgrid.local_wave_start.z()
                }};
                if (bsp_same) {
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
            complex_field_read(h, fields[i].location.c_str(), dst,
                               dst_strides[2], dst_strides[1], dst_strides[0]);
        }


        // If necessary, interpolate between B-spline bases.
        // Relies heavily on both bases being collocation-based.
        // This will change dramatically if we ever go the L_2 route.
        if (!bsp_same) {

            // Step 0: Obtain collocation points for new basis
            ArrayXr points(b.n());
            for (int k = 0; k < b.n(); ++k) points[k] = b.collocation_point(k);

            // Step 1: Affine transformation of points into old basis domain
            // Allows stretching of old range onto new range.
            const double newmin = b.collocation_point(0);
            const double newmax = b.collocation_point(b.n() - 1);
            const double oldmin = Fb->collocation_point(0);
            const double oldmax = Fb->collocation_point(Fb->n() - 1);
            // Only pay for floating point loss if strictly necessary
#pragma warning(push,disable:1572)
            if (oldmin != newmin || oldmax != newmax) {
#pragma warning(pop)
                points = (points - newmin)
                    * ((oldmax - oldmin) / (newmax - newmin)) + oldmin;
            }

            // Step 2: Evaluate old basis + coefficients at points into new
            const int jmax = numeric_cast<int>(field.shape()[1]);
            const int kmax = numeric_cast<int>(field.shape()[2]);
            for (int k = 0; k < kmax; ++k) {
                for (int j = 0; j < jmax; ++j) {
                    Fb->linear_combination(0, &(*tmp)[0][j][k], b.n(),
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
}

void save_collocation_values(
        const esio_handle h,
        const std::vector<field> &fields,
        contiguous_state<4,complex_t>& swave,
        const grid_specification& grid,
        const pencil_grid& dgrid,
        const bsplineop& cop)
{
    // Ensure state storage meets this routine's assumptions
    SUZERAIN_ENSURE(                  swave.shape()[0]  == fields.size());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[1]) == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[3]) == dgrid.local_wave_extent.z());

    // Convert coefficients into collocation point values
    // Transforms state from full-wave coefficients to full-physical points
    operator_tools otool(grid, dgrid, cop);
    for (size_t i = 0; i < swave.shape()[0]; ++i) {
        otool.bop_apply(0, 1, swave, i);
        otool.zero_dealiasing_modes(swave, i);
        dgrid.transform_wave_to_physical(
                reinterpret_cast<real_t *>(swave[i].origin()));
    }

    // Establish size of collective writes across all ranks
    esio_field_establish(h, grid.dN.y(), dgrid.local_physical_start.y(),
                                         dgrid.local_physical_extent.y(),
                            grid.dN.z(), dgrid.local_physical_start.z(),
                                         dgrid.local_physical_extent.z(),
                            grid.dN.x(), dgrid.local_physical_start.x(),
                                         dgrid.local_physical_extent.x());

    // Write each field in turn
    for (size_t i = 0; i < fields.size(); ++i) {
        std::string comment = fields[i].description;
        comment += " stored row-major YZX on the 3D rectilinear grid defined"
                   " by taking the outer product of arrays"
                   " /collocation_points_y, /collocation_points_z, and"
                   " /collocation_points_z";

        esio_field_write(h, fields[i].location.c_str(),
                reinterpret_cast<real_t *>(swave[i].origin()),
                0, 0, 0, comment.c_str());
    }
}

void load_collocation_values(
        const esio_handle h,
        const std::vector<field> &fields,
        contiguous_state<4,complex_t>& state,
        const grid_specification& grid,
        const pencil_grid& dgrid,
        const bsplineop& cop,
        const bspline& b)
{
    // Ensure state storage meets this routine's assumptions
    SUZERAIN_ENSURE(                  state.shape()[0]  == fields.size());
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[1]) == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());

    // This routine does no grid interpolation.  Yell loudly if necessary
    {
        // Check that restart file size matches runtime dealiased extents
        int cg, bg, ag;
        if (    ESIO_SUCCESS
            != esio_field_size(h, fields[0].location.c_str(), &cg, &bg, &ag)) {
            SUZERAIN_ERROR_VOID("Unable to find fields[0] size from restart",
                                SUZERAIN_EFAILED);
        }
        if (cg != grid.dN.y() || bg != grid.dN.z() || ag != grid.dN.x()) {
            ERROR0("Physical space restart fields have row-major YZX extents "
                   << "(" << cg << "," << bg << "," << ag << ")" << " but "
                   << "(" << grid.dN.y() << "," << grid.dN.z() << ","
                   << grid.dN.x() << ") are required");
            SUZERAIN_ERROR_VOID(
                    "Cannot interpolate during physical space restart",
                    SUZERAIN_EFAILED);
        }

        // Check that restart file specifies the same B-spline basis.
        // TODO Too restrictive?  Any floating point differences kill us.
        shared_ptr<bspline> Fb;
        shared_ptr<bsplineop> Fbop;
        support::load(h, Fb, Fbop);
        const double bsp_dist = b.distance_to(*Fb);
        const bool   bsp_same = bsp_dist < suzerain_bspline_distance_distinct;
        if (!bsp_same) {
            ERROR0("Physical restart has different wall-normal bases ("
                   << bsp_dist << ")");
            SUZERAIN_ERROR_VOID(
                    "Cannot interpolate during physical space restart",
                    SUZERAIN_EFAILED);
        }
    }

    // Prepare physical space view of the state storage
    physical_view<> sphys(dgrid, state);

    // Establish size of collective reads across all ranks
    esio_field_establish(h, grid.dN.y(), dgrid.local_physical_start.y(),
                                         dgrid.local_physical_extent.y(),
                            grid.dN.z(), dgrid.local_physical_start.z(),
                                         dgrid.local_physical_extent.z(),
                            grid.dN.x(), dgrid.local_physical_start.x(),
                                         dgrid.local_physical_extent.x());

    // Read each field from collocation values on disk into the same in memory
    for (size_t i = 0; i < fields.size(); ++i) {
        esio_field_read(h, fields[i].location.c_str(), &sphys(i,0), 0, 0, 0);
    }

    // Collectively convert physical state to wave space coefficients
    // Build FFT normalization constant into Y direction's mass matrix
    bsplineop_luz massluz(cop);
    const complex_t scale_factor = 1 / dgrid.chi();
    massluz.opform(1, &scale_factor, cop);
    massluz.factor();
    operator_tools otool(grid, dgrid, cop);
    for (size_t i = 0; i < state.shape()[0]; ++i) {
        dgrid.transform_physical_to_wave(&sphys.coeffRef(i, 0)); // X, Z
        otool.zero_dealiasing_modes(state, i);
        otool.bop_solve(massluz, state, i);                      // Y
    }
}

} // end namespace suzerain

} // end namespace support
