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
 * @copydoc quantities.hpp
 */

#include "quantities.hpp"

#include <esio/esio.h>
#include <esio/error.h>

#include <suzerain/blas_et_al.hpp>
#include <suzerain/coalescing_pool.hpp>
#include <suzerain/countof.h>
#include <suzerain/diffwave.hpp>
#include <suzerain/error.h>
#include <suzerain/grid_specification.hpp>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/operator_tools.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/rholut.hpp>
#include <suzerain/shared_range.hpp>
#include <suzerain/state.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/support.hpp>

#include "scenario_definition.hpp"

using boost::numeric_cast;

namespace suzerain {

namespace perfect {

static const char default_who[] = "quantities";

quantities::quantities()
    : t(std::numeric_limits<real_t>::quiet_NaN())
    , who(default_who)
{
    // NOP
}

quantities::quantities(real_t )
    : t(t)
    , who(default_who)
{
    // NOP
}

quantities::quantities(
        real_t t,
        quantities::storage_type::Index Ny)
    : t(t),
      storage(storage_type::Zero(Ny, storage_type::ColsAtCompileTime))
    , who(default_who)
{
    // NOP
}

// Helper for the quantities::save(...) implementation just below
class quantities_saver
{

public:

    quantities_saver(const std::string& who,
                     const esio_handle  esioh,
                     const std::string& prefix)
        : who(who), esioh(esioh), prefix(prefix) {}

    template< typename EigenArray >
    void operator()(const std::string& name, const EigenArray& dat) const
    {
        int procid;
        esio_handle_comm_rank(esioh, &procid);
        esio_field_establish(esioh,
                             1,          0, (procid == 0 ? 1          : 0),
                             dat.cols(), 0, (procid == 0 ? dat.cols() : 0),
                             dat.rows(), 0, (procid == 0 ? dat.rows() : 0));

        std::string key(prefix);
        key.append(name);
        esio_field_write(esioh, key.c_str(), dat.data(), 0, 0, 0,
            "Mean quantity sample stored using row-major indices (B-spline"
            " coefficient, tensor component, sample number) where the"
            " B-spline basis is defined by /Ny, /breakpoints_y, and /knots");
    }

private:

    std::string who;
    esio_handle esioh;
    std::string prefix;

};

void quantities::save(const esio_handle h) const
{
    if (this->storage.size()) {
        quantities_saver f(this->who, h, "bar_");
        this->foreach(f);
    } else {
        WARN0(who, "No mean quantity samples saved--"
                   " trivial storage needs detected");
    }
}

// Helper for the quantities::load(...) implementation just below
class quantities_loader
{

public:

    quantities_loader(const std::string& who,
                      const esio_handle esioh,
                      const std::string& prefix)
        : who(who), esioh(esioh), prefix(prefix) {}

    template< typename EigenArray >
    void operator()(const std::string& name, const EigenArray& dat_) const {

        // http://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html
        EigenArray& dat = const_cast<EigenArray&>(dat_);

        int procid;
        esio_handle_comm_rank(esioh, &procid);

        std::string key(prefix);
        key.append(name);

        int cglobal, bglobal, aglobal;
        if (   ESIO_SUCCESS == esio_field_size(esioh, key.c_str(),
                                               &cglobal, &bglobal, &aglobal)
            && (   EigenArray::ColsAtCompileTime == Dynamic
                || dat.cols() == bglobal)) {
            dat.resize(aglobal, bglobal);
            esio_field_establish(esioh,
                                 1,          0, (procid == 0 ? 1          : 0),
                                 dat.cols(), 0, (procid == 0 ? dat.cols() : 0),
                                 dat.rows(), 0, (procid == 0 ? dat.rows() : 0));
            esio_field_read(esioh, key.c_str(), dat.data(), 0, 0, 0);
        } else {
            WARN0(who, "Unable to load " << key
                  << " for nscalar = " << dat.cols());
            dat.fill(std::numeric_limits<
                     typename EigenArray::Scalar>::quiet_NaN());
        }
    }

private:

    std::string who;
    esio_handle esioh;
    std::string prefix;

};

bool quantities::load(const esio_handle h)
{
    // Were any quantities loaded from file?
    bool success = false;

    // Defensively NaN out all storage in the instance prior to load.
    this->storage.fill(std::numeric_limits<real_t>::quiet_NaN());

    int cglobal, bglobal, aglobal;
    if (ESIO_SUCCESS == esio_field_size(h, "bar_rho",
                                        &cglobal, &bglobal, &aglobal)) {
        this->storage.resize(aglobal, NoChange);
        quantities_loader f(this->who, h, "bar_");
        this->foreach(f);
        success = true;
    } else {
        WARN0(who, "No mean quantity samples loaded--"
                   " unable to anticipate storage needs");
    }

    return success;
}


// This looks like logic from nonlinear_operator.hpp but does not belong there.
// Reading through that file, especially the apply_operator implementation, is
// recommended before reviewing this logic.  This routine is definitely
// suboptimal but is expected to be invoked relatively infrequently.
quantities sample_quantities(
        const scenario_definition &scenario,
        const grid_specification &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        contiguous_state<4,complex_t> &swave,
        const real_t t)
{
    // We are only prepared to handle a fixed number of fields in this routine
    enum { state_count = 5 };

    // Shorthand
    const std::size_t Ny = swave.shape()[1];
    namespace acc = boost::accumulators;
    typedef contiguous_state<4,complex_t> state_type;

    // State enters method as coefficients in X, Y, and Z directions

    // We need auxiliary scalar-field storage.  Prepare logical indices using a
    // struct for scoping (e.g. aux::rho_y).  Ordering will match usage below.
    struct aux { enum {
        e_y,   e_x,   e_z,
        mx_y,  mx_x,  mx_z,
        my_y,  my_x,  my_z,
        mz_y,  mz_x,  mz_z,
        rho_y, rho_x, rho_z,
        count // Sentry
    }; };

    // Obtain the auxiliary storage (likely from a pool to avoid fragmenting).
    // We assume no garbage values in the memory will impact us (for speed).
    scoped_ptr<state_type> _auxw_ptr(
            support::allocate_padded_state<state_type>(
                aux::count, dgrid)); // RAII
    state_type &auxw = *_auxw_ptr;                                 // Shorthand

    // Ensure state storage meets this routine's assumptions

    // Sanity check incoming swave's and auxw's shape and contiguity
    SUZERAIN_ENSURE(                  swave.shape()[0]  == state_count);
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[1]) == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[3]) == dgrid.local_wave_extent.z());
    SUZERAIN_ENSURE((unsigned) swave.strides()[1] == 1u);
    SUZERAIN_ENSURE((unsigned) swave.strides()[2] == swave.shape()[1]);
    SUZERAIN_ENSURE((unsigned) swave.strides()[3] == swave.shape()[1]*swave.shape()[2]);
    SUZERAIN_ENSURE(std::equal(swave.shape() + 1, swave.shape() + 4,
                               auxw.shape() + 1));
    SUZERAIN_ENSURE(std::equal(swave.strides() + 1, swave.strides() + 4,
                               auxw.strides() + 1));

    // Rank-specific details accumulated in ret to be MPI_Reduce-d later
    quantities ret(t, Ny);

    // Obtain samples available in wave-space from mean conserved state.
    // These coefficients are inherently averaged across the X-Z plane.
    if (dgrid.has_zero_zero_modes()) {
        ret.rho_E()        = Map<VectorXc>(swave[ndx::e  ].origin(), Ny).real();
        ret.rho_u().col(0) = Map<VectorXc>(swave[ndx::mx ].origin(), Ny).real();
        ret.rho_u().col(1) = Map<VectorXc>(swave[ndx::my ].origin(), Ny).real();
        ret.rho_u().col(2) = Map<VectorXc>(swave[ndx::mz ].origin(), Ny).real();
        ret.rho()          = Map<VectorXc>(swave[ndx::rho].origin(), Ny).real();
    }

    // Obtain access to helper routines for differentiation
    operator_tools otool(grid, dgrid, cop);

    // Compute Y derivatives of total energy at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    otool.zero_dealiasing_modes(   swave, ndx::e);
    otool.bop_accumulate(1,    1., swave, ndx::e, 0., auxw, aux::e_y);
    otool.bop_apply     (0,    1., swave, ndx::e);

    // Compute X- and Z- derivatives of total energy at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    otool.diffwave_accumulate(1, 0, 1., swave, ndx::e, 0., auxw, aux::e_x);
    otool.diffwave_accumulate(0, 1, 1., swave, ndx::e, 0., auxw, aux::e_z);

    // Compute Y derivatives of X momentum at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    otool.zero_dealiasing_modes(   swave, ndx::mx);
    otool.bop_accumulate(1,    1., swave, ndx::mx, 0., auxw, aux::mx_y);
    otool.bop_apply     (0,    1., swave, ndx::mx);

    // Compute X- and Z- derivatives of X momentum at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    otool.diffwave_accumulate(1, 0, 1., swave, ndx::mx, 0., auxw, aux::mx_x);
    otool.diffwave_accumulate(0, 1, 1., swave, ndx::mx, 0., auxw, aux::mx_z);

    // Compute Y derivatives of Y momentum at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    otool.zero_dealiasing_modes(   swave, ndx::my);
    otool.bop_accumulate(1,    1., swave, ndx::my, 0., auxw, aux::my_y);
    otool.bop_apply     (0,    1., swave, ndx::my);

    // Compute X- and Z- derivatives of Y momentum at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    otool.diffwave_accumulate(1, 0, 1., swave, ndx::my, 0., auxw, aux::my_x);
    otool.diffwave_accumulate(0, 1, 1., swave, ndx::my, 0., auxw, aux::my_z);

    // Compute Y derivatives of Z momentum at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    otool.zero_dealiasing_modes(   swave, ndx::mz);
    otool.bop_accumulate(1,    1., swave, ndx::mz, 0., auxw, aux::mz_y);
    otool.bop_apply     (0,    1., swave, ndx::mz);

    // Compute X- and Z- derivatives of Z momentum at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    otool.diffwave_accumulate(1, 0, 1., swave, ndx::mz, 0., auxw, aux::mz_x);
    otool.diffwave_accumulate(0, 1, 1., swave, ndx::mz, 0., auxw, aux::mz_z);

    // Compute Y derivatives of density at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    otool.zero_dealiasing_modes(   swave, ndx::rho);
    otool.bop_accumulate(1,    1., swave, ndx::rho, 0., auxw, aux::rho_y);
    otool.bop_apply     (0,    1., swave, ndx::rho);

    // Compute X- and Z- derivatives of density at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    otool.diffwave_accumulate(1, 0, 1., swave, ndx::rho,  0., auxw, aux::rho_x);
    otool.diffwave_accumulate(0, 1, 1., swave, ndx::rho,  0., auxw, aux::rho_z);

    // Collectively convert swave and auxw to physical space using parallel
    // FFTs. In physical space, we'll employ views to reshape the 4D row-major
    // (F, Y, Z, X) with contiguous (Y, Z, X) into a 2D (F, Y*Z*X) layout where
    // we know F a priori.  Reducing the dimensionality encourages linear
    // access and eases indexing overhead.
    physical_view<aux::count> auxp(dgrid, auxw);
    physical_view<state_count>sphys(dgrid, swave);
    for (std::size_t i = 0; i < state_count; ++i) {
        dgrid.transform_wave_to_physical(&sphys.coeffRef(i,0));
    }
    for (std::size_t i = 0; i < aux::count; ++i) {
        dgrid.transform_wave_to_physical(&auxp.coeffRef(i,0));
    }

    // Physical space is traversed linearly using a single offset 'offset'.
    // The three loop structure is present to provide the global absolute
    // positions x(i), y(j), and z(k) where necessary.
    for (int offset = 0, j = dgrid.local_physical_start.y();
         j < dgrid.local_physical_end.y();
         ++j) {

        // Type used to accumulate mean quantities versus wall-normal position
        typedef acc::accumulator_set<
                real_t,
#if BOOST_VERSION >= 104700
                acc::stats<acc::tag::sum_kahan>
#else
                acc::stats<acc::tag::sum>
#endif
            > accumulator_type;

        // Accumulators for each mean quantity computed in physical space.
        // For example, quantity "foo" has accumulator "sum_foo".
        // Declared within Y loop so they are reset on each Y iteration.
#define DECLARE(r, data, tuple)                                              \
        accumulator_type BOOST_PP_CAT(sum_,BOOST_PP_TUPLE_ELEM(2, 0, tuple)) \
                [BOOST_PP_TUPLE_ELEM(2, 1, tuple)];
        BOOST_PP_SEQ_FOR_EACH(DECLARE,,
                SUZERAIN_PERFECT_QUANTITIES_PHYSICAL)
#undef DECLARE

        for (int k = dgrid.local_physical_start.z();
            k < dgrid.local_physical_end.z();
            ++k) {

            for (int i = dgrid.local_physical_start.x();
                i < dgrid.local_physical_end.x();
                ++i, /* NB */ ++offset) {

                // Unpack total energy-related quantities
                const real_t e(sphys(ndx::e, offset));
                const Vector3r grad_e(auxp(aux::e_x, offset),
                                      auxp(aux::e_y, offset),
                                      auxp(aux::e_z, offset));

                // Unpack momentum-related quantities
                const Vector3r m(sphys(ndx::mx, offset),
                                 sphys(ndx::my, offset),
                                 sphys(ndx::mz, offset));
                const real_t div_m = auxp(aux::mx_x, offset)
                                   + auxp(aux::my_y, offset)
                                   + auxp(aux::mz_z, offset);
                const Matrix3r grad_m;
                const_cast<Matrix3r&>(grad_m) <<
                                     auxp(aux::mx_x,  offset),
                                     auxp(aux::mx_y,  offset),
                                     auxp(aux::mx_z,  offset),
                                     auxp(aux::my_x,  offset),
                                     auxp(aux::my_y,  offset),
                                     auxp(aux::my_z,  offset),
                                     auxp(aux::mz_x,  offset),
                                     auxp(aux::mz_y,  offset),
                                     auxp(aux::mz_z,  offset);

                // Unpack density-related quantities
                const real_t rho(sphys(ndx::rho, offset));
                const Vector3r grad_rho(auxp(aux::rho_x, offset),
                                        auxp(aux::rho_y, offset),
                                        auxp(aux::rho_z, offset));

                // Compute local quantities based upon state.
                const Vector3r u   = rholut::u(
                                        rho, m);
                const real_t div_u = rholut::div_u(
                                        rho, grad_rho, m, div_m);
                const Matrix3r grad_u = rholut::grad_u(
                                        rho, grad_rho, m, grad_m);

                real_t p, T, mu, lambda;
                Vector3r grad_p, grad_T, grad_mu, grad_lambda;
                rholut::p_T_mu_lambda(
                    scenario.alpha, scenario.beta, scenario.gamma, scenario.Ma,
                    rho, grad_rho, m, grad_m, e, grad_e,
                    p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);

                const Matrix3r tau = rholut::tau(
                                        mu, lambda, div_u, grad_u);
                const Vector3r tau_u = tau * u;

                // Accumulate quantities into sum_XXX using function syntax.
                sum_E[0](e / rho);

                sum_T[0](T);

                sum_mu[0](mu);

                sum_nu[0](mu / rho);

                sum_u[0](u.x());
                sum_u[1](u.y());
                sum_u[2](u.z());

                sum_sym_grad_u[0]( grad_u(0,0)                   );
                sum_sym_grad_u[1]((grad_u(0,1) + grad_u(1,0)) / 2);
                sum_sym_grad_u[2]((grad_u(0,2) + grad_u(2,0)) / 2);
                sum_sym_grad_u[3]( grad_u(1,1)                   );
                sum_sym_grad_u[4]((grad_u(1,2) + grad_u(2,1)) / 2);
                sum_sym_grad_u[5]( grad_u(2,2)                   );

                sum_sym_rho_grad_u[0](rho *  grad_u(0,0)                   );
                sum_sym_rho_grad_u[1](rho * (grad_u(0,1) + grad_u(1,0)) / 2);
                sum_sym_rho_grad_u[2](rho * (grad_u(0,2) + grad_u(2,0)) / 2);
                sum_sym_rho_grad_u[3](rho *  grad_u(1,1)                   );
                sum_sym_rho_grad_u[4](rho * (grad_u(1,2) + grad_u(2,1)) / 2);
                sum_sym_rho_grad_u[5](rho *  grad_u(2,2)                   );

                sum_grad_T[0](grad_T.x());
                sum_grad_T[1](grad_T.y());
                sum_grad_T[2](grad_T.z());

                sum_rho_grad_T[0](rho * grad_T.x());
                sum_rho_grad_T[1](rho * grad_T.y());
                sum_rho_grad_T[2](rho * grad_T.z());

                sum_tau_colon_grad_u[0]((tau.transpose()*grad_u).trace());

                sum_tau[0](tau(0,0));
                sum_tau[1](tau(0,1));
                sum_tau[2](tau(0,2));
                sum_tau[3](tau(1,1));
                sum_tau[4](tau(1,2));
                sum_tau[5](tau(2,2));

                sum_tau_u[0](tau_u.x());
                sum_tau_u[1](tau_u.y());
                sum_tau_u[2](tau_u.z());

                sum_p_div_u[0](p*div_u);

                sum_rho_u_u[0](rho * u.x() * u.x());
                sum_rho_u_u[1](rho * u.x() * u.y());
                sum_rho_u_u[2](rho * u.x() * u.z());
                sum_rho_u_u[3](rho * u.y() * u.y());
                sum_rho_u_u[4](rho * u.y() * u.z());
                sum_rho_u_u[5](rho * u.z() * u.z());

                sum_rho_u_u_u[0](rho * u.x() * u.x() * u.x());
                sum_rho_u_u_u[1](rho * u.x() * u.x() * u.y());
                sum_rho_u_u_u[2](rho * u.x() * u.x() * u.z());
                sum_rho_u_u_u[3](rho * u.x() * u.y() * u.y());
                sum_rho_u_u_u[4](rho * u.x() * u.y() * u.z());
                sum_rho_u_u_u[5](rho * u.x() * u.z() * u.z());
                sum_rho_u_u_u[6](rho * u.y() * u.y() * u.y());
                sum_rho_u_u_u[7](rho * u.y() * u.y() * u.z());
                sum_rho_u_u_u[8](rho * u.y() * u.z() * u.z());
                sum_rho_u_u_u[9](rho * u.z() * u.z() * u.z());

                sum_rho_T_u[0](rho * T * u.x());
                sum_rho_T_u[1](rho * T * u.y());
                sum_rho_T_u[2](rho * T * u.z());

                sum_rho_mu[0](rho * mu);

                sum_mu_S[0](mu * ( grad_u(0,0)                    - div_u / 3));
                sum_mu_S[1](mu * ((grad_u(0,1) + grad_u(1,0)) / 2            ));
                sum_mu_S[2](mu * ((grad_u(0,2) + grad_u(2,0)) / 2            ));
                sum_mu_S[3](mu * ( grad_u(1,1)                    - div_u / 3));
                sum_mu_S[4](mu * ((grad_u(1,2) + grad_u(2,1)) / 2            ));
                sum_mu_S[5](mu * ( grad_u(2,2)                    - div_u / 3));

                sum_mu_div_u[0](mu * div_u);

                sum_mu_grad_T[0](mu * grad_T.x());
                sum_mu_grad_T[1](mu * grad_T.y());
                sum_mu_grad_T[2](mu * grad_T.z());

                // TODO Sum mean slow growth forcing contributions (Redmine #2496)

                sum_SrhoE[0](0);

                sum_Srhou[0](0);
                sum_Srhou[1](0);
                sum_Srhou[2](0);

                sum_Srho[0](0);

                sum_Srhou_dot_u[0](0);

            } // end X

        } // end Z

        // Move y-specific sums into MPI-reduction-ready storage for y(j) using
        // Eigen comma initialization syntax.  Yes, this is not stride 1.
#define EXTRACT_SUM(z, n, q) acc::sum(sum_##q[n])
#define MOVE_SUM_INTO_TMP(r, data, tuple)                                 \
        ret.BOOST_PP_TUPLE_ELEM(2, 0, tuple)().row(j) <<                  \
            BOOST_PP_ENUM(BOOST_PP_TUPLE_ELEM(2, 1, tuple),               \
                          EXTRACT_SUM, BOOST_PP_TUPLE_ELEM(2, 0, tuple));
        BOOST_PP_SEQ_FOR_EACH(MOVE_SUM_INTO_TMP,,
                SUZERAIN_PERFECT_QUANTITIES_PHYSICAL)
#undef EXTRACT_SUM
#undef MOVE_SUM_INTO_TMP

    } // end Y

    // Notice dgrid.rank_zero_zero_modes already contains "wave-sampled"
    // quantities while other ranks have zeros in those locations.

    // Reduce sums onto rank zero and then return garbage from non-zero ranks
    if (mpi::comm_rank(MPI_COMM_WORLD) == 0) {

        // Reduce operation requires no additional storage on rank-zero
        SUZERAIN_MPICHKR(MPI_Reduce(
                MPI_IN_PLACE, ret.storage.data(), ret.storage.size(),
                mpi::datatype<quantities::storage_type::Scalar>::value,
                MPI_SUM, /* root */ 0, MPI_COMM_WORLD));

    } else {

        // Reduce operation requires temporary storage on non-zero ranks
        ArrayXXr tmp;
        tmp.resizeLike(ret.storage);
        tmp.setZero();
        SUZERAIN_MPICHKR(MPI_Reduce(ret.storage.data(), tmp.data(), tmp.size(),
                mpi::datatype<quantities::storage_type::Scalar>::value,
                MPI_SUM, /* root */ 0, MPI_COMM_WORLD));

        // Force non-zero ranks contain all NaNs to help detect usage errors
        ret.storage.fill(std::numeric_limits<
                quantities::storage_type::Scalar>::quiet_NaN());

        // Return from all non-zero ranks
        return ret;

    }

    // Only rank zero reaches this logic because of return statement just above
    assert(mpi::comm_rank(MPI_COMM_WORLD) == 0);

    // Physical space sums, which are at collocation points, need to be
    // divided by the dealiased extents and converted to coefficients.
    const real_t scale_factor = 1 / dgrid.chi();
    bsplineop_lu scaled_mass(cop);
    scaled_mass.opform(1, &scale_factor, cop);
    scaled_mass.factor();
    scaled_mass.solve(quantities::nscalars::physical,
            ret.storage.middleCols<quantities::nscalars::physical>(
                quantities::nscalars::wave).data(),
            ret.storage.innerStride(), ret.storage.outerStride());

    // Fill with NaNs those samples which were not computed by this method
#define FILL(r, data, tuple)                                         \
    ret.BOOST_PP_TUPLE_ELEM(2, 0, tuple)().fill(std::numeric_limits< \
            quantities::storage_type::Scalar>::quiet_NaN());
    BOOST_PP_SEQ_FOR_EACH(FILL,,SUZERAIN_PERFECT_QUANTITIES_IMPLICIT)
#undef FILL

    return ret;
}

} // namespace perfect

} // namespace suzerain
