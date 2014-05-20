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

/** @file
 * @copydoc slowgrowth.hpp
 */

#include <suzerain/slowgrowth.hpp>

#include <largo/largo.h>

#include <suzerain/common.hpp>
#include <suzerain/baseflow.hpp>
#include <suzerain/error.h>
#include <suzerain/largo_formulation.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/operator_tools.hpp>
#include <suzerain/rholut.hpp>
#include <suzerain/specification_largo.hpp>
#include <suzerain/state.hpp>
#include <suzerain/timers.h>

// TODO Implementation is very perfect gas-specific

namespace suzerain {

slowgrowth::slowgrowth(
        const specification_largo& sg,
        const real_t code_Ma)
    : sg(sg)
    , inv_codeMa2(1 / (code_Ma * code_Ma))
    , meanrms(0)
    , meanrms_y(0)
{
}

void
slowgrowth::calculate_baseflow(
    const real_t y,
    largo_state &base,
    largo_state &dy,
    largo_state &dt,
    largo_state &dx,
    largo_state &src) const
{
    // Prepare any necessary baseflow information at the wall
    // assuming that baseflow state already matched code_Ma.
    if (sg.baseflow) {
        sg.baseflow->conserved(y, base.as_is(), dy.as_is(), dx.as_is());
        sg.baseflow->pressure (y, base.p,       dy.p,       dx.p);
    } else {
        base.zero();
        dy.zero();
        dx.zero();
    }

    // When a nontrivial inviscid base flow is present, compute
    // model-appropriate residual so that it might be eradicated.  The base
    // flow is intended to be stationary so any residual is a numerical
    // artifact and not part of the intended problem.
    if (base.trivial()) {

        dt.zero();
        src.zero();

    } else {

        // TODO Reference relevant equations from Largo model document
        //
        // The largo_prestep_baseflow(..., srcbase) parameter
        // requires semi-trivial, model-specific computations.
        // The differences lie in the treatment of slow derivatives.
        Vector3r grad_rho, grad_e, grad_p;
        Matrix3r grad_m;
        if (sg.formulation.is_strictly_temporal()) {

            // Notice &this->basewall == &base works because base was populated
            const real_t scale_u = basewall.u();

            dt.rho = scale_u * dx.rho;      // Modeled slow time
            dt.mx  = scale_u * dx.mx;       // derivative uses
            dt.my  = scale_u * dx.my;       // wall base flow
            dt.mz  = scale_u * dx.mz;
            dt.e   = scale_u * dx.e;

            grad_rho << 0, dy.rho, 0;       // Only account for
            grad_m   << 0, dy.mx,  0,       // wall-normal variation
                        0, dy.my,  0,       // Euler residual below.
                        0, dy.mz,  0;       // Notice spanwise trivial.
            grad_e   << 0, dy.e,   0;
            grad_p   << 0, dy.p,   0;

        } else {

            dt.zero();  // Trivial slow time derivative inherent to formulation

            grad_rho << dx.rho, dy.rho, 0;  // Account for wall-normal
            grad_m   << dx.mx,  dy.mx,  0,  // and slow streamwise
                        dx.my,  dy.my,  0,  // variation in Euler
                        dx.mz,  dy.mz,  0;  // residual below.
            grad_e   << dx.e,   dy.e,   0;  // Notice spanwise trivial.
            grad_p   << dx.p,   dy.p,   0;

        }

        // Compute spatial right hand side for inviscid Euler equations
        // agnostic of model-appropriate derivatives supplied above.
        // Base flow state from local variable "base" inlined below.
        const Vector3r m           (base.mx, base.my, base.mz);
        const real_t   div_m   =   grad_m.trace();
        const real_t   rhs_rho = - div_m;
        const Vector3r u       =   rholut::u(base.rho, m);
        const real_t   div_u   =   rholut::div_u(base.rho, grad_rho,
                                                 m,        div_m);
        const Vector3r rhs_m   = - rholut::div_u_outer_m(m, grad_m,
                                                         u, div_u)
                                 - grad_p * inv_codeMa2;
        const real_t   rhs_e   = - rholut::div_e_u(base.e, grad_e,
                                                   u,      div_u)
                                 - rholut::div_p_u(base.p, grad_p,
                                                   u,      div_u);

        // Compute largo_prestep_baseflow..., srcbase) for Largo
        src.rho = dt.rho - rhs_rho;
        src.mx  = dt.mx  - rhs_m.x();
        src.my  = dt.my  - rhs_m.y();
        src.mz  = dt.mz  - rhs_m.z();
        src.e   = dt.e   - rhs_e;

    }
}

void
slowgrowth::initialize(
        const std::size_t substep_index)
{
    // If necessary, perform initialization calls for Largo
    if (sg.formulation.enabled() && substep_index == 0) {

        SUZERAIN_TIMER_SCOPED("slowgrowth::initialize");

        // Initialize the slow growth workspace
        // Avoids debugging-related memory allocation if at all possible
        SUZERAIN_ENSURE(sg.gramp_mean.size());
        SUZERAIN_ENSURE(sg.gramp_mean.size() == sg.gramp_rms.size());
        if (SUZERAIN_UNLIKELY(sg.ignore_gramp_mean || sg.ignore_gramp_rms)) {
            std::vector<real_t> z(sg.gramp_mean.size(), 0);
            largo_init(sg.workspace,
                       sg.grdelta,
                       sg.ignore_gramp_mean ? &z[0] : &sg.gramp_mean[0],
                       sg.ignore_gramp_rms  ? &z[0] : &sg.gramp_rms [0]);
        } else {
            largo_init(sg.workspace,
                       sg.grdelta,
                       &sg.gramp_mean[0],
                       &sg.gramp_rms[0]);
        }

        // Prepare any necessary baseflow information at the wall
        // As a side-effect, populates this->basewall appropriately
        largo_state dy, dt, dx, src;
        calculate_baseflow(0.0, this->basewall, dy, dt, dx, src);

        // Present the baseflow information to Largo adjusting for code_Ma
        largo_init_wall_baseflow(sg.workspace,
                                    basewall.rescale(inv_codeMa2),
                                    dy      .rescale(inv_codeMa2),
                                    dt      .rescale(inv_codeMa2),
                                    dx      .rescale(inv_codeMa2),
                                    src     .rescale(inv_codeMa2));
    }
}

// TODO Avoid fluctuation computations when not required by formulation
void
slowgrowth::gather_wavexz(
        const operator_tools &otool,
        const contiguous_state<4,complex_t> &swave)
{
    if (sg.formulation.enabled()) {
        SUZERAIN_TIMER_SCOPED("slowgrowth::gather_wavexz");

        // With state being Fourier in X and Z but collocation in Y...
        // ...compute L^2_{xz} of state at each collocation point
        meanrms = compute_field_L2xz(swave, otool.grid, otool.dgrid);

        // ...and rescale to convert to root-mean-square (RMS) fluctuations
        // (mean L2 values are unused so also defensively NaN that storage)
        const real_t rms_adjust = 1 / sqrt(otool.grid.L.x()*otool.grid.L.z());
        for (size_t i = 0; i < meanrms.size(); ++i) {
            meanrms[i].mean.setConstant(
                    std::numeric_limits<real_t>::quiet_NaN());
            meanrms[i].fluctuating *= rms_adjust;
        }
    }
}

slowgrowth::physical_cons::~physical_cons()
{
}

void
slowgrowth::gather_physical_cons(
        const operator_tools &otool,
        const slowgrowth::physical_cons &data)
{
    if (sg.formulation.enabled()) {
        SUZERAIN_TIMER_SCOPED("slowgrowth::gather_physical_cons");

        // Slow growth requires mean conserved state at collocation points.
        // Abuse unused pieces within 'meanrms' to avoid more allocations.
        SUZERAIN_ENSURE(meanrms.size() == 5);
        meanrms[ndx::e  ].mean = data.rhoE();
        meanrms[ndx::mx ].mean = data.rhou();
        meanrms[ndx::my ].mean = data.rhov();
        meanrms[ndx::mz ].mean = data.rhow();
        meanrms[ndx::rho].mean = data.rho();

        // Slow growth requires mean pressure and its fluctuation profiles.
        // Abuse 'meanrms' to store the information after conserved state.
        // The RMS of fluctuating pressure computation can be numerically
        // noisy hence require a non-negative result prior to the sqrt.
        meanrms.resize(meanrms.size() + 1);
        meanrms.back().mean        = data.p();
        meanrms.back().fluctuating = (data.p2() - data.p().square())
                                     .max(0).sqrt();

        // Wall-normal derivatives of every mean RMS quantity are required
        meanrms_y.resize(meanrms.size());
        ArrayX2r tmp(otool.dgrid.global_physical_extent.y(), 2);
        std::vector<field_L2xz>::const_iterator src = meanrms.begin();
        std::vector<field_L2xz>::const_iterator end = meanrms.end();
        std::vector<field_L2xz>::      iterator dst = meanrms_y.begin();
        for (/*just above*/; src != end; ++src, ++dst) {
            tmp.col(0) = (*src).mean;
            tmp.col(1) = (*src).fluctuating;
            otool.masslu()->solve(tmp.cols(), tmp.data(),
                                  tmp.innerStride(), tmp.outerStride());
            (*dst).mean       .resizeLike(tmp.col(0));
            (*dst).fluctuating.resizeLike(tmp.col(1));
            otool.cop.accumulate(1, 1.0, tmp.col(0).data(), tmp.innerStride(),
                                    0.0, (*dst).mean.data(), 1);
            otool.cop.accumulate(1, 1.0, tmp.col(1).data(), tmp.innerStride(),
                                    0.0, (*dst).fluctuating.data(), 1);
        }
    }
}

slowgrowth::physical_rqq::~physical_rqq()
{
}

void
slowgrowth::gather_physical_rqq(
        const operator_tools &otool,
        const physical_rqq &data)
{
    if (sg.formulation.enabled()) {
        SUZERAIN_TIMER_SCOPED("slowgrowth::gather_physical_rqq");

        // Obtain "rqq" values for tensorially-consistent homogenization
        rqq.resize(otool.dgrid.global_physical_extent.y(), NoChange);
        rqq.col(ndx::rho) = data.rho  ();
        rqq.col(ndx::mx ) = data.rhouu();
        rqq.col(ndx::my ) = data.rhovv();
        rqq.col(ndx::mz ) = data.rhoww();
        rqq.col(ndx::e  ) = data.rhoEE();

        // Compute derivatives of these "rqq" values
        rqq_y = rqq;
        ArrayXr tmp;
        for (int i = 0; i < rqq_y.cols(); ++i) {
            tmp = rqq_y.col(i);
            otool.masslu()->solve(tmp.cols(), tmp.data(),
                                  tmp.innerStride(), tmp.outerStride());
            otool.cop.accumulate(1, 1.0, tmp.data(), tmp.innerStride(),
                                    0.0, rqq_y.col(i).data(), 1);
        }
    }
}

void
slowgrowth::inner_y(
        const int j,
        const real_t y_j)
{
    if (sg.formulation.enabled()) {
        SUZERAIN_TIMER_SCOPED("slowgrowth::inner_y");

        // Provide any baseflow-dependent information to Largo
        if (sg.baseflow) {
            largo_state base, dy, dt, dx, src;
            calculate_baseflow(y_j, base, dy, dt, dx, src);
            largo_prestep_baseflow(sg.workspace,
                                   base.rescale(inv_codeMa2),
                                   dy  .rescale(inv_codeMa2),
                                   dt  .rescale(inv_codeMa2),
                                   dx  .rescale(inv_codeMa2),
                                   src .rescale(inv_codeMa2));
        }

        // Repack Y-dependent profiles into a form consumable by Largo
        assert(meanrms.size() == 5 + 1); // State plus pressure
        largo_state mean  (meanrms  [ndx::e  ].mean[j],
                           meanrms  [ndx::mx ].mean[j],
                           meanrms  [ndx::my ].mean[j],
                           meanrms  [ndx::mz ].mean[j],
                           meanrms  [ndx::rho].mean[j],
                           meanrms  .back()   .mean[j]);
        largo_state mean_y(meanrms_y[ndx::e  ].mean[j],
                           meanrms_y[ndx::mx ].mean[j],
                           meanrms_y[ndx::my ].mean[j],
                           meanrms_y[ndx::mz ].mean[j],
                           meanrms_y[ndx::rho].mean[j],
                           meanrms_y.back()   .mean[j]);
        largo_state rms   (meanrms  [ndx::e  ].fluctuating[j],
                           meanrms  [ndx::mx ].fluctuating[j],
                           meanrms  [ndx::my ].fluctuating[j],
                           meanrms  [ndx::mz ].fluctuating[j],
                           meanrms  [ndx::rho].fluctuating[j],
                           meanrms  .back()   .fluctuating[j]);
        largo_state rms_y (meanrms_y[ndx::e  ].fluctuating[j],
                           meanrms_y[ndx::mx ].fluctuating[j],
                           meanrms_y[ndx::my ].fluctuating[j],
                           meanrms_y[ndx::mz ].fluctuating[j],
                           meanrms_y[ndx::rho].fluctuating[j],
                           meanrms_y.back()   .fluctuating[j]);
        largo_state mean_rqq  (rqq  (j, ndx::e  ),    // Notice pressure
                               rqq  (j, ndx::mx ),    // entry is NaN as
                               rqq  (j, ndx::my ),    // it is allegedly
                               rqq  (j, ndx::mz ),    // unused.  This
                               rqq  (j, ndx::rho),    // NaN makes sure.
                               std::numeric_limits<real_t>::quiet_NaN());
        largo_state mean_rqq_y(rqq_y(j, ndx::e  ),    // Ditto re: NaN
                               rqq_y(j, ndx::mx ),
                               rqq_y(j, ndx::my ),
                               rqq_y(j, ndx::mz ),
                               rqq_y(j, ndx::rho),
                               std::numeric_limits<real_t>::quiet_NaN());

        // If requested, have Largo ignore all fluctuations (for debugging)
        if (SUZERAIN_UNLIKELY(sg.ignore_fluctuations)) {
            // Hide RMS fluctuations from Largo
            rms.zero();
            rms_y.zero();

            // Have "rqq" quantities reflect only the mean profile
            mean_rqq.rho   = mean.rho;
            mean_rqq.mx    = mean.mx * mean.u();
            mean_rqq.my    = mean.my * mean.v();
            mean_rqq.mz    = mean.mz * mean.w();
            mean_rqq.e     = mean.e  * mean.E();

            // Have "rqq" derivatives reflect only the mean profile
            mean_rqq_y.rho = mean_y.rho;
            mean_rqq_y.mx  = mean.u()*(2*mean_y.mx - mean.u()*mean_y.rho);
            mean_rqq_y.mx  = mean.v()*(2*mean_y.my - mean.v()*mean_y.rho);
            mean_rqq_y.mx  = mean.w()*(2*mean_y.mz - mean.w()*mean_y.rho);
            mean_rqq_y.e   = mean.E()*(2*mean_y.e  - mean.E()*mean_y.rho);
        }

        // Present the baseflow information to Largo
        largo_prestep_seta_innery(
                    sg.workspace,
                    y_j,
                    mean      .rescale(inv_codeMa2            ),
                    rms       .rescale(inv_codeMa2            ),
                    mean_rqq  .rescale(inv_codeMa2*inv_codeMa2),
                    mean_y    .rescale(inv_codeMa2            ),
                    rms_y     .rescale(inv_codeMa2            ),
                    mean_rqq_y.rescale(inv_codeMa2*inv_codeMa2));
    }
}

void
slowgrowth::inner_xz(
        largo_state &local_state,
        largo_state &slowgrowth_forcing)
{
    if (sg.formulation.enabled()) {
        largo_prestep_seta_innerxz(sg.workspace,
                                   local_state.rescale(inv_codeMa2));
        largo_seta(sg.workspace, 0., 1.,
                   slowgrowth_forcing.rescale(inv_codeMa2));
    } else {
        slowgrowth_forcing.zero();
    }
}

} // namespace suzerain
