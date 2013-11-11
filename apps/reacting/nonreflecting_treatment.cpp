//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
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
 * @copydoc nonreflecting_treatment.hpp
 */

#include "nonreflecting_treatment.hpp"

#include <suzerain/bspline.hpp>
#include <suzerain/error.h>
#include <suzerain/grid_specification.hpp>
#include <suzerain/inorder.hpp>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/support/logging.hpp>

#include "nonlinear_operator_fwd.hpp"
#include "reacting.hpp"

#pragma float_control(precise, on)
#pragma fenv_access(on)
#pragma float_control(except, on)
#pragma fp_contract(off)
static inline
suzerain::real_t twopiover(const suzerain::real_t L)
{
    return 2*M_PI/L;
}
#pragma float_control(except, off)
#pragma fenv_access(off)
#pragma float_control(precise, off)
#pragma fp_contract(on)

namespace suzerain {

namespace reacting {

nonreflecting_treatment::nonreflecting_treatment(
        const grid_specification &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b,
        operator_common_block &common)
    : operator_base(grid, dgrid, cop, b)
    , common(common)
    , who("nonreflecting_treatment")
{
}

std::vector<real_t> nonreflecting_treatment::apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const lowstorage::method_interface<complex_t> &method,
            const std::size_t substep_index) const
{
    // Implementation approach:
    //   1) Preserve the wavenumber-dependent state at the upper boundary
    //   2) Invoke the wrapped nonlinear operator in the usual fashion
    //   3) Build the various matrices we need from the reference state.
    //   4) Prepare pre-computable products of the various matrices
    //   5) Modify the right hand side in a wave-number dependent fashion.

    // State enters method as coefficients in X, Y, and Z directions

    // Make a local, stride-1 copy of I times upper boundary state point values.
    // Notice that boundary coefficients are 1-1 with boundary point values.
    // That is, applying the mass matrix to the boundary is an ignorable NOP.

    // Number of state variables
    const size_t state_count = swave.shape()[0];
    const size_t Ns          = state_count - 4;

    // Name indices to make (hopefully) the code more legible
    // Local indices for primitive variables U=(rho_s(0,...,Ns-1),u,v,w,p)
    // indices of species densities start at 0 and run through Ns-1
    // index 0 corresponds to the diluter
    const int irhos0 = 0;
    const int iu     = Ns;
    const int iv     = Ns+1;
    const int iw     = Ns+2;
    const int ip     = Ns+3;

    // Local indices for conserved variables V=(rho,rho_s(1,...,Ns-1),rhoU,rhoV,rhoW,rhoE)
    // indices of species densities start at irho+1 and run through irho+Ns-1
    const int irho   = 0;
    const int irhoU  = Ns;
    const int irhoV  = Ns+1;
    const int irhoW  = Ns+2;
    const int irhoE  = Ns+3;

    // Indices for characteristic (ic*) variables
    const int icrhos0 = 0;
    const int icv     = Ns;
    const int icw     = Ns+1;
    const int icu_p   = Ns+2;    // u+a characteristic
    const int icu_m   = Ns+3;    // u-a characteristic

    MatrixXXc i_stash(state_count, swave.shape()[2] * swave.shape()[3]);
    {
        const complex_t imag_unit(0, 1);
        const int ku = boost::numeric_cast<int>(swave.shape()[0]);
        const int l  = boost::numeric_cast<int>(swave.shape()[1]) - 1; // Upper
        const int mu = boost::numeric_cast<int>(swave.shape()[2]);
        const int nu = boost::numeric_cast<int>(swave.shape()[3]);
        for (int k = 0; k < ku; ++k) {
            for (int n = 0; n < nu; ++n) {
                for (int m = 0; m < mu; ++m) {
                    i_stash(k, m + n * mu) = imag_unit * swave[k][l][m][n];
                }
            }
        }
    }

    // Prepare storage for chemistry sources to get added later
    common.chemsrcw.resize(Ns-1, swave.shape()[2] * swave.shape()[3]);
    common.chemsrcw.setZero();

    // Invoke the wrapped nonlinear operator (which may compute references!)
    const std::vector<real_t> retval = N->apply_operator(
            time, swave, method, substep_index);

    // State is now coefficients in X and Z directions
    // State is now collocation point values in Y direction

    // Prepare the rotation and its inverse that reorders from
    // ndx::{e, mx, my, mz, rho, rho_s} to {rho_s, rho, my, mz, mx, e}.
    // All remaining matrices/logic within the routine uses the latter order!
    // FIXME
    MatrixXXr RY(MatrixXXr::Zero(state_count, state_count));
    {
        for (unsigned int is_local = 1; is_local < Ns; ++is_local) {
            RY(irho+is_local, ndx::rho+is_local) = 1;
        }
        RY(irho,  ndx::rho) = 1;
        RY(irhoU, ndx::my ) = 1;
        RY(irhoV, ndx::mz ) = 1;
        RY(irhoW, ndx::mx ) = 1;
        RY(irhoE, ndx::e  ) = 1;
    }
    // FIXME
    MatrixXXr inv_RY(MatrixXXr::Zero(state_count, state_count));
    {
        inv_RY(ndx::e  , irhoE) = 1;
        inv_RY(ndx::mx , irhoW) = 1;
        inv_RY(ndx::my , irhoU) = 1;
        inv_RY(ndx::mz , irhoV) = 1;
        inv_RY(ndx::rho, irho ) = 1;
        for (unsigned int is_local = 1; is_local < Ns; ++is_local) {
            inv_RY(ndx::rho+is_local, irho+is_local) = 1;
        }
    }

    // Ny from grid
    const int Ny   = dgrid.global_wave_extent.y();

    // Retrieve upper boundary reference state from the common block.
    // These assignments are odd looking to permit others as documented.
    const real_t rho   = common.rho_ref;   // swave[ndx::rho][Ny-1][0][0].real(); //FIXME: Check if I get mean(rho) at the upper boundary
    const real_t u     = common.v_ref;     // common.ref_uy   ().tail<1>()[0]; // Ref u = v' per RY
    const real_t v     = common.w_ref;     // common.ref_uz   ().tail<1>()[0]; // Ref v = w' per RY
    const real_t w     = common.u_ref;     // common.ref_ux   ().tail<1>()[0]; // Ref w = u' per RY
    const real_t a     = common.a_ref;     // common.ref_a    ().tail<1>()[0];
    const real_t gamma = common.gamma_ref; // common.ref_gamma().tail<1>()[0];
    const real_t T     = common.T_ref;
    const real_t Cv    = common.Cv_ref;

    VectorXr etots;
    etots.resize(Ns);
    etots = common.etots_ref;

    // Update mean species densities
    // ... only at initial substep
    if (substep_index == 0 ) {
        rhos.resize(Ns);
        // Initialize diluter density to total density
        rhos(0) = rho;
        DEBUG0(who, "rho     = " << rho );
        for (unsigned int is_local = 1; is_local < Ns; ++is_local) {
            rhos(is_local) = rho * common.cs_ref(is_local);
            DEBUG0(who, "rhos(" << is_local << ") = " << rhos(is_local));

            // substract to compute diluter density
            rhos(0) -= rhos(is_local);
        }
        DEBUG0(who, "rhos(0) = " << rhos(0) );
    }


    // Prepare oft-used quantities derived from the reference state
    const real_t inv_rho = 1 / rho;
    const real_t u2      = u * u;
    const real_t v2      = v * v;
    const real_t w2      = w * w;
          real_t inv_a   = 1 / a;
          real_t a2      = a * a;
          real_t inv_a2  = 1 / a2;

    // Prepare oft-used scenario-related constants
    const real_t chi        = dgrid.chi();
    const real_t half       = static_cast<real_t>(1) / 2;
    const real_t inv_gamma  = 1 / gamma;
    const real_t gamma1     = gamma - 1;
    const real_t inv_gamma1 = 1 / gamma1;

    // Build the conserved -> primitive transformation and its inverse
    MatrixXXr S(MatrixXXr::Zero(state_count, state_count));
    {
        S(irhos0, irho) =   1;
        S(iu,     irho) = - u * inv_rho;
        S(iv,     irho) = - v * inv_rho;
        S(iw,     irho) = - w * inv_rho;
        S(ip,     irho) =   (half*(u2+v2+w2)-etots(0)+Cv*T)*gamma1;

        for (unsigned int is_local = 1; is_local < Ns; ++is_local) {
            S(irhos0+is_local, irho+is_local) =  1;
            S(irhos0,          irho+is_local) = -1;
            S(ip,              irho+is_local) =  (-etots(is_local)
                                                  +etots(0))*gamma1;
        }

        S(iu, irhoU) =   inv_rho;
        S(ip, irhoU) = - gamma1 * u;

        S(iv, irhoV) =   inv_rho;
        S(ip, irhoV) = - gamma1 * v;

        S(iw, irhoW) =   inv_rho;
        S(ip, irhoW) = - gamma1 * w;

        S(ip, irhoE) =   gamma1;
    }
    MatrixXXr inv_S(MatrixXXr::Zero(state_count, state_count));
    {
        for (unsigned int is_local = 0; is_local < Ns; ++is_local) {
            inv_S(irho,          irhos0+is_local) = 1;
            inv_S(irho+is_local, irhos0+is_local) = 1;
            inv_S(irhoU,         irhos0+is_local) = u ;
            inv_S(irhoV,         irhos0+is_local) = v ;
            inv_S(irhoW,         irhos0+is_local) = w ;
            inv_S(irhoE,         irhos0+is_local) = half * (u2+v2+w2)
                                       + etots(is_local)
                                       - Cv * T;
        }

        inv_S(irhoU, iu) =  rho;
        inv_S(irhoE, iu) =  rho * u;

        inv_S(irhoV, iv) =  rho;
        inv_S(irhoE, iv) =  rho * v;

        inv_S(irhoW, iw) =  rho;
        inv_S(irhoE, iw) =  rho * w;

        inv_S(irhoE, ip) =  inv_gamma1;
    }

    // Build the primitive -> characteristic transformation and its inverse
    MatrixXXr VL(MatrixXXr::Zero(state_count, state_count));
    {
        for (unsigned int is_local = 0; is_local < Ns; ++is_local) {
            VL(icrhos0+is_local, irhos0+is_local) = -a2;
        }

        VL(icu_p, iu) =   rho * a;
        VL(icu_m, iu) = - rho * a;

        VL(icv,   iv) =   rho * a;

        VL(icw,   iw) =   rho * a;

        for (unsigned int is_local = 0; is_local < Ns; ++is_local) {
            VL(icrhos0+is_local, ip) = rhos(is_local) * inv_rho;
        }
        VL(icu_p, ip) =   1;
        VL(icu_m, ip) =   1;
    }
    MatrixXXr inv_VL(MatrixXXr::Zero(state_count, state_count));
    {
        // First Ns rows
        for (unsigned int is_local = 0; is_local < Ns; ++is_local) {
            inv_VL(irhos0+is_local, icrhos0+is_local) = -         inv_a2;
            inv_VL(irhos0+is_local, icu_p           ) =
                               half * rhos(is_local) * inv_rho * inv_a2;
            inv_VL(irhos0+is_local, icu_m           ) =
                               half * rhos(is_local) * inv_rho * inv_a2;
        }

        // Last 4 rows
        inv_VL(iu, icu_p) =   half * inv_rho * inv_a;
        inv_VL(iu, icu_m) = - half * inv_rho * inv_a;

        inv_VL(iv, icv  ) =   inv_rho * inv_a;

        inv_VL(iw, icw  ) =   inv_rho * inv_a;

        inv_VL(ip, icu_p) =   half;
        inv_VL(ip, icu_m) =   half;
    }

    // The upper boundary is an inflow when the reference velocity is negative.
    // In the event of u == 0, treat the boundary like an outflow.
    const bool inflow = u < 0;

    // Build the in-vs-outflow characteristic-preserving projection
    // FIXME: Incoming waves for inflow or outflow
    //        correspond to SUBSONIC UPPER BOUNDARY
    //        (Generalize implementation for arbitrary application to
    //        lower or upper boundary and supersonic flow?)
    MatrixXXr PG(MatrixXXr::Zero(state_count, state_count));
    if (inflow) {
        for (unsigned int is_local = 0; is_local < Ns; ++is_local) {
            PG(irhos0+is_local, irhos0+is_local) = 1;
        }
        PG(Ns,   Ns  ) = 1;
        PG(Ns+1, Ns+1) = 1;
        PG(Ns+3, Ns+3) = 1;
    } else {
        PG(Ns+3, Ns+3) = 1;
    }

//     MatrixXXr BG(MatrixXXr::Zero(state_count, state_count));  // Medida's B^G_1
    // Declare this matrix for transverse wave approximation,
    // see perfect implementation/documentation, and the original
    // work of Giles and Medida

//     MatrixXXr CG(MatrixXXr::Zero(state_count, state_count));  // Medida's C^G_1
    // Declare this matrix for transverse wave approximation,
    // see perfect implementation/documentation, and the original
    // work of Giles and Medida

    // Prepare all necessary real-valued products of the above matrices
    const MatrixXXr VL_S_RY           = VL * S * RY;
    // FIXME: Uncomment if BG, CG are declared
//     const Matrix5r BG_VL_S_RY_by_chi = BG * VL_S_RY / chi;
//     const Matrix5r CG_VL_S_RY_by_chi = CG * VL_S_RY / chi;
    const MatrixXXr ImPG_VL_S_RY      = (MatrixXXr::Identity(state_count, state_count) - PG) * VL_S_RY;
    const MatrixXXr inv_VL_S_RY       = inv_RY * inv_S * inv_VL;

    // Wavenumber traversal modeled after those found in suzerain/diffwave.c
    // Ny previously declared
    const int Nx   = grid.N.x();
    const int dNx  = grid.dN.x();
    const int dkbx = dgrid.local_wave_start.x();
    const int dkex = dgrid.local_wave_end.x();
    const int Nz   = grid.N.z();
    const int dNz  = grid.dN.z();
    const int dkbz = dgrid.local_wave_start.z();
    const int dkez = dgrid.local_wave_end.z();
    const real_t twopioverLx = twopiover(grid.L.x());  // Weird looking...
    const real_t twopioverLz = twopiover(grid.L.z());  // ...for FP control

    // Traverse wavenumbers updating the RHS with the Giles boundary condition
    using suzerain::inorder::wavenumber;
    const int mu = boost::numeric_cast<int>(swave.shape()[2]);
    VectorXc N(state_count);
    VectorXc tmp(state_count);
    for (int n = dkbz; n < dkez; ++n) {
        const int wn = wavenumber(dNz, n);
        const real_t kn = twopioverLz*wn;

        for (int m = dkbx; m < dkex; ++m) {
            const int wm = wavenumber(dNx, m);
            const real_t km = twopioverLx*wm;

            // Unpack upper boundary RHS into a contiguous buffer
            for (int f = 0; f < N.size(); ++f) {
                N(f) = swave[f][Ny - 1][m - dkbx][n - dkbz];
            }

            // Modify the packed RHS per
            // "Implementation primarily within the nonlinear explicit operator"
            // The imaginary unit has already been included within i_stash.
            tmp.noalias()  = ImPG_VL_S_RY.cast<complex_t>() * N;
            // FIXME: Uncomment if BG, CG are declared
//             tmp.noalias() += kn
//                            * BG_VL_S_RY_by_chi.cast<complex_t>()
//                            * i_stash.col((m - dkbx) + mu*(n - dkbz));
//             tmp.noalias() += km
//                            * CG_VL_S_RY_by_chi.cast<complex_t>()
//                            * i_stash.col((m - dkbx) + mu*(n - dkbz));
            N.noalias()    = inv_VL_S_RY.cast<complex_t>() * tmp;

            // Pack new upper boundary RHS from the contiguous buffer
            for (int f = 0; f < N.size(); ++f) {
                swave[f][Ny - 1][m - dkbx][n - dkbz] = N(f);
            }

            // Add chemistry sources if the bc is an inflow
            if (inflow) {
                assert(Ns >= 1);
                for (unsigned int is = 0; is<Ns-1; ++is) {
                    int f = is + (N.size()-Ns+1);
                    swave[f][Ny - 1][m - dkbx][n - dkbz] +=
                        common.chemsrcw(is, (m - dkbx) + mu*(n - dkbz));
                }
            }
        }
    }

    return retval;
}

} // namespace reacting

} // namespace suzerain
