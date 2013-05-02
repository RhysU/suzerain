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

#include "nonlinear_operator_fwd.hpp"
#include "reacting.hpp"
// #include "scenario_definition.hpp"

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
    // NOP
}

std::vector<real_t> nonreflecting_treatment::apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
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
    // Local indices for primitive variables U=(rho_s(1,...,Ns),u,v,w,p)
    const int irhos0 = 0;       // indices of species densities start at 0 
                                // and run through Ns-1
    const int iu     = Ns;
    const int iv     = Ns+1;
    const int iw     = Ns+2;
    const int ip     = Ns+3;

    // Local indices for conserved variables V=(rho_s(1,...,Ns-1),rho,rhoU,rhoV,rhoW,rhoE)
    //const int irhos0 = 0;     // indices of species densities start at 0 
                                // and run through Ns-2
    const int irho   = Ns-1;
    const int irhoU  = Ns;
    const int irhoV  = Ns+1;
    const int irhoW  = Ns+2;
    const int irhoE  = Ns+3;

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

    // Invoke the wrapped nonlinear operator (which may compute references!)
    const std::vector<real_t> retval = N->apply_operator(
            time, swave, evmaxmag_real, evmaxmag_imag, substep_index);

    // State is now coefficients in X and Z directions
    // State is now collocation point values in Y direction

    // Prepare the rotation and its inverse that reorders from
    // ndx::{e, mx, my, mz, rho, rho_s} to {rho_s, rho, my, mz, mx, e}.
    // All remaining matrices/logic within the routine uses the latter order!
    // FIXME
    MatrixXXr RY(MatrixXXr::Zero(state_count, state_count));
    {
        for (unsigned int is_local = 0; is_local < Ns-1; ++is_local) {
            RY(irhos0+is_local, ndx::rho+1+is_local) = 1;
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
        inv_RY(ndx::mx , irhoV) = 1;
        inv_RY(ndx::my , irhoW) = 1;
        inv_RY(ndx::mz , irhoU) = 1;
        inv_RY(ndx::rho, irho ) = 1;
        for (unsigned int is_local = 0; is_local < Ns-1; ++is_local) {
            RY(ndx::rho+1+is_local, irhos0+is_local) = 1;
        }
    }

    // Ny from grid
    const int Ny   = dgrid.global_wave_extent.y();

    // Retrieve upper boundary reference state from the common block.
    // These assignments are odd looking to permit others as documented.
    const real_t rho   = swave[ndx::rho][Ny-1][0][0].real();
    const real_t u     = 0; // common.ref_uy   ().tail<1>()[0]; // Ref u = v' per RY
    const real_t v     = 0; // common.ref_uz   ().tail<1>()[0]; // Ref v = w' per RY
    const real_t w     = 0; // common.ref_ux   ().tail<1>()[0]; // Ref w = u' per RY
          real_t a     = common.ref_a    ().tail<1>()[0];
          real_t gamma = common.ref_gamma().tail<1>()[0];

    // Update mean species densities
    // ... only at initial substep
    if (substep_index == 0 ) {    
        rhos.resize(Ns);
        // Initialize diluter density to total density
        rhos(irhos0+Ns-1) = rho;
        for (unsigned int is_local = 0; is_local < Ns-1; ++is_local) {
            // density of species is_local
            rhos[irhos0+is_local] = 
                swave[ndx::rho+1+is_local][Ny-1][0][0].real();
            // substract for diluter density
            rhos[irhos0+Ns-1] -= rhos[irhos0+is_local]; 
        }
    }

    SUZERAIN_MPICHKR(MPI_Bcast(&rhos, Ns,
                mpi::datatype_of(rhos[0]),
                dgrid.rank_zero_zero_modes, MPI_COMM_WORLD));

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
        for (unsigned int is_local = 0; is_local < Ns-1; ++is_local) {
            S(irhos0+is_local, irhos0+is_local) =   1;
            S(irho,            irhos0+is_local) =  -1;
            S(ip,              irhos0+is_local) =   (-common.etots_upper(is_local+1)
                                                     +common.etots_upper(0))*gamma1;
        }

        S(irho, irho) =   1;
        S(iu,   irho) = - u * inv_rho;
        S(iv,   irho) = - v * inv_rho;
        S(iw,   irho) = - w * inv_rho;
        S(ip,   irho) =   (half*(u2+v2+w2)-common.etots_upper(0))*gamma1;

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
            inv_S(irhos0+is_local, irhos0+is_local) = 1;
            inv_S(irho,            irhos0+is_local) = 1;
            inv_S(irhoU,           irhos0+is_local) = u ;
            inv_S(irhoV,           irhos0+is_local) = v ;
            inv_S(irhoW,           irhos0+is_local) = w ;
            inv_S(irhoE,           irhos0+is_local) = half * (u2+v2+w2)
                                       + common.etots_upper((is_local+1)%Ns);
            // The energy for the diluter is in index 0 in etots, 
            // and it's index Ns-1 in inv_S.
            // Use the mod operation to shift indices
            // FIXME: make things more consistent with the rest of the code
        }

        inv_S(irho,  iu) =   rho;
        inv_S(irhoE, iu) =   rho * u;

        inv_S(irho,  iv) =   rho;
        inv_S(irhoE, iv) =   rho * v;

        inv_S(irho,  iw) =   rho;
        inv_S(irhoE, iw) =   rho * w;

        inv_S(irhoE, ip) =   inv_gamma1;
    }

    // Build the primitive -> characteristic transformation and its inverse
    MatrixXXr VL(MatrixXXr::Zero(state_count, state_count));
    {
        for (unsigned int is_local = 0; is_local < Ns; ++is_local) {
            VL(irhos0+is_local, irhos0+is_local) = -a2;
        }

        VL(iw, iu) =   rho * a;
        VL(ip, iu) = - rho * a;

        VL(iu, iv) =   rho * a;

        VL(iv, iw) =   rho * a;

        for (unsigned int is_local = 0; is_local < Ns; ++is_local) {
            VL(irhos0+is_local, ip) = -rhos[is_local]/rho;
        }
        VL(iw, ip) =   1;
        VL(ip, ip) =   1;
    }
    MatrixXXr inv_VL(MatrixXXr::Zero(state_count, state_count));
    {
        // First Ns rows
        for (unsigned int is_local = 0; is_local < Ns; ++is_local) {
            inv_VL(irhos0+is_local, irhos0+is_local) = 
                                    - rhos[is_local] * inv_rho * inv_a2;
            inv_VL(irhos0+is_local, Ns+2           ) =   
                                      half * rhos[is_local] * inv_rho * inv_a2;
            inv_VL(irhos0+is_local, Ns+3           ) =   
                                      half * rhos[is_local] * inv_rho * inv_a2;
        }

        // Last 4 rows
        inv_VL(Ns,   Ns+2) =   half * inv_rho * inv_a;
        inv_VL(Ns,   Ns+3) = - half * inv_rho * inv_a;
      
        inv_VL(Ns+1, Ns  ) =   inv_rho * inv_a;

        inv_VL(Ns+2, Ns+1) =   inv_rho * inv_a;

        inv_VL(Ns+3, Ns+2) =   half;
        inv_VL(Ns+3, Ns+3) =   half;
    }

    // The upper boundary is an inflow when the reference velocity is negative.
    // In the event of u == 0, treat the boundary like an outflow.
    const bool inflow = u < 0;

    // Build the in-vs-outflow characteristic-preserving projection
    MatrixXXr PG(MatrixXXr::Zero(state_count, state_count));
    if (inflow) {
        for (unsigned int is_local = 0; is_local < Ns; ++is_local) {
            PG(irhos0+is_local, irhos0+is_local) = 1;
        }
        PG(Ns,   Ns  ) = 1;
        PG(Ns+1, Ns+1) = 1;
        PG(Ns+2, Ns+2) = 1;
    } else {
        PG(Ns+3, Ns+3) = 1;
    }

    MatrixXXr BG(MatrixXXr::Zero(state_count, state_count));  // Medida's B^G_1
//     if (inflow) {
//         BG(1, 1) = v;
//         BG(3, 1) = half * (a - u);
//         BG(2, 2) = v;
//         BG(1, 3) = half * (a + u);
//         BG(3, 3) = v;
//         BG(1, 4) = half * (a - u);
//     } else {
//         BG(4, 1) = u;
//         BG(4, 4) = v;
//     }

    MatrixXXr CG(MatrixXXr::Zero(state_count, state_count));  // Medida's C^G_1
//     if (inflow) {
//         CG(1, 1) = w;
//         CG(2, 2) = w;
//         CG(3, 2) = half * (a - u);
//         CG(2, 3) = half * (a + u);
//         CG(3, 3) = w;
//         CG(2, 4) = half * (a - u);
//     } else {
//         CG(4, 2) = u;
//         CG(4, 4) = w;
//     }

    // Prepare all necessary real-valued products of the above matrices
    const MatrixXXr VL_S_RY           = VL * S * RY;
//     const Matrix5r BG_VL_S_RY_by_chi = BG * VL_S_RY / chi;
//     const Matrix5r CG_VL_S_RY_by_chi = CG * VL_S_RY / chi;
    const MatrixXXr ImPG_VL_S_RY      = (MatrixXXr::Identity(state_count, state_count) - PG) * VL_S_RY;
    const MatrixXXr inv_VL_S_RY       = inv_RY * inv_S * inv_VL;

    // Wavenumber traversal modeled after those found in suzerain/diffwave.c
//     const int Ny   = dgrid.global_wave_extent.y();
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
    for (int n = dkbz; n < dkez; ++n) {
        const int wn = wavenumber(dNz, n);
        const real_t kn = twopioverLz*wn;

        for (int m = dkbx; m < dkex; ++m) {
            const int wm = wavenumber(dNx, m);
            const real_t km = twopioverLx*wm;

            // Unpack upper boundary RHS into a contiguous buffer
            VectorXc N(state_count);
            for (int f = 0; f < N.size(); ++f) {
                N(f) = swave[f][Ny - 1][m - dkbx][n - dkbz];
            }

            // Modify the packed RHS per
            // "Implementation primarily within the nonlinear explicit operator"
            // The imaginary unit has already been included within i_stash.
            VectorXc tmp(state_count);
            tmp.noalias()  = ImPG_VL_S_RY.cast<complex_t>() * N;
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

        }
    }

    return retval;
}

} // namespace reacting

} // namespace suzerain
