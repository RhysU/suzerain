//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 Rhys Ulerich
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
// nonlinear.hpp: Building blocks for nonlinear Navier--Stokes operators
// $Id$

#ifndef NONLINEAR_HPP
#define NONLINEAR_HPP

#include "nonlinear_fwd.hpp"

#include <suzerain/error.h>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/rholut.hpp>
#include <suzerain/state.hpp>
#include <suzerain/timers.h>

#include "../support.hpp"

#pragma warning(disable:280 383 1572)

namespace suzerain { namespace reacting {

  // Indices for state fields
  struct ndx { enum {
      e,
      mx,
      my,
      mz,
      rho,
      species
    }; };

  // Indices for directions
  struct dir { enum {
      y,
      x,
      z,
      count
    }; };

  // Indices for auxiliary storage.
  //
  // NOTE: Could eliminate 1 field here by accumulating T derivatives
  // back into T storage.  Makes code more complex for minimal gain,
  // so I didn't do it.
  //
  struct aux { enum {
      kap     = 0,
      T       = 1,
      gT      = 2,
      e       = 2 + 1*dir::count,
      mx      = 2 + 2*dir::count,
      my      = 2 + 3*dir::count,
      mz      = 2 + 4*dir::count,
      rho     = 2 + 5*dir::count,
      species = 2 + 6*dir::count
    }; };


template<bool ZerothSubstep,
         linearize::type Linearize,
         class ConstitutiveLaws, // name here doesn't really fit.  Anything better??
         class ManufacturedSolution>
std::vector<real_t> applyNonlinearOperator(
            const suzerain::OperatorBase<real_t> &o,
            OperatorCommonBlock &common,
            const boost::shared_ptr<const ManufacturedSolution>& msoln,
            const real_t time,
            suzerain::ContiguousState<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            boost::shared_ptr<ConstitutiveLaws>& claws, // TODO: Decide if claws carries Ns (almost has to, right?)
            const unsigned int Ns )
{
    SUZERAIN_TIMER_SCOPED("applyNonlinearOperator");

    assert(Ns>0);

    // Shorthand
    typedef suzerain::ContiguousState<4,complex_t> state_type;
    using Eigen::Vector3r;
    using Eigen::Matrix3r;
    using Eigen::VectorXr;
    using Eigen::MatrixX3r;
    using std::size_t;

    // State enters method as coefficients in X, Y, and Z directions

    // Get total number of conserved state fields
    const int state_count = Ns+4; // = (Ns-1) species + 5 (rho, ru, rv, rw, rE)

    // Get total number of fields in aux storage
    const int aux_count = (aux::species      +   // everything *except* species derivatives
                           dir::count*(Ns-1) );  // spatial derivatives of species


    // Obtain the auxiliary storage (likely from a pool to avoid fragmenting).
    // We assume no garbage values in the memory will impact us (for speed).
    typename boost::scoped_ptr<state_type> _auxw_ptr(
            allocate_padded_state<state_type>(aux_count, o.dgrid)); // RAII
    state_type &auxw = *_auxw_ptr;                                  // Brevity

    // Sanity check incoming swave's and auxw's shape and contiguity
    assert(swave.shape()[0] == state_count);
    assert(swave.shape()[1] == (unsigned) o.dgrid.local_wave_extent.y());
    assert(swave.shape()[2] == (unsigned) o.dgrid.local_wave_extent.x());
    assert(swave.shape()[3] == (unsigned) o.dgrid.local_wave_extent.z());
    assert((unsigned) swave.strides()[1] == 1u);
    assert((unsigned) swave.strides()[2] == swave.shape()[1]);
    assert((unsigned) swave.strides()[3] == swave.shape()[1]*swave.shape()[2]);
    assert(std::equal(swave.shape() + 1, swave.shape() + 4,
                      auxw.shape() + 1));
    assert(std::equal(swave.strides() + 1, swave.strides() + 4,
                      auxw.strides() + 1));



    // Prepare common-block-like storage used to pass details from N to L.
    // Zeroing is done carefully as accumulated means and reference quantities
    // must survive from nonzero substep to substep while instant profiles do not.
    if (ZerothSubstep) common.setZero(/* Ny */ swave.shape()[1]);
    common.u().setZero();

    // Maintain stable time step values to return to the caller
    boost::array<real_t, 2> delta_t_candidates = {{
            std::numeric_limits<real_t>::max(),
            std::numeric_limits<real_t>::max()
    }};
    real_t &convective_delta_t = delta_t_candidates[0];
    real_t &diffusive_delta_t  = delta_t_candidates[1];


    //********************************************************************
    //
    // From coefficients, compute values at collocation points for
    // all state and derivatives of state *except* derivatives of
    // total energy.
    //
    // This is the first wave to physical pass.  There will be another
    // to get the temperature gradient.
    //

    // NOTE: The indexing here *assumes* that the total energy is
    // stored first.

    // Energy (no derivatives)
    o.zero_dealiasing_modes(swave, ndx::e);
    o.bop_apply     (0,    1, swave, ndx::e);

    // Everything else (all spatial derivatives)
    for (int var = ndx::mx; var<state_count; ++var) {

      // Indexing note: aux::mx is beginning of derivatives of
      // conserved state in aux ordering.  Things appearing ahead of
      // aux::mx are not used until later.

      // Compute Y derivatives of variable var at collocation points
      // Zero wavenumbers present only for dealiasing along the way
      o.zero_dealiasing_modes(swave, var);
      o.bop_accumulate(1,    1, swave, var, 0, auxw, aux::e + dir::count*var + dir::y);
      o.bop_apply     (0,    1, swave, var);

      // Compute X- and Z- derivatives of variable var at collocation points
      // Zeros wavenumbers present only for dealiasing in the target storage
      o.diffwave_accumulate(1, 0, 1, swave, var,  0, auxw, aux::e + dir::count*var + dir::x );
      o.diffwave_accumulate(0, 1, 1, swave, var,  0, auxw, aux::e + dir::count*var + dir::z );
    }

    // Collectively convert swave and auxw to physical space using parallel
    // FFTs. In physical space, we'll employ views to reshape the 4D row-major
    // (F, Y, Z, X) with contiguous (Y, Z, X) into a 2D (F, Y*Z*X) layout where
    // we know F a priori.  Reducing the dimensionality encourages linear
    // access and eases indexing overhead.

    typename support::physical_view<>::type auxp
        = support::physical_view<>::create(o.dgrid, auxw, aux_count);
    typename support::physical_view<>::type sphys
        = support::physical_view<>::create(o.dgrid, swave, state_count);
    for (size_t i = 0; i < state_count; ++i) {
      o.dgrid.transform_wave_to_physical(&sphys.coeffRef(i,0));
    }
    for (size_t i = aux::mx; i < aux_count; ++i) {
      o.dgrid.transform_wave_to_physical(&auxp.coeffRef(i,0));
    }

    //
    // Done with first wave to physical pass.
    //
    //********************************************************************



    // TODO: Removed scenario info and constants calc here.  Have to
    // refactor how scenario-like info comes in.  Is it through the
    // chemistry/constitutive law class?


    // TODO: Removed boost accumulator typedef here.  Add back when
    // necessary (should be when adding computation of reference
    // profiles for linearization).


    //********************************************************************
    //
    // Physical space is traversed linearly using a single offset 'offset'.
    // The three loop structure is present to provide the global absolute
    // positions x(i), y(j), and z(k) where necessary.
    //
    // Four traversals occur:
    // (1) Computing reference quantities and mean velocity OR mean quantities
    //     (depending on linearization and which substep is being performed).
    //
    // (2) Computing as much of the sources and fluxes as possible
    //     given *only* conserved state and conserved state gradients.
    //     In theory this is everything.  In practice, we don't do the
    //     heat flux in this pass b/c the temperature gradient isn't
    //     available.
    //
    // (3) (After going to wave space with T and back down to physical
    //     space with grad(T)), Compute heat flux contribution.
    //
    // (4) Computing any manufactured solution forcing (when enabled).
    //
    // The traversal pattern is best embodied by the last pass.  The
    // first passes use slightly more compact structure as only y(j)
    // must be known within them.  Study (4) and convince yourself
    // that the rest are equivalent but lack information on x(i)
    // and z(k).
    //
    // TODO: Combine some of these traversals (e.g., 1 and 2)


    // Traversal:
    // (1) Computing reference quantities and mean velocity OR mean velocity
    //     (depending on linearization and which substep is being performed).
    //

    // TODO: Write this code!  Removed block originally copied from
    // ../perfect/nonlinear.hpp for clarity/simplicity during
    // refactor.  Use that as a template when adding this capability
    // back.


    // Instantiate Eigen variable length arrays to store species
    // specific info
    VectorXr species(Ns); // species densities
    VectorXr Ds     (Ns); // mass diffusivities
    VectorXr hs     (Ns); // species enthalpies
    VectorXr om     (Ns); // reaction source terms
    VectorXr cs     (Ns); // species mass fractions

    // TODO: check eigen syntax here.  It appears that we need to pass
    // 3 to the ctor even though it knows this size at compile time.
    // Am I crazy?
    //
    MatrixX3r grad_species (Ns,3); // spatial derivatives of species densities
    MatrixX3r grad_cs      (Ns,3); // spatial derivatives of species mass fractions
    MatrixX3r sdiff        (Ns,3); // diffusive fluxes for species equations

    // Dereference the claws smart pointer outside the compute loop
    ConstitutiveLaws &cl = *claws;

    // Traversal:
    // (2) Compute most of sources and fluxes (everything but the heat flux)
    SUZERAIN_TIMER_BEGIN("nonlinear right hand sides");
    size_t offset = 0;
    for (int j = o.dgrid.local_physical_start.y();
         j < o.dgrid.local_physical_end.y();
         ++j) {


      // TODO: Removed lambda gets here.  Add back for time step calc.


      // Iterate across the j-th ZX plane
      const size_t last_zxoffset = offset
        + o.dgrid.local_physical_extent.z()
        * o.dgrid.local_physical_extent.x();
      for (; offset < last_zxoffset; ++offset) {

        // Unpack density-related quantities
        const real_t   rho         ( sphys(ndx::rho,    offset));

        const real_t  irho = 1.0/rho;

        const Vector3r grad_rho    (  auxp(aux::rho+dim::x,  offset),
                                      auxp(aux::rho+dim::y,  offset),
                                      auxp(aux::rho+dim::z,  offset));

        // Unpack momentum-related quantities
        const Vector3r m    ( sphys(ndx::mx, offset),
                              sphys(ndx::my, offset),
                              sphys(ndx::mz, offset));
        const real_t   div_m(  auxp(aux::mx+dim::x, offset)
                             + auxp(aux::my+dim::y, offset)
                             + auxp(aux::mz+dim::z, offset));
        const Matrix3r grad_m;
        const_cast<Matrix3r&>(grad_m) <<
                                         auxp(aux::mx+dim::x,  offset),
                                         auxp(aux::mx+dim::y,  offset),
                                         auxp(aux::mx+dim::z,  offset),
                                         auxp(aux::my+dim::x,  offset),
                                         auxp(aux::my+dim::y,  offset),
                                         auxp(aux::my+dim::z,  offset),
                                         auxp(aux::mz+dim::x,  offset),
                                         auxp(aux::mz+dim::y,  offset),
                                         auxp(aux::mz+dim::z,  offset);

        // Unpack total energy-related quantities
        const real_t e        (sphys(ndx::e,       offset));

        // Unpack species variables
        // NOTE: In species vector, idx 0 is the dilluter (the species
        // that is not explicitly part of the state vector)
        species(0) = rho;

        grad_species(0,0) = grad_rho(0);
        grad_species(0,1) = grad_rho(1);
        grad_species(0,2) = grad_rho(2);

        for (unsigned int s=1; s<Ns; ++s) {
          species(s) = sphys(ndx::species + s - 1, offset);

          grad_species(s,0) = auxp(aux::species + s - 1 + dim::x, offset);
          grad_species(s,1) = auxp(aux::species + s - 1 + dim::y, offset);
          grad_species(s,2) = auxp(aux::species + s - 1 + dim::z, offset);

          // dilluter density = rho_0 = rho - sum_{s=1}^{Ns-1} rho_s
          species(0)        -= species(s);

          grad_species(0,0) -= grad_species(s,0);
          grad_species(0,1) -= grad_species(s,1);
          grad_species(0,2) -= grad_species(s,2);
        }

        // Compute mass fractions and mass fraction gradients
        for (unsigned int s=0; s<Ns; ++s) {

          cs(s) = irho * species(s);

          grad_cs(s,0) = irho * (grad_species(s,0) - cs(s) * grad_rho(0));
          grad_cs(s,1) = irho * (grad_species(s,1) - cs(s) * grad_rho(1));
          grad_cs(s,2) = irho * (grad_species(s,2) - cs(s) * grad_rho(2));

        }


        // Compute velocity-related quantities
        const Vector3r u          = suzerain::rholut::u(rho, m);
        const real_t div_u        = suzerain::rholut::div_u(rho, grad_rho, m, div_m);
        const Matrix3r grad_u     = suzerain::rholut::grad_u(rho, grad_rho, m, grad_m);


        // Compute temperature, pressure, mass diffusivities,
        // viscosity, thermal conductivity, species enthalpies, and
        // reaction source terms
        real_t T, p, mu, kap;
        cl.evaluate(e, m.data(), rho, species.data(), cs.data(),
                    T, p, Ds.data(), mu, kap, hs.data(), om.data());

        // TODO: Compute bulk viscosity.
        // TODO: Get alpha here
        const reat_t lam = (alpha - 2.0/3.0)*mu;

        // Compute quantities related to the viscous stress tensor
        const Matrix3r tau     = suzerain::rholut::tau(mu, lam, div_u, grad_u);


        // Place to sum species fluxes from Fick's model
        Vector3r sdifftot (0.0, 0.0, 0.0);

        // Compute Fick's model contribution to diffusive fluxes and sum
        for (unsigned int s=0; s<Ns; ++s) {

          sdiff(s,0) = rho * Ds(s) * grad_cs(s,0);
          sdiff(s,1) = rho * Ds(s) * grad_cs(s,1);
          sdiff(s,2) = rho * Ds(s) * grad_cs(s,2);

          sdifftot(0) += sdiff(s,0);
          sdifftot(1) += sdiff(s,1);
          sdifftot(2) += sdiff(s,2);

        }

        // Subtract off cs*sdifftot to get SCEBD fluxes
        for (unsigned int s=0; s<Ns; ++s) {

          sdiff(s,0) -= cs(s) * sdifftot(0);
          sdiff(s,1) -= cs(s) * sdifftot(1);
          sdiff(s,2) -= cs(s) * sdifftot(2);

        }


        // Source terms get accumulated into state storage
        //
        // NOTE: Sign correct for source appearing on the RHS---i.e.,
        // U_t + div(F) = S.
        //
        // TODO: Add slow growth terms.
        sphys(ndx::e  , offset) = 0.0;
        sphys(ndx::mx , offset) = 0.0;
        sphys(ndx::my , offset) = 0.0;
        sphys(ndx::mz , offset) = 0.0;
        sphys(ndx::rho, offset) = 0.0;

        for (unsigned int s=1; s<Ns; ++s) {
          sphys(ndx::rho+s, offset) = om(s);
        }


        // Fluxes get accumulated into auxp
        //
        // NOTE: Sign correct for fluxes appearing on the LHS---i.e.,
        // U_t + div(F) = S(U).
        //
        // TODO: "Eigenify" these calcs where appropriate

        Vector3r vwork = tau*u; // TODO: check that this does what I think

        //----------------------------------------------------------------------
        // ENERGY                    = u     * (rho*H) - viscous work ...
        auxp(aux::e +dir::x, offset) = u.x() * (e + p) - vwork.x() ;
        auxp(aux::e +dir::y, offset) = u.y() * (e + p) - vwork.y() ;
        auxp(aux::e +dir::z, offset) = u.z() * (e + p) - vwork.z() ;

        // ... - species enthalpy term
        if (Ns>1) {
          // NOTE: If Ns=1, we should have sdiff(0,*) = 0.0.  Thus,
          // this loop would be entered but shouldn't do anything.
          for (unsigned int s=0; s<Ns; ++s) {
            auxp(aux::e +dir::x, offset) -= sdiff(s, 0) * hs(s);
            auxp(aux::e +dir::y, offset) -= sdiff(s, 1) * hs(s);
            auxp(aux::e +dir::z, offset) -= sdiff(s, 2) * hs(s);
          }
        }
        // NOTE: rest of energy flux (i.e., the heat flux) is accumulated below.
        //----------------------------------------------------------------------

        //----------------------------------------------------------------------
        // MOMENTUM                  = convection     - viscous    + pressure
        auxp(aux::mx+dir::x, offset) = u.x() * m.x()  -  tau(0,0)  +  p ;
        auxp(aux::mx+dir::y, offset) = u.x() * m.y()  -  tau(1,0)       ;
        auxp(aux::mx+dir::z, offset) = u.x() * m.z()  -  tau(2,0)       ;

        auxp(aux::my+dir::x, offset) = u.y() * m.x()  -  tau(0,1)       ;
        auxp(aux::my+dir::y, offset) = u.y() * m.y()  -  tau(1,1)  +  p ;
        auxp(aux::my+dir::z, offset) = u.y() * m.z()  -  tau(2,1)       ;

        auxp(aux::mz+dir::x, offset) = u.z() * m.x()  -  tau(0,2)       ;
        auxp(aux::mz+dir::y, offset) = u.z() * m.y()  -  tau(1,2)       ;
        auxp(aux::mz+dir::z, offset) = u.z() * m.z()  -  tau(2,2)  +  p ;
        //----------------------------------------------------------------------

        //----------------------------------------------------------------------
        // mass                       = mass flux
        auxp(aux::rho+dir::x, offset) = m.x();
        auxp(aux::rho+dir::y, offset) = m.y();
        auxp(aux::rho+dir::z, offset) = m.z();
        //----------------------------------------------------------------------

        //----------------------------------------------------------------------
        // species
        for (unsigned int s=0; s<Ns-1; ++s) {
          // NOTE: species(0) is the species that is not explicitly carried!

          //                                = convection     - diffusion
          auxp(aux::species+dir::x, offset) = cs(s+1)*m.x()  - sdiff(s+1,0);
          auxp(aux::species+dir::y, offset) = cs(s+1)*m.y()  - sdiff(s+1,1);
          auxp(aux::species+dir::z, offset) = cs(s+1)*m.z()  - sdiff(s+1,2);

        }
        //----------------------------------------------------------------------

        // Finally, put temperature and thermal conductivity data
        // into auxp for later use
        auxp(aux::kap, offset) = kap;
        auxp(aux::T  , offset) = T;


        // TODO: Add in time step handling once Rhys has it sorted out
        // in perfect gas case.

        // TODO: Generalize time step handling for reacting case.

      } // end X // end Z

    } // end Y
    SUZERAIN_TIMER_END("nonlinear right hand sides");

    //------------------------------------------------------------------
    //
    // At this point, we have done everything except the heat flux.
    // To do this, we will
    //
    // (1) Take temperature from physical to wave space.
    // (2) Differentiate temperature in wave space and bring grad(T)
    //     back to physical space.
    //
    // This is the second transform stage.


    // (1) Temperature from physical to wave space
    o.dgrid.transform_physical_to_wave(&auxp.coeffRef(aux::T,0));


    // TODO: Apply inverse mass matrix to get to pure coefficient space

    // ...and zero wavenumbers present only for dealiasing
    o.zero_dealiasing_modes(auxw, aux::T);


    // (2) Get derivatives.
    //
    // Compute Y derivatives of variable var at collocation points
    o.bop_accumulate(1,    1, auxw, aux::T, 0, auxw, aux::gT + dir::y);

    // Compute X- and Z- derivatives of variable var at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    o.diffwave_accumulate(1, 0, 1, auxw, aux::T,  0, auxw, aux::gT + dir::x );
    o.diffwave_accumulate(0, 1, 1, auxw, aux::T,  0, auxw, aux::gT + dir::z );

    // FFTs to get to physical space
    o.dgrid.transform_wave_to_physical(&auxp.coeffRef(aux::gT+dir::y,0));
    o.dgrid.transform_wave_to_physical(&auxp.coeffRef(aux::gT+dir::x,0));
    o.dgrid.transform_wave_to_physical(&auxp.coeffRef(aux::gT+dir::z,0));

    // Now have temperature gradient at collocation points.
    //-----------------------------------------------------------------


    // Traversal:
    // (3) (After going back to wave space with T and back down the
    //      grad(T)), Compute heat flux contribution.
    SUZERAIN_TIMER_BEGIN("nonlinear right hand sides");
    size_t offset = 0;
    for (int j = o.dgrid.local_physical_start.y();
         j < o.dgrid.local_physical_end.y();
         ++j) {

      // Iterate across the j-th ZX plane
      const size_t last_zxoffset = offset
        + o.dgrid.local_physical_extent.z()
        * o.dgrid.local_physical_extent.x();
      for (; offset < last_zxoffset; ++offset) {

        // Thermal conductivity
        const real_t   kap       ( auxp(aux::kap     ,   offset));

        // Extract grad(T)
        const Vector3r grad_T    ( auxp(aux::gT+dim::x,  offset),
                                   auxp(aux::gT+dim::y,  offset),
                                   auxp(aux::gT+dim::z,  offset));

        // Accumulate into energy fluxes
        //
        // NOTE: Sign correct for fluxes appearing on the LHS---i.e.,
        // U_t + div(F) = S(U).
        auxp(aux::e +dir::x, offset) -= kap*grad_T.x();
        auxp(aux::e +dir::y, offset) -= kap*grad_T.y();
        auxp(aux::e +dir::z, offset) -= kap*grad_T.z();


      } // end x,z
    } // end y

    // Traversal:
    // (4) Computing any manufactured solution forcing (when enabled).
    // Isolating this pass allows skipping the work when unnecessary
    if (msoln) {
      SUZERAIN_TIMER_BEGIN("manufactured forcing");

      // Dereference the msoln smart pointer outside the compute loop
      const support::manufactured_solution &ms = *msoln;

      offset = 0;
      for (int j = o.dgrid.local_physical_start.y();
           j < o.dgrid.local_physical_end.y();
           ++j) {

        const real_t y = o.y(j);

        for (int k = o.dgrid.local_physical_start.z();
             k < o.dgrid.local_physical_end.z();
             ++k) {

          const real_t z = o.z(k);

          for (int i = o.dgrid.local_physical_start.x();
               i < o.dgrid.local_physical_end.x();
               ++i, /* NB */ ++offset) {

            const real_t x = o.x(i);

            real_t Q_rho, Q_rhou, Q_rhov, Q_rhow, Q_rhoe;
            ms.Q_conservative(x, y, z, time,
                              Q_rho, Q_rhou, Q_rhov, Q_rhow, Q_rhoe);

            sphys(ndx::rho, offset) += Q_rho;
            sphys(ndx::mx , offset) += Q_rhou;
            sphys(ndx::my , offset) += Q_rhov;
            sphys(ndx::mz , offset) += Q_rhow;
            sphys(ndx::e  , offset) += Q_rhoe;

          } // end X

        } // end Z

      } // end Y

      SUZERAIN_TIMER_END("manufactured forcing");
    } // end msoln


    // TODO: Anything related to BCs that must be done in physical
    // space.


    // At this point we have accumulated sources (into sphys/swave)
    // and fluxes (into auxp/auxw) in physical space.  Now we need to
    // get back to wave space to finish accumulating the RHS.


    // Collectively convert source to wave space using parallel FFTs
    for (size_t i = 0; i < state_count; ++i) {
      o.dgrid.transform_physical_to_wave(&sphys.coeffRef(i,0));
    }

    // Collectively convert fluxes to wave space using parallel FFTs
    for (size_t i = aux::e; i < aux_count; ++i) {
      o.dgrid.transform_physical_to_wave(&auxp.coeffRef(i,0));
    }


    // Now we have sources and fluxes in wave space.  More
    // specifically, we have X,Z coefficients and Y collocation pts.
    // The last piece is to apply the divergence to the fluxes and
    // accumulate into the sources.

    // Note that the divergence is never formed explicitly.  We
    // accumulate the Y, X, and Z derivatives into the source
    // separately.
    {

      suzerain::bsplineop_lu boplu(o.bop);

      // FIXME: Need form mass matrix call
      // FIXME: Move these ops outside this function

      for (size_t i = 0; i < state_count; ++i) {

        // Apply inverse mass matrix to Y flux to get to pure
        // coefficient representation
        o.bop_solve(boplu, auxw, aux::e + dir::count*i + dir::y);

        //o.diffwave_apply(0, 0, 1, auxw, aux::e + dir::count*i + dir::y);

        // Accumulate Y derivative of Y flux (at collocation pts in Y
        // and coefficients in X,Z) into source.  Note that alpha = -1
        // b/c we need to subtract divergence from source
        o.bop_accumulate(1,   -1, auxw , aux::e + dir::count*i + dir::y,
                               1, swave, i );

        // Accumulate X and Z derivatives of X and Z fluxes into source

        // alpha = -1 b/c we need to subtract divergence from source
        o.diffwave_accumulate(1, 0, -1, auxw , aux::e + dir::count*i + dir::x,
                                     1, swave, i );

        // alpha = -1 b/c we need to subtract divergence from source
        o.diffwave_accumulate(0, 1, -1, auxw , aux::e + dir::count*i + dir::z,
                                     1, swave, i );

        // and zero wavenumbers present only for dealiasing to
        // prevent "leakage" of dealiasing modes to other routines.
        o.zero_dealiasing_modes(swave, i);

      } // end for

    } // end accumulate

    // Return the stable time step criteria separately on each rank.  The time
    // stepping logic must perform the Allreduce.  Delegating the Allreduce
    // responsibility allows reducing additional info with minimal overhead.
    return std::vector<real_t>(delta_t_candidates.begin(),
                               delta_t_candidates.end());

    // State leaves method as coefficients in X and Z directions
    // State leaves method as collocation point values in Y direction

} // end applyNonlinearOperator

} /* namespace reacting */ } /* namespace suzerain */

#endif  /* NONLINEAR_HPP */
