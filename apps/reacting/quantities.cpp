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
 * @copydoc quantities.hpp
 */

#include "quantities.hpp"

#include <esio/esio.h>
#include <esio/error.h>

#include <suzerain/blas_et_al.hpp>
#include <suzerain/coalescing_pool.hpp>
#include <suzerain/common.hpp>
#include <suzerain/countof.h>
#include <suzerain/diffwave.hpp>
#include <suzerain/error.h>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/operator_tools.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/rholut.hpp>
#include <suzerain/shared_range.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/support.hpp>

#include "reacting_ndx.hpp"


using boost::numeric_cast;

namespace suzerain {

namespace reacting {

static const char default_who[] = "quantities";

quantities_base::quantities_base()
    : t(std::numeric_limits<real_t>::quiet_NaN())
    , who(default_who)
{
}

quantities_base::quantities_base(real_t t)
    : t(t)
    , who(default_who)
{
}

quantities_base::quantities_base(
        real_t t,
        quantities_base::storage_type::Index Ny)
    : t(t),
      storage(storage_type::Zero(Ny, storage_type::ColsAtCompileTime))
    , who(default_who)
{
}

// Helper for the quantities_base::save(...) implementation just below
class quantities_saver
{

public:

    quantities_saver(const std::string& who,
                     const esio_handle  esioh,
                     const std::string& prefix)
        : who(who), esioh(esioh), prefix(prefix) {}

    template< typename EigenArray >
    bool operator()(const std::string& name, const EigenArray& dat) const
    {
        int procid;
        esio_handle_comm_rank(esioh, &procid);
        esio_field_establish(esioh,
                             1,          0, (procid == 0 ? 1          : 0),
                             dat.cols(), 0, (procid == 0 ? dat.cols() : 0),
                             dat.rows(), 0, (procid == 0 ? dat.rows() : 0));

        std::string key(prefix);
        key.append(name);
        return ESIO_SUCCESS == esio_field_write(
            esioh, key.c_str(), dat.data(), 0, 0, 0,
            "Mean quantity sample stored using row-major indices (B-spline"
            " coefficient, tensor component, sample number) where the"
            " B-spline basis is defined by /Ny, /breakpoints_y, and /knots");
    }

private:

    std::string who;
    esio_handle esioh;
    std::string prefix;

};

bool quantities_base::save(const esio_handle h) const
{
    if (this->storage.size()) {
        quantities_saver f(this->who, h, "bar_");
        return this->foreach(f);
    } else {
        WARN0(who, "No mean quantity samples saved--"
                   " trivial storage needs detected");
        return true; // Warning occurs but behavior was "successful".
    }
}

// Helper for the quantities_base::load(...) implementation just below
class quantities_loader
{

public:

    quantities_loader(const std::string& who,
                      const esio_handle esioh,
                      const std::string& prefix)
        : who(who), esioh(esioh), prefix(prefix) {}

    template< typename EigenArray >
    bool operator()(const std::string& name, const EigenArray& dat_) const {

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
            return ESIO_SUCCESS == esio_field_read(
                    esioh, key.c_str(), dat.data(), 0, 0, 0);
        } else {
            WARN0(who, "Unable to load " << key
                  << " for nscalar = " << dat.cols());
            dat.fill(std::numeric_limits<
                     typename EigenArray::Scalar>::quiet_NaN());
            return false;
        }
    }

private:

    std::string who;
    esio_handle esioh;
    std::string prefix;

};

bool quantities_base::load(const esio_handle h)
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
        success = this->foreach(f);
    } else {
        WARN0(who, "No mean quantity samples loaded--"
                   " unable to anticipate storage needs");
    }

    return success;
}

quantities::quantities()
    : quantities_base()
    , Ns(0)
{
}

quantities::quantities(real_t t)
    : quantities_base(t)
    , Ns(0)
{
}

quantities::quantities(
        real_t t,
        quantities::storage_type::Index Ny,
        quantities::storage_type::Index Ns)
    : quantities_base(t, Ny)
    , Ns(Ns)
{
    species_storage.setZero(Ny, 12*Ns);
}


bool quantities::save(const esio_handle h) const
{
    bool retval = super::save(h);

    DEBUG0("Saving species statistics.");
    if (this->species_storage.size()) {
        quantities_saver f("quantities", h, "bar_");
        retval &= f("rho_s",             this->rho_s(0, Ns));
        retval &= f("cs_s" ,             this->cs_s (0, Ns));
        retval &= f("om_s" ,             this->om_s (0, Ns));
        retval &= f("rho_s_u"          , this->rho_s_u          ());
        retval &= f("rho_Ds_grad_cs"   , this->rho_Ds_grad_cs   ());
        retval &= f("rho_Ds_grad_cs_hs", this->rho_Ds_grad_cs_hs());
    } else {
        WARN0("quantities", "No mean quantity samples saved--"
                            " trivial storage needs detected");
    }

    return retval;
}

bool quantities::load(const esio_handle h)
{
    bool retval = super::load(h);

    // Defensively NaN out all storage in the instance prior to load.
    this->species_storage.fill(std::numeric_limits<real_t>::quiet_NaN());

    DEBUG0("Loading species statistics.");

    int cglobal, bglobal, aglobal;
    if (ESIO_SUCCESS == esio_field_size(h, "bar_rho_s",
                                        &cglobal, &bglobal, &aglobal)) {
        this->species_storage.resize(aglobal, 12*bglobal);

        if (this->Ns==0)
            { this->Ns = bglobal; }
        else
            { assert(static_cast<int>(this->Ns)==bglobal); }

        quantities_loader f("quantities", h, "bar_");
        retval &= f("rho_s", this->rho_s(0, Ns));
        retval &= f("cs_s" , this->cs_s (0, Ns));
        retval &= f("om_s" , this->om_s (0, Ns));
        retval &= f("rho_s_u"          , this->rho_s_u          ());
        retval &= f("rho_Ds_grad_cs"   , this->rho_Ds_grad_cs   ());
        retval &= f("rho_Ds_grad_cs_hs", this->rho_Ds_grad_cs_hs());
    } else {
        WARN0("quantities", "No mean quantity samples loaded--"
                            " unable to anticipate storage needs");
    }

    return retval;
}


// This looks like logic from operator_nonlinear.hpp but does not belong there.
// Reading through that file, especially the apply_operator implementation, is
// recommended before reviewing this logic.  This routine is definitely
// suboptimal but is expected to be invoked relatively infrequently.
quantities sample_quantities(
        const antioch_constitutive& cmods,
        const specification_grid &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        contiguous_state<4,complex_t> &swave,
        const real_t t)
{
    // FIXME: This function supports reacting in that all calculations
    // here use appropriate constitutive laws, etc.  However, not all
    // stats we want, particularly species-related quantities, are
    // currently computed.
    if (cmods.Ns()>1) {
        WARN0("sample_quantitites",
              "sample_quantities does not fully support reacting flow yet!");
    }

    const size_t Ns = cmods.Ns();
    const size_t state_count = Ns+4;

    // Shorthand
    const std::size_t Ny = swave.shape()[1];
    namespace acc = boost::accumulators;
    typedef contiguous_state<4,complex_t> state_type;

    // State enters method as coefficients in X, Y, and Z directions

    // We need auxiliary scalar-field storage.  Prepare logical indices using a
    // struct for scoping (e.g. aux::rho_y).  Ordering will match usage below.

    // Indices for directions
    struct dir { enum {
            y,
            x,
            z,
            count
        }; };

    // Indices for auxiliary storage for state derivatives
    // and heat flux related quantities
    struct aux { enum {
            T       = 0,
            gT      = 1,
            e       = 1 + 1*dir::count,
            mx      = 1 + 2*dir::count,
            my      = 1 + 3*dir::count,
            mz      = 1 + 4*dir::count,
            rho     = 1 + 5*dir::count,
            species = 1 + 6*dir::count
        }; };

    // Get total number of fields in aux storage
    const size_t aux_count = (aux::species      +
                              dir::count*(Ns-1) );

    // Obtain the auxiliary storage (likely from a pool to avoid fragmenting).
    // We assume no garbage values in the memory will impact us (for speed).
    scoped_ptr<state_type> _auxw_ptr(
            support::allocate_padded_state<state_type>(
                aux_count, dgrid));                                // RAII
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

    // Rank-specific details accumulated in ret to be MPI_Allreduce later
    //quantities ret(t, Ny);
    quantities ret(t, Ny, Ns);

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

    for (size_t var = 0; var<state_count; ++var) {

      // Compute Y derivatives of variable var at collocation points
      // Zero wavenumbers present only for dealiasing along the way
      otool.zero_dealiasing_modes(swave, var);
      otool.bop_accumulate(1,  1, swave, var,
                               0, auxw , aux::e + dir::count*var + dir::y);
      otool.bop_apply     (0,  1, swave, var);

      // Compute X- and Z- derivatives of variable var at collocation points
      // Zeros wavenumbers present only for dealiasing in the target storage
      otool.diffwave_accumulate(1, 0, 1, swave, var,
                                      0, auxw , aux::e + dir::count*var + dir::x );
      otool.diffwave_accumulate(0, 1, 1, swave, var,
                                      0, auxw , aux::e + dir::count*var + dir::z );

    }

    // Collectively convert swave and auxw to physical space using parallel
    // FFTs. In physical space, we'll employ views to reshape the 4D row-major
    // (F, Y, Z, X) with contiguous (Y, Z, X) into a 2D (F, Y*Z*X) layout where
    // we know F a priori.  Reducing the dimensionality encourages linear
    // access and eases indexing overhead.
    physical_view<> auxp(dgrid, auxw);
    physical_view<>sphys(dgrid, swave);
    for (std::size_t i = 0; i < state_count; ++i) {
        dgrid.transform_wave_to_physical(&sphys.coeffRef(i,0));
    }
    for (std::size_t i = 0; i < aux_count; ++i) {
        dgrid.transform_wave_to_physical(&auxp.coeffRef(i,0));
    }



    // Before traversal, instantiate Eigen variable length arrays
    // to store species specific info
    VectorXr species(Ns); // species densities
    VectorXr Ds     (Ns); // mass diffusivities
    VectorXr hs     (Ns); // species enthalpies
    VectorXr oms    (Ns); // reaction source terms
    VectorXr cs     (Ns); // species mass fractions
    VectorXr etots  (Ns); // species internal energies

    Matrix3Xr grad_species (3,Ns); // spatial derivatives of species densities
    Matrix3Xr grad_cs      (3,Ns); // spatial derivatives of species mass fracs
    Matrix3Xr sdiff        (3,Ns); // diffusive fluxes for species equations


    //------------------------------------------------------------------
    // Adapted from nonlinear, to get temperature gradient and heat flux
    //------------------------------------------------------------------
    //
    // (1) Compute temperature and store the field
    // (2) Take temperature from physical to wave space.
    // (3) Differentiate temperature in wave space and bring grad(T)
    //     back to physical space.
    // (4) (After going back to wave space with T and back down the
    //      grad(T)), Compute heat flux along wit all the other
    //      variables for sampling.

    // (1) Traversal to compute temperature
    real_t Tguess = -1;
    for (int offset = 0, j = dgrid.local_physical_start.y();
         j < dgrid.local_physical_end.y();
         ++j) {

        real_t T = -1;
        for (int k = dgrid.local_physical_start.z();
            k < dgrid.local_physical_end.z();
            ++k) {

            for (int i = dgrid.local_physical_start.x();
                i < dgrid.local_physical_end.x();
                ++i, /* NB */ ++offset) {

                // Unpack total energy, momentum, and density
                const real_t   e  (sphys(ndx::e,   offset));
                const Vector3r m  (sphys(ndx::mx,  offset),
                        sphys(ndx::my,  offset),
                        sphys(ndx::mz,  offset));
                const real_t   rho(sphys(ndx::rho, offset));

                const real_t irho = 1.0/rho;

                // Unpack species variables
                // NOTE: In species vector, idx 0 is the dilluter (the species
                // that is not explicitly part of the state vector)
                // Compute mass fractions
                species(0) = rho;
                for (unsigned int s=1; s<Ns; ++s) {
                    species(s) = sphys(ndx::species + s - 1, offset);
                    cs(s) = irho * species(s);

                    // dilluter density = rho_0 = rho - sum_{s=1}^{Ns-1} rho_s
                    species(0)        -= species(s);
                }
                cs(0) = irho * species(0);

                // Compute temperature from internal energy
                // (assuming thermal equilibrium)
                const real_t re_internal = e -
                    0.5*irho*(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
                T = cmods.sm_thermo->T_from_e_tot(irho*re_internal, cs, Tguess);
                Tguess = T;

                auxp(aux::T  , offset) = T;
            } // end x
        } // end z
    } // end y

    // (2a) Temperature from physical to wave space (still collocation in y)
    otool.dgrid.transform_physical_to_wave(&auxp.coeffRef(aux::T,0));

    // (2b) Apply inverse mass matrix to get to pure coefficient space
    otool.bop_solve(*otool.massluz(), auxw, aux::T);

    // ... and zero wavenumbers present only for dealiasing
    // while simultaneously normalizing FFT
    otool.diffwave_apply(0, 0, otool.dgrid.chi(), auxw, aux::T);

    // (3) Get derivatives.
    //
    // Compute Y derivatives of variable var at collocation points
    otool.bop_accumulate(1,    1, auxw, aux::T, 0, auxw, aux::gT + dir::y);

    // Get to collocation points in preparation for computing x and z
    // derivatives
    otool.bop_apply     (0,    1, auxw, aux::T);

    // Compute X- and Z- derivatives of variable var at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    otool.diffwave_accumulate(1, 0, 1, auxw, aux::T, 0, auxw, aux::gT+dir::x);
    otool.diffwave_accumulate(0, 1, 1, auxw, aux::T, 0, auxw, aux::gT+dir::z);

    // FFTs to get to physical space
    otool.dgrid.transform_wave_to_physical(&auxp.coeffRef(aux::gT+dir::y,0));
    otool.dgrid.transform_wave_to_physical(&auxp.coeffRef(aux::gT+dir::x,0));
    otool.dgrid.transform_wave_to_physical(&auxp.coeffRef(aux::gT+dir::z,0));

    // Bring T to physical space as well to use the exact value as guess
    // I'm not sure if doing this is less expensive though ...
    otool.dgrid.transform_wave_to_physical(&auxp.coeffRef(aux::T,0));

    // Now have temperature gradient at collocation points.
    //-----------------------------------------------------------------


    // (4) Physical space is traversed linearly using a single offset 'offset'.
    // The three loop structure is present to provide the global absolute
    // positions x(i), y(j), and z(k) where necessary.
    for (int offset = 0, j = dgrid.local_physical_start.y();
         j < dgrid.local_physical_end.y();
         ++j) {

        // Type used to accumulate mean quantities versus wall-normal position
        typedef acc::accumulator_set<
                real_t, acc::stats<acc::tag::sum_kahan>
            > accumulator_type;

        // Accumulators for each mean quantity computed in physical space.
        // For example, quantity "foo" has accumulator "sum_foo".
        // Declared within Y loop so they are reset on each Y iteration.
#define DECLARE(r, data, tuple)                                         \
        accumulator_type BOOST_PP_CAT(sum_,BOOST_PP_TUPLE_ELEM(2, 0, tuple)) \
                [BOOST_PP_TUPLE_ELEM(2, 1, tuple)];
        BOOST_PP_SEQ_FOR_EACH(DECLARE,,
                SUZERAIN_REACTING_FLOW_QUANTITIES_PHYSICAL)
#undef DECLARE

        // Vectors of species quantity accumulators
        std::vector<accumulator_type> sum_rho_s  (Ns);
        std::vector<accumulator_type> sum_cs_s   (Ns);
        std::vector<accumulator_type> sum_om_s   (Ns);
        std::vector<accumulator_type> sum_rho_s_u          (Ns*dir::count);
        std::vector<accumulator_type> sum_rho_Ds_grad_cs   (Ns*dir::count);
        std::vector<accumulator_type> sum_rho_Ds_grad_cs_hs(Ns*dir::count);

        for (int k = dgrid.local_physical_start.z();
            k < dgrid.local_physical_end.z();
            ++k) {

            for (int i = dgrid.local_physical_start.x();
                i < dgrid.local_physical_end.x();
                ++i, /* NB */ ++offset) {

                // Unpack total energy-related quantities
                const real_t e(sphys(ndx::e, offset));
                const Vector3r grad_e(auxp(aux::e+dir::x, offset),
                                      auxp(aux::e+dir::y, offset),
                                      auxp(aux::e+dir::z, offset));

                // Unpack momentum-related quantities
                const Vector3r m(sphys(ndx::mx, offset),
                                 sphys(ndx::my, offset),
                                 sphys(ndx::mz, offset));
                const real_t div_m = (  auxp(aux::mx+dir::x, offset)
                                      + auxp(aux::my+dir::y, offset)
                                      + auxp(aux::mz+dir::z, offset));

                const Matrix3r grad_m;
                const_cast<Matrix3r&>(grad_m) <<
                    auxp(aux::mx+dir::x,  offset),
                    auxp(aux::mx+dir::y,  offset),
                    auxp(aux::mx+dir::z,  offset),
                    auxp(aux::my+dir::x,  offset),
                    auxp(aux::my+dir::y,  offset),
                    auxp(aux::my+dir::z,  offset),
                    auxp(aux::mz+dir::x,  offset),
                    auxp(aux::mz+dir::y,  offset),
                    auxp(aux::mz+dir::z,  offset);

                // Unpack density-related quantities
                const real_t  rho(sphys(ndx::rho, offset));
                const real_t irho = 1.0/rho;

                const Vector3r grad_rho(auxp(aux::rho+dir::x, offset),
                                        auxp(aux::rho+dir::y, offset),
                                        auxp(aux::rho+dir::z, offset));


                // Unpack species variables
                // NOTE: In species vector, idx 0 is the dilluter (the species
                // that is not explicitly part of the state vector)
                species(0) = rho;

                grad_species(0,0) = grad_rho(0);
                grad_species(1,0) = grad_rho(1);
                grad_species(2,0) = grad_rho(2);

                for (unsigned int s=1; s<Ns; ++s) {
                    species(s) = sphys(ndx::species + s - 1, offset);

                    unsigned int si = aux::species + (s-1)*dir::count;
                    grad_species(0,s) = auxp(si + dir::x, offset);
                    grad_species(1,s) = auxp(si + dir::y, offset);
                    grad_species(2,s) = auxp(si + dir::z, offset);

                    // dilluter density = rho_0 = rho - sum_{s=1}^{Ns-1} rho_s
                    species(0)        -= species(s);

                    grad_species(0,0) -= grad_species(0,s);
                    grad_species(1,0) -= grad_species(1,s);
                    grad_species(2,0) -= grad_species(2,s);
                }

                // Compute mass fractions and mass fraction gradients
                for (unsigned int s=0; s<Ns; ++s) {

                    cs(s) = irho * species(s);

                    grad_cs(0,s) = irho*(grad_species(0,s) - cs(s)*grad_rho(0));
                    grad_cs(1,s) = irho*(grad_species(1,s) - cs(s)*grad_rho(1));
                    grad_cs(2,s) = irho*(grad_species(2,s) - cs(s)*grad_rho(2));

                }

                // Compute velocity-related quantities
                const Vector3r u     = suzerain::rholut::u(rho, m);

                const real_t div_u   =
                    suzerain::rholut::div_u(rho, grad_rho, m, div_m);

                const Matrix3r grad_u=
                    suzerain::rholut::grad_u(rho, grad_rho, m, grad_m);

                // ... this is vorticity
                const Vector3r om       (grad_u(2,1) - grad_u(1,2),
                                         grad_u(0,2) - grad_u(2,0),
                                         grad_u(1,0) - grad_u(0,1));

                // Compute temperature, pressure, mass diffusivities,
                // viscosity, thermal conductivity, species enthalpies, and
                // reaction source terms
                real_t T, p, mu, kap, a, Cv, Cp;
                Tguess = auxp(aux::T, offset);
                cmods.evaluate(e, m, rho, species, cs, Tguess,
                               T, p, Ds, mu, kap, hs, oms, a, Cv, Cp);

                // Extract grad(T)
                const Vector3r grad_T ( auxp(aux::gT+dir::x, offset),
                                        auxp(aux::gT+dir::y, offset),
                                        auxp(aux::gT+dir::z, offset));

                const real_t lam = (cmods.alpha - 2.0/3.0)*mu;

                // Compute quantities related to the viscous stress tensor
                const Matrix3r tau = suzerain::rholut::tau(mu, lam, div_u, grad_u);
                const Vector3r tau_u = tau * u;

                // Compute Mach number
                const real_t M = sqrt(u.x() * u.x() +
                                      u.y() * u.y() +
                                      u.z() * u.z()) / a;

                // Accumulate quantities into sum_XXX using function syntax.
                for (unsigned int s=0; s<Ns; ++s) {
                    sum_rho_s[s](species[s]);
                    sum_cs_s [s](cs[s]);
                    sum_om_s [s](oms[s]);

                    sum_rho_s_u[s*dir::count+0](species[s] * u.x());
                    sum_rho_s_u[s*dir::count+1](species[s] * u.y());
                    sum_rho_s_u[s*dir::count+2](species[s] * u.z());

                    sum_rho_Ds_grad_cs[s*dir::count+0]
                      (rho * Ds[s] * grad_cs(0,s));
                    sum_rho_Ds_grad_cs[s*dir::count+1]
                      (rho * Ds[s] * grad_cs(1,s));
                    sum_rho_Ds_grad_cs[s*dir::count+2]
                      (rho * Ds[s] * grad_cs(2,s));

                    sum_rho_Ds_grad_cs_hs[s*dir::count+0]
                      (rho * Ds[s] * grad_cs(0,s) * hs[s]);
                    sum_rho_Ds_grad_cs_hs[s*dir::count+1]
                      (rho * Ds[s] * grad_cs(1,s) * hs[s]);
                    sum_rho_Ds_grad_cs_hs[s*dir::count+2]
                      (rho * Ds[s] * grad_cs(2,s) * hs[s]);
                }

                sum_E[0](e / rho);

                sum_T[0](T);

                sum_p[0](p);

                sum_a[0](a);

                sum_M[0](M);

                sum_Cv[0](Cv);

                sum_Cp[0](Cp);

                sum_mu[0](mu);

                sum_nu[0](mu / rho);

                sum_kappa[0](kap);

                // NOTE: D0 is meaningful alone only in the case of
                // constant Lewis number
                sum_D0[0](Ds[0]);

                sum_rho_D0[0](rho * Ds[0]);

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

                sum_rho_E_u[0](e * u.x());
                sum_rho_E_u[1](e * u.y());
                sum_rho_E_u[2](e * u.z());

                sum_p_u[0](p * u.x());
                sum_p_u[1](p * u.y());
                sum_p_u[2](p * u.z());

                sum_rho_T[0](rho * T);

                sum_rho_mu[0](rho * mu);

                sum_mu_S[0](mu * ( grad_u(0,0)                    - div_u / 3));
                sum_mu_S[1](mu * ((grad_u(0,1) + grad_u(1,0)) / 2            ));
                sum_mu_S[2](mu * ((grad_u(0,2) + grad_u(2,0)) / 2            ));
                sum_mu_S[3](mu * ( grad_u(1,1)                    - div_u / 3));
                sum_mu_S[4](mu * ((grad_u(1,2) + grad_u(2,1)) / 2            ));
                sum_mu_S[5](mu * ( grad_u(2,2)                    - div_u / 3));

                sum_mu_div_u[0](mu * div_u);

                sum_kappa_grad_T[0](kap * grad_T.x());
                sum_kappa_grad_T[1](kap * grad_T.y());
                sum_kappa_grad_T[2](kap * grad_T.z());

                sum_rho_rho[0](rho * rho);

                sum_T_T[0](T * T);

                sum_p_p[0](p * p);

                sum_a_a[0](a * a);

                sum_M_M[0](M * M);

                sum_mu_mu[0](mu * mu);

                sum_u_u[0](u.x() * u.x());
                sum_u_u[1](u.x() * u.y());
                sum_u_u[2](u.x() * u.z());
                sum_u_u[3](u.y() * u.y());
                sum_u_u[4](u.y() * u.z());
                sum_u_u[5](u.z() * u.z());

                sum_om[0](om.x());
                sum_om[1](om.y());
                sum_om[2](om.z());

                sum_om_om[0](om.x() * om.x());
                sum_om_om[1](om.x() * om.y());
                sum_om_om[2](om.x() * om.z());
                sum_om_om[3](om.y() * om.y());
                sum_om_om[4](om.y() * om.z());
                sum_om_om[5](om.z() * om.z());

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
        // Eigen comma initialization syntax.  Yes, this is not stride 1. ...

        // ... flow quantities...
#define EXTRACT_SUM(z, n, q) acc::sum(sum_##q[n])
#define MOVE_SUM_INTO_TMP(r, data, tuple)                                 \
        ret.BOOST_PP_TUPLE_ELEM(2, 0, tuple)().row(j) <<                  \
            BOOST_PP_ENUM(BOOST_PP_TUPLE_ELEM(2, 1, tuple),               \
                          EXTRACT_SUM, BOOST_PP_TUPLE_ELEM(2, 0, tuple));
        BOOST_PP_SEQ_FOR_EACH(MOVE_SUM_INTO_TMP,,
                SUZERAIN_REACTING_FLOW_QUANTITIES_PHYSICAL)
#undef EXTRACT_SUM
#undef MOVE_SUM_INTO_TMP

        // ... species quantities
        for (unsigned int s=0; s<Ns; ++s) {
            ret.rho_s(s)[j] = acc::sum(sum_rho_s[s]);
            ret.cs_s (s)[j] = acc::sum(sum_cs_s [s]);
            ret.om_s (s)[j] = acc::sum(sum_om_s [s]);

            ret.rho_s_u(s,0)[j] = acc::sum(sum_rho_s_u[s*dir::count+0]);
            ret.rho_s_u(s,1)[j] = acc::sum(sum_rho_s_u[s*dir::count+1]);
            ret.rho_s_u(s,2)[j] = acc::sum(sum_rho_s_u[s*dir::count+2]);

            ret.rho_Ds_grad_cs(s,0)[j] =
              acc::sum(sum_rho_Ds_grad_cs[s*dir::count+0]);
            ret.rho_Ds_grad_cs(s,1)[j] =
              acc::sum(sum_rho_Ds_grad_cs[s*dir::count+1]);
            ret.rho_Ds_grad_cs(s,2)[j] =
              acc::sum(sum_rho_Ds_grad_cs[s*dir::count+2]);

            ret.rho_Ds_grad_cs_hs(s,0)[j] =
              acc::sum(sum_rho_Ds_grad_cs_hs[s*dir::count+0]);
            ret.rho_Ds_grad_cs_hs(s,1)[j] =
              acc::sum(sum_rho_Ds_grad_cs_hs[s*dir::count+1]);
            ret.rho_Ds_grad_cs_hs(s,2)[j] =
              acc::sum(sum_rho_Ds_grad_cs_hs[s*dir::count+2]);
        }

    } // end Y

    // Notice dgrid.rank_zero_zero_modes already contains "wave-sampled"
    // quantities while other ranks have zeros in those locations.

    // Allreduce to obtain global sums on every rank
    SUZERAIN_MPICHKR(MPI_Allreduce(
            MPI_IN_PLACE, ret.storage.data(), ret.storage.size(),
            mpi::datatype<quantities::storage_type::Scalar>::value,
            MPI_SUM, MPI_COMM_WORLD));
    SUZERAIN_MPICHKR(MPI_Allreduce(
            MPI_IN_PLACE, ret.species_storage.data(),
                          ret.species_storage.size(),
            mpi::datatype<quantities::storage_type::Scalar>::value,
            MPI_SUM, MPI_COMM_WORLD));

    // Physical space sums, which are at collocation points, need to be
    // divided by the dealiased extents and converted to coefficients.
    const real_t scale_factor = 1 / dgrid.chi();
    bsplineop_lu scaled_mass(cop);
    scaled_mass.opform(1, &scale_factor, cop);
    scaled_mass.factor();

    // Convert flow and species data on collocation points to coefficients
    scaled_mass.solve(quantities::nscalars::physical,
            ret.storage.middleCols<quantities::nscalars::physical>(
                quantities::start::physical).data(),
            ret.storage.innerStride(), ret.storage.outerStride());
    scaled_mass.solve(12*Ns, ret.species_storage.data(),
                      ret.species_storage.innerStride(),
                      ret.species_storage.outerStride());

    // Fill with NaNs those samples that were not computed by this method
#define FILL(r, data, tuple)                                         \
    ret.BOOST_PP_TUPLE_ELEM(2, 0, tuple)().fill(std::numeric_limits< \
            quantities::storage_type::Scalar>::quiet_NaN());
    BOOST_PP_SEQ_FOR_EACH(FILL,,SUZERAIN_REACTING_FLOW_QUANTITIES_IMPLICIT)
#undef FILL

    return ret;
}

} // namespace reacting

} // namespace suzerain
