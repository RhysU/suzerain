//--------------------------------------------------------------------------
//
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
// cantera_interface.hpp: Interface to cantera library for chemistry, etc
// $Id$

#ifndef CANTERA_INTERFACE_HPP
#define CANTERA_INTERFACE_HPP

// No HAVE_CANTERA #ifdef as any HAVE_CANTERA checks occur
// well before anyone tries to compile this file.

#include <cassert>
#include <cstddef>
#include <string>

#include <cantera/Cantera.h>
#include <cantera/IdealGasMix.h>
#include <cantera/transport.h>

namespace suzerain { namespace reacting {

class canteraInterface {

public:

    // Ctor: require cantera input file and desired mixture construct
    canteraInterface (const std::string& cantera_chem_file,
                      const std::string& cantera_mixture_name);

    ~canteraInterface ();

    // Pass in conserved state and get back everything in one shot
    void evaluate (const double  e,
                   const double* m,
                   const double  rho,
                   const double* species,
                   const double* cs,
                   double& T,
                   double& p,
                   double* Ds,
                   double& mu,
                   double& kap,
                   double* hs,
                   double* om);

    std::size_t Ns() { return _cgas.nSpecies(); }

private:

    // TODO: Is "species" actually rho_s? c_s?  "species" is wildly ambiguous.
    void setCanteraState (const double  e,
                          const double* m,
                          const double  rho,
                          const double* species);

    void TPhs (double& T,
               double& p,
               double* hs);

    void transport (double* Ds,
                    double& mu,
                    double& kap);

    void enthalpies (double* hs);

    void sources (double *om);

    // TODO Address member ownership semantics.
    //
    // Either canteraInterface shouldn't be copyable or assignable or a
    // copy/assignment operator combo should be declared or these resources
    // could be shared amongst multiple instances using smart pointers.
    //
    // Odd to just have naked pointer members otherwise.

    // Cantera ideal gas object
    Cantera::IdealGasMix* _cgas; // any namespace issue here??

    // Cantera transport object
    Cantera::Transport* _ctrans;

};

// Public interface

// Ctor instantiates cantera ideal gas and cantera transport
canteraInterface::canteraInterface (const std::string& cantera_chem_file,
                                    const std::string& cantera_mixture_name)
    :
    _cgas(new Cantera::IdealGasMix(cantera_chem_file, cantera_mixture_name)),
    _ctrans(Cantera::newTransportMgr("Pecos", _cgas))
{
    // NOP
}

// Dtor deletes cantera ideal gas and cantera transport
canteraInterface::~canteraInterface ()
{
    delete _cgas;
    delete _ctrans;
}

// Evaluate: takes state and gives back everything we need from
// cantera, including temp, pres, transport props, enthalpies, and
// reaction rates
void canteraInterface::evaluate (const double  e,
                                 const double* m,
                                 const double  rho,
                                 const double* species,
                                 const double* cs,
                                 double& T,
                                 double& p,
                                 double* Ds,
                                 double& mu,
                                 double& kap,
                                 double* hs,
                                 double* om)
{
    setCanteraState(e, m, rho, species);
    TPhs(T, p, hs);
    transport(Ds, mu, kap);
    sources(om);
}


// Private methods

// Set cantera state given input info
void canteraInterface::setCanteraState (const double  e,
                                        const double* m,
                                        const double  rho,
                                        const double* species,
                                        const double* cs)
{
    assert(m);
    assert(species);
    assert(cs);

    double irho = 1.0/rho;
    double einternal = e - 0.5*irho*(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);

    // Set cantera state by setting mass fractions, internal energy
    // (per unit mass), and specific volume
    _cgas->setMassFractions(cs);
    _cgas->setState_UV(e, irho);
}


// Get temperature and pressure
void canteraInterface::TPhs (double& T,
                             double& p,
                             double* hs)
{
    T = _cgas->temperature ();
    p = _cgas->pressure ();

    // TODO: Do we have to do this?  Cantera can't give us hs directly?
    //
    // update: we might be able to: enthalpy_mole() looks like it might be correct
    //         are we trying to pull down the molar enthalpy? Or is this the enthalpy of formation?
    //
    // http://cantera.github.com/docs/doxygen/html/classCantera_1_1IdealGasPhase.html#a6dd87c68aea566830f1e6af8b2412d63

    _cgas->getEnthalpy_RT(hs);

    for (int s=0; s<Ns; ++s)
        hs[s] *= Cantera::GasConstant * T / _cgas->molecularWeight(s);
}

// Get transport properties
void canteraInterface::transport (double* Ds,
                                  double& mu,
                                  double& kap)
{
    // TODO: This right diffusivity call?  Check w/ Nick re: mass flux
    // vs. mole flux resolution.
    // UPDATE: Looks good to me (Nick)

    // the call below returns coefficients for calculating
    // the diffusive mass fluxes
    _ctrans->getMixDiffCoeffsMass (Ds);
    mu  = _ctrans->viscosity ();
    kap = _ctrans->thermalConductivity ();
}

// Get reaction rates
void canteraInterface::sources (double* om)
{
    // Get sources in kmol/m^3/s
    _cgas->getNetProductionRates(om);

    // TODO: Do we have to do this?  Cantera can't give us om in kg/m^3/s directly

    // Convert to kg/m^3/s
    for (int s=0; s<Ns; ++s)
        om[s] *= _cgas->molecularWeight(s);

}

} /* namespace reacting */ } /* namespace suzerain */

#endif
