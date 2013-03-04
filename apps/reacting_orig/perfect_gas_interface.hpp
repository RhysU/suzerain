//--------------------------------------------------------------------------
//
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

#ifndef PERFECT_GAS_INTERFACE_HPP
#define PERFECT_GAS_INTERFACE_HPP

/** @file
 * Interface to perfect gas transport and thermodynamics.
 */

#include <cstddef>

namespace suzerain {

namespace reacting {

class perfectGasInterface {

public:

    // Ctor
    perfectGasInterface (const double Cp,
                         const double Cv,
                         const double Pr,
                         const double T0,
                         const double mu0,
                         const double beta);

    // Dtor
    ~perfectGasInterface ();

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

    // Number of species
    std::size_t Ns() { return _Ns; }

private:

    // Number of species
    std::size_t _Ns;

    // Thermo constants
    double _Cp;
    double _Cv;
    double _R;
    double _gam;
    double _gmi;

    // Prandtl number
    double _Pr;

    // Power law viscosity
    double _T0;
    double _mu0;
    double _beta;

};

// Public interface

// Ctor
perfectGasInterface::perfectGasInterface (const double Cp,
                                          const double Cv,
                                          const double Pr,
                                          const double T0,
                                          const double mu0,
                                          const double beta)
    :
    _Ns(1u),
    _Cp(Cp),
    _Cv(Cv),
    _R(Cp - Cv),
    _gam(Cp/Cv),
    _gmi(_gam-1.0),
    _Pr(Pr),
    _T0(T0),
    _mu0(mu0),
    _beta(beta)
{
    // NOP
}

// Dtor
perfectGasInterface::~perfectGasInterface ()
{
    // NOP
}

// Evaluate: takes state and gives back everything we need from
// Cantera, including temp, pres, transport props, enthalpies, and
// reaction rates
void perfectGasInterface::evaluate (const double  e,
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

    // NOTE: Really only have e, m, rho input and T, p, mu, kap
    // output.  Everything else is used for the multispecies case.

    double irho = 1.0/rho;

    p = _gmi * (e - 0.5*irho*(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]));

    T = irho * p / _R;

    mu = _mu0 * pow( T/_T0, _beta );

    kap = mu * _Cp / _Pr;

}

} // namespace reacting

} // namespace suzerain

#endif