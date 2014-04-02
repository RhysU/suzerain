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

#ifndef SUZERAIN_SAMPLES_H
#define SUZERAIN_SAMPLES_H

/** @file
 * <a
 * href="www.boost.org/doc/libs/release/libs/preprocessor/">Boost.Preprocessor</a>
 * definitions driving \ref suzerain::samples. Isolated to ease use, to simplify
 * re-use, and to facilitate debugging.
 */

/**
 * A Boost.Preprocessor sequence of tuples of scalar and vector quantities
 * generally sampled in wave space.  The second element of the tuple is a
 * Boost.Preprocessor sequence of tuples naming and describing each scalar
 * component.  The number of scalar components is the length of that sequence.
 * These are part of the data found in \ref suzerain::samples.
 */
#define SUZERAIN_SAMPLES_WAVE                                                           \
      (( rho,   (( rho,    "Reynolds-averaged density"                               )) \
    ))(( rho_u, (( rho_u,  "Reynolds-averaged streamwise momentum"                   )) \
                (( rho_v,  "Reynolds-averaged wall-normal momentum"                  )) \
                (( rho_w,  "Reynolds-averaged spanwise momentum"                     )) \
    ))(( rho_E, (( rho_E,  "Reynolds-averaged total (intrinsic plus kinetic) energy" )) \
    ))

/**
 * A Boost.Preprocessor sequence of tuples of scalar, vector, and tensor
 * quantities generally sampled in physical space.
 * \copydetails SUZERAIN_SAMPLES_WAVE
 */
#define SUZERAIN_SAMPLES_PHYSICAL                                                                                                              \
      (( E,                (( E,                "Reynolds-averaged total (intrinsic plus kinetic) energy per unit mass"                     )) \
    ))(( p,                (( p,                "Reynolds-averaged pressure"                                                                )) \
    ))(( T,                (( T,                "Reynolds-averaged temperature"                                                             )) \
    ))(( a,                (( a,                "Reynolds-averaged speed of sound"                                                          )) \
    ))(( h0,               (( h0,               "Reynolds-averaged stagnation enthalpy"                                                     )) \
    ))(( H0,               (( H0,               "Reynolds-averaged stagnation enthalpy per unit mass"                                       )) \
    ))(( ke,               (( ke,               "Reynolds-averaged kinetic energy per unit mass"                                            )) \
    ))(( mu,               (( mu,               "Reynolds-averaged dynamic viscosity"                                                       )) \
    ))(( nu,               (( nu,               "Reynolds-averaged kinematic viscosity"                                                     )) \
    ))(( u,                (( u,                "Reynolds-averaged streamwise velocity"                                                     )) \
                           (( v,                "Reynolds-averaged wall-normal velocity"                                                    )) \
                           (( w,                "Reynolds-averaged spanwise velocity"                                                       )) \
    ))(( sym_grad_u,       (( symxx_grad_u,     "Symmetric part (x,x)-component of Reynolds-averaged velocity gradient"                     )) \
                           (( symxy_grad_u,     "Symmetric part (x,y)-component of Reynolds-averaged velocity gradient"                     )) \
                           (( symxz_grad_u,     "Symmetric part (x,z)-component of Reynolds-averaged velocity gradient"                     )) \
                           (( symyy_grad_u,     "Symmetric part (y,y)-component of Reynolds-averaged velocity gradient"                     )) \
                           (( symyz_grad_u,     "Symmetric part (y,z)-component of Reynolds-averaged velocity gradient"                     )) \
                           (( symzz_grad_u,     "Symmetric part (z,z)-component of Reynolds-averaged velocity gradient"                     )) \
    ))(( sym_rho_grad_u,   (( symxx_rho_grad_u, "Symmetric part (x,x)-component of Reynolds-averaged density times velocity gradient"       )) \
                           (( symxy_rho_grad_u, "Symmetric part (x,y)-component of Reynolds-averaged density times velocity gradient"       )) \
                           (( symxz_rho_grad_u, "Symmetric part (x,z)-component of Reynolds-averaged density times velocity gradient"       )) \
                           (( symyy_rho_grad_u, "Symmetric part (y,y)-component of Reynolds-averaged density times velocity gradient"       )) \
                           (( symyz_rho_grad_u, "Symmetric part (y,z)-component of Reynolds-averaged density times velocity gradient"       )) \
                           (( symzz_rho_grad_u, "Symmetric part (z,z)-component of Reynolds-averaged density times velocity gradient"       )) \
    ))(( grad_T,           (( gradx_T,          "Reynolds-averaged x-component of temperature gradient"                                     )) \
                           (( grady_T,          "Reynolds-averaged y-component of temperature gradient"                                     )) \
                           (( gradz_T,          "Reynolds-averaged z-component of temperature gradient"                                     )) \
    ))(( rho_grad_T,       (( rho_gradx_T,      "Reynolds-averaged x-component of density times temperature gradient"                       )) \
                           (( rho_grady_T,      "Reynolds-averaged y-component of density times temperature gradient"                       )) \
                           (( rho_gradz_T,      "Reynolds-averaged z-component of density times temperature gradient"                       )) \
    ))(( tau_colon_grad_u, (( tau_colon_grad_u, "Reynolds-averaged contraction of the viscous stress tensor against the velocity gradient"  )) \
    ))(( tau,              (( tauxx,            "Reynolds-averaged (x,x)-component of the viscous stress tensor"                            )) \
                           (( tauxy,            "Reynolds-averaged (x,y)-component of the viscous stress tensor"                            )) \
                           (( tauxz,            "Reynolds-averaged (x,z)-component of the viscous stress tensor"                            )) \
                           (( tauyy,            "Reynolds-averaged (y,y)-component of the viscous stress tensor"                            )) \
                           (( tauyz,            "Reynolds-averaged (y,z)-component of the viscous stress tensor"                            )) \
                           (( tauzz,            "Reynolds-averaged (z,z)-component of the viscous stress tensor"                            )) \
    ))(( tau_u,            (( tauux,            "Reynolds-averaged x-component of the viscous stress tensor applied to the velocity"        )) \
                           (( tauuy,            "Reynolds-averaged y-component of the viscous stress tensor applied to the velocity"        )) \
                           (( tauuz,            "Reynolds-averaged z-component of the viscous stress tensor applied to the velocity"        )) \
    ))(( p_div_u,          (( p_div_u,          "Reynolds-averaged pressure times divergence of the velocity"                               )) \
    ))(( rho_u_u,          (( rho_u_u,          "Reynolds-averaged (x,x)-component of the momentum times the velocity"                      )) \
                           (( rho_u_v,          "Reynolds-averaged (x,y)-component of the momentum times the velocity"                      )) \
                           (( rho_u_w,          "Reynolds-averaged (x,z)-component of the momentum times the velocity"                      )) \
                           (( rho_v_v,          "Reynolds-averaged (y,y)-component of the momentum times the velocity"                      )) \
                           (( rho_v_w,          "Reynolds-averaged (y,z)-component of the momentum times the velocity"                      )) \
                           (( rho_w_w,          "Reynolds-averaged (z,z)-component of the momentum times the velocity"                      )) \
    ))(( rho_u_u_u,        (( rho_u_u_u,        "Reynolds-averaged (x,x,x)-component of the momentum times the velocity times the velocity" )) \
                           (( rho_u_u_v,        "Reynolds-averaged (x,x,y)-component of the momentum times the velocity times the velocity" )) \
                           (( rho_u_u_w,        "Reynolds-averaged (x,x,z)-component of the momentum times the velocity times the velocity" )) \
                           (( rho_u_v_v,        "Reynolds-averaged (x,y,y)-component of the momentum times the velocity times the velocity" )) \
                           (( rho_u_v_w,        "Reynolds-averaged (x,y,z)-component of the momentum times the velocity times the velocity" )) \
                           (( rho_u_w_w,        "Reynolds-averaged (x,z,z)-component of the momentum times the velocity times the velocity" )) \
                           (( rho_v_v_v,        "Reynolds-averaged (y,y,y)-component of the momentum times the velocity times the velocity" )) \
                           (( rho_v_v_w,        "Reynolds-averaged (y,y,z)-component of the momentum times the velocity times the velocity" )) \
                           (( rho_v_w_w,        "Reynolds-averaged (y,z,z)-component of the momentum times the velocity times the velocity" )) \
                           (( rho_w_w_w,        "Reynolds-averaged (z,z,z)-component of the momentum times the velocity times the velocity" )) \
    ))(( rho_T_u,          (( rho_T_u,          "Reynolds-averaged x-component of the temperature times the momentum"                       )) \
                           (( rho_T_v,          "Reynolds-averaged y-component of the temperature times the momentum"                       )) \
                           (( rho_T_w,          "Reynolds-averaged z-component of the temperature times the momentum"                       )) \
    ))(( rho_mu,           (( rho_mu,           "Reynolds-averaged dynamic viscosity times the density"                                     )) \
    ))(( mu_S,             (( mu_Sxx,           "Reynolds-averaged (x,x)-component of the deviatoric portion of the strain rate"            )) \
                           (( mu_Sxy,           "Reynolds-averaged (x,y)-component of the deviatoric portion of the strain rate"            )) \
                           (( mu_Sxz,           "Reynolds-averaged (x,z)-component of the deviatoric portion of the strain rate"            )) \
                           (( mu_Syy,           "Reynolds-averaged (y,y)-component of the deviatoric portion of the strain rate"            )) \
                           (( mu_Syz,           "Reynolds-averaged (y,z)-component of the deviatoric portion of the strain rate"            )) \
                           (( mu_Szz,           "Reynolds-averaged (z,z)-component of the deviatoric portion of the strain rate"            )) \
    ))(( mu_div_u,         (( mu_div_u,         "Reynolds-averaged dynamic viscosity times divergence of the velocity"                      )) \
    ))(( mu_grad_T,        (( mu_gradx_T,       "Reynolds-averaged x-component of dynamic viscosity times the temperature gradient"         )) \
                           (( mu_grady_T,       "Reynolds-averaged y-component of dynamic viscosity times the temperature gradient"         )) \
                           (( mu_gradz_T,       "Reynolds-averaged z-component of dynamic viscosity times the temperature gradient"         )) \
    ))

/**
 * A Boost.Preprocessor sequence of tuples of scalar and vector quantities
 * generally computed during implicit forcing.
 * \copydetails SUZERAIN_SAMPLES_WAVE
 */
#define SUZERAIN_SAMPLES_IMPLICIT                                                                                                    \
      (( SrhoE,        (( SrhoE,       "Reynolds-averaged total energy contributions due to slow growth forcing"                   )) \
    ))(( Srhou,        (( Srhou,       "Reynolds-averaged streamwise momentum contributions due to slow growth forcing"            )) \
                       (( Srhov,       "Reynolds-averaged wall-normal momentum contributions due to slow growth forcing"           )) \
                       (( Srhow,       "Reynolds-averaged spanwise momentum contributions due to slow growth forcing"              )) \
    ))(( Srho,         (( Srho,        "Reynolds-averaged mass contributions due to slow growth forcing"                           )) \
    ))(( Srhou_dot_u,  (( Srhou_dot_u, "Reynolds-averaged energy contribution due to slow growth forcing work"                     )) \
    ))(( f,            (( fx,          "Reynolds-averaged x-component of the momentum forcing"                                     )) \
                       (( fy,          "Reynolds-averaged y-component of the momentum forcing"                                     )) \
                       (( fz,          "Reynolds-averaged z-component of the momentum forcing"                                     )) \
    ))(( f_dot_u,      (( f_dot_u,     "Reynolds-averaged energy contribution due to momentum forcing work"                        )) \
    ))(( qb,           (( qb,          "Reynolds-averaged volumetric energy forcing"                                               )) \
    ))(( CrhoE,        (( CrhoE,       "Reynolds-averaged total energy contributions due to integral constraints"                  )) \
    ))(( Crhou,        (( Crhou,       "Reynolds-averaged streamwise momentum contributions due to integral constraints"           )) \
                       (( Crhov,       "Reynolds-averaged wall-normal momentum contributions due to integral constraints"          )) \
                       (( Crhow,       "Reynolds-averaged spanwise momentum contributions due to integral constraints"             )) \
    ))(( Crho,         (( Crho,        "Reynolds-averaged mass contributions due to integral constraints"                          )) \
    ))(( Crhou_dot_u,  (( Crhou_dot_u, "Reynolds-averaged energy contribution due to work by integral constraints"                 )) \
    ))(( C2rhoE,       (( C2rhoE,       "Reynolds-averaged squared total energy contributions due to integral constraints"         )) \
    ))(( C2rhou,       (( C2rhou,       "Reynolds-averaged squared streamwise momentum contributions due to integral constraints"  )) \
                       (( C2rhov,       "Reynolds-averaged squared wall-normal momentum contributions due to integral constraints" )) \
                       (( C2rhow,       "Reynolds-averaged squared spanwise momentum contributions due to integral constraints"    )) \
    ))(( C2rho,        (( C2rho,        "Reynolds-averaged squared mass contributions due to integral constraints"                 )) \
    ))(( C2rhou_dot_u, (( C2rhou_dot_u, "Reynolds-averaged squared energy contribution due to work by integral constraints"        )) \
    ))

/**
 * A Boost.Preprocessor sequence of tuples of all sampled quantities.
 * \copydetails SUZERAIN_SAMPLES_WAVE
 */
#define SUZERAIN_SAMPLES      \
    SUZERAIN_SAMPLES_WAVE     \
    SUZERAIN_SAMPLES_PHYSICAL \
    SUZERAIN_SAMPLES_IMPLICIT


/**
 * A Boost.Preprocessor sequence of (name, description) tuples for all scalar
 * components comprising SUZERAIN_SAMPLES.
 *
 * Parameter \c pre is a prefix to be prepended onto each component's name.
 */
#define SUZERAIN_SAMPLES_COMPONENTS(pre)                                              \
    BOOST_PP_SEQ_FOR_EACH(SUZERAIN_SAMPLES_COMPONENTS_HELPER2,pre,                    \
        BOOST_PP_SEQ_FOR_EACH(SUZERAIN_SAMPLES_COMPONENTS_HELPER1,,SUZERAIN_SAMPLES))

#ifndef SUZERAIN_PARSED_BY_DOXYGEN

#define SUZERAIN_SAMPLES_COMPONENTS_HELPER1(r, pre, elem) \
    BOOST_PP_TUPLE_ELEM(2, 1, elem)

#define SUZERAIN_SAMPLES_COMPONENTS_HELPER2(r, pre, tuple) \
    ((BOOST_PP_CAT(pre, BOOST_PP_TUPLE_ELEM(2, 0, tuple)), BOOST_PP_TUPLE_ELEM(2, 1, tuple)))

#endif /* SUZERAIN_PARSED_BY_DOXYGEN */


/**
 * An Boost.Preprocessor-like iteration construct invoking <code>
 *     macro(quantity, component, offset, description)
 * </code> for each component present in \c collection
 * (e.g. \ref SUZERAIN_SAMPLES_PHYSICAL).
 */
#define SUZERAIN_SAMPLES_COMPONENTS_FOR_EACH(macro, collection) \
    BOOST_PP_SEQ_FOR_EACH(SUZERAIN_SAMPLES_COMPONENTS_FOR_EACH_HELPER1, macro, collection)

#ifndef SUZERAIN_PARSED_BY_DOXYGEN

#define SUZERAIN_SAMPLES_COMPONENTS_FOR_EACH_HELPER1(s, macro, quantity_components)  \
    BOOST_PP_SEQ_FOR_EACH_I(SUZERAIN_SAMPLES_COMPONENTS_FOR_EACH_HELPER2,            \
                            (macro, BOOST_PP_TUPLE_ELEM(2, 0, quantity_components)), \
                            BOOST_PP_TUPLE_ELEM(2, 1, quantity_components))

#define SUZERAIN_SAMPLES_COMPONENTS_FOR_EACH_HELPER2(r, macro_quantity, i, component_description) \
        BOOST_PP_TUPLE_ELEM(2, 0, macro_quantity)(                                                \
            BOOST_PP_TUPLE_ELEM(2, 1, macro_quantity),                                            \
            BOOST_PP_TUPLE_ELEM(2, 0, component_description),                                     \
            i,                                                                                    \
            BOOST_PP_TUPLE_ELEM(2, 1, component_description)                                      \
        )                                                                                         \

#endif /* SUZERAIN_PARSED_BY_DOXYGEN */


#endif // SUZERAIN_SAMPLES_H
