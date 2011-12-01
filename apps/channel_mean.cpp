//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2010 The PECOS Development Team
//
// Please see http://pecos.ices.utexas.edu for more information.
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
//
// channel_mean.cpp: Dump mean statistics for one or more restart files
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// FIXME Support multiple sample collections per file
// FIXME Add help-like option listing names and descriptions
// FIXME Allow excluding particular time ranges in the output

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <Eigen/Core>
#include <esio/esio.h>
#include <suzerain/math.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/program_options.hpp>

#include "logging.hpp"
#include "precision.hpp"
#include "channel.hpp"

// Introduce shorthand for common names
using boost::make_shared;
using boost::math::constants::pi;
using boost::numeric_cast;
using boost::shared_ptr;
using std::auto_ptr;
using std::numeric_limits;

// Provided by channel_mean_svnrev.{c,h} to speed recompilation
extern "C" const char revstr[];

#pragma warning(disable:383 1572)

/** Details on the sampled and computed quantities */
namespace quantity {

/** A Boost.Preprocessor sequence of tuples of grid-related details */
#define SEQ_GRID                                     \
    ((t, "Simulation time"))                         \
    ((y, "Wall-normal collocation point locations"))

/** A Boost.Preprocessor sequence of tuples of directly sampled quantities.  */
#define SEQ_SAMPLED                                                                                                       \
    ((bar_rho,              "Reynolds-averaged density"))                                                                 \
    ((bar_rho_u,            "Reynolds-averaged streamwise momentum"))                                                     \
    ((bar_rho_v,            "Reynolds-averaged wall-normal momentum"))                                                    \
    ((bar_rho_w,            "Reynolds-averaged spanwise momentum"))                                                       \
    ((bar_rho_E,            "Reynolds-averaged total (intrinsic plus kinetic) energy"))                                   \
    ((bar_mu,               "Reynolds-averaged dynamic viscosity"))                                                       \
    ((bar_u,                "Reynolds-averaged streamwise velocity"))                                                     \
    ((bar_v,                "Reynolds-averaged wall-normal velocity"))                                                    \
    ((bar_w,                "Reynolds-averaged spanwise velocity"))                                                       \
    ((bar_symxx_rho_grad_u, "Symmetric part (x,x)-component of Reynolds-averaged density times velocity"))                \
    ((bar_symxy_rho_grad_u, "Symmetric part (x,y)-component of Reynolds-averaged density times velocity"))                \
    ((bar_symxz_rho_grad_u, "Symmetric part (x,z)-component of Reynolds-averaged density times velocity"))                \
    ((bar_symyy_rho_grad_u, "Symmetric part (y,y)-component of Reynolds-averaged density times velocity"))                \
    ((bar_symyz_rho_grad_u, "Symmetric part (y,z)-component of Reynolds-averaged density times velocity"))                \
    ((bar_symzz_rho_grad_u, "Symmetric part (z,z)-component of Reynolds-averaged density times velocity"))                \
    ((bar_rho_gradx_T,      "Reynolds-averaged x-component of density times temperature gradient"))                       \
    ((bar_rho_grady_T,      "Reynolds-averaged y-component of density times temperature gradient"))                       \
    ((bar_rho_gradz_T,      "Reynolds-averaged z-component of density times temperature gradient"))                       \
    ((bar_tau_colon_grad_u, "Reynolds-averaged contraction of the viscous stress tensor against the velocity gradient"))  \
    ((bar_tauxx,            "Reynolds-averaged (x,x)-component of the viscous stress tensor"))                            \
    ((bar_tauxy,            "Reynolds-averaged (x,y)-component of the viscous stress tensor"))                            \
    ((bar_tauxz,            "Reynolds-averaged (x,z)-component of the viscous stress tensor"))                            \
    ((bar_tauyy,            "Reynolds-averaged (y,y)-component of the viscous stress tensor"))                            \
    ((bar_tauyz,            "Reynolds-averaged (y,z)-component of the viscous stress tensor"))                            \
    ((bar_tauzz,            "Reynolds-averaged (z,z)-component of the viscous stress tensor"))                            \
    ((bar_tauux,            "Reynolds-averaged x-component of the viscous stress tensor applied to the velocity"))        \
    ((bar_tauuy,            "Reynolds-averaged y-component of the viscous stress tensor applied to the velocity"))        \
    ((bar_tauuz,            "Reynolds-averaged z-component of the viscous stress tensor applied to the velocity"))        \
    ((bar_p_div_u,          "Reynolds-averaged pressure times divergence of the velocity"))                               \
    ((bar_rho_u_u,          "Reynolds-averaged (x,x)-component of the momentum times the velocity"))                      \
    ((bar_rho_u_v,          "Reynolds-averaged (x,y)-component of the momentum times the velocity"))                      \
    ((bar_rho_u_w,          "Reynolds-averaged (x,z)-component of the momentum times the velocity"))                      \
    ((bar_rho_v_v,          "Reynolds-averaged (y,y)-component of the momentum times the velocity"))                      \
    ((bar_rho_v_w,          "Reynolds-averaged (y,z)-component of the momentum times the velocity"))                      \
    ((bar_rho_w_w,          "Reynolds-averaged (z,z)-component of the momentum times the velocity"))                      \
    ((bar_rho_u_u_u,        "Reynolds-averaged (x,x,x)-component of the momentum times the velocity times the velocity")) \
    ((bar_rho_u_u_v,        "Reynolds-averaged (x,x,y)-component of the momentum times the velocity times the velocity")) \
    ((bar_rho_u_u_w,        "Reynolds-averaged (x,x,z)-component of the momentum times the velocity times the velocity")) \
    ((bar_rho_u_v_v,        "Reynolds-averaged (x,y,y)-component of the momentum times the velocity times the velocity")) \
    ((bar_rho_u_v_w,        "Reynolds-averaged (x,y,z)-component of the momentum times the velocity times the velocity")) \
    ((bar_rho_u_w_w,        "Reynolds-averaged (x,z,z)-component of the momentum times the velocity times the velocity")) \
    ((bar_rho_v_v_v,        "Reynolds-averaged (y,y,y)-component of the momentum times the velocity times the velocity")) \
    ((bar_rho_v_v_w,        "Reynolds-averaged (y,y,z)-component of the momentum times the velocity times the velocity")) \
    ((bar_rho_v_w_w,        "Reynolds-averaged (y,z,z)-component of the momentum times the velocity times the velocity")) \
    ((bar_rho_w_w_w,        "Reynolds-averaged (z,z,z)-component of the momentum times the velocity times the velocity")) \
    ((bar_rho_T_u,          "Reynolds-averaged x-component of the temperature times the velocity"))                       \
    ((bar_rho_T_v,          "Reynolds-averaged y-component of the temperature times the velocity"))                       \
    ((bar_rho_T_w,          "Reynolds-averaged z-component of the temperature times the velocity"))                       \
    ((bar_mu_Sxx,           "Reynolds-averaged (x,x)-component of the deviatoric portion of the strain rate"))            \
    ((bar_mu_Sxy,           "Reynolds-averaged (x,y)-component of the deviatoric portion of the strain rate"))            \
    ((bar_mu_Sxz,           "Reynolds-averaged (x,z)-component of the deviatoric portion of the strain rate"))            \
    ((bar_mu_Syy,           "Reynolds-averaged (y,y)-component of the deviatoric portion of the strain rate"))            \
    ((bar_mu_Syz,           "Reynolds-averaged (y,z)-component of the deviatoric portion of the strain rate"))            \
    ((bar_mu_Szz,           "Reynolds-averaged (z,z)-component of the deviatoric portion of the strain rate"))            \
    ((bar_mu_div_u,         "Reynolds-averaged dynamic viscosity times divergence of the velocity"))                      \
    ((bar_mu_gradx_T,       "Reynolds-averaged x-component of dynamic viscosity times the temperature gradient"))         \
    ((bar_mu_grady_T,       "Reynolds-averaged y-component of dynamic viscosity times the temperature gradient"))         \
    ((bar_mu_gradz_T,       "Reynolds-averaged z-component of dynamic viscosity times the temperature gradient"))         \
    ((bar_fx,               "Reynolds-averaged x-component of the momentum forcing"))                                     \
    ((bar_fy,               "Reynolds-averaged y-component of the momentum forcing"))                                     \
    ((bar_fz,               "Reynolds-averaged z-component of the momentum forcing"))                                     \
    ((bar_qb,               "Reynolds-averaged volumetric energy forcing"))                                               \
    ((bar_f_dot_u,          "Reynolds-averaged energy contribution due to momentum forcing work"))

/** A Boost.Preprocessor sequence of tuples of directly sampled quantities.  */
#define SEQ_DERIVED                                                                                                                                                 \
    ((tilde_u,                "Favre-averaged streamwise velocity"))                                                                                                \
    ((tilde_v,                "Favre-averaged wall-normal velocity"))                                                                                               \
    ((tilde_w,                "Favre-averaged spanwise velocity"))                                                                                                  \
    ((tilde_E,                "Favre-averaged total (intrinsic plus kinetic) energy"))                                                                              \
    ((tilde_u_u,              "Favre-averaged (x,x)-component of the velocity times itself"))                                                                       \
    ((tilde_u_v,              "Favre-averaged (x,y)-component of the velocity times itself"))                                                                       \
    ((tilde_u_w,              "Favre-averaged (x,z)-component of the velocity times itself"))                                                                       \
    ((tilde_v_v,              "Favre-averaged (y,y)-component of the velocity times itself"))                                                                       \
    ((tilde_v_w,              "Favre-averaged (y,z)-component of the velocity times itself"))                                                                       \
    ((tilde_w_w,              "Favre-averaged (z,z)-component of the velocity times itself"))                                                                       \
    ((tilde_upp_upp,          "Favre-averaged (x,x)-component of the fluctuating velocity times itself"))                                                           \
    ((tilde_upp_vpp,          "Favre-averaged (x,y)-component of the fluctuating velocity times itself"))                                                           \
    ((tilde_upp_wpp,          "Favre-averaged (x,z)-component of the fluctuating velocity times itself"))                                                           \
    ((tilde_vpp_vpp,          "Favre-averaged (y,y)-component of the fluctuating velocity times itself"))                                                           \
    ((tilde_vpp_wpp,          "Favre-averaged (y,z)-component of the fluctuating velocity times itself"))                                                           \
    ((tilde_wpp_wpp,          "Favre-averaged (z,z)-component of the fluctuating velocity times itself"))                                                           \
    ((tilde_k,                "Turbulent kinetic energy density"))                                                                                                  \
    ((tilde_T,                "Favre-averaged temperature"))                                                                                                        \
    ((tilde_H,                "Favre-averaged total enthalpy"))                                                                                                     \
    ((bar_p,                  "Reynolds-averaged pressure"))                                                                                                        \
    ((bar_tau_colon_grad_upp, "Reynolds-averaged contraction of the viscous stress tensor against the fluctuating velocity gradient"))                              \
    ((tilde_epsilon,          "Turbulent kinetic energy dissipation rate density"))                                                                                 \
    ((bar_upp,                "Reynolds-averaged Favre-fluctuating streamwise velocity"))                                                                           \
    ((bar_vpp,                "Reynolds-averaged Favre-fluctuating wall-normal velocity"))                                                                          \
    ((bar_wpp,                "Reynolds-averaged Favre-fluctuating spanwise velocity"))                                                                             \
    ((bar_f_dot_upp,          "Reynolds-averaged momentum forcing dotted with fluctuating velocity"))                                                               \
    ((bar_tauuppx,            "Reynolds-averaged x-component of the viscous stress tensor times the Favre-fluctuating velocity"))                                   \
    ((bar_tauuppy,            "Reynolds-averaged y-component of the viscous stress tensor times the Favre-fluctuating velocity"))                                   \
    ((bar_tauuppz,            "Reynolds-averaged z-component of the viscous stress tensor times the Favre-fluctuating velocity"))                                   \
    ((bar_pp_div_upp,         "Reynolds-averaged Reynolds-fluctuating pressure times the divergence of the Favre-fluctuating velocity"))                            \
    ((tilde_u_u_u,            "Favre-averaged (x,x,x)-component of velocity times itself times itself"))                                                            \
    ((tilde_u_u_v,            "Favre-averaged (x,x,y)-component of velocity times itself times itself"))                                                            \
    ((tilde_u_u_w,            "Favre-averaged (x,x,z)-component of velocity times itself times itself"))                                                            \
    ((tilde_u_v_v,            "Favre-averaged (x,y,y)-component of velocity times itself times itself"))                                                            \
    ((tilde_u_v_w,            "Favre-averaged (x,y,z)-component of velocity times itself times itself"))                                                            \
    ((tilde_u_w_w,            "Favre-averaged (x,z,z)-component of velocity times itself times itself"))                                                            \
    ((tilde_v_v_v,            "Favre-averaged (y,y,y)-component of velocity times itself times itself"))                                                            \
    ((tilde_v_v_w,            "Favre-averaged (y,y,z)-component of velocity times itself times itself"))                                                            \
    ((tilde_v_w_w,            "Favre-averaged (y,z,z)-component of velocity times itself times itself"))                                                            \
    ((tilde_w_w_w,            "Favre-averaged (z,z,z)-component of velocity times itself times itself"))                                                            \
    ((tilde_u2u,              "Favre-averaged x-component of the square of the velocity times the streamwise velocity"))                                            \
    ((tilde_u2v,              "Favre-averaged y-component of the square of the velocity times the wall-normal velocity"))                                           \
    ((tilde_u2w,              "Favre-averaged z-component of the square of the velocity times the spanwise velocity"))                                              \
    ((tilde_upp2upp,          "Favre-averaged x-component of the square of the fluctuating velocity times the velocity"))                                           \
    ((tilde_upp2vpp,          "Favre-averaged y-component of the square of the fluctuating velocity times the velocity"))                                           \
    ((tilde_upp2wpp,          "Favre-averaged z-component of the square of the fluctuating velocity times the velocity"))                                           \
    ((tilde_T_u,              "Favre-averaged temperature times the streamwise velocity"))                                                                          \
    ((tilde_T_v,              "Favre-averaged temperature times the wall-normal velocity"))                                                                         \
    ((tilde_T_w,              "Favre-averaged temperature times the spanwise velocity"))                                                                            \
    ((tilde_Tpp_upp,          "Favre-averaged fluctuating temperature times the fluctuating streamwise velocity"))                                                  \
    ((tilde_Tpp_vpp,          "Favre-averaged fluctuating temperature times the fluctuating wall-normal velocity"))                                                 \
    ((tilde_Tpp_wpp,          "Favre-averaged fluctuating temperature times the fluctuating spanwise velocity"))                                                    \
    ((tilde_nu,               "Favre-averaged kinematic viscosity"))                                                                                                \
    ((tilde_symxx_grad_u,     "Symmetric part (x,x)-component of Favre-averaged velocity gradient"))                                                                \
    ((tilde_symxy_grad_u,     "Symmetric part (x,y)-component of Favre-averaged velocity gradient"))                                                                \
    ((tilde_symxz_grad_u,     "Symmetric part (x,z)-component of Favre-averaged velocity gradient"))                                                                \
    ((tilde_symyy_grad_u,     "Symmetric part (y,y)-component of Favre-averaged velocity gradient"))                                                                \
    ((tilde_symyz_grad_u,     "Symmetric part (y,z)-component of Favre-averaged velocity gradient"))                                                                \
    ((tilde_symzz_grad_u,     "Symmetric part (z,z)-component of Favre-averaged velocity gradient"))                                                                \
    ((tilde_Sxx,              "Favre-averaged (x,x)-component of the deviatoric portion of the strain rate"))                                                       \
    ((tilde_Sxy,              "Favre-averaged (x,y)-component of the deviatoric portion of the strain rate"))                                                       \
    ((tilde_Sxz,              "Favre-averaged (x,z)-component of the deviatoric portion of the strain rate"))                                                       \
    ((tilde_Syy,              "Favre-averaged (y,y)-component of the deviatoric portion of the strain rate"))                                                       \
    ((tilde_Syz,              "Favre-averaged (y,z)-component of the deviatoric portion of the strain rate"))                                                       \
    ((tilde_Szz,              "Favre-averaged (z,z)-component of the deviatoric portion of the strain rate"))                                                       \
    ((tilde_nupp_Sppxx,       "Favre-averaged (x,x)-component of the fluctuating kinematic viscosity times the deviatoric portion of the fluctuating strain rate")) \
    ((tilde_nupp_Sppxy,       "Favre-averaged (x,y)-component of the fluctuating kinematic viscosity times the deviatoric portion of the fluctuating strain rate")) \
    ((tilde_nupp_Sppxz,       "Favre-averaged (x,z)-component of the fluctuating kinematic viscosity times the deviatoric portion of the fluctuating strain rate")) \
    ((tilde_nupp_Sppyy,       "Favre-averaged (y,y)-component of the fluctuating kinematic viscosity times the deviatoric portion of the fluctuating strain rate")) \
    ((tilde_nupp_Sppyz,       "Favre-averaged (y,z)-component of the fluctuating kinematic viscosity times the deviatoric portion of the fluctuating strain rate")) \
    ((tilde_nupp_Sppzz,       "Favre-averaged (z,z)-component of the fluctuating kinematic viscosity times the deviatoric portion of the fluctuating strain rate")) \
    ((tilde_nupp_div_upp,     "Favre-averaged fluctuating kinematic viscosity times the divergence of the fluctuating velocity"))                                   \
    ((tilde_nupp_gradxTpp,    "Favre-averaged x-component of the kinematic viscosity times the fluctuating temperature gradient"))                                  \
    ((tilde_nupp_gradyTpp,    "Favre-averaged y-component of the kinematic viscosity times the fluctuating temperature gradient"))                                  \
    ((tilde_nupp_gradzTpp,    "Favre-averaged z-component of the kinematic viscosity times the fluctuating temperature gradient"))

/** A Boost.Preprocessor sequence of tuples of locally computed quantities */
#define SEQ_LOCALS                                                                     \
    ((local_a,    "Local speed of sound formed via sqrt(tilde_T)"))                    \
    ((local_Ma,   "Local Mach number formed via bar_u / local_a"))                     \
    ((local_Mat,  "Local turbulent Mach number formed via sqrt(2*tilde_k) / local_a")) \
    ((local_Re,   "Local Reynolds number formed via bar_u L / tilde_nu"))

/**
 * A Boost.Preprocessor sequence of tuples of stationary, time-invariant
 * equation residuals.  In the unsteady case, these are nothing but the time
 * derivative of the state variables.
 */
#define SEQ_RESIDUAL                                                                                         \
    ((bar_rho__t,   "Residual of stationary Favre-averaged density equation"))                               \
    ((bar_rho_u__t, "Residual of stationary Favre-averaged streamwise momentum equation"))                   \
    ((bar_rho_v__t, "Residual of stationary Favre-averaged wall-normal momentum equation"))                  \
    ((bar_rho_w__t, "Residual of stationary Favre-averaged spanwise momentum equation"))                     \
    ((bar_rho_E__t, "Residual of stationary Favre-averaged total (intrinsic plus kinetic) energy equation")) \
    ((bar_rho_k__t, "Residual of stationary Favre-averaged turbulent kinetic energy equation"))

// Helpers for working with sequences of tuples
#define NAME(tuple)  BOOST_PP_TUPLE_ELEM(2,0,tuple)
#define SNAME(tuple) BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(2,0,tuple))
#define DESC(tuple)  BOOST_PP_TUPLE_ELEM(2,1,tuple)

// Prepare sequences of first wall-normal derivatives
#define TRANSFORM_Y(r, data, tuple)                                                 \
    (BOOST_PP_CAT(NAME(tuple),__y),  "Wall-normal first derivative of "DESC(tuple))
#define SEQ_SAMPLED_Y BOOST_PP_SEQ_TRANSFORM(TRANSFORM_Y,,SEQ_SAMPLED)
#define SEQ_DERIVED_Y BOOST_PP_SEQ_TRANSFORM(TRANSFORM_Y,,SEQ_DERIVED)

// Prepare sequences of second wall-normal derivatives
#define TRANSFORM_YY(r, data, tuple)                                                 \
    (BOOST_PP_CAT(NAME(tuple),__yy), "Wall-normal second derivative of "DESC(tuple))
#define SEQ_SAMPLED_YY BOOST_PP_SEQ_TRANSFORM(TRANSFORM_YY,,SEQ_SAMPLED)
#define SEQ_DERIVED_YY BOOST_PP_SEQ_TRANSFORM(TRANSFORM_YY,,SEQ_DERIVED)

// Building a Boost.Preprocessor sequence of all data of interest
//   #define SEQ_ALL
//       SEQ_GRID SEQ_SAMPLED    SEQ_DERIVED    SEQ_LOCALS SEQ_RESIDUAL
//                SEQ_SAMPLED_Y  SEQ_DERIVED_Y
//                SEQ_SAMPLED_YY SEQ_DERIVED_YY
// appears to be possible as the sequence has more than 256 elements.
// Instead, we have to invoke on each component sequence in turn.

    /** Number of scalar quantities processed as a wall-normal function */
    static const std::size_t count = BOOST_PP_SEQ_SIZE(SEQ_GRID)
                                   + BOOST_PP_SEQ_SIZE(SEQ_SAMPLED)
                                   + BOOST_PP_SEQ_SIZE(SEQ_DERIVED)
                                   + BOOST_PP_SEQ_SIZE(SEQ_LOCALS)
                                   + BOOST_PP_SEQ_SIZE(SEQ_RESIDUAL)
                                   + BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y)
                                   + BOOST_PP_SEQ_SIZE(SEQ_DERIVED_Y)
                                   + BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_YY)
                                   + BOOST_PP_SEQ_SIZE(SEQ_DERIVED_YY);

#define OP(s, data, tuple) NAME(tuple)
    /** Provides named index constants for each quantity */
    enum index {
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_GRID))      ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_SAMPLED))   ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_DERIVED))   ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_LOCALS))    ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_RESIDUAL))  ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_SAMPLED_Y)) ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_DERIVED_Y)) ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_SAMPLED_YY)),
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_DERIVED_YY))
    };
#undef OP

#define OP(s, data, tuple) SNAME(tuple)
    /** Provides names indexed on \ref index */
    static const char * name[count] = {
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_GRID))      ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_SAMPLED))   ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_DERIVED))   ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_LOCALS))    ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_RESIDUAL))  ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_SAMPLED_Y)) ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_DERIVED_Y)) ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_SAMPLED_YY)),
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_DERIVED_YY))
    };
#undef OP

#define OP(s, data, tuple) DESC(tuple)
    /** Provides human-readable descriptions indexed on \ref index */
    static const char * desc[count] = {
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_GRID))      ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_SAMPLED))   ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_DERIVED))   ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_LOCALS))    ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_RESIDUAL))  ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_SAMPLED_Y)) ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_DERIVED_Y)) ,
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_SAMPLED_YY)),
        BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(OP, , SEQ_DERIVED_YY))
    };
#undef OP

    /** Type used to store all quantities in a single contiguous region */
    typedef Eigen::Array<real_t, Eigen::Dynamic, count> storage_type;

    /**
     * Map type used to manage and sort samples across a time series.
     */
    typedef boost::ptr_map<real_t, storage_type> storage_map_type;

    /** Output names in a manner suitable for \ref iofmt-ed column header */
    static void write_names(std::ostream &out)
    {
        for (size_t i = 0; i < quantity::count; ++i) {  // Headings
            out << std::setw(numeric_limits<real_t>::digits10 + 11)
                << quantity::name[i];
            if (i < quantity::count - 1) out << ", ";
        }
        out << std::endl;
    }

    /** Used for formatting output data to match \ref quantity::write_names. */
    static const Eigen::IOFormat iofmt(
            Eigen::FullPrecision, 0, ",     ", "\n", "    ");

} // namespace quantity


/**
 * Compute all quantities from namespace \ref quantity using the sample
 * collections present in \c filename using the wall-normal discretization from
 * \c filename.
 *
 * @param filename To be loaded.
 * @param i_b      If <tt>!!i_b</tt> on entry, after computation interpolate
 *                 the results onto the collocation points given by \c i_b.
 *                 Otherwise, perform no additional interpolation and
 *                 update \c i_b with the basis in \c filename.
 * @param i_bop    Handled identically to \c i_b.
 * @param i_boplu  Handled identically to \c i_b.
 *
 * @return A map of quantities keyed on the nondimensional simulation time.
 */
static quantity::storage_map_type process(
        const std::string& filename,
        shared_ptr<suzerain::bspline>& i_b,
        shared_ptr<suzerain::bsplineop>& i_bop,
        shared_ptr<suzerain::bsplineop_lu>& i_boplu);

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);                         // Initialize MPI
    atexit((void (*) ()) MPI_Finalize);             // Finalize MPI at exit
    logging::initialize(MPI_COMM_WORLD,             // Initialize logging
                        channel::log4cxx_config_console);

    // Process incoming arguments
    std::vector<std::string> restart_files;
    std::string outfile;
    bool use_stdout = false;
    bool use_hdf5   = false;
    {
        suzerain::ProgramOptions options(
                "Suzerain-based channel mean quantity computations",
                "RESTART-OR-SAMPLE-HDF5-FILE...",
"Invocable in three distinct ways:\n"
"\t1) channel_mean               INFILE.h5 ...\n"
"\t2) channel_mean -s            INFILE.h5 ...\n"
"\t3) channel_mean -o OUTFILE.h5 INFILE.h5 ...\n"
"\n"
"The first way processes each INFILE.h5 in turn outputting a corresponding "
"INFILE.mean containing a comma-separated table of means from the first "
"sample collection in the file.  The second way sends the data from all "
"sample collections to standard output sorted according to the simulation "
"time with a blank line separating adjacent times.  The third way outputs "
"a single HDF5 file called OUTFILE.h5 containing all sample collections. "
                , revstr);
        options.add_options()
            ("stdout,s", "Write results to standard output?")
            ("outfile,o", boost::program_options::value(&outfile),
                         "Write results to an HDF5 output file")
            ;
        restart_files = options.process(argc, argv);
        switch (options.verbose()) {
            case 0:                   break;
            case 1:  DEBUG0_ENABLE(); break;
            default: TRACE0_ENABLE(); break;
        }
        switch (options.verbose_all()) {
            case 0:                   break;
            case 1:  DEBUG_ENABLE();  break;
            default: TRACE_ENABLE();  break;
        }
        use_stdout = options.variables().count("stdout");
        use_hdf5   = options.variables().count("outfile");
    }

    // Ensure that we're running in a single processor environment
    if (suzerain::mpi::comm_size(MPI_COMM_WORLD) > 1) {
        FATAL(argv[0] << " only intended to run on single rank");
        return EXIT_FAILURE;
    }

    // Processing differs slightly when done file-by-file versus
    // aggregated across multiple files...
    if (!use_hdf5 && !use_stdout) {

        BOOST_FOREACH(const std::string& filename, restart_files) {

            // Numerics details reset on each new file
            shared_ptr<suzerain::bspline> b;
            shared_ptr<suzerain::bsplineop> bop;
            shared_ptr<suzerain::bsplineop_lu> boplu;

            // Load data from filename
            quantity::storage_map_type data = process(filename, b, bop, boplu);

            // Save quantities to `basename filename .h5`.mean
            static const char suffix[] = ".h5";
            const std::size_t suffix_len = sizeof(suffix) - 1;
            std::string outname;
            if (filename.rfind(suffix) == filename.length() - suffix_len) {
                outname = filename.substr(
                        0, filename.length() - suffix_len) + ".mean";
            } else {
                outname = filename + ".mean";
            }
            DEBUG("Saving nondimensional quantities to " << outname);

            // Write header followed by data values separated by blanks
            std::ofstream ofs(outname.c_str());
            quantity::write_names(ofs);
            BOOST_FOREACH(quantity::storage_map_type::value_type i, data) {
                ofs << i->second->format(quantity::iofmt) << std::endl
                    << std::endl;
            }
            ofs.close();
        }

    } else {

        // A single map of data is stored across all files.  Because the map
        // key is the simulation time, we automatically get a well-ordered,
        // unique set of data across all files.
        quantity::storage_map_type pool;

        // Numerics details preserved across multiple files
        // Notice that the first file determines the grid.
        shared_ptr<suzerain::bspline> b;
        shared_ptr<suzerain::bsplineop> bop;
        shared_ptr<suzerain::bsplineop_lu> boplu;

        BOOST_FOREACH(const std::string& filename, restart_files) {

            // Load data from filename
            quantity::storage_map_type data = process(filename, b, bop, boplu);

            // Transfer data into larger pool (which erases it from data)
            pool.transfer(data);

            // Warn on any duplicate values which were not transfered
            BOOST_FOREACH(quantity::storage_map_type::value_type i, data) {
                WARN0("Duplicate sample time "
                      << i->first << " from " << filename << " ignored");
            }

        }

        if (use_stdout) {
            // Write header followed by data values separated by blanks
            quantity::write_names(std::cout);
            BOOST_FOREACH(quantity::storage_map_type::value_type i, pool) {
                std::cout << i->second->format(quantity::iofmt) << std::endl
                          << std::endl;
            }
        }

        if (use_hdf5) {
            FATAL("HDF5 output not yet implemented"); // FIXME Implement
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

static quantity::storage_map_type process(
        const std::string& filename,
        shared_ptr<suzerain::bspline>& i_b,
        shared_ptr<suzerain::bsplineop>& i_bop,
        shared_ptr<suzerain::bsplineop_lu>& i_boplu)
{
    using quantity::storage_type;
    using quantity::storage_map_type;

    quantity::storage_map_type retval;

    // Create a file-specific ESIO handle using RAII
    shared_ptr<boost::remove_pointer<esio_handle>::type> h(
            esio_handle_initialize(MPI_COMM_WORLD), esio_handle_finalize);

    DEBUG("Loading file " << filename);
    esio_file_open(h.get(), filename.c_str(), 0 /* read-only */);

    // Load time, scenario, grid, and B-spline details from file
    real_t time;
    suzerain::problem::ScenarioDefinition<real_t> scenario;
    suzerain::problem::GridDefinition grid;
    shared_ptr<suzerain::bspline> b;
    shared_ptr<suzerain::bsplineop> bop;
    channel::load_time(h.get(), time);
    channel::load(h.get(), scenario);
    channel::load(h.get(), grid);
    channel::load(h.get(), b, bop);
    assert(b->n() == grid.N.y());

    // Compute factorized mass matrix
    shared_ptr<suzerain::bsplineop_lu> boplu
        = make_shared<suzerain::bsplineop_lu>(*bop.get());
    boplu->form_mass(*bop.get());

    // Load samples as coefficients
    auto_ptr<channel::mean> m(new channel::mean(time, b->n()));
    channel::load(h.get(), *m.get());
    if (m->t >= 0) {
        DEBUG0("Successfully loaded sample collection from " << filename);
    } else {
        WARN0("No valid sample collection found in " << filename);
        return retval;
    }

    // Convert samples into collocation point values in s
    auto_ptr<storage_type> s(new storage_type(b->n(),
                (storage_type::Index) storage_type::ColsAtCompileTime));
#ifndef NDEBUG
    s->fill(numeric_limits<real_t>::quiet_NaN());  // ++paranoia
#endif

#define ACCUMULATE(coeff_name, coeff_col, point_name)                 \
    bop->accumulate(0, 1.0, m->coeff_name().col(coeff_col).data(), 1, \
                    0.0, s->col(quantity::point_name).data(),   1)
    ACCUMULATE(rho,                      0, bar_rho               );
    ACCUMULATE(rhou,                     0, bar_rho_u             );
    ACCUMULATE(rhou,                     1, bar_rho_v             );
    ACCUMULATE(rhou,                     2, bar_rho_w             );
    ACCUMULATE(rhoe,                     0, bar_rho_E             );
    ACCUMULATE(mu,                       0, bar_mu                );
    ACCUMULATE(u,                        0, bar_u                 );
    ACCUMULATE(u,                        1, bar_v                 );
    ACCUMULATE(u,                        2, bar_w                 );
    ACCUMULATE(sym_rho_grad_u,           0, bar_symxx_rho_grad_u  );
    ACCUMULATE(sym_rho_grad_u,           1, bar_symxy_rho_grad_u  );
    ACCUMULATE(sym_rho_grad_u,           2, bar_symxz_rho_grad_u  );
    ACCUMULATE(sym_rho_grad_u,           3, bar_symyy_rho_grad_u  );
    ACCUMULATE(sym_rho_grad_u,           4, bar_symyz_rho_grad_u  );
    ACCUMULATE(sym_rho_grad_u,           5, bar_symzz_rho_grad_u  );
    ACCUMULATE(rho_grad_T,               0, bar_rho_gradx_T       );
    ACCUMULATE(rho_grad_T,               1, bar_rho_grady_T       );
    ACCUMULATE(rho_grad_T,               2, bar_rho_gradz_T       );
    ACCUMULATE(tau_colon_grad_u,         0, bar_tau_colon_grad_u  );
    ACCUMULATE(tau,                      0, bar_tauxx             );
    ACCUMULATE(tau,                      1, bar_tauxy             );
    ACCUMULATE(tau,                      2, bar_tauxz             );
    ACCUMULATE(tau,                      3, bar_tauyy             );
    ACCUMULATE(tau,                      4, bar_tauyz             );
    ACCUMULATE(tau,                      5, bar_tauzz             );
    ACCUMULATE(tau_u,                    0, bar_tauux             );
    ACCUMULATE(tau_u,                    1, bar_tauuy             );
    ACCUMULATE(tau_u,                    2, bar_tauuz             );
    ACCUMULATE(p_div_u,                  0, bar_p_div_u           );
    ACCUMULATE(rho_u_otimes_u,           0, bar_rho_u_u           );
    ACCUMULATE(rho_u_otimes_u,           1, bar_rho_u_v           );
    ACCUMULATE(rho_u_otimes_u,           2, bar_rho_u_w           );
    ACCUMULATE(rho_u_otimes_u,           3, bar_rho_v_v           );
    ACCUMULATE(rho_u_otimes_u,           4, bar_rho_v_w           );
    ACCUMULATE(rho_u_otimes_u,           5, bar_rho_w_w           );
    ACCUMULATE(rho_u_otimes_u_otimes_u,  0, bar_rho_u_u_u         );
    ACCUMULATE(rho_u_otimes_u_otimes_u,  1, bar_rho_u_u_v         );
    ACCUMULATE(rho_u_otimes_u_otimes_u,  2, bar_rho_u_u_w         );
    ACCUMULATE(rho_u_otimes_u_otimes_u,  3, bar_rho_u_v_v         );
    ACCUMULATE(rho_u_otimes_u_otimes_u,  4, bar_rho_u_v_w         );
    ACCUMULATE(rho_u_otimes_u_otimes_u,  5, bar_rho_u_w_w         );
    ACCUMULATE(rho_u_otimes_u_otimes_u,  6, bar_rho_v_v_v         );
    ACCUMULATE(rho_u_otimes_u_otimes_u,  7, bar_rho_v_v_w         );
    ACCUMULATE(rho_u_otimes_u_otimes_u,  8, bar_rho_v_w_w         );
    ACCUMULATE(rho_u_otimes_u_otimes_u,  9, bar_rho_w_w_w         );
    ACCUMULATE(rho_T_u,                  0, bar_rho_T_u           );
    ACCUMULATE(rho_T_u,                  1, bar_rho_T_v           );
    ACCUMULATE(rho_T_u,                  2, bar_rho_T_w           );
    ACCUMULATE(mu_S,                     0, bar_mu_Sxx            );
    ACCUMULATE(mu_S,                     1, bar_mu_Sxy            );
    ACCUMULATE(mu_S,                     2, bar_mu_Sxz            );
    ACCUMULATE(mu_S,                     3, bar_mu_Syy            );
    ACCUMULATE(mu_S,                     4, bar_mu_Syz            );
    ACCUMULATE(mu_S,                     5, bar_mu_Szz            );
    ACCUMULATE(mu_div_u,                 0, bar_mu_div_u          );
    ACCUMULATE(mu_grad_T,                0, bar_mu_gradx_T        );
    ACCUMULATE(mu_grad_T,                1, bar_mu_grady_T        );
    ACCUMULATE(mu_grad_T,                2, bar_mu_gradz_T        );
    ACCUMULATE(f,                        0, bar_fx                );
    ACCUMULATE(f,                        1, bar_fy                );
    ACCUMULATE(f,                        2, bar_fz                );
    ACCUMULATE(qb,                       0, bar_qb                );
    ACCUMULATE(f_dot_u,                  0, bar_f_dot_u           );
#undef ACCUMULATE

    // Store time and collocation points into s.
    // Not strictly necessary, but very useful for textual output
    // and as a sanity check of any later grid projection.
    s->col(quantity::t).fill(m->t);
    for (int i = 0; i < b->n(); ++i)
        s->col(quantity::y)[i] = b->collocation_point(i);

    // Free coefficient-related resources
    m.reset();

    // Introduce shorthand for constants
    const real_t Ma    = scenario.Ma;
    const real_t Re    = scenario.Re;    SUZERAIN_UNUSED(Re); // FIXME
    const real_t Pr    = scenario.Pr;    SUZERAIN_UNUSED(Pr); // FIXME
    const real_t gamma = scenario.gamma;

    // Shorthand for referring to a particular column
#define C(name) s->col(quantity::name)

    // Shorthand for computing derivatives within a particular __y, __yy
#define D(name)                                      \
    C(name##__y) = C(name);                          \
    boplu->solve(1, C(name##__y).data(), 1, b->n()); \
    C(name##__yy) = C(name##__y);                    \
    bop->apply(1, 1.0, C(name##__y).data(),  1);     \
    bop->apply(2, 1.0, C(name##__yy).data(), 1)

    // Computations following "Sampling logistics" in writeups
    C(tilde_u) = C(bar_rho_u)/C(bar_rho);
    C(tilde_v) = C(bar_rho_v)/C(bar_rho);
    C(tilde_w) = C(bar_rho_w)/C(bar_rho);
    C(tilde_E) = C(bar_rho_E)/C(bar_rho);
    C(tilde_u_u) = C(bar_rho_u_u)/C(bar_rho);
    C(tilde_u_v) = C(bar_rho_u_v)/C(bar_rho);
    C(tilde_u_w) = C(bar_rho_u_w)/C(bar_rho);
    C(tilde_v_v) = C(bar_rho_v_v)/C(bar_rho);
    C(tilde_v_w) = C(bar_rho_v_w)/C(bar_rho);
    C(tilde_w_w) = C(bar_rho_w_w)/C(bar_rho);
    C(tilde_upp_upp) = C(tilde_u_u) - C(tilde_u).square();
    C(tilde_upp_vpp) = C(tilde_u_v) - C(tilde_u)*C(tilde_v);
    C(tilde_upp_wpp) = C(tilde_u_w) - C(tilde_u)*C(tilde_w);
    C(tilde_vpp_vpp) = C(tilde_v_v) - C(tilde_v).square();
    C(tilde_vpp_wpp) = C(tilde_v_w) - C(tilde_v)*C(tilde_w);
    C(tilde_wpp_wpp) = C(tilde_w_w) - C(tilde_w).square();
    C(tilde_k) = (C(tilde_upp_upp) + C(tilde_vpp_vpp) + C(tilde_wpp_wpp)) / 2;
    C(tilde_T) = (gamma*(gamma-1))*C(tilde_E) - Ma*Ma*(
                  (  C(tilde_u).square()
                   + C(tilde_v).square()
                   + C(tilde_w).square() ) / 2
               +  C(tilde_k)
               );
    C(tilde_H) = C(tilde_E) + C(tilde_T)/gamma;
    C(bar_p)   = C(bar_rho) * C(tilde_T)/gamma;
    D(tilde_u);  // Form derivatives
    D(tilde_v);  // Form derivatives
    D(tilde_w);  // Form derivatives
    C(bar_tau_colon_grad_upp) = C(bar_tau_colon_grad_u)
                              - C(bar_tauxy)*C(tilde_u__y)
                              - C(bar_tauyy)*C(tilde_v__y)
                              - C(bar_tauyz)*C(tilde_w__y);
    C(tilde_epsilon) = C(bar_tau_colon_grad_upp)/C(bar_rho);
    C(bar_upp) = C(bar_u) - C(tilde_u);
    C(bar_vpp) = C(bar_v) - C(tilde_v);
    C(bar_wpp) = C(bar_w) - C(tilde_w);
    C(bar_f_dot_upp) = C(bar_f_dot_u)
                     - C(bar_fx) * C(tilde_u)
                     - C(bar_fy) * C(tilde_v)
                     - C(bar_fz) * C(tilde_w);
    C(bar_tauuppx) = C(bar_tauux)
                   - C(bar_tauxx)*C(tilde_u)
                   - C(bar_tauxy)*C(tilde_v)
                   - C(bar_tauxz)*C(tilde_w);
    C(bar_tauuppy) = C(bar_tauuy)
                   - C(bar_tauxy)*C(tilde_u)
                   - C(bar_tauyy)*C(tilde_v)
                   - C(bar_tauyz)*C(tilde_w);
    C(bar_tauuppz) = C(bar_tauuz)
                   - C(bar_tauxz)*C(tilde_u)
                   - C(bar_tauxy)*C(tilde_v)
                   - C(bar_tauxz)*C(tilde_w);
    D(bar_v);  // Form derivatives
    C(bar_pp_div_upp) = C(bar_p_div_u) - C(bar_p)*C(bar_v__y);
    C(tilde_u_u_u) = C(bar_rho_u_u_u)/C(bar_rho);
    C(tilde_u_u_v) = C(bar_rho_u_u_v)/C(bar_rho);
    C(tilde_u_u_w) = C(bar_rho_u_u_w)/C(bar_rho);
    C(tilde_u_v_v) = C(bar_rho_u_v_v)/C(bar_rho);
    C(tilde_u_v_w) = C(bar_rho_u_v_w)/C(bar_rho);
    C(tilde_u_w_w) = C(bar_rho_u_w_w)/C(bar_rho);
    C(tilde_v_v_v) = C(bar_rho_v_v_v)/C(bar_rho);
    C(tilde_v_v_w) = C(bar_rho_v_v_w)/C(bar_rho);
    C(tilde_v_w_w) = C(bar_rho_v_w_w)/C(bar_rho);
    C(tilde_w_w_w) = C(bar_rho_w_w_w)/C(bar_rho);
    C(tilde_u2u) = C(tilde_u_u_u) + C(tilde_u_v_v) + C(tilde_u_w_w);
    C(tilde_u2v) = C(tilde_u_u_v) + C(tilde_v_v_v) + C(tilde_v_w_w);
    C(tilde_u2w) = C(tilde_u_u_w) + C(tilde_v_v_w) + C(tilde_w_w_w);
    C(tilde_upp2upp) = C(tilde_u2u)
                     - (C(tilde_u_u) + C(tilde_v_v) + C(tilde_w_w))*C(tilde_u)
                     - 2*(C(tilde_u_u)*C(tilde_u) + C(tilde_u_v)*C(tilde_v) + C(tilde_u_w)*C(tilde_w))
                     + 2*(C(tilde_u).square() + C(tilde_v).square() + C(tilde_w).square())*C(tilde_u);
    C(tilde_upp2vpp) = C(tilde_u2v)
                     - (C(tilde_u_u) + C(tilde_v_v) + C(tilde_w_w))*C(tilde_v)
                     - 2*(C(tilde_u_v)*C(tilde_u) + C(tilde_v_v)*C(tilde_v) + C(tilde_v_w)*C(tilde_w))
                     + 2*(C(tilde_u).square() + C(tilde_v).square() + C(tilde_w).square())*C(tilde_v);
    C(tilde_upp2wpp) = C(tilde_u2w)
                     - (C(tilde_u_u) + C(tilde_v_v) + C(tilde_w_w))*C(tilde_w)
                     - 2*(C(tilde_u_w)*C(tilde_u) + C(tilde_v_w)*C(tilde_v) + C(tilde_w_w)*C(tilde_w))
                     + 2*(C(tilde_u).square() + C(tilde_v).square() + C(tilde_w).square())*C(tilde_w);
    C(tilde_T_u) = C(bar_rho_T_u)/C(bar_rho);
    C(tilde_T_v) = C(bar_rho_T_v)/C(bar_rho);
    C(tilde_T_w) = C(bar_rho_T_w)/C(bar_rho);
    C(tilde_Tpp_upp) = C(tilde_T_u) - C(tilde_T)*C(tilde_u);
    C(tilde_Tpp_vpp) = C(tilde_T_v) - C(tilde_T)*C(tilde_v);
    C(tilde_Tpp_wpp) = C(tilde_T_w) - C(tilde_T)*C(tilde_w);
    C(tilde_nu) = C(bar_mu)/C(bar_rho);
    C(tilde_symxx_grad_u) = C(bar_symxx_rho_grad_u)/C(bar_rho);
    C(tilde_symxy_grad_u) = C(bar_symxy_rho_grad_u)/C(bar_rho);
    C(tilde_symxz_grad_u) = C(bar_symxz_rho_grad_u)/C(bar_rho);
    C(tilde_symyy_grad_u) = C(bar_symyy_rho_grad_u)/C(bar_rho);
    C(tilde_symyz_grad_u) = C(bar_symyz_rho_grad_u)/C(bar_rho);
    C(tilde_symzz_grad_u) = C(bar_symzz_rho_grad_u)/C(bar_rho);
    C(tilde_Sxx) = C(tilde_symxx_grad_u)
                 - (C(tilde_symxx_grad_u) + C(tilde_symyy_grad_u) + C(tilde_symzz_grad_u)) / 3;
    C(tilde_Sxy) = C(tilde_symxy_grad_u);
    C(tilde_Sxz) = C(tilde_symxz_grad_u);
    C(tilde_Syy) = C(tilde_symyy_grad_u)
                 - (C(tilde_symxx_grad_u) + C(tilde_symyy_grad_u) + C(tilde_symzz_grad_u)) / 3;
    C(tilde_Syz) = C(tilde_symyz_grad_u);
    C(tilde_Szz) = C(tilde_symzz_grad_u)
                 - (C(tilde_symxx_grad_u) + C(tilde_symyy_grad_u) + C(tilde_symzz_grad_u)) / 3;
    C(tilde_nupp_Sppxx) = C(bar_mu_Sxx)/C(bar_rho) - C(tilde_nu)*C(tilde_Sxx);
    C(tilde_nupp_Sppxy) = C(bar_mu_Sxy)/C(bar_rho) - C(tilde_nu)*C(tilde_Sxy);
    C(tilde_nupp_Sppxz) = C(bar_mu_Sxz)/C(bar_rho) - C(tilde_nu)*C(tilde_Sxz);
    C(tilde_nupp_Sppyy) = C(bar_mu_Syy)/C(bar_rho) - C(tilde_nu)*C(tilde_Syy);
    C(tilde_nupp_Sppyz) = C(bar_mu_Syz)/C(bar_rho) - C(tilde_nu)*C(tilde_Syz);
    C(tilde_nupp_Sppzz) = C(bar_mu_Szz)/C(bar_rho) - C(tilde_nu)*C(tilde_Szz);
    C(tilde_nupp_div_upp) = C(bar_mu_div_u)/C(bar_rho)
                          - C(tilde_nu)*(C(tilde_symxx_grad_u) + C(tilde_symyy_grad_u) + C(tilde_symzz_grad_u));
    C(tilde_nupp_gradxTpp) = (C(bar_mu_gradx_T) - C(tilde_nu)*C(bar_rho_gradx_T))/C(bar_rho);
    C(tilde_nupp_gradyTpp) = (C(bar_mu_grady_T) - C(tilde_nu)*C(bar_rho_grady_T))/C(bar_rho);
    C(tilde_nupp_gradzTpp) = (C(bar_mu_gradz_T) - C(tilde_nu)*C(bar_rho_gradz_T))/C(bar_rho);

    // Computations of local quantities (see descriptions for definitions)
    C(local_a)   = C(tilde_T).sqrt();
    C(local_Ma)  = C(bar_u) / C(local_a);
    C(local_Mat) = (std::sqrt(real_t(2))*C(tilde_k).sqrt()) / C(local_a);
    C(local_Re)  = C(bar_u) /* L = 1 */ / C(tilde_nu);

    // Differentiate SAMPLED
    // Uses that bar_rho{,__y,__yy} is the first entry in SAMPLED{,_Y,_YY}
    s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y)>(quantity::bar_rho__y)
        = s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED)>(quantity::bar_rho);
    boplu->solve(BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y),
            s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y)>(quantity::bar_rho__y).data(),
            1, b->n());
    s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_YY)>(quantity::bar_rho__yy)
        = s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y)>(quantity::bar_rho__y);
    bop->apply(1, BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y), 1.0,
            s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y)>(quantity::bar_rho__y).data(),
            1, b->n());
    bop->apply(2, BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_Y), 1.0,
            s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_SAMPLED_YY)>(quantity::bar_rho__yy).data(),
            1, b->n());

    // Differentiate DERIVED
    // Uses that tilde_u{,__y,__yy} is the first entry in DERIVED{,_Y,_YY}
    s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_DERIVED_Y)>(quantity::tilde_u__y)
        = s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_DERIVED)>(quantity::tilde_u);
    boplu->solve(BOOST_PP_SEQ_SIZE(SEQ_DERIVED_Y),
            s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_DERIVED_Y)>(quantity::tilde_u__y).data(),
            1, b->n());
    s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_DERIVED_YY)>(quantity::tilde_u__yy)
        = s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_DERIVED_Y)>(quantity::tilde_u__y);
    bop->apply(1, BOOST_PP_SEQ_SIZE(SEQ_DERIVED_Y), 1.0,
            s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_DERIVED_Y)>(quantity::tilde_u__y).data(),
            1, b->n());
    bop->apply(2, BOOST_PP_SEQ_SIZE(SEQ_DERIVED_Y), 1.0,
            s->middleCols<BOOST_PP_SEQ_SIZE(SEQ_DERIVED_YY)>(quantity::tilde_u__yy).data(),
            1, b->n());

    // Computation of Favre-averaged equation residuals following writeup
    // FIXME Implement the stationary equation residuals
    C(bar_rho__t).fill(numeric_limits<real_t>::quiet_NaN());
    C(bar_rho_u__t).fill(numeric_limits<real_t>::quiet_NaN());
    C(bar_rho_v__t).fill(numeric_limits<real_t>::quiet_NaN());
    C(bar_rho_w__t).fill(numeric_limits<real_t>::quiet_NaN());
    C(bar_rho_E__t).fill(numeric_limits<real_t>::quiet_NaN());
    C(bar_rho_k__t).fill(numeric_limits<real_t>::quiet_NaN());

#undef C

    // Use b and friends if i_b was not supplied by the caller
    // This mutates the shared_ptrs provided by the caller
    if (!i_b) {
        i_b     = b;
        i_bop   = bop;
        i_boplu = boplu;
    }

    const real_t bsplines_dist = channel::distance(*b, *i_b);
    if (bsplines_dist <= channel::bsplines_distinct_distance) {

        // Results match target numerics to within acceptable tolerance.
        retval.insert(time, s);

    } else {

        // Results do not match target numerics.
        // Must project onto target collocation points.

        // Convert all results in s to coefficients
        boplu->solve(quantity::count, s->data(), 1, b->n());

        // Obtain target collocation points
        Eigen::ArrayXr buf(i_b->n());
        for (int i = 0; i < i_b->n(); ++i) buf[i] = i_b->collocation_point(i);

        // Evaluate coefficients onto the target collocation points
        auto_ptr<storage_type> r(new storage_type(i_b->n(),
                    (storage_type::Index) storage_type::ColsAtCompileTime));
        for (std::size_t i = 0; i < quantity::count; ++i) {
            b->linear_combination(0, s->col(i).data(),
                                  buf.size(), buf.data(), r->col(i).data());
        }

        // Notice that quantity::t, being a constant, and quantity::y, being a
        // linear, should have been converted to the target collocation points
        // without more than epsilon-like floating point loss.

        retval.insert(time, r);

    }

    return retval;
}
