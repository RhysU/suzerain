//--------------------------------------------------------------------------
//
// Copyright (C) 2011,2012 Rhys Ulerich
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at http://mozilla.org/MPL/2.0/.
//
// This is one of two files implementing a manufactured solution
// nondimensionalized for Suzerain based on a dimensional solution from
// @inproceedings{Ulerich2012,
//   title={A transient manufactured solution for the compressible
//          {N}avier–-{S}tokes equations with a power law viscosity},
//   author={Rhys Ulerich and Kemelli C. Estacio-Hiroms
//           and Nicholas Malaya and Robert D. Moser},
//   booktitle={{10th World Congress on Computational Mechanics}},
//   address={{S}\~{a}o {P}aulo, {B}razil},
//   month={{J}uly},
//   year={2012}
// }
//
//--------------------------------------------------------------------------

// Using this file will require #include-ing <cmath>, <limits>, <sstream>, and
// <string>.  Those headers are not #include-d here to allow #include-ing this
// file inside any namespace, as is the consistent use of ::std instead of
// merely std.  Running Doxygen on this file will require setting
// EXTRA_PACKAGES = "amsmath accents".

#ifndef NSCTPL_RHOLUT_FWD_HPP
#define NSCTPL_RHOLUT_FWD_HPP

/** @file
 * Provides a manufactured solution for a nondimensional compressible
 * Navier-Stokes formulation using a reference density \f$\rho_0\f$, length
 * \f$l_0\f$, velocity \f$u_0\f$, and temperature \f$T_0\f$.
 */

/**
 * Provides a manufactured solution for the nondimensional compressible
 * Navier-Stokes formulation using a reference density \f$\rho_0\f$, length
 * \f$l_0\f$, velocity \f$u_0\f$, and temperature \f$T_0\f$.  The exact
 * formulation used is
 * \f{align*}{
 *   \frac{\partial}{\partial{}t}\rho
 * &=
 *   - \nabla\cdot\rho{}\vec{u}
 *   + Q_{\rho}
 *   \\
 *   \frac{\partial{}}{\partial{}t}\rho{}\vec{u}
 * &=
 *   - \nabla\cdot(\vec{u}\otimes{}\rho{}\vec{u})
 *   - \frac{1}{\mbox{Ma}^{2}} \nabla{} p
 *   + \frac{1}{\mbox{Re}} \nabla\cdot{} \accentset{\leftrightarrow}{\tau}
 *   + \vec{Q}_{\rho{}u}
 *   \\
 *   \frac{\partial}{\partial{}t} \rho{}e
 * &=
 *   - \nabla\cdot{}\rho{}e\vec{u}
 *   - \nabla\cdot{} p \vec{u}
 *   - \nabla\cdot{} \vec{q}
 *   + \frac{\mbox{Ma}^2}{\mbox{Re}}
 *     \nabla\cdot{} \accentset{\leftrightarrow}{\tau} \vec{u}
 *   + Q_{\rho{}e}
 * \f}
 * aided by the auxiliary relations
 * \f{align*}{
 *   p &=   \left(\gamma-1\right)\left(\rho{}e
 *           - \mbox{Ma}^{2} \rho\frac{\vec{u}\cdot{}\vec{u}}{2} \right)
 *   &
 *   T &= \gamma\frac{p}{\rho}
 *   \\
 *   \mu &= T^{\beta}
 *   &
 *   \lambda &= \left(\alpha-\frac{2}{3}\right) \mu
 *   \\
 *   \accentset{\leftrightarrow}{\tau} &=
 *        \mu \left( \nabla{}\vec{u} + {\nabla{}\vec{u}}^{\mathsf{T}} \right)
 *      + \lambda \left( \nabla\cdot{}\vec{u} \right) I
 *   &
 *   \vec{q} &= - \frac{1}{\mbox{Re}\mbox{Pr}\left(\gamma-1\right)} \mu \nabla{} T
 * \f}
 * where the nondimensional quantities
 * \f{align*}{
 *   \mbox{Re} &= \frac{\rho_{0}u_{0}l_{0}}{\mu_{0}}
 *   &
 *   \mbox{Ma} &= \frac{u_{0}}{a_{0}}
 *   &
 *   \mbox{Pr} &= \frac{\mu_{0}C_{p}}{\kappa_{0}}
 *   &
 *   \gamma &= \frac{C_{p}}{C_{v}}
 * \f}
 * are the constant Reynolds number, Mach number, and Prandtl number, and ratio
 * of specific heats, respectively.
 */
namespace nsctpl_rholut {

/**
 * Class providing operations for evaluating an analytical solution at a
 * given location and time.  The solution is of the form
 * \verbatim
         a_0                                         *cos(f_0 *t + g_0 )
       + a_x  * cos(b_x *x + c_x )                   *cos(f_x *t + g_x )
       + a_xy * cos(b_xy*x + c_xy)*cos(d_xy*y + e_xy)*cos(f_xy*t + g_xy)
       + a_xz * cos(b_xz*x + c_xz)*cos(d_xz*z + e_xz)*cos(f_xz*t + g_xz)
       + a_y  * cos(b_y *y + c_y )                   *cos(f_y *t + g_y )
       + a_yz * cos(b_yz*y + c_yz)*cos(d_yz*z + e_yz)*cos(f_yz*t + g_yz)
       + a_z  * cos(b_z *z + c_z )                   *cos(f_z *t + g_z )
   \endverbatim
 * The member method names are non-traditional but permitted by the language
 * standards and well-aligned with mathematical notation.  Location and
 * time arguments are templated to allow extended precision intermediate
 * computations (followed by truncation) when possible.
 *
 * The foreach_parameter() member allows invoking an operation on each solution
 * parameter to aid parameter registration, initialization, or output.
 * Parameters may have an infix name added for use with foreach_parameter().
 * For example, setting the name 'phi' will cause foreach_parameter() to report
 * names like 'a_phix'.
 */
template <typename Scalar>
class primitive {

public:

    typedef Scalar scalar_type;  //!< Scalar type employed

     //! X macro per http://drdobbs.com/blogs/cpp/228700289 which will
     //! apply a macro on all solution parameter prefix/suffix pairs.
#define NSCTPL_RHOLUT_FOR_EACH_PRIMITIVE_PARAMETER(apply)                      \
    apply(a_,0)                                                                \
    apply(a_,x) apply(a_,xy) apply(a_,xz) apply(a_,y) apply(a_,yz) apply(a_,z) \
    apply(b_,x) apply(b_,xy) apply(b_,xz) apply(b_,y) apply(b_,yz) apply(b_,z) \
    apply(c_,x) apply(c_,xy) apply(c_,xz) apply(c_,y) apply(c_,yz) apply(c_,z) \
                apply(d_,xy) apply(d_,xz)             apply(d_,yz)             \
                apply(e_,xy) apply(e_,xz)             apply(e_,yz)             \
    apply(f_,0)                                                                \
    apply(f_,x) apply(f_,xy) apply(f_,xz) apply(f_,y) apply(f_,yz) apply(f_,z) \
    apply(g_,0)                                                                \
    apply(g_,x) apply(g_,xy) apply(g_,xz) apply(g_,y) apply(g_,yz) apply(g_,z)

    // Declare all solution parameters as members, e.g. 'a_xy'
#define NSCTPL_RHOLUT_APPLY(pre,suf) Scalar pre##suf;
    NSCTPL_RHOLUT_FOR_EACH_PRIMITIVE_PARAMETER(NSCTPL_RHOLUT_APPLY)
#undef NSCTPL_RHOLUT_APPLY

    //! The name used infix in foreach_parameter' names.  For example,
    //! name = 'phi' implies parameters names like 'a_phix'.
    ::std::string name;

    const Scalar* Lx;           //!< Domain extent in x direction
    const Scalar* Ly;           //!< Domain extent in y direction
    const Scalar* Lz;           //!< Domain extent in z direction

    //! Construct an instance using \c name in the reported parameter names.
    //! The domain sizes are referenced from some external Lx, Ly, and Lz.
    //! All parameters set to zero at construction time.
    explicit primitive(const ::std::string &name,
                       const Scalar& Lx,
                       const Scalar& Ly,
                       const Scalar& Lz)
#define NSCTPL_RHOLUT_APPLY(pre,suf) pre##suf(0),
        : NSCTPL_RHOLUT_FOR_EACH_PRIMITIVE_PARAMETER(NSCTPL_RHOLUT_APPLY)  // ,
#undef NSCTPL_RHOLUT_APPLY
          name(name), Lx(&Lx), Ly(&Ly), Lz(&Lz)
    {}

    //! Construct an instance using \c name in the reported parameter names.
    //! Domain sizes Lx, Ly, and Lz \b must be set prior to invoking any
    //! member methods.  All parameters set to zero at construction time.
    explicit primitive(const ::std::string &name = "")
#define NSCTPL_RHOLUT_APPLY(pre,suf) pre##suf(0),
        : NSCTPL_RHOLUT_FOR_EACH_PRIMITIVE_PARAMETER(NSCTPL_RHOLUT_APPLY)  // ,
#undef NSCTPL_RHOLUT_APPLY
          name(name), Lx(NULL), Ly(NULL), Lz(NULL)
    {}

#define NSCTPL_RHOLUT_APPLY_STRINGIFY(s) #s
#define NSCTPL_RHOLUT_APPLY(pre,suf)              \
        os.clear(); os.str("");                   \
        os << NSCTPL_RHOLUT_APPLY_STRINGIFY(pre)  \
           << this->name                          \
           << NSCTPL_RHOLUT_APPLY_STRINGIFY(suf); \
        f(os.str(), this->pre##suf);

    //! Invoke the binary function f on each parameter name and its value.
    template <typename BinaryFunction>
    void foreach_parameter(BinaryFunction f) const {
        ::std::ostringstream os;
        NSCTPL_RHOLUT_FOR_EACH_PRIMITIVE_PARAMETER(NSCTPL_RHOLUT_APPLY)
    }

    //! Invoke the binary function f on each parameter name and its value.
    template <typename BinaryFunction>
    void foreach_parameter(BinaryFunction f) {
        ::std::ostringstream os;
        NSCTPL_RHOLUT_FOR_EACH_PRIMITIVE_PARAMETER(NSCTPL_RHOLUT_APPLY)
    }

#undef NSCTPL_RHOLUT_APPLY
#undef NSCTPL_RHOLUT_APPLY_STRINGIFY

    //! Evaluate the solution
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar operator()(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to time
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _t(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to \c x
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _x(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's second derivative with respect to \c x
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _xx(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to \c x and \c y
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _xy(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to \c x and \c z
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _xz(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to \c y
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _y(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's second derivative with respect to \c y
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _yy(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to \c y and \c z
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _yz(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to \c z
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _z(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's second derivative with respect to \c z
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _zz(T1 x, T2 y, T3 z, T4 t) const;

private:

    static const Scalar twopi;  //!< \f$2\pi\f$ to \c Scalar precision

}; // end class

// It is handy to template manufactured_solution (just below) on the
// primitive functional forms chosen for rho, u, v, w, and T.  However
// SWIG, even version 2.0.3, has trouble with template template parameters:
// https://sourceforge.net/tracker/?func=detail&atid=101645&aid=1861407&group_id=1645
// Therefore only employ the template template parameter when SWIG is not used.

/**
 * Template that, given a primitive function for \c rho, \c u, \c v, \c w, and
 * \c T along with a floating point type, evaluates a manufactured solution and
 * the required forcing for the transient, compressible Navier--Stokes
 * equations with a power law viscosity.  Location and time arguments are
 * templated to allow extended precision intermediate computations (followed by
 * truncation) when possible.
 */
#ifndef SWIG
template <typename Scalar,
          int IndexBase = 0,
          template <typename> class Primitive = primitive>
#else
template <typename Scalar, int IndexBase = 0>
#endif /* SWIG */
class manufactured_solution {

public:

    typedef Scalar scalar_type;                //!< Scalar type employed
    static const int index_base;               //!< Indexing base for gradients
#ifndef SWIG
    typedef Primitive<Scalar> primitive_type;  //!< Primitive solution type
#else
    typedef primitive<Scalar> primitive_type;  //!< Primitive solution type
#endif /* SWIG */

    //! Scenario parameters
    //!@{
    Scalar alpha;            //!< Ratio of bulk to dynamic viscosity
    Scalar beta;             //!< Temperature power law exponent
    Scalar gamma;            //!< Constant ratio of specific heats
    Scalar Ma;               //!< Mach number
    Scalar Pr;               //!< Prandtl number
    Scalar Re;               //!< Reynolds number
    //!@}

    //! Analytic solutions (which contain additional parameters)
    //!@{
    primitive_type rho;      //!< Analytic solution for rho
    primitive_type u;        //!< Analytic solution for u
    primitive_type v;        //!< Analytic solution for v
    primitive_type w;        //!< Analytic solution for w
    primitive_type T;        //!< Analytic solution for T
    //!@}

    //! Domain extents
    //!@{
    Scalar Lx;               //!< Domain extent in x direction
    Scalar Ly;               //!< Domain extent in y direction
    Scalar Lz;               //!< Domain extent in z direction
    //!@}

    /*! Default constructor */
    manufactured_solution()
        : alpha(0), beta(0), gamma(0), Ma(0), Pr(0), Re(0),
          rho("rho", Lx, Ly, Lz),
          u  ("u"  , Lx, Ly, Lz),
          v  ("v"  , Lx, Ly, Lz),
          w  ("w"  , Lx, Ly, Lz),
          T  ("T"  , Lx, Ly, Lz),
          Lx(1), Ly(1), Lz(1)
    {
    }

    /*! Invoke the binary function f on each parameter name
     *  and its constant value. */
    template <typename BinaryFunction>
    void foreach_parameter(BinaryFunction f) const {
        f(::std::string("alpha"), alpha);
        f(::std::string("beta"),  beta );
        f(::std::string("gamma"), gamma);
        f(::std::string("Ma"),    Ma);
        f(::std::string("Pr"),    Pr);
        f(::std::string("Re"),    Re);
        rho.foreach_parameter(f);
        u.foreach_parameter(f);
        v.foreach_parameter(f);
        w.foreach_parameter(f);
        T.foreach_parameter(f);
        f(::std::string("Lx"), Lx);
        f(::std::string("Ly"), Ly);
        f(::std::string("Lz"), Lz);
    }

    /*! Invoke the binary function f on each parameter name
     * and its mutable value. */
    template <typename BinaryFunction>
    void foreach_parameter(BinaryFunction f) {
        f(::std::string("alpha"), alpha);
        f(::std::string("beta"),  beta );
        f(::std::string("gamma"), gamma);
        f(::std::string("Ma"),    Ma);
        f(::std::string("Pr"),    Pr);
        f(::std::string("Re"),    Re);
        rho.foreach_parameter(f);
        u.foreach_parameter(f);
        v.foreach_parameter(f);
        w.foreach_parameter(f);
        T.foreach_parameter(f);
        f(::std::string("Lx"), Lx);
        f(::std::string("Ly"), Ly);
        f(::std::string("Lz"), Lz);
    }

    // Analytically determined quantities
    // Note that Primitive members can be used directly.
    // For example, T(x,y,z,t) or T._xx(x,y,z,t)

    /*! Compute component \c index of \f$\nabla\rho\f$. */
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_rho (T1 x, T2 y, T3 z, T4 t, int index) const;

    /*! Compute component \c index of \f$\nabla{}u\f$. */
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_u   (T1 x, T2 y, T3 z, T4 t, int index) const;

    /*! Compute component \c index of \f$\nabla{}v\f$. */
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_v   (T1 x, T2 y, T3 z, T4 t, int index) const;

    /*! Compute component \c index of \f$\nabla{}w\f$. */
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_w   (T1 x, T2 y, T3 z, T4 t, int index) const;

    /*! Compute component \c index of \f$\nabla{}T\f$. */
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_T   (T1 x, T2 y, T3 z, T4 t, int index) const;

    // Quantities built from the analytical solutions
    // TODO Build up q_u, q_v, q_w, q_e, q_T, q_p

    /*! Compute \f$e\f$. */
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar e(T1 x, T2 y, T3 z, T4 t) const;

    /*! Compute \f$p\f$. */
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar p(T1 x, T2 y, T3 z, T4 t) const;

    /*! Compute \f$\mu\f$. */
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar mu(T1 x, T2 y, T3 z, T4 t) const;

    /*! Compute \f$\rho{}u\f$. */
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar rhou(T1 x, T2 y, T3 z, T4 t) const;

    /*! Compute \f$\rho{}v\f$. */
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar rhov(T1 x, T2 y, T3 z, T4 t) const;

    /*! Compute \f$\rho{}w\f$. */
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar rhow(T1 x, T2 y, T3 z, T4 t) const;

    /*! Compute \f$\rho{}e\f$. */
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar rhoe(T1 x, T2 y, T3 z, T4 t) const;

    /*! Compute component \c index of \f$\nabla{}e\f$. */
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_e(T1 x, T2 y, T3 z, T4 t, int index) const;

    /*! Compute component \c index of \f$\nabla{}p\f$. */
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_p(T1 x, T2 y, T3 z, T4 t, int index) const;

    /*! Compute component \c index of \f$\nabla\mu\f$. */
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_mu(T1 x, T2 y, T3 z, T4 t, int index) const;

    /*! Compute forcing \f$Q_{\rho}\f$ */
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar Q_rho(T1 x, T2 y, T3 z, T4 t) const;

    /*! Compute forcing \f$Q_{\rho{}u}\f$ */
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar Q_rhou(T1 x, T2 y, T3 z, T4 t) const;

    /*! Compute forcing \f$Q_{\rho{}v}\f$ */
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar Q_rhov(T1 x, T2 y, T3 z, T4 t) const;

    /*! Compute forcing \f$Q_{\rho{}w}\f$ */
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar Q_rhow(T1 x, T2 y, T3 z, T4 t) const;

    /*! Compute forcing \f$Q_{\rho{}e}\f$ */
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar Q_rhoe(T1 x, T2 y, T3 z, T4 t) const;

    /**
     * A single method for simultaneously computing forcing for all five
     * equations in a conservative formulation at a given location.  That is,
     * computing Q_rho(), Q_rhou(), Q_rhov(), Q_rhow(), and Q_rhoe().  Faster
     * than using those five methods separately because it eliminates redundant
     * computations.
     *
     * @param[in] x X coordinate at which to evaluate forcing
     * @param[in] y Y coordinate at which to evaluate forcing
     * @param[in] z Z coordinate at which to evaluate forcing
     * @param[in] t Time coordinate at which to evaluate forcing
     * @param[out] Q_rho  Equivalent to <tt>Q_rho(x,y,z,t)</tt>
     * @param[out] Q_rhou Equivalent to <tt>Q_rhou(x,y,z,t)</tt>
     * @param[out] Q_rhov Equivalent to <tt>Q_rhov(x,y,z,t)</tt>
     * @param[out] Q_rhow Equivalent to <tt>Q_rhow(x,y,z,t)</tt>
     * @param[out] Q_rhoe Equivalent to <tt>Q_rhoe(x,y,z,t)</tt>
     */
    template <typename T1, typename T2, typename T3, typename T4,
              typename Result>
    void Q_conservative(T1 x, T2 y, T3 z, T4 t,
                        Result& Q_rho, Result& Q_rhou,
                        Result& Q_rhov, Result& Q_rhow, Result& Q_rhoe) const;

}; // end class


// Utilities take generic templated types to allow reuse by
// any variants satisfying the necessary APIs.

/** Zero all of an instance's parameters using foreach_parameter */
template <class T> void zero(T& t);

/** Set recommended isothermal channel problem parameters per write up */
template <class T> void isothermal_channel(T& t);

/** Set recommended isothermal flat plate problem parameters per write up */
template <class T> void isothermal_flat_plate(T& t);

} // end namespace nsctpl_rholut

#endif /* NSCTPL_RHOLUT_FWD_HPP */
