//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#ifndef SUZERAIN_EXPRPARSE_HPP
#define SUZERAIN_EXPRPARSE_HPP

/** @file
 * Constant arithmetic expression evaluation via Boost Spirit.
 *
 * To reduce compilation times, only a \c float and \c double instantiation of
 * the parsing methods is exposed for \c std::string and <tt>const char *</tt>
 * data.
 */

#include <string>

namespace suzerain {

/**
 * Parse and evaluate a constant-valued arithmetic expression.  Parsing obeys a
 * case- and whitespace- insensitive, C-like grammar and operator precedence.
 * The logic supports many unary operators (<tt>+</tt> <tt>-</tt>), binary
 * operators (<tt>+</tt> <tt>-</tt> <tt>*</tt> <tt>/</tt> <tt>**</tt>), unary
 * functions (<tt>abs</tt> <tt>acos</tt> <tt>asin</tt> <tt>atan</tt>
 * <tt>ceil</tt> <tt>cos</tt> <tt>cosh</tt> <tt>exp</tt> <tt>floor</tt>
 * <tt>log</tt> <tt>log10</tt> <tt>sin</tt> <tt>sinh</tt> <tt>sqrt</tt>
 * <tt>tan</tt> <tt>tanh</tt>), binary functions (<tt>atan2</tt> <tt>max</tt>
 * <tt>min</tt> <tt>pow</tt>), and constants (<tt>digits</tt> <tt>digits10</tt>
 * <tt>e</tt> <tt>pi</tt> <tt>inf</tt> <tt>nan</tt>).
 *
 * @param[in]  s    Constant arithmetic expression to be parsed.
 * @param[out] v    Floating point value found from evaluating \c s.
 * @param[in]  name Optional name of what is being parsed for error reporting.
 *
 * @throws std::invalid_argument if the given input cannot be \em completely
 *         consumed by the underlying parser.  \c v is not modified.
 */
void exprparse(const char *s, float& v, const char *name = NULL);

/** @copydoc exprparse(const char *, float&, const char *) */
void exprparse(const char *s, double& v, const char *name = NULL);

/** @copydoc exprparse(const char *, float&, const char *) */
void exprparse(const std::string& s, float& v, const char *name = NULL);

/** @copydoc exprparse(const char *, float&, const char *) */
void exprparse(const std::string& s, double& v, const char *name = NULL);

/**
 * Parse and evaluate a constant-valued arithmetic expression.  Identical to
 * calling <tt>exprparse(const std::string&, FPT&, const char*)</tt> where the
 * result is returned by value rather than being a mutable argument.
 *
 * @tparam FPT  Type of floating point to evaluate and return.
 *              An exprparse() overload must exist.
 * @param  s    Constant arithmetic expression to be parsed.
 * @param  name Optional name of what is being parsed for error reporting.
 *
 * @return Floating point value found from evaluating \c s.
 */
template<typename FPT>
FPT exprparse(const std::string& s, const char *name)
{
    FPT t;
    exprparse(s, t, name);
    return t;
}

/**
 * Parse and evaluate a constant-valued arithmetic expression.  Identical to
 * calling <tt>exprparse(const char*, FPT&, const char*)</tt> where the result
 * is returned by value rather than being a mutable argument.
 *
 * @tparam FPT  Type of floating point to evaluate and return.
 *              An exprparse() overload must exist.
 * @param  s    Constant arithmetic expression to be parsed.
 * @param  name Optional name of what is being parsed for error reporting.
 *
 * @return Floating point value found from evaluating \c s.
 */
template<typename FPT>
FPT exprparse(const char *s, const char *name)
{
    FPT t;
    exprparse(s, t, name);
    return t;
}

/**
 * Parse and evaluate a constant-valued arithmetic expression.  Identical to
 * calling <tt>exprparse(const std::string&, FPT&)</tt> where the result is
 * returned by value rather than being a mutable argument.
 *
 * @tparam FPT  Type of floating point to evaluate and return.
 *              An exprparse() overload must exist.
 * @param  s    Constant arithmetic expression to be parsed.
 *
 * @return Floating point value found from evaluating \c s.
 */
template<typename FPT>
FPT exprparse(const std::string& s)
{
    FPT t;
    exprparse(s, t);
    return t;
}

/**
 * Parse and evaluate a constant-valued arithmetic expression.  Identical to
 * calling <tt>exprparse(const char*, FPT&)</tt> where the result is returned
 * by value rather than being a mutable argument.
 *
 * @tparam FPT  Type of floating point to evaluate and return.
 *              An exprparse() overload must exist.
 * @param  s    Constant arithmetic expression to be parsed.
 *
 * @return Floating point value found from evaluating \c s.
 */
template<typename FPT>
FPT exprparse(const char *s)
{
    FPT t;
    exprparse(s, t);
    return t;
}

/**
 * Parse a range with default minimum and maximum values.
 *
 * That is, parses "min:max", "min:[defaultmax]", or "[defaultmin]:max" into
 * valmin, \c valmax where \c absmin <= \c valmin <= \c valmax <= \c absmax is
 * enforced with the outer two inequalities being considered a validation
 * failure.
 *
 * @throws std::invalid_argument if the given input cannot be \em completely
 *         consumed by the underlying parser.
 */
void exprparse_range(const char *s,
                     float *valmin,
                     float *valmax,
                     const float defaultmin,
                     const float defaultmax,
                     const float absmin,
                     const float absmax,
                     const char *name = NULL);

/** @copydoc exprparse(const char *, float*, float*, const float, const float, const float, const float, const char *) */
void exprparse_range(const char *s,
                     double *valmin,
                     double *valmax,
                     const double defaultmin,
                     const double defaultmax,
                     const double absmin,
                     const double absmax,
                     const char *name = NULL);

/** @copydoc exprparse(const char *, float*, float*, const float, const float, const float, const float, const char *) */
void exprparse_range(const std::string& s,
                     float *valmin,
                     float *valmax,
                     const float defaultmin,
                     const float defaultmax,
                     const float absmin,
                     const float absmax,
                     const char *name = NULL);

/** @copydoc exprparse(const char *, float*, float*, const float, const float, const float, const float, const char *) */
void exprparse_range(const std::string& s,
                     double *valmin,
                     double *valmax,
                     const double defaultmin,
                     const double defaultmax,
                     const double absmin,
                     const double absmax,
                     const char *name = NULL);

/**
 * @copydoc exprparse(const char *, float*, float*, const float, const float, const float, const float, const char *)
 *
 * A templated overload easing selecting one of the concrete implementations.
 */
template <typename S, typename T>
void exprparse_range(S s,
                     T *valmin,
                     T *valmax,
                     const T defaultmin,
                     const T defaultmax,
                     const T absmin,
                     const T absmax,
                     const char *name = NULL)
{
    return exprparse_range(
            s, valmin, valmax, defaultmin, defaultmax, absmin, absmax, name);
}

} // end namespace suzerain

#endif // SUZERAIN_EXPRPARSE_HPP
