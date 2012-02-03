//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// exprparse.hpp: constant arithmetic expression evaluation
// $Id$

#ifndef __SUZERAIN_EXPRPARSE_HPP
#define __SUZERAIN_EXPRPARSE_HPP

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
 * <tt>tan</tt> <tt>tanh</tt>), binary functions (<tt>pow</tt> <tt>atan2</tt>),
 * and constants (<tt>pi</tt> <tt>inf</tt> <tt>nan</tt>).
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
 * Parse an evaluate a constant-valued arithmetic expression.  Identical to
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
 * Parse an evaluate a constant-valued arithmetic expression.  Identical to
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
 * Parse an evaluate a constant-valued arithmetic expression.  Identical to
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
 * Parse an evaluate a constant-valued arithmetic expression.  Identical to
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

} // end namespace suzerain

#endif // __SUZERAIN_EXPRPARSE_HPP
