/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2011 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * exprparse.hpp: constant arithmetic expression evaluation
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

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
 * @warning If a Boost Spirit 2.1 or later is not available, the logic falls
 * back onto boost::lexical_cast and only floating point literals may be
 * parsed.
 *
 * @param[in]  s    Constant arithmetic expression to be parsed.
 * @param[out] v    Floating point value found from evaluating \c s.
 * @param[in]  name Optional name of what is being parsed for error reporting.
 *
 * @throws std::invalid_argument if the given input cannot be \i completely
 *         consumed by the underlying parser.  \c v is not modified.
 */
void exprparse(const char *s, float& v, const char *name = NULL);

/** @copydoc exprparse(const char *, float&) */
void exprparse(const char *s, double& v, const char *name = NULL);

/** @copydoc exprparse(const char *, float&) */
void exprparse(const char *s, long double& v, const char *name = NULL);

/** @copydoc exprparse(const char *, float&) */
void exprparse(const std::string& s, float& v, const char *name = NULL);

/** @copydoc exprparse(const char *, float&) */
void exprparse(const std::string& s, double& v, const char *name = NULL);

/** @copydoc exprparse(const char *, float&) */
void exprparse(const std::string& s, long double& v, const char *name = NULL);

} // end namespace suzerain

#endif // __SUZERAIN_EXPRPARSE_HPP
