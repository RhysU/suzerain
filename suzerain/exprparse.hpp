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
 * @internal To mitigate long compilation times associated with Boost Spirit.Qi
 * parsers, only a small number of explicit template instantiations are
 * provided.  These are compiled once and may be linked quickly.
 *
 * @tparam FPT        Floating point precision to return.  Must be one of
 *                    <tt>float</tt>, <tt>double</tt>, or <tt>long double</tt>.
 * @tparam StringType String-like type to parse.  Must be either
 *                    <tt>char *</tt>, or <tt>std::string&</tt>.
 * @param s           Constant arithmetic expression to be parsed.
 *
 * @return The floating point result found from parsing and evaluating the
 *         input \c s.
 *
 * @throws std::invalid_argument if the given input cannot be \i completely
 *         consumed by the underlying parser.
 */
template<typename FPT, typename StringType> FPT exprparse(const StringType s);

template<> float       exprparse(const char *s);
template<> double      exprparse(const char *s);
template<> long double exprparse(const char *s);
template<> float       exprparse(const std::string& s);
template<> double      exprparse(const std::string& s);
template<> long double exprparse(const std::string& s);

} // end namespace suzerain

#endif // __SUZERAIN_EXPRPARSE_HPP
