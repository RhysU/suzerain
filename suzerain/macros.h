/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Suzerain is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * macros.h: miscellaneous utility macros
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef PECOS_SUZERAIN_MACROS_H
#define PECOS_SUZERAIN_MACROS_H

/* Macro follows LOG4CXX_UNLIKELY */
#if !defined(SUZERAIN_UNLIKELY)
#if __GNUC__ >= 3
/**
Provides optimization hint to the compiler
to optimize for the expression being false.
@param expr boolean expression.
@returns value of expression.
*/
#define SUZERAIN_UNLIKELY(expr) __builtin_expect(expr, 0)
#else
/**
Provides optimization hint to the compiler
to optimize for the expression being false.
@param expr boolean expression.
@returns value of expression.
**/
#define SUZERAIN_UNLIKELY(expr) expr
#endif
#endif

#endif // PECOS_SUZERAIN_MACROS_H
