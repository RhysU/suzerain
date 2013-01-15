/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2012, 2013 Rhys Ulerich
 * Copyright (C) 2012, 2013 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *--------------------------------------------------------------------------
 */

/** @file
 * @copydoc kahan.h
 */

#include <suzerain/kahan.h>

// --------------------------------------------------
// Generate routines for various floating point types
// --------------------------------------------------

#define NAME(infix,affix) suzerain_ ## infix ## affix
#define SCALAR  double
#include "kahan.def"

#define NAME(infix,affix) suzerain_ ## infix ## f ## affix
#define SCALAR  float
#include "kahan.def"

#define NAME(infix,affix) suzerain_ ## infix ## z ## affix
#define SCALAR  complex_double
#include "kahan.def"

#define NAME(infix,affix) suzerain_ ## infix ## c ## affix
#define SCALAR  complex_float
#include "kahan.def"
