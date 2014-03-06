//--------------------------------------------------------------------------
//
// Copyright (C) 2014 Rhys Ulerich
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

/** \file
 * \copydoc helm.h
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/helm.h>

// C99 extern declaration for static inline function
extern
void
helm_reset(struct helm_state * const h);

// C99 extern declaration for static inline function
extern
void
helm_approach(struct helm_state * const h);

// C99 extern declaration for static inline function
extern
double
helm_steady(struct helm_state * const h,
            const double dt,
            const double r,
            const double u,
            const double v,
            const double y);
