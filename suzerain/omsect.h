/*
 * Copyright (C) 2013 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef OMSECT_H

/** @file
 * Provides \ref omsect() adapted from https://github.com/RhysU/intersection.
 */

#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Given points \c a, \c b, \c x, and \c y under the assumption <code>a <=
 * b</code>, do intervals \c ab and \c xy overlap and if so what is the
 * order-matching interval?  Order-matching means an interval with endpoints
 * (p, q) such that p < q whenever x <= y and p > q whenever x > y.
 *
 * \param[in ] a Left endpoint on the first interval
 * \param[in ] b Right endpoint on the first interval
 * \param[in ] x Endpoint on the second interval
 * \param[in ] y Another on the second interval
 * \param[out] l if method returns \c true, set to be the lower
 *               endpoint of the order-matching intersection.
 * \param[out] h if method returns \c true, set to be the lower
 *               endpoint of the order-matching intersection.
 *
 * \return True if the segments intersect.  False otherwise.
 */
static inline
int omsect(double a, double b, double x, double y, double *l, double *u)
{
    assert(a<b);
    int D=b<x, E=b<y, C=a<y, B=a<x;
    int ret = (D&~E)|(~D&E)|(C&~(D&E))|(~C&(D|E))|(B&~(C&D&E))|(~B&(C|D|E));
    if (ret) {
        int F=x<y;
        *l  = a*(~B&C&~D&F);                       /* alower */
        *u  = a*(B&~C&~E&~F);                      /* aupper */
        *l += b*(B&D&~E&~F);                       /* blower */
        *u += b*(C&~D&E&F);                        /* bupper */
        *l += x*(B&(~D&((~E&~F)|(C&(~E|F)))));     /* xlower */
        *u += y*(C&(~E&((B&((~D)|(~F)))|(~D&F)))); /* yupper */
    }
    return ret;
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
