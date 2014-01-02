/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2012-2014 Rhys Ulerich
 * Copyright (C) 2012-2014 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
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
 */

#include <suzerain/parcel.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_test.h>

#include <suzerain/common.h>

// Random number generation for all methods managed in main()
static gsl_rng * rng = NULL;

static
void
test_one()
{
    for (int j = 0; j < 50; ++j) {

        // Generate random size and alignment requirements
        const size_t sc = 1 + gsl_rng_uniform_int(rng, 10U);
        const size_t sd = 5 + gsl_rng_uniform_int(rng, 21U);
        const size_t si = 2 + gsl_rng_uniform_int(rng, 10U);
        const size_t sf = 7 + gsl_rng_uniform_int(rng, 31U);
        const size_t ad = 1 << (3 + gsl_rng_uniform_int(rng,  2U));
        const size_t af = 1 << (2 + gsl_rng_uniform_int(rng,  2U));

        // Malloc and parcel out memory according to those needs
        SUZERAIN_PARCEL(p, malloc, ((c, char,   sc, const,     0))
                                   ((d, double, sd, restrict, ad))
                                   ((i, int,    si,         ,  0))
                                   ((f, float,  sd, restrict, af)));

        // Ensure the alignments came back as requested
        const unsigned long int ac = __alignof__(char);
        const unsigned long int ai = __alignof__(int);
        gsl_test(!p                          , "j=%d,sc=%zu,sd=%zu,si=%zu,sf=%zu:        alloc p  OK", j, sc, sd, si, sf    );
        gsl_test(p_parcel < sc + sd + si + sf, "j=%d,sc=%zu,sd=%zu,si=%zu,sf=%zu:        p_parcel OK", j, sc, sd, si, sf    );
        gsl_test(!((uintptr_t)c & ~(ac - 1)) , "j=%d,sc=%zu,sd=%zu,si=%zu,sf=%zu,ac=%zu: align c  OK", j, sc, sd, si, sf, ac);
        gsl_test(!((uintptr_t)d & ~(ad - 1)) , "j=%d,sc=%zu,sd=%zu,si=%zu,sf=%zu,ad=%zu: align d  OK", j, sc, sd, si, sf, ad);
        gsl_test(!((uintptr_t)i & ~(ai - 1)) , "j=%d,sc=%zu,sd=%zu,si=%zu,sf=%zu,ai=%zu: align d  OK", j, sc, sd, si, sf, ai);
        gsl_test(!((uintptr_t)f & ~(af - 1)) , "j=%d,sc=%zu,sd=%zu,si=%zu,sf=%zu,af=%zu: align f  OK", j, sc, sd, si, sf, af);

        // Ensure the pointers came back with increasing addresses
        gsl_test((uintptr_t)p > (uintptr_t)c, "j=%d,sc=%zu,sd=%zu,si=%zu:        p <= c   OK", j, sc, sd, si, ac);
        gsl_test((uintptr_t)c > (uintptr_t)d, "j=%d,sc=%zu,sd=%zu,si=%zu:        c <= d   OK", j, sc, sd, si, ac);
        gsl_test((uintptr_t)d > (uintptr_t)i, "j=%d,sc=%zu,sd=%zu,si=%zu:        d <= i   OK", j, sc, sd, si, ac);
        gsl_test((uintptr_t)i > (uintptr_t)f, "j=%d,sc=%zu,sd=%zu,si=%zu:        i <= f   OK", j, sc, sd, si, ac);

        // Write data through the pointers to check for segfaults
        for (size_t k = 0; k < sc; ++k) c[k] = 'y';
        for (size_t k = 0; k < sd; ++k) d[k] = 123.456;
        for (size_t k = 0; k < si; ++k) i[k] = 7;
        for (size_t k = 0; k < sd; ++k) f[k] = 789;

        // Free all the memory in one swoop
        free(p);

        // Ensure a second usage in the same scope is without name conflicts
        SUZERAIN_PARCEL(qp, malloc, ((qc, char,   (sc-1), const restrict, 0))
                                    ((qd, double, (sd-2),               , 0))
                                    ((qi, int,    (si-1),               , 0)));
        (void) qi;
        free(qp);
    }
}

int
main(int argc, char **argv)
{
    SUZERAIN_UNUSED(argc);
    SUZERAIN_UNUSED(argv);

    rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, 0);

    test_one();

    gsl_rng_free(rng);

    exit(gsl_test_summary());
}
