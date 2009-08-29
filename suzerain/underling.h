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
 * underling.h: A parallel, three dimensional FFT library atop MPI
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef PECOS_SUZERAIN_UNDERLING_H
#define PECOS_SUZERAIN_UNDERLING_H

#include <stdlib.h>
#ifdef __cplusplus
#include <iosfwd>
#endif
#include <suzerain/error.h>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif
__BEGIN_DECLS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

typedef enum underling_state {
    UNDERLING_STATE_UNINITIALIZED  = 0,
    UNDERLING_STATE_PHYSICAL       = 1,
    UNDERLING_STATE_WAVE           = 2,
    UNDERLING_STATE_NOTTRANSFORMED = 4
} underling_state;

typedef struct underling_dimension underling_dimension;
struct underling_dimension {
    int size;
    int stride;
    int global_size;
    int global_start;
    double dealias_by;
    underling_state state;
    underling_dimension *next_r2c;
    underling_dimension *next_c2r;
};

typedef struct underling_stage {
    underling_dimension *dim;
} underling_stage;

typedef struct underling_workspace {
    int                  ndim;
    int                  nstage;
    underling_stage     *stage;
} underling_workspace;

typedef double         underling_real;

typedef underling_real underling_complex[2];

underling_workspace *
underling_workspace_alloc(const int ndim, const int nstage);

void
underling_workspace_free(underling_workspace * const w);

int
underling_prepare_physical_size(underling_workspace * const w,
                                const int *physical_size);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
__END_DECLS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifdef __cplusplus
inline
bool operator==(const underling_dimension& a, const underling_dimension& b)
{
    return (a.size == b.size)
        && (a.stride == b.stride)
        && (a.global_size == b.global_size)
        && (a.global_start == b.global_start)
        && (a.dealias_by == b.dealias_by)
        && (a.state == b.state)
        && (a.next_r2c == b.next_r2c)
        && (a.next_c2r == b.next_c2r);
}

inline
std::ostream& operator<<(std::ostream& os, const underling_dimension& ud)
{
    return os << "("
              << "size=" << ud.size << ","
              << "stride=" << ud.stride << ","
              << "global_size=" << ud.global_size << ","
              << "global_start=" << ud.global_start << ","
              << "dealias_by=" << ud.dealias_by << ","
              << "state=" << ud.state << ","
              << "next_r2c=" << ud.next_r2c << ","
              << "next_c2r=" << ud.next_c2r
              << ")";
}

#endif /* __cplusplus */

#endif // PECOS_SUZERAIN_UNDERLING_H
