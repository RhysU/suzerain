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

#ifndef SUZERAIN_PARCEL_H
#define SUZERAIN_PARCEL_H

/** @file
 * Parcels out pieces of monolithically <tt>malloc</tt>-ed memory
 */

#include <inttypes.h>
#include <stddef.h>
#include <boost/preprocessor.hpp>

#ifdef __cplusplus
# error "SUZERAIN_PARCEL implementation cannot be used in C++ code"
#endif

// TODO Employ malloc guarantees of alignment of __alignof__(double)
// Doing so on the first decl will allow it to have no padding when
// minalign < __alignof__(double).

/**
 * Perform a single malloc that will be parceled out into multiple logical
 * regions in a type- and alignment-safe manner.  Macro arguments are
 * <pre>
 *   ptr      - pointer name (check for NULL and free appropriately)
 *   malloc   - malloc(3)-like function used for allocation
 *   decl_seq - a Boost.Preprocessor seq of "decl"s.
 * </pre>
 *
 * A "decl" is a tuple (varname, datatype, varcount, varqual, minalign)
 * where
 * <pre>
 *   varname  - pointer variable name (e.g. my_double_ptr)
 *   datatype - variable type (e.g. 'double')
 *   varcount - number of elements of type datatype to allocate (e.g. 1, 5)
 *   varqual  - qualifiers to add to pointer variable (e.g. const, restrict)
 *   minalign - minimum allowable alignment of data in bytes
 * </pre>
 *
 * One sample invocation is
 * <pre>
 *   SUZERAIN_PARCEL(p, malloc, ((c, char,   3, const,     0))
 *                              ((d, double, 5, restrict, 16))
 *                              ((i, int,    1,         ,  0)));
 * </pre>
 * where the extra parenthesis are crucial to make a sequence of decls.  This
 * invocation allocates a single buffer of sufficient size to hold a char[3], a
 * double[5], and a single int.  The buffer starts at p and is of size p_parcel
 * bytes.  c is a const pointer to a naturally aligned location for the
 * char[3], d points to a 16-byte aligned, restrict-ed location for the
 * double[5], and i points to a naturally aligned location for a single int.
 * The three pointers have increasing addresses.  The buffer spans [p, p + sz).
 * malloc(3) was used for the allocation but any other unary allocation
 * function may be supplied.  The entire work buffer should be deallocated
 * using free(p).
 *
 * The packing is suboptimal and the memory footprint slightly larger compared
 * to using a struct with a single flexible array member.  That is, either
 * <pre>
 *   struct { char c[3], int i, double d[]; } s;
 * </pre>
 * or
 * <pre>
 *   struct { double d[5], int i, char c[]; } s;
 * </pre>
 * is more efficient from a memory footprint perspective.  This is because we
 * must assume a worst-case allocation and cannot pack members.  The upside of
 * SUZERAIN_PARCEL is that multiple arrays may be allocated via one malloc(3)
 * call in a way that satisfies their alignment needs.
 *
 * Ideas lifted verbatim from http://stackoverflow.com/questions/227897/.
 * Attempting to find a cleaner solution via
 * http://stackoverflow.com/questions/8839566/.  These macros generate
 * code which the compiler should trivially optimize into the ground.
 */
#define SUZERAIN_PARCEL(ptr, malloc, decl_seq)                                     \
    BOOST_PP_SEQ_FOR_EACH(SUZERAIN_PARCEL_DATASIZE,,decl_seq)                      \
    BOOST_PP_SEQ_FOR_EACH(SUZERAIN_PARCEL_ALIGNMENT,,decl_seq)                     \
    BOOST_PP_SEQ_FOR_EACH(SUZERAIN_PARCEL_PAD,,decl_seq)                           \
    BOOST_PP_SEQ_FOR_EACH(SUZERAIN_PARCEL_MASK,,decl_seq)                          \
    const size_t BOOST_PP_CAT(ptr,_parcel)                                         \
        = BOOST_PP_SEQ_FOR_EACH(SUZERAIN_PARCEL_SIZESUM,,decl_seq);                \
    void * const ptr = malloc(BOOST_PP_CAT(ptr,_parcel))                           \
    BOOST_PP_TUPLE_ELEM(3, 2, BOOST_PP_SEQ_FOLD_LEFT(SUZERAIN_PARCEL_OP,           \
        ((char*) ptr                                                               \
         BOOST_PP_COMMA() 0                                                        \
         BOOST_PP_COMMA()), decl_seq))

#ifndef SUZERAIN_PARSED_BY_DOXYGEN

#define SUZERAIN_PARCEL_APPLY(macro, tuple_args) macro tuple_args

#define SUZERAIN_PARCEL_DATASIZE(r,data,elem)                                      \
    SUZERAIN_PARCEL_APPLY(SUZERAIN_PARCEL_DATASIZE_,elem)
#define SUZERAIN_PARCEL_DATASIZE_(varname, datatype, varcount, varqual, minalign)  \
    const size_t BOOST_PP_CAT(_PARCEL_datasize_,varname)                           \
        = varcount*sizeof(datatype);

#define SUZERAIN_PARCEL_ALIGNMENT(r,data,elem)                                     \
    SUZERAIN_PARCEL_APPLY(SUZERAIN_PARCEL_ALIGNMENT_,elem)
#define SUZERAIN_PARCEL_ALIGNMENT_(varname, datatype, varcount, varqual, minalign) \
    const size_t BOOST_PP_CAT(_PARCEL_alignment_,varname)                          \
        = __alignof__(datatype) > minalign ? __alignof__(datatype) : minalign;     \
    assert(!(BOOST_PP_CAT(_PARCEL_alignment_,varname)                              \
           & (BOOST_PP_CAT(_PARCEL_alignment_,varname) - 1)));

#define SUZERAIN_PARCEL_PAD(r,data,elem)                                           \
    SUZERAIN_PARCEL_APPLY(SUZERAIN_PARCEL_PAD_,elem)
#define SUZERAIN_PARCEL_PAD_(varname, datatype, varcount, varqual, minalign)       \
    const size_t BOOST_PP_CAT(_PARCEL_pad_,varname)                                \
        = BOOST_PP_CAT(_PARCEL_alignment_,varname) - 1;

#define SUZERAIN_PARCEL_MASK(r,data,elem)                                          \
    SUZERAIN_PARCEL_APPLY(SUZERAIN_PARCEL_MASK_,elem)
#define SUZERAIN_PARCEL_MASK_(varname, datatype, varcount, varqual, minalign)      \
    const uintptr_t BOOST_PP_CAT(_PARCEL_mask_,varname)                            \
        = ~(uintptr_t)(BOOST_PP_CAT(_PARCEL_alignment_,varname) - 1);

#define SUZERAIN_PARCEL_SIZESUM(r,data,elem)                                       \
    SUZERAIN_PARCEL_APPLY(SUZERAIN_PARCEL_SIZESUM_,elem)
#define SUZERAIN_PARCEL_SIZESUM_(varname, datatype, varcount, varqual, minalign)   \
    + (  BOOST_PP_CAT(_PARCEL_datasize_,varname)                                   \
       + BOOST_PP_CAT(_PARCEL_pad_,varname))

#define SUZERAIN_PARCEL_OP(s, state, elem) (                                       \
        BOOST_PP_TUPLE_ELEM(5, 0, elem) BOOST_PP_COMMA()                           \
        BOOST_PP_TUPLE_ELEM(5, 2, elem) BOOST_PP_COMMA()                           \
        BOOST_PP_TUPLE_ELEM(3, 2, state) SUZERAIN_PARCEL_OP_(                      \
            BOOST_PP_TUPLE_ELEM(5, 0, elem),                                       \
            BOOST_PP_TUPLE_ELEM(5, 1, elem),                                       \
            BOOST_PP_TUPLE_ELEM(5, 2, elem),                                       \
            BOOST_PP_TUPLE_ELEM(5, 3, elem),                                       \
            BOOST_PP_TUPLE_ELEM(3, 0, state),                                      \
            BOOST_PP_TUPLE_ELEM(3, 1, state)                                       \
        )                                                                          \
    )
#define SUZERAIN_PARCEL_OP_(varname, datatype, varcount, varqual, prev, prevcount) \
    ; datatype * varqual varname = (void *)                                        \
    ((((uintptr_t)(prev + prevcount)+BOOST_PP_CAT(_PARCEL_pad_,varname))           \
      & BOOST_PP_CAT(_PARCEL_mask_,varname)))

#endif // SUZERAIN_PARSED_BY_DOXYGEN

#endif // SUZERAIN_PARCEL_H
