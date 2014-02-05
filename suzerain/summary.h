//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Suzerain is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------

#ifndef SUZERAIN_SUMMARY_H
#define SUZERAIN_SUMMARY_H

/** @file
 * <a
 * href="www.boost.org/doc/libs/release/libs/preprocessor/">Boost.Preprocessor</a>
 * definitions driving \ref suzerain::summary. Isolated to ease use, to simplify
 * re-use, and to facilitate debugging.
 *
 * An all-encompassing sequence like
 * <pre>
 * #define SUZERAIN_SUMMARY        \
 *     SUZERAIN_SUMMARY_GRID       \
 *     SUZERAIN_SUMMARY_SAMPLED    \
 *     SUZERAIN_SUMMARY_SAMPLED_Y  \
 *     SUZERAIN_SUMMARY_SAMPLED_YY
 * </pre>
 * is deliberately not defined.  It would contain more than BOOST_PP_LIMIT_SEQ
 * elements and be utterly useless for work with Boost.Preprocessor.  However,
 * building blocks called \ref SUZERAIN_SUMMARY_FOR_EACH and \ref
 * SUZERAIN_SUMMARY_ENUM_TRANSFORM are provided which circumvents this limit.
 */
#include <suzerain/samples.h>

/**
 * A Boost.Preprocessor sequence of tuples of grid-related details.
 * Each tuple contains a name and a description.
 */
#define SUZERAIN_SUMMARY_GRID                                                                         \
    ((t,            "Simulation time"                                                              )) \
    ((y,            "Wall-normal collocation point locations"                                      )) \
    ((bulk_weights, "Take dot product of these weights against any quantity to find the bulk value"))

/**
 * A Boost.Preprocessor sequence of tuples of directly-sampled quantity
 * components all prefixed by "bar_".  Each tuple contains a name and a
 * description. Automatically build from \ref SUZERAIN_SAMPLES in \ref
 * samples.h.
 *
 * \warning This will break gloriously if more than 255 elements arise.
 */
#define SUZERAIN_SUMMARY_SAMPLED      \
    SUZERAIN_SAMPLES_COMPONENTS(bar_)

/**
 * A Boost.Preprocessor sequence of tuples of directly-sampled quantity
 * components all prefixed by "bar_" and suffixed by "__y".
 * \copydetails SUZERAIN_SUMMARY_SAMPLED
 */
#define SUZERAIN_SUMMARY_SAMPLED_Y \
    BOOST_PP_SEQ_TRANSFORM(SUZERAIN_SUMMARY_SAMPLED_Y_TRANSFORM,,SUZERAIN_SUMMARY_SAMPLED)

#ifndef SUZERAIN_PARSED_BY_DOXYGEN

#define SUZERAIN_SUMMARY_SAMPLED_TRANSFORM_Y(r, data, tuple)           \
    (BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM(2,0,tuple),__y),                 \
     "Wall-normal first derivative of "BOOST_PP_TUPLE_ELEM(2,1,tuple))

#endif /* SUZERAIN_PARSED_BY_DOXYGEN */

/**
 * A Boost.Preprocessor sequence of tuples of directly-sampled quantity
 * components all prefixed by "bar_" and suffixed by "__yy".
 * \copydetails SUZERAIN_SUMMARY_SAMPLED
 */
#define SUZERAIN_SUMMARY_SAMPLED_YY \
    BOOST_PP_SEQ_TRANSFORM(SUZERAIN_SUMMARY_SAMPLED_YY_TRANSFORM,,SUZERAIN_SUMMARY_SAMPLED)

#ifndef SUZERAIN_PARSED_BY_DOXYGEN

#define SUZERAIN_SUMMARY_SAMPLED_TRANSFORM_YY(r, data, tuple)           \
    (BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM(2,0,tuple),__yy),                 \
     "Wall-normal second derivative of "BOOST_PP_TUPLE_ELEM(2,1,tuple))

#endif /* SUZERAIN_PARSED_BY_DOXYGEN */

/**
 * An Boost.Preprocessor-like iteration construct invoking <code>
 *     macro(data, name, description, offset)
 * </code> for each component present in \ref SUZERAIN_SUMMARY.
 *
 * On account of BOOST_PP_ADD capping out at BOOST_PP_LIMIT_MAG, beware \c
 * offset will be a constant, integer-only arithmetic expression wrapped in
 * parenthesis and not a single integer token.
 */
#define SUZERAIN_SUMMARY_FOR_EACH(macro, data)               \
    BOOST_PP_SEQ_FOR_EACH(                                   \
            SUZERAIN_SUMMARY_FOR_EACH_HELPER,                \
            (macro,                                          \
              0,                                             \
            data),                                           \
            SUZERAIN_SUMMARY_GRID)                           \
    BOOST_PP_SEQ_FOR_EACH(                                   \
            SUZERAIN_SUMMARY_FOR_EACH_HELPER,                \
            (macro,                                          \
              BOOST_PP_SEQ_SIZE(SUZERAIN_SUMMARY_GRID),      \
            data),                                           \
            SUZERAIN_SUMMARY_SAMPLED)                        \
    BOOST_PP_SEQ_FOR_EACH(                                   \
            SUZERAIN_SUMMARY_FOR_EACH_HELPER,                \
            (macro,                                          \
              BOOST_PP_SEQ_SIZE(SUZERAIN_SUMMARY_GRID)       \
            + BOOST_PP_SEQ_SIZE(SUZERAIN_SUMMARY_SAMPLED),   \
            data),                                           \
            SUZERAIN_SUMMARY_SAMPLED_Y)                      \
    BOOST_PP_SEQ_FOR_EACH(                                   \
            SUZERAIN_SUMMARY_FOR_EACH_HELPER,                \
            (macro,                                          \
              BOOST_PP_SEQ_SIZE(SUZERAIN_SUMMARY_GRID)       \
            + BOOST_PP_SEQ_SIZE(SUZERAIN_SUMMARY_SAMPLED),   \
            + BOOST_PP_SEQ_SIZE(SUZERAIN_SUMMARY_SAMPLED_Y), \
            data),                                           \
            SUZERAIN_SUMMARY_SAMPLED_YY)

#ifndef SUZERAIN_PARSED_BY_DOXYGEN

#define SUZERAIN_SUMMARY_FOR_EACH_HELPER(r, macro_base_data, i, name_description) \
    BOOST_PP_TUPLE_ELEM(3, 0, macro_base_data)(                                   \
        BOOST_PP_TUPLE_ELEM(3, 2, macro_base_data),                               \
        BOOST_PP_TUPLE_ELEM(2, 0, name_description),                              \
        BOOST_PP_TUPLE_ELEM(2, 1, name_description),                              \
        (BOOST_PP_TUPLE_ELEM(3, 1, macro_base_data)+i)                            \
    )                                                                             \

#endif /* SUZERAIN_PARSED_BY_DOXYGEN */

/**
 * Transform component according to <code>
 *     op(s, data, name_description_tuple)
 * </code> and enumerate the results. This macro effectively combines
 * BOOST_PP_SEQ_TRANSFORM and BOOST_PP_SEQ_ENUM in a way avoiding
 * BOOST_PP_LIMIT_SEQ.
 */
SUZERAIN_SUMMARY_ENUM_TRANSFORM(op, data)                                          \
    BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(op,data,SUZERAIN_SUMMARY_GRID)),      \
    BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(op,data,SUZERAIN_SUMMARY_SAMPLED)),   \
    BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(op,data,SUZERAIN_SUMMARY_SAMPLED_Y)), \
    BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(op,data,SUZERAIN_SUMMARY_SAMPLED_YY))

#endif /* SUZERAIN_SUMMARY_H */
