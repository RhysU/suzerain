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
 * exprgrammar.hpp: constant arithmetic expression evaluation via Boost Spirit
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_EXPRGRAMMAR_HPP
#define __SUZERAIN_EXPRGRAMMAR_HPP

//////////////////////////////////////////////////////////////////////
// Is Spirit new enough?  Compilation errors are astounding otherwise.
//////////////////////////////////////////////////////////////////////
#include <boost/version.hpp>
#if BOOST_VERSION < 104100
# error "Boost 1.41 or newer is required for exprgrammar.hpp functionality"
#endif

//////////////////////////////////////////////////////////////////
// Implementation uses X macros to reduce boilerplate.
// See http://drdobbs.com/blogs/cpp/228700289 for X macro details.
//////////////////////////////////////////////////////////////////

#include <cmath>
#include <limits>

#include <boost/math/constants/constants.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/qi.hpp>

namespace suzerain {

/**
 * Provides a Boost Spirit.Qi-based constant arithmetic expression evaluation
 * utilities.  More details on this parsing implementation can be found at
 * http://agentzlerich.blogspot.com/2011/06/using-boost-spirit-21-to-evaluate.html
 *
 * This code is expensive to compile and a little bit tricky when it comes to
 * handling errors.  Most users will find the pre-compiled function overloads
 * suzerain::exprparse() declared in exprparse.hpp to be sufficient and easier
 * to use.
 */
namespace exprgrammar {

#define SUZERAIN_EXPRGRAMMAR_FORALL_UNARY_FUNCTIONS(apply) \
    apply(abs)                                             \
    apply(acos)                                            \
    apply(asin)                                            \
    apply(atan)                                            \
    apply(ceil)                                            \
    apply(cos)                                             \
    apply(cosh)                                            \
    apply(exp)                                             \
    apply(floor)                                           \
    apply(log)                                             \
    apply(log10)                                           \
    apply(sin)                                             \
    apply(sinh)                                            \
    apply(sqrt)                                            \
    apply(tan)                                             \
    apply(tanh)

#define SUZERAIN_EXPRGRAMMAR_FORALL_BINARY_FUNCTIONS(apply) \
    apply(pow)                                              \
    apply(atan2)


namespace detail {

#define DECLARE_UNARY_FUNCTOR(name)                          \
struct name ## _                                             \
{                                                            \
    template <typename X> struct result { typedef X type; }; \
                                                             \
    template <typename X> X operator()(X x) const            \
    {                                                        \
        using namespace std;                                 \
        return name(x);                                      \
    }                                                        \
};
SUZERAIN_EXPRGRAMMAR_FORALL_UNARY_FUNCTIONS(DECLARE_UNARY_FUNCTOR)
#undef DECLARE_UNARY_FUNCTOR

#define DECLARE_BINARY_FUNCTOR(name)                                     \
struct name ## _                                                         \
{                                                                        \
    template <typename X, typename Y> struct result { typedef X type; }; \
                                                                         \
    template <typename X, typename Y> X operator()(X x, Y y) const       \
    {                                                                    \
        using namespace std;                                             \
        return name(x, y);                                               \
    }                                                                    \
};
SUZERAIN_EXPRGRAMMAR_FORALL_BINARY_FUNCTIONS(DECLARE_BINARY_FUNCTOR)
#undef DECLARE_BINARY_FUNCTOR

} // end namespace detail

/** A Boost Spirit.Qi-based constant arithmetic expression grammar. */
template <typename FPT, typename Iterator>
struct grammar
    : boost::spirit::qi::grammar<
            Iterator, FPT(), boost::spirit::ascii::space_type
        >
{
    boost::spirit::qi::rule<
            Iterator, FPT(), boost::spirit::ascii::space_type
        > expression, term, factor, primary;

    grammar() : grammar::base_type(expression)
    {
        namespace constants = boost::math::constants;
#define DECLARE_PHOENIX_FUNCTION(name)                    \
        boost::phoenix::function<detail::name ## _> name;
SUZERAIN_EXPRGRAMMAR_FORALL_UNARY_FUNCTIONS(DECLARE_PHOENIX_FUNCTION)
SUZERAIN_EXPRGRAMMAR_FORALL_BINARY_FUNCTIONS(DECLARE_PHOENIX_FUNCTION)
#undef DECLARE_PHOENIX_FUNCTION

        using boost::spirit::qi::real_parser;
        using boost::spirit::qi::real_policies;
        real_parser<FPT,real_policies<FPT> > real;

        using boost::spirit::qi::_1;
        using boost::spirit::qi::_2;
        using boost::spirit::qi::no_case;
        using boost::spirit::qi::_val;

        expression =
            term                   [_val =  _1]
            >> *(  ('+' >> term    [_val += _1])
                |  ('-' >> term    [_val -= _1])
                )
            ;

        term =
            factor                 [_val =  _1]
            >> *(  ('*' >> factor  [_val *= _1])
                |  ('/' >> factor  [_val /= _1])
                )
            ;

        factor =
            primary                [_val =  _1          ]
            >> *(  ("**" >> factor [_val = pow(_val, _1)])
                )
            ;

#define STRINGIFY(s) #s
#define UNARY_FUNCTION_RULE(name)                                             \
    | no_case[STRINGIFY(name)] >> '(' >> expression [_val =  name(_1)] >> ')'
#define BINARY_FUNCTION_RULE(name)                                            \
    | ( no_case[STRINGIFY(name)]>> '('                                        \
        >> expression >> ',' >> expression >> ')' )[_val =  name(_1,_2)]

        primary =
            real                   [_val =  _1]
            |   '(' >> expression  [_val =  _1] >> ')'
            |   ('-' >> primary    [_val = -_1])
            |   ('+' >> primary    [_val =  _1])
            |   no_case["pi"]      [_val =  constants::pi<FPT>()]
                SUZERAIN_EXPRGRAMMAR_FORALL_UNARY_FUNCTIONS(UNARY_FUNCTION_RULE)
                SUZERAIN_EXPRGRAMMAR_FORALL_BINARY_FUNCTIONS(BINARY_FUNCTION_RULE)
            ;

#undef UNARY_RULE
#undef BINARY_RULE
#undef STRINGIFY

    }
};

#undef SUZERAIN_EXPRGRAMMAR_FORALL_UNARY_FUNCTIONS
#undef SUZERAIN_EXPRGRAMMAR_FORALL_BINARY_FUNCTIONS

/**
 * Parse a constant arithmetic expression in [<tt>iter</tt>,<tt>end</tt>)
 * according to a grammar instance \c g using
 * <tt>boost::spirit::qi::phrase_parse</tt>.
 */
template <typename FPT, typename Iterator>
bool parse(Iterator &iter,
           Iterator end,
           const grammar<FPT,Iterator> &g,
           FPT &result)
{
    return boost::spirit::qi::phrase_parse(
                iter, end, g, boost::spirit::ascii::space, result);
}

/**
 * Parse a constant arithmetic expression in [<tt>iter</tt>,<tt>end</tt>) using
 * <tt>boost::spirit::qi::phrase_parse</tt>.
 */
template <typename FPT, typename Iterator>
bool parse(Iterator &iter,
           Iterator end,
           FPT &result)
{
    // Idiom: construct on first use using local static
    static const grammar<FPT,Iterator> g;

    return parse(iter, end, g, result);
}

} // end namespace exprgrammar

} // end namespace suzerain

#endif // __SUZERAIN_EXPRGRAMMAR_HPP
