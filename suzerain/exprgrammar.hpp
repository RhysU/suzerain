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

#include <cmath>
#include <limits>
#include <iterator>

#include <boost/spirit/version.hpp>
#if !defined(SPIRIT_VERSION) || SPIRIT_VERSION < 0x2010
#error "At least Spirit version 2.1 required"
#endif
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

namespace { // anonymous

struct lazy_pow_
{
    template <typename X, typename Y>
    struct result { typedef X type; };

    template <typename X, typename Y>
    X operator()(X x, Y y) const
    {
        return std::pow(x, y);
    }
};

struct lazy_ufunc_
{
    template <typename F, typename A1>
    struct result { typedef A1 type; };

    template <typename F, typename A1>
    A1 operator()(F f, A1 a1) const
    {
        return f(a1);
    }
};

struct lazy_bfunc_
{
    template <typename F, typename A1, typename A2>
    struct result { typedef A1 type; };

    template <typename F, typename A1, typename A2>
    A1 operator()(F f, A1 a1, A2 a2) const
    {
        return f(a1, a2);
    }
};

} // end namespace anonymous

/** A Boost Spirit.Qi-based constant arithmetic expression grammar. */
template <typename FPT, typename Iterator>
struct grammar
    : boost::spirit::qi::grammar<
            Iterator, FPT(), boost::spirit::ascii::space_type
        >
{

    // symbol table for constants like "pi"
    struct constant_
        : boost::spirit::qi::symbols<
                typename std::iterator_traits<Iterator>::value_type,
                FPT
            >
    {
        constant_()
        {
            this->add
                ("pi", boost::math::constants::pi<FPT>()  )
            ;
        }
    } constant;

    // symbol table for unary functions like "abs"
    struct ufunc_
        : boost::spirit::qi::symbols<
                typename std::iterator_traits<Iterator>::value_type,
                FPT (*)(FPT)
            >
    {
        ufunc_()
        {
            this->add
                ("abs"   , (FPT (*)(FPT)) std::abs  )
                ("acos"  , (FPT (*)(FPT)) std::acos )
                ("asin"  , (FPT (*)(FPT)) std::asin )
                ("atan"  , (FPT (*)(FPT)) std::atan )
                ("ceil"  , (FPT (*)(FPT)) std::ceil )
                ("cos"   , (FPT (*)(FPT)) std::cos  )
                ("cosh"  , (FPT (*)(FPT)) std::cosh )
                ("exp"   , (FPT (*)(FPT)) std::exp  )
                ("floor" , (FPT (*)(FPT)) std::floor)
                ("log"   , (FPT (*)(FPT)) std::log  )
                ("log10" , (FPT (*)(FPT)) std::log10)
                ("sin"   , (FPT (*)(FPT)) std::sin  )
                ("sinh"  , (FPT (*)(FPT)) std::sinh )
                ("sqrt"  , (FPT (*)(FPT)) std::sqrt )
                ("tan"   , (FPT (*)(FPT)) std::tan  )
                ("tanh"  , (FPT (*)(FPT)) std::tanh )
            ;
        }
    } ufunc;

    // symbol table for binary functions like "pow"
    struct bfunc_
        : boost::spirit::qi::symbols<
                typename std::iterator_traits<Iterator>::value_type,
                FPT (*)(FPT, FPT)
            >
    {
        bfunc_()
        {
            this->add
                ("pow"  , (FPT (*)(FPT, FPT)) std::pow  )
                ("atan2", (FPT (*)(FPT, FPT)) std::atan2)
            ;
        }
    } bfunc;

    boost::spirit::qi::rule<
            Iterator, FPT(), boost::spirit::ascii::space_type
        > expression, term, factor, primary;

    grammar() : grammar::base_type(expression)
    {
        using boost::spirit::qi::real_parser;
        using boost::spirit::qi::real_policies;
        real_parser<FPT,real_policies<FPT> > real;

        using boost::spirit::qi::_1;
        using boost::spirit::qi::_2;
        using boost::spirit::qi::_3;
        using boost::spirit::qi::no_case;
        using boost::spirit::qi::_val;

        boost::phoenix::function<lazy_pow_>   lazy_pow;
        boost::phoenix::function<lazy_ufunc_> lazy_ufunc;
        boost::phoenix::function<lazy_bfunc_> lazy_bfunc;

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
            primary                [_val =  _1]
            >> *(  ("**" >> factor [_val = lazy_pow(_val, _1)])
                )
            ;

        primary =
            real                   [_val =  _1]
            |   '(' >> expression  [_val =  _1] >> ')'
            |   ('-' >> primary    [_val = -_1])
            |   ('+' >> primary    [_val =  _1])
            |   no_case[constant]  [_val =  _1]
            |   (no_case[ufunc] >> '(' >> expression >> ')')
                                   [_val = lazy_ufunc(_1, _2)]
            |   (no_case[bfunc] >> '(' >> expression >> ','
                                       >> expression >> ')')
                                   [_val = lazy_bfunc(_1, _2, _3)]
            ;

    }
};

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
