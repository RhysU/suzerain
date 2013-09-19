//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#ifndef SUZERAIN_EXPRGRAMMAR_HPP
#define SUZERAIN_EXPRGRAMMAR_HPP

/** @file
 * A grammar for constant arithmetic expression evaluation via Boost Spirit
 */

#include <algorithm>
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
 * \ref exprparse() declared in exprparse.hpp to be sufficient and much easier
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

template <class T>
T max_by_value ( const T a, const T b ) {
    return std::max(a, b);
}

template <class T>
T min_by_value ( const T a, const T b ) {
    return std::min(a, b);
}

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
                ("digits",   std::numeric_limits<FPT>::digits    )
                ("digits10", std::numeric_limits<FPT>::digits10  )
                ("e" ,       boost::math::constants::e<FPT>()    )
                ("epsilon",  std::numeric_limits<FPT>::epsilon() )
                ("pi",       boost::math::constants::pi<FPT>()   )
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
                ("abs"  , static_cast<FPT (*)(FPT)>(&std::abs  ))
                ("acos" , static_cast<FPT (*)(FPT)>(&std::acos ))
                ("asin" , static_cast<FPT (*)(FPT)>(&std::asin ))
                ("atan" , static_cast<FPT (*)(FPT)>(&std::atan ))
                ("ceil" , static_cast<FPT (*)(FPT)>(&std::ceil ))
                ("cos"  , static_cast<FPT (*)(FPT)>(&std::cos  ))
                ("cosh" , static_cast<FPT (*)(FPT)>(&std::cosh ))
                ("exp"  , static_cast<FPT (*)(FPT)>(&std::exp  ))
                ("floor", static_cast<FPT (*)(FPT)>(&std::floor))
                ("log"  , static_cast<FPT (*)(FPT)>(&std::log  ))
                ("log10", static_cast<FPT (*)(FPT)>(&std::log10))
                ("sin"  , static_cast<FPT (*)(FPT)>(&std::sin  ))
                ("sinh" , static_cast<FPT (*)(FPT)>(&std::sinh ))
                ("sqrt" , static_cast<FPT (*)(FPT)>(&std::sqrt ))
                ("tan"  , static_cast<FPT (*)(FPT)>(&std::tan  ))
                ("tanh" , static_cast<FPT (*)(FPT)>(&std::tanh ))
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
                ("atan2", static_cast<FPT (*)(FPT, FPT)>(&std::atan2  ))
                ("max"  , static_cast<FPT (*)(FPT, FPT)>(&max_by_value))
                ("min"  , static_cast<FPT (*)(FPT, FPT)>(&min_by_value))
                ("pow"  , static_cast<FPT (*)(FPT, FPT)>(&std::pow    ))
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
            |   (no_case[ufunc] >> '(' >> expression >> ')')
                                   [_val = lazy_ufunc(_1, _2)]
            |   (no_case[bfunc] >> '(' >> expression >> ','
                                       >> expression >> ')')
                                   [_val = lazy_bfunc(_1, _2, _3)]
            |   no_case[constant]  [_val =  _1]
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

#endif // SUZERAIN_EXPRGRAMMAR_HPP
