#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <boost/version.hpp>
#include <suzerain/exprparse.hpp>
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>

#pragma warning(disable:1572)

using suzerain::exprparse;
const double close_enough = std::numeric_limits<double>::epsilon();

#define EXPRTEST(casename, expr, expected)                  \
BOOST_AUTO_TEST_CASE( casename )                            \
{                                                           \
    double result1;                                         \
    exprparse(expr, result1);                               \
    BOOST_CHECK_CLOSE(result1, (expected), close_enough);   \
                                                            \
    double result2;                                         \
    exprparse(std::string(expr), result2);                  \
    BOOST_CHECK_CLOSE(result2, (expected), close_enough);   \
                                                            \
    BOOST_CHECK_EQUAL(result1, result2);                    \
}

EXPRTEST(literal1, "1.234", 1.234)
EXPRTEST(literal2, "4.2e2", 420)
EXPRTEST(literal3, "5e-01", 0.5)
EXPRTEST(literal4, "-3",    -3)

#if BOOST_VERSION >= 104100 // Spirit.Qi 2.1+ available

EXPRTEST(literal5, "pi", boost::math::constants::pi<double>())

EXPRTEST(basicop1, " 2 +\t3\n",   5)       // Whitespace ignored
EXPRTEST(basicop2, " 2 -\t3\n",  -1)
EXPRTEST(basicop3, " 2 *\t3\n",   6)
EXPRTEST(basicop4, " 2 /\t3\n",   2./3.)   // Double division
EXPRTEST(basicop5, " 2 ** 3\n",  8)

EXPRTEST(pemdas1,  "2*3+4*5",   2*3+4*5)
EXPRTEST(pemdas2,  "2*(3+4)*5", 2*(3+4)*5)
EXPRTEST(pemdas3,  "2**3+4",    std::pow(2,3)+4)
EXPRTEST(pemdas4,  "2**2**-3",  std::pow(2.,std::pow(2.,-3.)))

EXPRTEST(unary1,   "-(2)",      -2)
EXPRTEST(unary2,   "-(-2)",      2)
EXPRTEST(unary3,   "+(-2)",     -2)
EXPRTEST(unary4,   "+(+2)",      2)

EXPRTEST(func_abs,   "abs (-1.0)", std::abs(-1.0))
EXPRTEST(func_acos,  "acos( 1.0)", std::acos(1.0))
EXPRTEST(func_asin,  "asin( 1.0)", std::asin(1.0))
EXPRTEST(func_atan,  "atan( 1.0)", std::atan(1.0))
EXPRTEST(func_ceil,  "ceil( 0.5)", std::ceil(0.5))
EXPRTEST(func_cos,   "cos ( 1.0)", std::cos(1.0))
EXPRTEST(func_cosh,  "cosh( 1.0)", std::cosh(1.0))
EXPRTEST(func_exp,   "exp ( 1.0)", std::exp(1.0))
EXPRTEST(func_floor, "floor(0.5)", std::floor(0.5))
EXPRTEST(func_log,   "log ( 1.0)", std::log(1.0))
EXPRTEST(func_log10, "log10(0.5)", std::log10(0.5))
EXPRTEST(func_sin,   "sin ( 1.0)", std::sin(1.0))
EXPRTEST(func_sinh,  "sinh( 1.0)", std::sinh(1.0))
EXPRTEST(func_sqrt,  "sqrt( 1.0)", std::sqrt(1.0))
EXPRTEST(func_tan,   "tan ( 1.0)", std::tan(1.0))
EXPRTEST(func_tanh,  "tanh( 1.0)", std::tanh(1.0))

EXPRTEST(func_pow,   "pow  (2.0, 3.0)", std::pow(2.0,3.0))
EXPRTEST(func_atan2, "atan2(2.0, 3.0)", std::atan2(2.0,3.0))

BOOST_AUTO_TEST_CASE( parse_inf )
{
    double result1;
    exprparse("inf", result1);
    BOOST_CHECK((boost::math::isinf)(result1));

    double result2;
    exprparse(std::string("inf"), result2);
    BOOST_CHECK((boost::math::isinf)(result2));
}

BOOST_AUTO_TEST_CASE( parse_infinity )
{
    double result1;
    exprparse("infinity", result1);
    BOOST_CHECK((boost::math::isinf)(result1));

    double result2;
    exprparse(std::string("infinity"), result2);
    BOOST_CHECK((boost::math::isinf)(result2));
}

BOOST_AUTO_TEST_CASE( parse_nan )
{
    double result1;
    exprparse("nan", result1);
    BOOST_CHECK((boost::math::isnan)(result1));

    double result2;
    exprparse(std::string("nan"), result2);
    BOOST_CHECK((boost::math::isnan)(result2));
}

#endif // BOOST_VERSION
