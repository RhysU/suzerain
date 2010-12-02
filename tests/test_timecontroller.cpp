#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>
#include <suzerain/timecontroller.hpp>
#include "test_tools.hpp"

class TestTimeController
    : public suzerain::timestepper::AbstractTimeController<double>
{
    double stepTime(double max_dt) const { return max_dt / 10; }
};

BOOST_AUTO_TEST_CASE( nop )
{
    // FIXME Stubbed
}
