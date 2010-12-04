#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>
#include <suzerain/timecontroller.hpp>
#include "test_tools.hpp"

using suzerain::timestepper::TimeController;

// Subclass TimeController to make the logic a little easier to test
// Mainly done to save use from having to specify the stepper over and over
template <typename Integer = unsigned int>
class TestTimeController : public TimeController<double, Integer>
{
public:
    TestTimeController(double initial_t = 0,
                       double min_dt = 1e-8,
                       double max_dt = std::numeric_limits<double>::max())
        : TimeController<double,Integer>(
                boost::bind(&TestTimeController::stepTime, this, _1),
                initial_t, min_dt, max_dt),
          actual_dt(std::numeric_limits<double>::max()) {};

    double actual_dt;

private:
    double stepTime(double max_dt) const
    {
        return std::min(actual_dt, max_dt);
    }
};

// Helper class used to count callbacks
// Noncopyable as we want to ensure we're handling references
template<typename FPT = double, typename Integer = unsigned long>
struct Callback : public boost::noncopyable
{
    bool retval;
    Integer count;
    FPT last_t;
    Integer last_nt;

    Callback(bool retval = true, Integer count = 0)
        : retval(retval), count(count), last_t(0), last_nt(0) {};

    bool operator()(FPT t, Integer nt)
    {
        ++count; last_t  = t; last_nt = nt; return retval;
    };
};

// Helpful constants used within some test cases
static const double       forever_t
    = std::numeric_limits<double>::max();
static const unsigned int forever_nt
    = std::numeric_limits<unsigned int>::max();

// Simple tests of the outer loop logic
BOOST_AUTO_TEST_CASE( basic )
{
    TestTimeController<> tc1;
    BOOST_CHECK_EQUAL(0,  tc1.current_t());
    BOOST_CHECK_EQUAL(0u, tc1.current_nt());

    TestTimeController<> tc2(5);
    BOOST_CHECK_EQUAL(5,  tc2.current_t());
    BOOST_CHECK_EQUAL(0u, tc2.current_nt());

    TestTimeController<> tc3(5, 1e-6);
    BOOST_CHECK_EQUAL(5,    tc3.current_t());
    BOOST_CHECK_EQUAL(0u,   tc3.current_nt());
    BOOST_CHECK_EQUAL(1e-6, tc3.min_dt());
    tc3.min_dt(1e-7);
    BOOST_CHECK_EQUAL(1e-7, tc3.min_dt());

    TestTimeController<> tc4(5, 1e-6, 1e7);
    BOOST_CHECK_EQUAL(5,    tc4.current_t());
    BOOST_CHECK_EQUAL(0u,   tc4.current_nt());
    BOOST_CHECK_EQUAL(1e-6, tc4.min_dt());
    BOOST_CHECK_EQUAL(1e+7, tc4.max_dt());
    tc4.max_dt(1e+8);
    BOOST_CHECK_EQUAL(1e+8, tc4.max_dt());
}

// Simple tests of the outer loop logic
BOOST_AUTO_TEST_CASE( no_callbacks )
{
    TestTimeController<> tc(0, 1e-8, 10);

    // No advance: final_t == current_t
    BOOST_REQUIRE(tc.advance(0.0, 100u));
    BOOST_REQUIRE_EQUAL(0.0, tc.current_t());
    BOOST_REQUIRE_EQUAL(0u,  tc.current_nt());

    // No advance: final_nt == current_nt
    BOOST_REQUIRE(tc.advance(100.0, 0u + tc.current_nt()));
    BOOST_REQUIRE_EQUAL(0.0, tc.current_t());
    BOOST_REQUIRE_EQUAL(0u,  tc.current_nt());

    // Advance halted by number of allowable time steps
    BOOST_REQUIRE(tc.advance(100.0, 9u + tc.current_nt()));
    BOOST_REQUIRE_EQUAL(90.0, tc.current_t());
    BOOST_REQUIRE_EQUAL(9u,   tc.current_nt());

    // Advance halted by final time
    BOOST_REQUIRE(tc.advance(100.0, 1000u + tc.current_nt()));
    BOOST_REQUIRE_EQUAL(100.0, tc.current_t());
    BOOST_REQUIRE_EQUAL(10u,   tc.current_nt());

    // No advance: final_t < current_t
    BOOST_REQUIRE(tc.advance(100.0 - 1.0, 100u + tc.current_nt()));
    BOOST_REQUIRE_EQUAL(100.0, tc.current_t());
    BOOST_REQUIRE_EQUAL(10u,   tc.current_nt());

    // No advance: final_nt < current_nt
    BOOST_REQUIRE(tc.advance(200.0, tc.current_nt() - 1));
    BOOST_REQUIRE_EQUAL(100.0, tc.current_t());
    BOOST_REQUIRE_EQUAL(10u,   tc.current_nt());

    // Advance halted by final time
    tc.max_dt(1);
    BOOST_REQUIRE(tc.advance(200.0, 1000u + tc.current_nt()));
    BOOST_REQUIRE_EQUAL(200.0, tc.current_t());
    BOOST_REQUIRE_EQUAL(110u,  tc.current_nt());
}

BOOST_AUTO_TEST_CASE( simple_periodic_callback_nt )
{
    Callback<> cb;

    TestTimeController<> tc(0, 1e-8, 1);
    tc.add_periodic_callback(forever_t, 10u, boost::ref(cb));

    BOOST_REQUIRE(tc.advance(9.0, forever_nt));
    BOOST_REQUIRE_EQUAL(9.0, tc.current_t());
    BOOST_CHECK_EQUAL(cb.count, 0u);

    BOOST_REQUIRE(tc.advance(15.0, forever_nt));
    BOOST_REQUIRE_EQUAL(15.0, tc.current_t());
    BOOST_CHECK_EQUAL(cb.count, 1u);
    BOOST_CHECK_EQUAL(cb.last_t, 10.0);
    BOOST_CHECK_EQUAL(cb.last_nt, 10u);

    BOOST_REQUIRE(tc.advance(25.0, forever_nt));
    BOOST_REQUIRE_EQUAL(25.0, tc.current_t());
    BOOST_CHECK_EQUAL(cb.count, 2u);
    BOOST_CHECK_EQUAL(cb.last_t, 20.0);
    BOOST_CHECK_EQUAL(cb.last_nt, 20u);
}

BOOST_AUTO_TEST_CASE( two_simultaneous_callbacks_nt )
{
    Callback<> cb1, cb2;

    TestTimeController<> tc(0, 1e-8, 1);
    tc.add_periodic_callback(forever_t, 10u, boost::ref(cb1));
    tc.add_periodic_callback(forever_t, 10u, boost::ref(cb2));

    BOOST_REQUIRE(tc.advance(15.0, forever_nt));
    BOOST_REQUIRE_EQUAL(15.0, tc.current_t());
    BOOST_CHECK_EQUAL(cb1.count, 1u);
    BOOST_CHECK_EQUAL(cb1.last_t, 10.0);
    BOOST_CHECK_EQUAL(cb1.last_nt, 10u);
    BOOST_CHECK_EQUAL(cb2.count, 1u);
    BOOST_CHECK_EQUAL(cb2.last_t, 10.0);
    BOOST_CHECK_EQUAL(cb2.last_nt, 10u);

    BOOST_REQUIRE(tc.advance(25.0, forever_nt));
    BOOST_REQUIRE_EQUAL(25.0, tc.current_t());
    BOOST_CHECK_EQUAL(cb1.count, 2u);
    BOOST_CHECK_EQUAL(cb1.last_t, 20.0);
    BOOST_CHECK_EQUAL(cb1.last_nt, 20u);
    BOOST_CHECK_EQUAL(cb2.count, 2u);
    BOOST_CHECK_EQUAL(cb2.last_t, 20.0);
    BOOST_CHECK_EQUAL(cb2.last_nt, 20u);
}

BOOST_AUTO_TEST_CASE( two_different_callbacks_nt )
{
    Callback<> cb1, cb2;

    TestTimeController<> tc(0, 1e-8, 1);
    tc.add_periodic_callback(forever_t, 10u, boost::ref(cb1));
    tc.add_periodic_callback(forever_t,  5u, boost::ref(cb2));

    BOOST_REQUIRE(tc.advance(15.0, forever_nt));
    BOOST_REQUIRE_EQUAL(15.0, tc.current_t());
    BOOST_CHECK_EQUAL(cb1.count, 1u);
    BOOST_CHECK_EQUAL(cb1.last_t, 10.0);
    BOOST_CHECK_EQUAL(cb1.last_nt, 10u);
    BOOST_CHECK_EQUAL(cb2.count, 3u);
    BOOST_CHECK_EQUAL(cb2.last_t, 15.0);
    BOOST_CHECK_EQUAL(cb2.last_nt, 15u);

    BOOST_REQUIRE(tc.advance(25.0, forever_nt));
    BOOST_REQUIRE_EQUAL(25.0, tc.current_t());
    BOOST_CHECK_EQUAL(cb1.count, 2u);
    BOOST_CHECK_EQUAL(cb1.last_t, 20.0);
    BOOST_CHECK_EQUAL(cb1.last_nt, 20u);
    BOOST_CHECK_EQUAL(cb2.count, 5u);
    BOOST_CHECK_EQUAL(cb2.last_t, 25.0);
    BOOST_CHECK_EQUAL(cb2.last_nt, 25u);
}

BOOST_AUTO_TEST_CASE( simple_periodic_callback_t )
{
    Callback<> cb;

    TestTimeController<> tc(0, 1e-8, 1);
    tc.add_periodic_callback(10.0, forever_nt, boost::ref(cb));

    BOOST_REQUIRE(tc.advance(forever_t, 9u));
    BOOST_REQUIRE_EQUAL(9.0, tc.current_t());
    BOOST_CHECK_EQUAL(cb.count, 0u);

    BOOST_REQUIRE(tc.advance(forever_t, 15u));
    BOOST_REQUIRE_EQUAL(15.0, tc.current_t());
    BOOST_CHECK_EQUAL(cb.count, 1u);
    BOOST_CHECK_EQUAL(cb.last_t, 10.0);
    BOOST_CHECK_EQUAL(cb.last_nt, 10u);

    BOOST_REQUIRE(tc.advance(forever_t, 25u));
    BOOST_REQUIRE_EQUAL(25.0, tc.current_t());
    BOOST_CHECK_EQUAL(cb.count, 2u);
    BOOST_CHECK_EQUAL(cb.last_t, 20.0);
    BOOST_CHECK_EQUAL(cb.last_nt, 20u);
}

BOOST_AUTO_TEST_CASE( two_simultaneous_callbacks_t )
{
    Callback<> cb1, cb2;

    TestTimeController<> tc(0, 1e-8, 1);
    tc.add_periodic_callback(10.0, forever_nt, boost::ref(cb1));
    tc.add_periodic_callback(10.0, forever_nt, boost::ref(cb2));

    BOOST_REQUIRE(tc.advance(forever_t, 15u));
    BOOST_REQUIRE_EQUAL(15.0, tc.current_t());
    BOOST_CHECK_EQUAL(cb1.count, 1u);
    BOOST_CHECK_EQUAL(cb1.last_t, 10.0);
    BOOST_CHECK_EQUAL(cb1.last_nt, 10u);
    BOOST_CHECK_EQUAL(cb2.count, 1u);
    BOOST_CHECK_EQUAL(cb2.last_t, 10.0);
    BOOST_CHECK_EQUAL(cb2.last_nt, 10u);

    BOOST_REQUIRE(tc.advance(forever_t, 25u));
    BOOST_REQUIRE_EQUAL(25.0, tc.current_t());
    BOOST_CHECK_EQUAL(cb1.count, 2u);
    BOOST_CHECK_EQUAL(cb1.last_t, 20.0);
    BOOST_CHECK_EQUAL(cb1.last_nt, 20u);
    BOOST_CHECK_EQUAL(cb2.count, 2u);
    BOOST_CHECK_EQUAL(cb2.last_t, 20.0);
    BOOST_CHECK_EQUAL(cb2.last_nt, 20u);
}

BOOST_AUTO_TEST_CASE( two_different_callbacks_t )
{
    Callback<> cb1, cb2;

    TestTimeController<> tc(0, 1e-8, 1);
    tc.add_periodic_callback(10.0, forever_nt, boost::ref(cb1));
    tc.add_periodic_callback( 5.0, forever_nt, boost::ref(cb2));

    BOOST_REQUIRE(tc.advance(forever_t, 15u));
    BOOST_REQUIRE_EQUAL(15.0, tc.current_t());
    BOOST_CHECK_EQUAL(cb1.count, 1u);
    BOOST_CHECK_EQUAL(cb1.last_t, 10.0);
    BOOST_CHECK_EQUAL(cb1.last_nt, 10u);
    BOOST_CHECK_EQUAL(cb2.count, 3u);
    BOOST_CHECK_EQUAL(cb2.last_t, 15.0);
    BOOST_CHECK_EQUAL(cb2.last_nt, 15u);

    BOOST_REQUIRE(tc.advance(forever_t, 25u));
    BOOST_REQUIRE_EQUAL(25.0, tc.current_t());
    BOOST_CHECK_EQUAL(cb1.count, 2u);
    BOOST_CHECK_EQUAL(cb1.last_t, 20.0);
    BOOST_CHECK_EQUAL(cb1.last_nt, 20u);
    BOOST_CHECK_EQUAL(cb2.count, 5u);
    BOOST_CHECK_EQUAL(cb2.last_t, 25.0);
    BOOST_CHECK_EQUAL(cb2.last_nt, 25u);
}

BOOST_AUTO_TEST_CASE( mixed_periodic_callback_criteria )
{
    Callback<> cb;

    TestTimeController<> tc(0); // Different from other cases
    tc.add_periodic_callback(10.0, 5u, boost::ref(cb));

    // Deliberately trip the nt-based criteria
    tc.max_dt(0.5);
    BOOST_REQUIRE(tc.advance(3.0, forever_nt));
    BOOST_REQUIRE_EQUAL(3.0, tc.current_t());
    BOOST_CHECK_EQUAL(cb.count, 1u);
    BOOST_CHECK_EQUAL(cb.last_t, 2.5);
    BOOST_CHECK_EQUAL(cb.last_nt, 5u);

    // Deliberately trip the t-based criteria
    tc.max_dt(5.0);
    BOOST_REQUIRE(tc.advance(13.0, forever_nt));
    BOOST_REQUIRE_EQUAL(13.0, tc.current_t());
    BOOST_REQUIRE_EQUAL(tc.current_nt(), 9u);
    BOOST_CHECK_EQUAL(cb.count, 2u);
    BOOST_CHECK_EQUAL(cb.last_t, 12.5);
    BOOST_CHECK_EQUAL(cb.last_nt, 8u);
}

BOOST_AUTO_TEST_CASE( callback_initiated_abort )
{
    Callback<> cb(false); // Callback returns false when invoked

    TestTimeController<> tc(0, 1e-8, 1);
    tc.add_periodic_callback(forever_t, 10u, boost::ref(cb));

    // Successful advance prior to callback
    BOOST_REQUIRE(tc.advance(8.0, forever_nt));
    BOOST_REQUIRE_EQUAL(8.0, tc.current_t());
    BOOST_CHECK_EQUAL(cb.count, 0u);

    // Callback returns false and stops advance
    BOOST_REQUIRE(false == tc.advance(16.0, forever_nt));
    BOOST_REQUIRE_EQUAL(10.0, tc.current_t());
    BOOST_CHECK_EQUAL(cb.count, 1u);
}

BOOST_AUTO_TEST_CASE( physics_initiated_abort )
{
    TestTimeController<> tc(0.0, 1 /* min_dt */, 2 /* max_dt */);

    // Successful advance for step smaller than min_dt
    // when the controller demands it.
    BOOST_REQUIRE(tc.advance(0.5, forever_nt));
    BOOST_REQUIRE_EQUAL(0.5, tc.current_t());

    // Force "physics" to return a dt < dt_min
    tc.actual_dt = 0.5;

    // Unsuccessful advance because physics dt < dt_min
    BOOST_REQUIRE(!tc.advance(forever_t, forever_nt));
    BOOST_REQUIRE_EQUAL(1.0, tc.current_t());
}

// Sane behavior when every_values should cause floating point overflow?
BOOST_AUTO_TEST_CASE( every_values_too_large )
{
    TestTimeController<> tc(0.0, 1e-8, std::numeric_limits<double>::max());
    BOOST_REQUIRE(tc.advance(1.0, forever_nt));
    BOOST_REQUIRE_EQUAL(1.0, tc.current_t());
    BOOST_REQUIRE_EQUAL(1u,  tc.current_nt());

    // Now every_XX = forever should cause a numeric overflow
    // Code should silently coerce these to maximum representable value
    Callback<> cb;
    tc.add_periodic_callback(forever_t, forever_nt, boost::ref(cb));

    // Advance to the end of days in a single step
    BOOST_REQUIRE(tc.advance(forever_t, forever_nt));
    BOOST_REQUIRE_EQUAL(forever_t,  tc.current_t());
    BOOST_REQUIRE_EQUAL(2u, tc.current_nt());

    // Check that the callback was only hit at the very end
    BOOST_CHECK_EQUAL(cb.count, 1u);
    BOOST_CHECK_EQUAL(cb.last_t,  forever_t);
    BOOST_CHECK_EQUAL(cb.last_nt, 2u);
}

// Sane behavior when we hit the maximum number of time steps?
BOOST_AUTO_TEST_CASE( largest_possible_time_step_count )
{
    // Use an artificially small step counting type to avoid being here all day
    const unsigned char forever_nt = std::numeric_limits<unsigned char>::max();

    // Set up a controller which will definitely exhaust forever_nt;
    const double epsilon = std::numeric_limits<double>::epsilon();
    TestTimeController<unsigned char> tc(
            0.0, epsilon, (1.0 / forever_nt) - epsilon);
    BOOST_REQUIRE_LE(tc.min_dt(), tc.max_dt());

    // Advance to the end of days in many steps
    BOOST_REQUIRE(tc.advance(forever_t, forever_nt));
    BOOST_REQUIRE_EQUAL(forever_nt, tc.current_nt());

    // There's nothing out there, really
    BOOST_REQUIRE(tc.advance(forever_t, forever_nt));
    BOOST_REQUIRE_EQUAL(forever_nt, tc.current_nt());
}
