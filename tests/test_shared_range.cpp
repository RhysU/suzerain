#include <suzerain/shared_range.hpp>

#include <algorithm>
#include <iterator>
#include <numeric>
#include <memory>
#include <boost/array.hpp>
#include <boost/concept/assert.hpp>
#include <boost/concept_check.hpp>
#include <boost/iterator/iterator_concepts.hpp>
#include <boost/range/concepts.hpp>
#include <boost/range.hpp>

#define BOOST_TEST_MODULE test_shared_range
#include <boost/test/included/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>

// Explicit instantiation to flush out compilation errors
template class suzerain::shared_range<int>;

// "Type under test"
typedef suzerain::shared_range<int> tut;

// Basic type requirements
BOOST_CONCEPT_ASSERT((boost::Assignable<tut>));
BOOST_CONCEPT_ASSERT((boost::SGIAssignable<tut>));
BOOST_CONCEPT_ASSERT((boost::DefaultConstructible<tut>));
BOOST_CONCEPT_ASSERT((boost::CopyConstructible<tut>));
BOOST_CONCEPT_ASSERT((boost::EqualityComparable<tut>));
BOOST_CONCEPT_ASSERT((boost::LessThanComparable<tut>));

// Does implementation satisfy the desired range concepts?
BOOST_CONCEPT_ASSERT((boost::SinglePassRangeConcept<tut>));
BOOST_CONCEPT_ASSERT((boost::ForwardRangeConcept<tut>));
BOOST_CONCEPT_ASSERT((boost::BidirectionalRangeConcept<tut>));
BOOST_CONCEPT_ASSERT((boost::RandomAccessRangeConcept<tut>));

// Do implementation's iterators satisfy the desired concepts?
BOOST_CONCEPT_ASSERT((boost_concepts::ReadableIteratorConcept<
    boost::range_iterator<tut>::type>));
BOOST_CONCEPT_ASSERT((boost_concepts::WritableIteratorConcept<
    boost::range_iterator<tut>::type>));
BOOST_CONCEPT_ASSERT((boost_concepts::SwappableIteratorConcept<
    boost::range_iterator<tut>::type>));
BOOST_CONCEPT_ASSERT((boost_concepts::LvalueIteratorConcept<
    boost::range_iterator<tut>::type>));

////////////////////////
// Executable test cases
////////////////////////

// Yes, it is funky that tests use non-range algorithms.
// It is because we want tests to compile on pre-1.42 boost.

BOOST_AUTO_TEST_CASE( default_construction ) // constructor1
{
    tut r; // Many iterator_range methods fail asserts as r.is_singular()
}

BOOST_AUTO_TEST_CASE( constructor2_and_basic_methods )
{
    static int data[] = { 1, 2, 3, 4, 5 };
    static const boost::iterator_range<int*> testrange(data, data + 5);

    // Instantiate and populate
    tut r(new int[5], 5);
    std::fill(r.begin(), r.end(), 1);
    std::partial_sum(r.begin(), r.end(), r.begin());

    // Member methods
    BOOST_CHECK_EQUAL(5, std::distance(r.begin(), r.end()));
    BOOST_CHECK(!!r);
    BOOST_CHECK(r.equal(r));
    BOOST_CHECK(!r.equal(testrange));
    BOOST_CHECK(!testrange.equal(r));
    BOOST_CHECK_EQUAL(r.front(), 1);
    BOOST_CHECK_EQUAL(r.back(),  5);
    BOOST_CHECK_EQUAL(r.empty(), false);
    BOOST_CHECK_EQUAL(r[0],      1);
    BOOST_CHECK_EQUAL(r[1],      2);
    BOOST_CHECK_EQUAL(r[2],      3);
    BOOST_CHECK_EQUAL(r[3],      4);
    BOOST_CHECK_EQUAL(r[4],      5);
    BOOST_CHECK_EQUAL(r(0),      1);
    BOOST_CHECK_EQUAL(r(1),      2);
    BOOST_CHECK_EQUAL(r(2),      3);
    BOOST_CHECK_EQUAL(r(3),      4);
    BOOST_CHECK_EQUAL(r(4),      5);
    BOOST_CHECK_EQUAL(r.size(),  5);
    {
        boost::test_tools::output_test_stream output;
        output << r;
        BOOST_CHECK(output.is_equal("12345"));
    }

    // Check comparison operators
    BOOST_CHECK((r == r));
    BOOST_CHECK((r == testrange));
    BOOST_CHECK((testrange == r));
    BOOST_CHECK(!(r < r));
    data[2] = 2; // Mutate test data
    BOOST_CHECK((testrange < r));
    BOOST_CHECK(!(r < testrange));
    BOOST_CHECK((r != testrange));
    BOOST_CHECK((testrange != r));

    // Mutate via advance_begin
    r.advance_begin(2);
    BOOST_CHECK_EQUAL(3, std::distance(r.begin(), r.end()));
    BOOST_CHECK_EQUAL(r.front(), 3);
    BOOST_CHECK_EQUAL(r[0],      3);
    BOOST_CHECK_EQUAL(r.back(),  5);
    BOOST_CHECK_EQUAL(r[2],      5);
    BOOST_CHECK_EQUAL(r.size(),  3);
    {
        boost::test_tools::output_test_stream output;
        output << r;
        BOOST_CHECK(output.is_equal("345"));
    }

    // Mutate via advance_end
    r.advance_end(-2);
    BOOST_CHECK_EQUAL(1, std::distance(r.begin(), r.end()));
    BOOST_CHECK_EQUAL(r.front(), 3);
    BOOST_CHECK_EQUAL(r[0],      3);
    BOOST_CHECK_EQUAL(r.back(),  3);
    BOOST_CHECK_EQUAL(r.size(),  1);
    {
        boost::test_tools::output_test_stream output;
        output << r;
        BOOST_CHECK(output.is_equal("3"));
    }

    // Mutate into zero-length range
    r.advance_begin(1);
    BOOST_CHECK_EQUAL(0, std::distance(r.begin(), r.end()));
    BOOST_CHECK_EQUAL(r.empty(), true);
    BOOST_CHECK_EQUAL(r.size(),  0);
    {
        boost::test_tools::output_test_stream output;
        output << r;
        BOOST_CHECK(output.is_empty());
    }
}

// constructor3 is trivially different from constructor2

BOOST_AUTO_TEST_CASE( constructor4_and_abbreviated_methods )
{
    static int data[] = { 1, 2, 3, 4, 5 };
    static const boost::iterator_range<int*> testrange(data, data + 5);

    int *x = new int[5];
    tut r(x, x+5);
    std::fill(r.begin(), r.end(), 1);
    std::partial_sum(r.begin(), r.end(), r.begin());

    // Small number of member methods
    BOOST_CHECK_EQUAL(5, std::distance(r.begin(), r.end()));
    BOOST_CHECK(!!r);
    BOOST_CHECK(r.equal(r));
    BOOST_CHECK(!r.equal(testrange));
    BOOST_CHECK(!testrange.equal(r));
    BOOST_CHECK_EQUAL(r.front(),  1);
    BOOST_CHECK_EQUAL(r.back(),   5);
    BOOST_CHECK_EQUAL(r.empty(),  false);
    BOOST_CHECK_EQUAL(r.unique(), true);
    BOOST_CHECK_EQUAL(r[0],       1);
    BOOST_CHECK_EQUAL(r[1],       2);
    BOOST_CHECK_EQUAL(r[2],       3);
    BOOST_CHECK_EQUAL(r[3],       4);
    BOOST_CHECK_EQUAL(r[4],       5);
    BOOST_CHECK_EQUAL(r.size(),   5);

    // Small number of comparison operators
    BOOST_CHECK((r == r));
    BOOST_CHECK((r == testrange));
    BOOST_CHECK((testrange == r));
}

// constructor5 is trivially different from constructor4

BOOST_AUTO_TEST_CASE( construct_from_shared_array ) // constructor6
{
    static int data[] = { 1, 2, 3, 4, 5 };
    static const boost::iterator_range<int*> testrange(data, data + 5);

    boost::shared_array<int> a(new int[5]);
    tut r(a, 5);
    std::fill(r.begin(), r.end(), 1);
    std::partial_sum(r.begin(), r.end(), r.begin());

    // Small number of member methods
    BOOST_CHECK_EQUAL(5, std::distance(r.begin(), r.end()));
    BOOST_CHECK(!!r);
    BOOST_CHECK(r.equal(r));
    BOOST_CHECK(!r.equal(testrange));
    BOOST_CHECK(!testrange.equal(r));
    BOOST_CHECK_EQUAL(r.front(),  1);
    BOOST_CHECK_EQUAL(r.back(),   5);
    BOOST_CHECK_EQUAL(r.empty(),  false);
    BOOST_CHECK_EQUAL(r.unique(), false);
    BOOST_CHECK_EQUAL(r[0],       1);
    BOOST_CHECK_EQUAL(r[1],       2);
    BOOST_CHECK_EQUAL(r[2],       3);
    BOOST_CHECK_EQUAL(r[3],       4);
    BOOST_CHECK_EQUAL(r[4],       5);
    BOOST_CHECK_EQUAL(r.size(),   5);
}

BOOST_AUTO_TEST_CASE( zero_size_construction ) // constructors 2 and 4
{
    tut r(new int[0], static_cast<std::size_t>(0));
    BOOST_CHECK_EQUAL(r.size(),   0);
    BOOST_CHECK_EQUAL(r.empty(),  true);
    BOOST_CHECK_EQUAL(r.unique(), true);
    BOOST_CHECK(!r);

    int *x = new int[0];
    tut s(x, x);
    BOOST_CHECK_EQUAL(s.size(),   0);
    BOOST_CHECK_EQUAL(s.empty(),  true);
    BOOST_CHECK_EQUAL(s.unique(), true);
    BOOST_CHECK(!s);

    // Swap seems to tickle edge cases in iterator range
    r.swap(s);
    BOOST_CHECK(r.begin() == x);
    BOOST_CHECK(r.end()   == x);
}

BOOST_AUTO_TEST_CASE( copy_like_construction )
{
    static int data[] = { 1, 2, 3, 4, 5 };
    static const boost::iterator_range<int*> testrange(data, data + 5);

    // Construction does not manage resources
    tut r(testrange);
    BOOST_CHECK_EQUAL(r.begin(),  testrange.begin());
    BOOST_CHECK_EQUAL(r.end(),    testrange.end());
    BOOST_CHECK_EQUAL(r.unique(), false);
    BOOST_CHECK(!!r);
}

BOOST_AUTO_TEST_CASE( copy_construction )
{
    tut r(new int[5], 5);
    std::fill(r.begin(), r.end(), 1);
    std::partial_sum(r.begin(), r.end(), r.begin());
    BOOST_CHECK_EQUAL(r.unique(), true);

    // Copy construct (with shared_ptr-like semantics)
    tut s(r);
    BOOST_CHECK_EQUAL(s.begin(), r.begin());
    BOOST_CHECK_EQUAL(s.end(),   r.end());
    BOOST_CHECK_EQUAL(r.unique(), false);
    BOOST_CHECK_EQUAL(s.unique(), false);

    // Mutate begin/end of r and see s remains unaffected
    r.advance_begin(2);
    r.advance_end(-2);
    BOOST_CHECK_EQUAL(r[0], 3);
    BOOST_CHECK_EQUAL(s[0], 1);
    BOOST_CHECK_EQUAL(s[1], 2);
    BOOST_CHECK_EQUAL(s[2], 3);
    BOOST_CHECK_EQUAL(s[3], 4);
    BOOST_CHECK_EQUAL(s[4], 5);

    // Mutate data from r and see it reflected in s
    r[0] = 7;
    BOOST_CHECK_EQUAL(r[0], 7);
    BOOST_CHECK_EQUAL(s[2], 7);

    // Reset r and see that s becomes unique
    r.reset();
    BOOST_CHECK_EQUAL(r.empty(),  true);
    BOOST_CHECK_EQUAL(s.empty(),  false);
    BOOST_CHECK_EQUAL(s.unique(), true);
}

BOOST_AUTO_TEST_CASE( swap_both_nonsingular )  // Cases for pre-1.43
{
    tut r(new int[5], 5);
    std::fill(r.begin(), r.end(), 1);
    std::partial_sum(r.begin(), r.end(), r.begin());

    tut s(new int[5], 5);
    std::fill(s.begin(), s.end(), 2);
    std::partial_sum(s.begin(), s.end(), s.begin());

    for (int i = 0; i < 5; ++i) BOOST_CHECK_EQUAL(r[i],   (i+1));
    for (int i = 0; i < 5; ++i) BOOST_CHECK_EQUAL(s[i], 2*(i+1));
    swap(r, s);
    for (int i = 0; i < 5; ++i) BOOST_CHECK_EQUAL(r[i], 2*(i+1));
    for (int i = 0; i < 5; ++i) BOOST_CHECK_EQUAL(s[i],   (i+1));

    BOOST_CHECK(!!r);
    BOOST_CHECK(!!s);
    s.reset();
    BOOST_CHECK(!!r);
    BOOST_CHECK(!s);
}

BOOST_AUTO_TEST_CASE( swap_lhs_singular )  // Cases for pre-1.43
{
    tut r;

    tut s(new int[5], 5);
    std::fill(s.begin(), s.end(), 2);
    std::partial_sum(s.begin(), s.end(), s.begin());

    BOOST_CHECK(!r);
    BOOST_CHECK(!!s);
    for (int i = 0; i < 5; ++i) BOOST_CHECK_EQUAL(s[i], 2*(i+1));
    swap(r, s);
    for (int i = 0; i < 5; ++i) BOOST_CHECK_EQUAL(r[i], 2*(i+1));
    BOOST_CHECK(!!r);
    BOOST_CHECK(!s);
}

BOOST_AUTO_TEST_CASE( swap_rhs_singular )  // Cases for pre-1.43
{
    tut r(new int[5], 5);
    std::fill(r.begin(), r.end(), 2);
    std::partial_sum(r.begin(), r.end(), r.begin());

    tut s;

    BOOST_CHECK(!!r);
    BOOST_CHECK(!s);
    for (int i = 0; i < 5; ++i) BOOST_CHECK_EQUAL(r[i], 2*(i+1));
    swap(r, s);
    for (int i = 0; i < 5; ++i) BOOST_CHECK_EQUAL(s[i], 2*(i+1));
    BOOST_CHECK(!r);
    BOOST_CHECK(!!s);
}

BOOST_AUTO_TEST_CASE( swap_both_singular )  // Cases for pre-1.43
{
    tut r, s;
    BOOST_CHECK(!r);
    BOOST_CHECK(!s);
    swap(r, s);
    BOOST_CHECK(!r);
    BOOST_CHECK(!s);
}

BOOST_AUTO_TEST_CASE( copy_assignment )
{
    tut r(new int[5], 5);
    std::fill(r.begin(), r.end(), 1);
    std::partial_sum(r.begin(), r.end(), r.begin());
    BOOST_CHECK_EQUAL(r.unique(), true);

    tut s;
    s = r;
    BOOST_CHECK_EQUAL(r.unique(), false);
    BOOST_CHECK_EQUAL(s.unique(), false);
    BOOST_CHECK_EQUAL(s.begin(), r.begin());
    BOOST_CHECK_EQUAL(s.end(),   r.end());
}

BOOST_AUTO_TEST_CASE( assignment )
{
    static int data[] = { 1, 2, 3, 4, 5 };
    static const boost::iterator_range<int*> testrange(data, data + 5);

    tut r;

    r = testrange;
    BOOST_CHECK_EQUAL(r.unique(), false);
    BOOST_CHECK_EQUAL(r.begin(),  testrange.begin());
    BOOST_CHECK_EQUAL(r.end(),    testrange.end());
}

BOOST_AUTO_TEST_CASE( make_helper )
{
    using suzerain::make_shared_range;

    tut r = make_shared_range<int>(5);
    BOOST_CHECK(!!r);
    BOOST_CHECK_EQUAL(r.empty(),  false);
    BOOST_CHECK_EQUAL(r.unique(), true);
    BOOST_CHECK_EQUAL(r.size(),   5);

    r = make_shared_range<int>(3);
    BOOST_CHECK(!!r);
    BOOST_CHECK_EQUAL(r.empty(),  false);
    BOOST_CHECK_EQUAL(r.unique(), true);
    BOOST_CHECK_EQUAL(r.size(),   3);
}

BOOST_AUTO_TEST_CASE( allocate_helper )
{
    using suzerain::allocate_shared_range;

    tut r = allocate_shared_range(std::allocator<int>(), 5); // Default value
    BOOST_CHECK(!!r);
    BOOST_CHECK_EQUAL(r.empty(),  false);
    BOOST_CHECK_EQUAL(r.unique(), true);
    BOOST_CHECK_EQUAL(r.size(),   5);
    BOOST_CHECK_EQUAL(r[0], 0);
    BOOST_CHECK_EQUAL(r[1], 0);
    BOOST_CHECK_EQUAL(r[2], 0);
    BOOST_CHECK_EQUAL(r[3], 0);
    BOOST_CHECK_EQUAL(r[4], 0);

    r = allocate_shared_range(std::allocator<int>(), 5, 1); // Specified value
    BOOST_CHECK(!!r);
    BOOST_CHECK_EQUAL(r.empty(),  false);
    BOOST_CHECK_EQUAL(r.unique(), true);
    BOOST_CHECK_EQUAL(r.size(),   5);
    BOOST_CHECK_EQUAL(r[0], 1);
    BOOST_CHECK_EQUAL(r[1], 1);
    BOOST_CHECK_EQUAL(r[2], 1);
    BOOST_CHECK_EQUAL(r[3], 1);
    BOOST_CHECK_EQUAL(r[4], 1);
}

BOOST_AUTO_TEST_CASE( use_clone )
{
    using suzerain::clone_shared_range;

    tut r(new int[3], 3);
    std::fill(r.begin(), r.end(), 1);
    std::partial_sum(r.begin(), r.end(), r.begin());

    tut s = clone_shared_range(r);
    BOOST_CHECK_EQUAL(s[0],     1);
    BOOST_CHECK_EQUAL(s[1],     2);
    BOOST_CHECK_EQUAL(s[2],     3);
    BOOST_CHECK_EQUAL(s.size(), 3);

    BOOST_CHECK_NE(r.begin(), s.begin());
    BOOST_CHECK_NE(r.end(),   s.end());
    BOOST_CHECK_EQUAL(r.unique(), true);
    BOOST_CHECK_EQUAL(s.unique(), true);
}

BOOST_AUTO_TEST_CASE( clone_with_custom_allocator )
{
    using suzerain::clone_shared_range;

    tut r(new int[3], 3);
    std::fill(r.begin(), r.end(), 1);
    std::partial_sum(r.begin(), r.end(), r.begin());

    tut s = clone_shared_range(std::allocator<int>(), r);
    BOOST_CHECK_EQUAL(s[0],     1);
    BOOST_CHECK_EQUAL(s[1],     2);
    BOOST_CHECK_EQUAL(s[2],     3);
    BOOST_CHECK_EQUAL(s.size(), 3);

    BOOST_CHECK_NE(r.begin(), s.begin());
    BOOST_CHECK_NE(r.end(),   s.end());
    BOOST_CHECK_EQUAL(r.unique(), true);
    BOOST_CHECK_EQUAL(s.unique(), true);
}
