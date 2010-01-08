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
 * common.hpp: C++ common definitions, utility macros, and inline functions
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_COMMON_HPP
#define __SUZERAIN_COMMON_HPP

// Include all of the C common material
#include <suzerain/common.h>

// Required standard C++ functionality used throughout Suzerain
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <ios>
#include <iostream>
#include <iterator>
#include <limits>
#include <new>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <utility>
#include <valarray>
#include <vector>

// Include Boost functionality used throughout Suzerain
#ifdef SUZERAIN_HAVE_BOOST
#include <boost/array.hpp>
#include <boost/concept/assert.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/function.hpp>
#include <boost/integer_traits.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/multi_array.hpp>
#include <boost/noncopyable.hpp>
#include <boost/numeric/conversion/converter.hpp>
#include <boost/preprocessor/comparison/greater.hpp>
#include <boost/preprocessor/comparison/less.hpp>
#include <boost/preprocessor/control/iif.hpp>
#include <boost/preprocessor/facilities/empty.hpp>
#include <boost/preprocessor/seq/elem.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/seq/size.hpp>
#include <boost/preprocessor/seq/to_tuple.hpp>
#include <boost/program_options.hpp>
#include <boost/scoped_array.hpp>
#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/static_assert.hpp>
#include <boost/swap.hpp>
#include <boost/test/utils/nullstream.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility.hpp>
#endif // SUZERAIN_HAVE_BOOST

// Include other functionality used throughout Suzerain
#ifdef SUZERAIN_HAVE_EIGEN
#include <Eigen/Core>
#endif // SUZERAIN_HAVE_EIGEN

#endif // __SUZERAIN_COMMON_HPP
