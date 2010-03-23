//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2010 The PECOS Development Team
//
// Please see http://pecos.ices.utexas.edu for more information.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------
//
// boost_multi_array_hack.hpp: Hacks used to access Boost.MultiArray details
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef __SUZERAIN_BOOST_MULTI_ARRAY_HACK_HPP
#define __SUZERAIN_BOOST_MULTI_ARRAY_HACK_HPP

// No #include <suzerain/common.hpp> because common.hpp includes us.

#ifndef DOXYGEN_SHOULD_SKIP_THIS

// Message from http://lists.boost.org/boost-users/2010/03/57634.php
// --------------------------------------------------------------------
//
// I'd like to obtain an N-dimensional MultiArray implementation given a
// base pointer, N extents, and N strides. My use case requires padding
// the stride in one dimension in a way seemingly not obtainable taking
// views of the usual boost::multi_array_ref.
//
// Code like the following would be ideal
//
//     boost::array<std::size_t,3> extents = { 2, 3, 4 };
//     boost::array<std::size_t,3> strides = { 1, 2, 7 }; // Note 7 not 6
//     boost::scoped_array<int> raw(new int[extents[2]*strides[2]]);
//
//     using boost::detail::multi_array::multi_array_view;
//     multi_array_view<int,3> a(raw.get(), extents, strides);
//
// except that the appropriate constructor in multi_array_view is private
// (boost/multi_array/view.hpp:442).
// boost::detail::multi_array::sub_array would also be ideal if it's
// constructor was accessible(boost/multi_array/subarray.hpp:370). I can
// make either accessible by #defining BOOST_NO_MEMBER_TEMPLATE_FRIENDS,
// but that's evil.
//
// Am I missing something in Boost.MultiArray? Or is there no publicly
// accessible way to provide a custom stride list?
//
// --------------------------------------------------------------------

// Until a better answer comes back on that thread, or until I cook
// up another way to accomplish what I need, the following is an
// ugly-but-effective hack to obtain access to the necessary constructors.

#ifdef BOOST_MULTI_ARRAY_VIEW_RG071301_HPP
    #error "boost/multi_array/view.hpp included before BOOST_NO_MEMBER_TEMPLATE_FRIENDS hack"
#else
    #ifdef BOOST_NO_MEMBER_TEMPLATE_FRIENDS
        #include <boost/multi_array.hpp>
    #else
        #define BOOST_NO_MEMBER_TEMPLATE_FRIENDS
        #include <boost/multi_array.hpp>
        #undef BOOST_NO_MEMBER_TEMPLATE_FRIENDS
    #endif
#endif
#ifndef BOOST_MULTI_ARRAY_VIEW_RG071301_HPP
    #error "boost/multi_array/view.hpp version incompatible with BOOST_NO_MEMBER_TEMPLATE_FRIENDS hack"
#endif

#endif // DOXYGEN_SHOULD_SKIP_THIS

#endif // __SUZERAIN_BOOST_MULTI_ARRAY_HACK_HPP
