//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2010 The PECOS Development Team
//
// Please see http://pecos.ices.utexas.edu for more information.
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
//
// storage.hpp: Marker classes indicating storage ordering
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef __SUZERAIN_STORAGE_HPP
#define __SUZERAIN_STORAGE_HPP

#include <suzerain/common.hpp>
#include <suzerain/iterator.hpp>
#include <suzerain/mpl.hpp>

/** @file
 * Provides marker types indicating state element type and element storage
 * configurations.
 **/

namespace suzerain
{

/**
 * Provides marker types indicating state element type and element storage
 * configurations.
 **/
namespace storage
{

// FIXME Document
template< typename Element, typename StorageOrderSequence >
class General {
private:
    typedef typename suzerain::mpl::sequence_array<StorageOrderSequence>
            sequence_array_type;

public:
    typedef Element element_type;
    typedef StorageOrderSequence storage_order_sequence;
    static const std::size_t dimensionality = sequence_array_type::static_size;

    typedef boost::general_storage_order<dimensionality>
            element_storage_order_type;

    static element_storage_order_type element_storage_order() {
        element_storage_order_type result(
                sequence_array_type().begin(),
                suzerain::iterator::make_infinite_constant(true));
        return result;
    }
};

// FIXME Document
template< typename Element >
struct Interleaved
    : public General<Element, boost::mpl::vector_c<std::size_t,0,1,2> > {};

// FIXME Document
template< typename Element >
struct NonInterleaved
    : public General<Element, boost::mpl::vector_c<std::size_t,1,2,0> > {};

} // namespace storage

} // namespace suzerain

#endif // __SUZERAIN_STORAGE_HPP
