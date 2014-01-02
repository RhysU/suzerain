//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
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

#ifndef SUZERAIN_ITERATOR_HPP
#define SUZERAIN_ITERATOR_HPP

/** @file
 * Iterator-related utilities.
 */

#include <suzerain/common.hpp>

namespace suzerain {

/**
 * Iterator-related utilities.
 */
namespace iterator {

/** An input_iterator that always returns a fixed constant. */
template<typename T>
class infinite_constant
    : public std::iterator<std::input_iterator_tag, const T>
{
public:
    infinite_constant(const T& t) : t_(t) {}

    infinite_constant(const infinite_constant& other) : t_(other.t_) {}

    infinite_constant& operator++() { return *this; }
    infinite_constant& operator++(int) { return *this; }

    bool operator==(const infinite_constant& rhs) const { return t_ == rhs.t_; }
    bool operator!=(const infinite_constant& rhs) const { return t_ == rhs.t_; }

    const T& operator*() const { return t_; }

private:
    const T t_;
};

template< typename charT, typename traits, typename T >
::std::basic_ostream<charT,traits>& operator<<(
        ::std::basic_ostream<charT,traits> &os,
        const infinite_constant<T> &ic)
{
    os << *ic;
    return os;
}

template<typename T>
infinite_constant<T> make_infinite_constant(const T& t)
{
    return infinite_constant<T>(t);
}

} // namespace iterator

} // namespace suzerain

#endif // SUZERAIN_ITERATOR_HPP
