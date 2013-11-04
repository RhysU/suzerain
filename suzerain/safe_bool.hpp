/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2013 Rhys Ulerich
 * Copyright (C) 2013 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
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
 */

#ifndef SUZERAIN_SAFE_BOOL_HPP
#define SUZERAIN_SAFE_BOOL_HPP

/** @file
 * Helpers for encapsulating <a
 * href="http://www.artima.com/cppsource/safebool.html">The Safe Bool Idiom</a>
 * as described by Bjorn Karlsson.  The solution has additionally been
 * modified per <a href="http://stackoverflow.com/questions/3657494/">a
 * question on Stack Overflow</a> to address compilation issues.
 */

namespace suzerain
{

/**
 * <a href="http://www.artima.com/cppsource/safebool3.html">The Safe Bool
 * Idiom</a> interface.
 */
class safe_bool_base
{
protected:
    typedef void (safe_bool_base::*bool_type)() const;
    void this_type_does_not_support_comparisons() const {}

    safe_bool_base() {}
    safe_bool_base(const safe_bool_base&) {}
    safe_bool_base& operator=(const safe_bool_base&)
    {
        return *this;
    }
    ~safe_bool_base() {}
};

/**
 * A Safe Bool Idiom implementation for <a
 * href="http://www.artima.com/cppsource/safebool3.html">non-virtual use via
 * the Curiously Recurring Template Pattern</a>.
 *
 * @see The \ref safe_bool_base interface for the required contract.
 */
template < typename T = void > class safe_bool : public safe_bool_base
{
public:
    operator bool_type() const
    {
        return (static_cast<const T*>(this))->boolean_test()
               ? &safe_bool<T>::this_type_does_not_support_comparisons : 0;
    }
protected:
    ~safe_bool() {}
};

/**
 * A Safe Bool Idiom implementation for <a
 * href="http://www.artima.com/cppsource/safebool3.html">virtual use</a>.
 *
 * @see The \ref safe_bool_base interface for the required contract.
 */
template<> class safe_bool<void> : public safe_bool_base
{
public:
    operator bool_type() const
    {
        return boolean_test() == true ?
               &safe_bool<void>::this_type_does_not_support_comparisons : 0;
    }
protected:
    virtual bool boolean_test() const = 0;
    virtual ~safe_bool() {}
};

template <typename T, typename U>
void operator==(const safe_bool<T>& lhs, const safe_bool<U>& rhs)
{
    lhs.this_type_does_not_support_comparisons();
    return false;
}

template <typename T, typename U>
void operator!=(const safe_bool<T>& lhs, const safe_bool<U>& rhs)
{
    lhs.this_type_does_not_support_comparisons();
    return false;
}

}

#endif /* SUZERAIN_SAFE_BOOL_HPP */
