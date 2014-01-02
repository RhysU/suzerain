//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
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

#ifndef SUZERAIN_CONTIGUOUS_MEMORY_HPP
#define SUZERAIN_CONTIGUOUS_MEMORY_HPP

#include <suzerain/common.hpp>

/** @file
 * Provides the contiguous_memory class.
 */

namespace suzerain
{

/**
 * At construction, allocates a contiguous block of memory containing
 * <tt>Element</tt>s allocated using <tt>Allocator</tt>.  At destruction the
 * memory is freed.  This seemingly goofy class is useful because it can be
 * inherited from as a way to workaround base class versus member construction
 * time semantics.
 *
 * @tparam Element   Type of elements to allocate.
 * @tparam Allocator Allocator to use for memory operations.
 */
template<
    typename Element,
    typename Allocator = std::allocator<Element>
>
class contiguous_memory
    : private Allocator // Permits empty base class optimization
{
public:
    /**
     * Create an block of memory suitable for \c count Elements.
     *
     * @param count Number of Elements for which to allocate memory.
     * @param dummy Used to disambiguate some overloads.  Ignore it.
     */
    template< typename I >
    contiguous_memory(
        I count,
        typename boost::enable_if<boost::is_integral<I> >::type* dummy = 0)
        : Allocator(),
          pbegin_(Allocator::allocate(boost::numeric_cast<
                      typename Allocator::size_type>(count))),
          pend_(pbegin_ + count) {
        SUZERAIN_UNUSED(dummy);
    }

    /** Create a block of memory which is a copy of \c other. */
    contiguous_memory(const contiguous_memory &other)
        : Allocator(other),
          pbegin_(Allocator::allocate(boost::numeric_cast<
                      typename Allocator::size_type>(
                      std::distance(other.pbegin_,other.pend_)))),
          pend_(pbegin_ + std::distance(other.pbegin_,other.pend_))
    {
        std::memcpy(pbegin_, other.pbegin_,
            std::distance(other.pbegin_, other.pend_) * sizeof(Element));
    }

    /**
     * Destruct the instance and free the associated memory
     */
    ~contiguous_memory() {
        Allocator::deallocate(pbegin_, std::distance(pbegin_,pend_));
    }

    /** @return a pointer to the first element. */
    typename Allocator::pointer memory_begin() const {
        return pbegin_;
    }

    /** @return a pointer to the last-plus-one element. */
    typename Allocator::pointer memory_end() const {
        return pend_;
    }

    /**
     * Swap the \em ownership of \c this and \c other's memory.
     *
     * This does \em not move any elements between the regions.  It does work
     * even when the two instances contain different numbers of elements.
     */
    void memory_swap(contiguous_memory& other) throw () {
        std::swap(pbegin_, other.pbegin_);
        std::swap(pend_, other.pend_);
    }

    /**
     * Assign *this from \c other
     * @return \c *this.
     */
    contiguous_memory& operator= (contiguous_memory other)
    {
        // FIXME Not appropriate given that instances are expected to be huge

        // Implementation uses copy-and-swap idiom written for copy elision.
        // See https://secure.wikimedia.org/wikibooks/en/wiki/More_C%2B%2B_Idioms/Copy-and-swap
        other.memory_swap (*this);
        return *this;
    }

private:
    typename Allocator::pointer pbegin_;
    typename Allocator::pointer pend_;
};

} // namespace suzerain

#endif // SUZERAIN_CONTIGUOUS_MEMORY_HPP
