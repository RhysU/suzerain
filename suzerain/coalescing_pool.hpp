//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#ifndef SUZERAIN_COALESCING_POOL_HPP
#define SUZERAIN_COALESCING_POOL_HPP

/** @file
 * Memory pools for contiguous state storage.
 */

#include <algorithm>
#include <cstring>
#include <functional>
#include <iosfwd>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <boost/intrusive/link_mode.hpp>
#include <boost/intrusive/set_hook.hpp>
#include <boost/intrusive/set.hpp>
#include <boost/next_prior.hpp>
#include <boost/ptr_container/ptr_set.hpp>
#include <boost/utility.hpp>

namespace suzerain {

/**
 * Default policies used within coalescing_pool.  Template specializations may
 * be used to globally change the default policy for a specific type.
 *
 * @see The template coalescing_pool for how such a policy is used.
 */
template <typename T>
struct coalescing_pool_policy
{
    /** Should all newly allocated arenas be zeroed? */
    static const bool zero_new_arenas
            = true;

    /**
     * What Boost Intrusive linking mode should be used for free list entries?
     * Valid options are \c normal_link, \c safe_link, or \c auto_unlink.
     *
     * @see Boost Intrusive for more details.
     */
    static const ::boost::intrusive::link_mode_type link_mode_type
            = boost::intrusive::normal_link;
};

/**
 * A free list-based memory pool designed for memory-constrained situations
 * where large, contiguous blocks are repeatedly acquired and released.  The
 * pool grows _as slowly as possible_ and in a way that optimizes serving
 * large, contiguous requests relative to the memory resources consumed.  This
 * policy of minimal growth and maximal contiguity makes it differ from many
 * other pool implementations.  Its simple implementation uses the Boost
 * Intrusive and Pointer Container libraries.
 *
 * Memory is returned to the operating system when the pool is either destroyed
 * or explicitly drained.  Thread safety, alignment, and object construction
 * and destruction issues are left to the user.  Block sizes below \c
 * minimum_blocksize may be padded to provide enough working space for free
 * list headers.
 *
 * @tparam T              Data type used by the pool.  All size_type arguments
 *                        work in units of <tt>sizeof(T)</tt>.
 * @tparam ArenaAllocator Allocator used to obtain memory regions from the
 *                        operating system, termed "arenas".
 * @tparam Policy         Policies driving the pool implementation.
 */
template<
    typename T,
    class ArenaAllocator = ::std::allocator<T>,
    class Policy = coalescing_pool_policy<T>
>
class coalescing_pool : public ::boost::noncopyable
{
private:

    // Forward declaration
    class fl_entry;

public:

    // Forward declaration
    class blocks;

    /** @{ */
    typedef ArenaAllocator                           allocator_type;
    typedef typename allocator_type::size_type       size_type;
    typedef typename allocator_type::difference_type difference_type;
    typedef typename allocator_type::value_type      value_type;
    typedef typename allocator_type::pointer         pointer;
    typedef typename allocator_type::const_pointer   const_pointer;
    typedef typename allocator_type::pointer         iterator;
    typedef typename allocator_type::const_pointer   const_iterator;
    typedef ::std::reverse_iterator<iterator>        reverse_iterator;
    typedef ::std::reverse_iterator<const_iterator>  const_reverse_iterator;
    /** @} */

    enum {
        /**
         * Minimum \c blocksize that can be accommodated by the pool.
         * Smaller sizes will be implicitly rounded up to this value.
         */
        minimum_blocksize = (sizeof(fl_entry) - 1)/sizeof(T) + 1
    };

    /**
     * Default construct a pool with a zero-sized blocks.
     * One must call blocksize(size_type) before using the pool.
     */
    coalescing_pool() : blocksize_(0) {}

    /**
     * Construct a pool handling <tt>blocksize</tt>-sized blocks where each
     * block is at least <tt>blocksize * sizeof(T)</tt> bytes.
     *
     * @param blocksize       Size of each individual block.
     * @param initial_nblocks Number of blocks to anticipate needing.
     */
    explicit coalescing_pool(size_type blocksize,
                             size_type initial_nblocks = 0)
        : blocksize_(::std::max((size_type) minimum_blocksize, blocksize)) {
        if (initial_nblocks) anticipate(initial_nblocks);
    }

    /**
     * Destroy the pool instance.
     * All blocks acquired from the pool will be implicitly released.
     */
    ~coalescing_pool() {}

    /**
     * Memory-related operations
     * @{
     */

    /**
     * Anticipate a request for at most \c nblocks contiguous blocks.  If the
     * pool could not immediately accommodate a call to
     * <tt>acquire(nblocks)</tt>, more memory will be allocated.  Users may
     * greatly improve the space and time efficiency by using this method.
     *
     * An out-of-memory condition will either result in a std::bad_malloc if
     * ArenaAllocator throws std::bad_malloc or a \c blocks instance that
     * evaluates to false in a boolean context.
     *
     * @param nblocks Number of blocks to anticipate needing.
     */
    void anticipate(size_type nblocks);

    /**
     * Acquire \c nblocks contiguous blocks from the pool.  Each individual
     * block will be of size blocksize().  The blocks may (and often should)
     * subsequently be released back to the pool using release();
     *
     * An out-of-memory condition will either result in a std::bad_malloc if
     * ArenaAllocator throws std::bad_malloc or a \c blocks instance that
     * evaluates to false in a boolean context.
     *
     * @param nblocks Number of contiguous blocks requested.
     *
     * @return The newly acquired \c blocks.
     *         The contents of the \c blocks are undefined.
     */
    blocks acquire(size_type nblocks);

    /**
     * Return previously allocated blocks to the pool.
     *
     * @param b Blocks to be returned.
     *
     * @throws std::invalid_argument if \c b did not come from this pool.
     */
    void release(blocks& b);

    /**
     * Drain all unused memory from the pool.  If
     * <tt>ArenaAllocator::deallocate</tt> supports it, the memory will be
     * returned to the operating system.
     *
     * @return The number of unused blocks drained from the pool.
     */
    size_type drain();

    /** @} */

    /**
     * Status queries
     * @{
     */

    /**
     * Obtain the block size currently used by the pool instance.
     * It will be at least \c minimum_blocksize.
     *
     * @return the block size.
     */
    size_type blocksize() const { return blocksize_; };

    /**
     * Specify a new block size for the pool instance.  The pool must not have
     * any blocks in use.
     *
     * @param newsize New block size where each block will be at least
     *                <tt>newsize * sizeof(T)</tt> bytes.
     *
     * @return the new block size which will be at least \c minimum_blocksize.
     * @throw std::logic_error if the pool is in use.
     */
    size_type blocksize(size_type newsize);

    /**
     * Are any blocks in the pool currently in use?
     *
     * @return True if any acquired blocks have not been released.
     */
    bool any_blocks_used();

    /**
     * Query the number of blocks from the pool that are in use.
     *
     * @return the number of blocks in use.
     */
    size_type nblocks();

    /**
     * Query the number of used blocks, free blocks, and the largest contiguous
     * number of blocks available in the pool.
     *
     * @param used   Number of blocks in use.  This value is also returned.
     * @param free   Number of free blocks in the pool.
     * @param contig Largest contiguous number of blocks available without
     *               forcing a memory allocation.
     *
     * @return the number of blocks in use.
     */
    template< typename U, typename V, typename W>
    size_type nblocks(U& used, V& free, W& contig);

    /**
     * Query the number of distinct memory arenas that have been obtained
     * from the operating system via ArenaAllocator.
     *
     * @return The number of memory arenas in use.
     */
    size_type narenas() { return arenas.size(); }

    /** @} */

    /**
     * A contiguous sequence of blocks acquired from a coalescing_pool.  Each
     * instance logically represents the memory spanned by those contiguous
     * blocks.
     */
    class blocks
    {
    public:

        /** Produces an empty, unused block. */
        blocks()                       : b_(0), e_(0)   {}

        /** Produces a block spanning the range <tt>[b,e)</tt>. */
        blocks(pointer b, pointer e)   : b_(b), e_(e)   {}

        /** Produces a block spanning the range <tt>[p,p+n)</tt>. */
        blocks(pointer p, size_type n) : b_(p), e_(p+n) {}

        /** Does this instance wholly contain the given block? */
        bool contains(const blocks& o) const {
            return    ::std::distance(b_, o.b_) >= 0
                   && ::std::distance(o.e_, e_) >= 0;
        }

        /** How many elements of type \c T are contained in this block? */
        size_type nelems() const { return ::std::distance(b_, e_); }

        /** Set the memory contained by this instance to zero. */
        void zero() { ::std::memset(b_, 0, nelems()*sizeof(T)); }

        /** Is this instance equivalent to another? */
        bool operator==(const blocks& o) const {
            return b_ == o.b_ && e_ == o.e_;
        }

        /** Is this instance not equivalent to another? */
        bool operator!=(const blocks& o) const {
            return !(*this == o);
        }

        /**
         * Causes a valid instance to evaluate to \c true in a \c bool context.
         */
        operator void *() { return b_; }

        /**
         * Causes a valid instance to evaluate to \c true in a \c bool context.
         */
        operator const void *() const { return b_; }

        /** @{ */
        iterator               begin()        { return b_; }
        const_iterator         begin()  const { return b_; }
        iterator               end()          { return e_; }
        const_iterator         end()    const { return e_; }
        reverse_iterator       rbegin()       { return reverse_iterator(e_); }
        const_reverse_iterator rbegin() const { return reverse_iterator(e_); }
        reverse_iterator       rend()         { return reverse_iterator(b_); }
        const_reverse_iterator rend()   const { return reverse_iterator(b_); }
        /** @} */

        /** Output the range spanned by this instance. */
        template< typename CharT, class Traits >
        friend ::std::basic_ostream<CharT,Traits>& operator<<(
            ::std::basic_ostream<CharT,Traits>& os, const blocks& b) {
            return os << '[' << b.b_ << ',' << b.e_ << ')';
        }

    protected:

        friend class coalescing_pool;

        /** Beginning of the spanned memory range. */
        pointer b_;

        /** The end of the spanned memory range. */
        pointer e_;
    };

    /**
     * A "scoped" blocks instance to facilitate the resource acquisition is
     * initialization pattern.  Blocks are acquired from a coalescing_pool at
     * construction time and are automatically released at destruction.
     */
    class scoped_blocks : public ::boost::noncopyable, private blocks
    {
    public:

        /** Produces an empty, unused scoped_block. */
        explicit scoped_blocks(coalescing_pool &p)
            : blocks(), p_(&p) {}

        /**
         * Acquire \c nblock contiguous blocks from pool \c p.
         *
         * An out-of-memory condition will either result in a std::bad_malloc
         * if ArenaAllocator throws std::bad_malloc or a \c blocks instance
         * that evaluates to false in a boolean context.
         *
         * @param p       Pool from which to request memory.
         * @param nblocks Number of contiguous blocks requested.
         */
        scoped_blocks(coalescing_pool &p, size_type nblocks)
            : blocks(p.acquire(nblocks)), p_(&p) {}

        /**
         * Destroy the instance.
         * All blocks acquired from the pool will be implicitly released.
         */
        ~scoped_blocks() { release(); }

        /**
         * Explicitly release any memory held by this instance.
         */
        void release() { p_->release(*this); }

        /**
         * Release any currently held memory back to the pool and acquire \c
         * nblocks.  Beware that \c nblocks refers to the blocksize()
         * <i>currently in use</i> by the pool provided at construction time.
         *
         * @param nblocks Number of contiguous blocks requested.
         */
        void reset(size_type nblocks) {
            release();
            *reinterpret_cast<blocks*>(this) = p_->acquire(nblocks);
        }

        /** Swap the contents of this instance with \c o. */
        void swap(scoped_blocks &o) {
            using ::std::swap;
            swap(this->b_, o.b_);
            swap(this->e_, o.e_);
            swap(this->p_, o.p_);
        }

        /** @copydoc blocks::contains */
        using blocks::contains;

        /** @copydoc blocks::nelems */
        using blocks::nelems;

        /** @copydoc blocks::zero */
        using blocks::zero;

        /** Obtain pool \c blocksize pool associated with this instance. */
        size_type blocksize() const { return p_->blocksize(); };

        /** @copydoc blocks::operator== */
        bool operator==(const scoped_blocks& o) const {
            return p_ == o.p_ &&   *reinterpret_cast<const blocks*>(this)
                                 == reinterpret_cast<const blocks&>(o);
        }

        /** @copydoc blocks::operator!= */
        bool operator!=(const scoped_blocks& o) const {
            return !(*this == o);
        }

        /** @copydoc blocks::operator void* */
        using blocks::operator void*;

        /** @{ */
        using blocks::begin;
        using blocks::end;
        using blocks::rbegin;
        using blocks::rend;
        /** @} */

        /** Output the range spanned by this instance. */
        template< typename CharT, class Traits >
        friend ::std::basic_ostream<CharT,Traits>& operator<<(
            ::std::basic_ostream<CharT,Traits>& os, const scoped_blocks& b) {
            return os << reinterpret_cast<const blocks&>(b);
        }

    protected:

        /** The pool associated with this instance. */
        coalescing_pool *p_;
    };

private:

    size_type blocksize_;

    struct compare_locations
        : public ::std::binary_function<const blocks&, const blocks&, bool>
    {
        bool operator() (const blocks& l, const blocks& r) const {
            return l.b_ != r.b_ ? l.b_ < r.b_ : l.e_ < r.e_;
        }

        bool operator() (const blocks& l, const_pointer p) const {
            return l.b_ < p;
        }

        bool operator() (const_pointer p, const blocks& r) const {
            return p < r.b_;
        }
    };

    struct compare_extents
        : public ::std::binary_function<const blocks&, const blocks&, bool>
    {
        bool operator() (const blocks& l, const blocks& r) const {
            const size_type ln = l.nelems(), rn = r.nelems();
            return ln != rn ? ln < rn : l.b_ < r.b_;
        }

        bool operator() (const blocks& l, size_type nelems) const {
            return l.e_ < l.b_ + nelems;
        }

        bool operator() (size_type nelems, const blocks& r) const {
            return r.b_ + nelems < r.e_;
        }
    };

    class fl_entry
        : public blocks,
          public ::boost::intrusive::set_base_hook<
                ::boost::intrusive::link_mode<Policy::link_mode_type>
          >
    {
    public:
        fl_entry(pointer b, pointer e)   : blocks(b, e) {}
        fl_entry(pointer p, size_type n) : blocks(p, n) {}
    };

    typedef typename ::boost::intrusive::set<
            fl_entry, ::boost::intrusive::compare<compare_locations>
        > fl_type;

    class arena : public  blocks,
                  public  ::boost::noncopyable,
                  private ArenaAllocator
    {
    public:
        arena(size_type nelems, const_pointer hint = 0)
            : blocks(this->allocate(nelems, hint), nelems) {
            if (Policy::zero_new_arenas) this->zero();
            fl.insert(*(new (this->b_) fl_entry(this->b_, this->e_)));
        }

        ~arena() { fl.clear(); this->deallocate(this->b_, this->nelems()); }

        fl_type fl;
    };

    typedef ::boost::ptr_set<arena,compare_extents> arenas_type;

    arenas_type arenas;

    void anticipate_(size_type nblocks,
                     typename arenas_type::iterator& i,
                     typename fl_type::iterator& j);
};

template< typename T, class ArenaAllocator, class Policy >
typename coalescing_pool<T,ArenaAllocator,Policy>::size_type
coalescing_pool<T,ArenaAllocator,Policy>::blocksize(
        size_type newsize)
{
    if (blocksize_ == newsize) return blocksize_;   // Short circuit on NOP

    if (any_blocks_used()) throw ::std::logic_error(
          ::std::string(__PRETTY_FUNCTION__)
        + " cannot set new block size when pool is actively in use");

    newsize = ::std::max((size_type) minimum_blocksize, newsize);  // Constrain

    // Purge existing arenas unless newsize evenly divides blocksize_
    if (blocksize_ % newsize) arenas.clear();

    return (blocksize_ = newsize);
}

template< typename T, class ArenaAllocator, class Policy >
typename coalescing_pool<T,ArenaAllocator,Policy>::size_type
coalescing_pool<T,ArenaAllocator,Policy>::nblocks()
{
    if (!blocksize_) throw ::std::logic_error(::std::string(__PRETTY_FUNCTION__)
            + " requires blocksize to be previously specified");

    size_type total = 0;
    typename arenas_type::iterator i, ie;
    for (i = arenas.begin(), ie = arenas.end(); i != ie; ++i) {
        total += i->nelems() / blocksize_;
    }
    return total;
}

template< typename T, class ArenaAllocator, class Policy >
template< typename U, typename V, typename W>
typename coalescing_pool<T,ArenaAllocator,Policy>::size_type
coalescing_pool<T,ArenaAllocator,Policy>::nblocks(
        U& used, V& free, W& contig)
{
    if (!blocksize_) throw ::std::logic_error(::std::string(__PRETTY_FUNCTION__)
            + " requires blocksize to be previously specified");

    size_type total = 0;
    free = contig = 0;

    typename arenas_type::iterator i, ie;
    for (i = arenas.begin(), ie = arenas.end(); i != ie; ++i) {
        total += i->nelems() / blocksize_;
        typename fl_type::iterator j, je;
        for (j = i->fl.begin(), je = i->fl.end(); j != je; ++j) {
            const size_type n = j->nelems() / blocksize_;
            free += n;
            contig = ::std::max(contig, static_cast<W>(n));
        }
    }

    used = total - free;
    return total;
}

template< typename T, class ArenaAllocator, class Policy >
bool coalescing_pool<T,ArenaAllocator,Policy>::any_blocks_used()
{
    typename arenas_type::iterator i, ie;
    for (i = arenas.begin(), ie = arenas.end(); i != ie; ++i) {
        if (i->fl.empty() || (*i->fl.begin()) != *i)
            return true;
    }

    return false;
}

template< typename T, class ArenaAllocator, class Policy >
typename coalescing_pool<T,ArenaAllocator,Policy>::size_type
coalescing_pool<T,ArenaAllocator,Policy>::drain()
{
    size_type n = 0;

    typename arenas_type::iterator i, ie;
    for (i = arenas.begin(), ie = arenas.end(); i != ie;) {
        if (i->fl.empty() || (*i->fl.begin()) != *i) {
            ++i;
        } else {
            n += i->nelems() / blocksize_;
            arenas.erase(i++);
        }
    }

    return n;
}

template< typename T, class ArenaAllocator, class Policy >
void coalescing_pool<T,ArenaAllocator,Policy>::anticipate(
        coalescing_pool::size_type nblocks)
{
    typename coalescing_pool::arenas_type::iterator i;
    typename coalescing_pool::fl_type::iterator j;
    return anticipate_(nblocks, i, j);
}

template< typename T, class ArenaAllocator, class Policy >
void coalescing_pool<T,ArenaAllocator,Policy>::anticipate_(
        coalescing_pool::size_type nblocks,
        typename coalescing_pool::arenas_type::iterator& i,
        typename coalescing_pool::fl_type::iterator& j)
{
    if (!blocksize_) throw ::std::logic_error(::std::string(__PRETTY_FUNCTION__)
            + " requires blocksize to be previously specified");

    if (!nblocks) return;  // Short circuit on degenerate request

    const size_type nelems = nblocks * blocksize_;

    // Search for first arena containing nblocks contiguous, unused blocks
    const typename arenas_type::iterator ie  = arenas.end();
    const typename arenas_type::iterator ilb = ::std::lower_bound(
            arenas.begin(), ie, nelems, compare_extents());
    for (i = ilb; i != ie; ++i) {
        const typename fl_type::iterator je = i->fl.end();
        for (j = i->fl.begin(); j != je; ++j) {
            if (j->e_ >= j->b_ + nelems)
                return;
        }
    }

    // No existing arena satisfies request; coalesce while allocating new arena
    i = arenas.insert(ilb, new arena(::std::max(blocksize_*drain(), nelems)));
    j = i->fl.begin();
}

template< typename T, class ArenaAllocator, class Policy >
typename coalescing_pool<T,ArenaAllocator,Policy>::blocks
coalescing_pool<T,ArenaAllocator,Policy>::acquire(
        coalescing_pool::size_type nblocks)
{
    if (!blocksize_) throw ::std::logic_error(::std::string(__PRETTY_FUNCTION__)
            + " requires blocksize to be previously specified");

    if (!nblocks) return blocks();  // Short circuit on degenerate request

    // Find or allocate the arena/free list from which we'll acquire nblocks
    typename coalescing_pool::arenas_type::iterator i;
    typename coalescing_pool::fl_type::iterator j;
    anticipate_(nblocks, /* mutated */ i, /* mutated */ j);
    fl_type &fl = i->fl;

    // Hold out for a perfect fit in arena while tracking a best fit candidate
    const size_type nelems = nblocks * blocksize_;
    const typename fl_type::iterator je = fl.end();
    typename fl_type::iterator b = je;
    for (; j != je; ++j) {
        const size_type jn = j->nelems();
        if (jn == nelems) {
            fl.erase(j);
            return *j;
        } else if (jn > nelems && (b == je || jn < (*b).nelems())) {
            b = j;
        }
    }

    if (b == je) throw ::std::runtime_error(::std::string(__PRETTY_FUNCTION__)
            + " sanity failure: anticipate_ did not return a usable arena");

    // Split head off *b, record tail in fl, and used head to satisfy request
    blocks head((*b).b_, nelems);
    fl.insert(b, *(new (head.e_) fl_entry(head.e_, (*b).nelems() - nelems)));
    fl.erase(b);
    return head;
}

template< typename T, class ArenaAllocator, class Policy >
void coalescing_pool<T,ArenaAllocator,Policy>::release(
        coalescing_pool::blocks& b)
{
    using ::std::invalid_argument;
    using ::std::string;

    if (!b) return;  // Short circuit on release of inactive block

    const typename arenas_type::iterator end = arenas.end();
    typename arenas_type::iterator i = ::std::lower_bound(
            arenas.begin(), end, b.nelems(), compare_extents());
    for (; i != end; ++i) {
        if (i->contains(b)) {

            // Prepare the fl_entry header in block b and add to arena fl
            // Double release detection necessary since copies of b may exist
            ::std::pair<typename fl_type::iterator,bool> j
                = i->fl.insert(*(new (b.b_) fl_entry(b.b_, b.e_)));
            if (!j.second) throw invalid_argument(string(__PRETTY_FUNCTION__)
                    + " detected invalid argument or double release");

            // Mark b as having been released
            b.b_ = b.e_ = 0;

            // Coalescing requires only one lookahead or lookbehind
            // because we perform it on every release invocation

            // Possibly coalesce the next fl entry into j
            typename fl_type::iterator n = ::boost::next(j.first);
            if (n != i->fl.end() && j.first->e_ == n->b_) {
                j.first->e_ = n->e_;
                i->fl.erase(n);
            }

            // Possibly coalesce j into the prior fl entry
            typename fl_type::iterator p = ::boost::prior(j.first); // Invalid?
            if (j.first != i->fl.begin() && p->e_ == j.first->b_) { // Valid!
                p->e_ = j.first->e_;
                i->fl.erase(j.first);
            }

            return;
        }
    }

    throw invalid_argument(string(__PRETTY_FUNCTION__)
            + " called on blocks from a different pool");
}

/** A free swap function delegating to the member swap function. */
template< typename T, class ArenaAllocator, class Policy >
void swap(
        typename coalescing_pool<T,ArenaAllocator,Policy>::scoped_blocks& a,
        typename coalescing_pool<T,ArenaAllocator,Policy>::scoped_blocks& b)
{
    a.swap(b);
}


} // namespace suzerain

#endif /* SUZERAIN_COALESCING_POOL_HPP */
