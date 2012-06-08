/***********************************************************************\
 *
 * File:           RngStream.cpp for multiple streams of Random Numbers
 * Language:       C++ (ISO 1998)
 * Copyright:      Pierre L'Ecuyer, University of Montreal
 * Notice:         This code can be used freely for personal, academic,
 *                 or non-commercial purposes. For commercial purposes,
 *                 please contact P. L'Ecuyer at: lecuyer@iro.umontreal.ca
 * Date:           14 August 2001
 *
 * $Id$
 *
\***********************************************************************/

/** @file
 * Declares <tt>RngStream</tt> by L'Ecuyer <i>et al</i>, 2002.
 */

#ifndef SUZERAIN_RNGSTREAM_HPP
#define SUZERAIN_RNGSTREAM_HPP

#include <iosfwd>
#include <string>

namespace suzerain {

/**
 * <tt>RngStream</tt> by L'Ecuyer, Simard, Chen, and Kelton  as presented in
 * "An object-oriented random-number package with many long streams and
 * substreams", Operations Research (2002), pages 1073-1075.  An <a
 * href="http://www.iro.umontreal.ca/~lecuyer/myftp/streams00/c++/streams4.pdf">
 * extended version of the paper</a> can be found online.
 *
 * The <a
 * href="http://www.iro.umontreal.ca/~lecuyer/myftp/streams00/c++/">original
 * source</a> has been modified in three ways:
 * \li It has been moved into the ::suzerain namespace for encapsulation.
 * \li Doxygen per section 3 of the extended paper has been added.
 * \li The WriteState() and WriteStateFull() methods have been modified
 *     to take a <tt>std::ostream</tt> on which to output information.
 * \li A RandN01() method has been added for drawing from the standard normal.
 *
 * The original authors retain all ownership and copyrights, of course.
 */
class RngStream
{
public:

    /**
     * This constructor creates a new stream with (optional) descriptor \c
     * name.  It initializes its seed \f$I_g\f$ , and sets \f$B_g\f$ and
     * \f$C_g\f$ to \f$I_g\f$ . It also sets its \c anti and \c incPrec
     * switches to \c false. The seed \f$I_g\f$ is equal to the initial
     * seed of the package if this is the first stream created; otherwise it is
     * Z steps ahead of the seed of the most recently created stream.
     */
    RngStream (const char *name = "");

    /**
     * Sets the initial seed \f$s_0\f$ of the package to the six integers in
     * the vector seed. The first 3 integers in the seed must all be less
     * than \f$m_1\f$ = 4294967087, and not all 0; and the last 3 integers must
     * all be less than \f$m_2\f$ = 4294944443, and not all 0. If this method
     * is not called, the default initial seed is (12345, 12345, 12345, 12345,
     * 12345, 12345).
     *
     * @return \c false for invalid seeds, and \c true otherwise.
     */
    static bool SetPackageSeed (const unsigned long seed[6]);

    /**
     * Reinitializes the stream to its initial state:
     * \f$C_g\f$ and \f$B_g\f$ are set to \f$I_g\f$.
     */
    void ResetStartStream ();

    /**
     * Reinitializes the stream to the beginning of its current substream:
     * \f$C_g\f$ is set to \f$B_g\f$.
     */
    void ResetStartSubstream ();

    /**
     * Reinitializes the stream to the beginning of its next substream:
     * \f$N_g\f$ is computed, and \f$C_g\f$ and \f$B_g\f$ are set to \f$N_g\f$.
     */
    void ResetNextSubstream ();

    /**
     * If \c a = \c true, the stream will start generating antithetic variates,
     * i.e., 1 − U instead of U , until this method is called again with
     * \c a = \c false.
     */
    void SetAntithetic (bool a);

    /**
     * After calling this method with \c incp = \c true, each call to the
     * generator (direct or indirect) for this stream will return a uniform
     * random number with more bits of resolution (53 bits if machine
     * follows IEEE 754 standard) instead of 32 bits, and will advance the
     * state of the stream by 2 steps instead of 1. More precisely, if \c s
     * is a stream of the class \c RngStream, in the non- antithetic case,
     * the instruction "<code>u = s.RandU01()</code>" will be equivalent
     * to "<code>u = (s.RandU01() + s.RandU01() * fact) % 1.0</code>"
     * where the constant fact is equal to \f$2^{−24}\f$ . This also
     * applies when calling \c RandU01() indirectly (e.g., via \c RandInt(),
     * etc.). By default, or if this method is called again with \c incp
     * = \c false, each call to \c RandU01() for this stream advances the
     * state by 1 step and returns a number with 32 bits of resolution.
     */
    void IncreasedPrecis (bool incp);

    /*
     * Sets the initial seed \f$I_g\f$ of the stream to the vector
     * \c seed. The vector \c seed should contain valid seed values as
     * described in \c SetPackageSeed(). The state of the stream is then reset
     * to this initial seed. The states and seeds of the other streams are
     * not modified. As a result, after calling this method, the initial
     * seeds of the streams are no longer spaced \f$Z\f$ values apart. We
     * discourage the use of this method; proper use of the Reset* methods
     * is preferable.
     *
     * @return \c false for invalid seeds, and \c true otherwise.
     */
    bool SetSeed (const unsigned long seed[6]);

    /**
     * Advances the state by \f$n\f$ steps (see below for the meaning
     * of \f$n\f$), without modifying the states of other streams or the
     * values of \f$B_g\f$ and \f$I_g\f$ in the current object. If \f$e >
     * 0\f$, then \f$n = 2e + c\f$; if \f$e < 0\f$, then \f$n = −2−e +
     * c\f$; and if \f$e = 0\f$, then \f$n = c\f$. Note: c is allowed to
     * take negative values.  We discourage the use of this method.
     */
    void AdvanceState (long e, long c);

    /**
     * Returns in <tt>seed[0..5]</tt> the current state \f$C_g\f$ of this
     * stream. This is convenient if we want to save the state for
     * subsequent use.
     */
    void GetState (unsigned long seed[6]) const;

    /**
     * Writes (to standard output) the current state \f$C_g\f$ of this stream.
     */
    void WriteState (std::ostream &out = std::cout) const;

    /**
     * Writes (to standard output) the value of all the internal variables of
     * this stream: \c name, \c anti, \c incPrec, \f$I_g\f$, \f$B_g\f$,
     * \f$C_g\f$.
     */
    void WriteStateFull (std::ostream &out = std::cout) const;

    /**
     * Normally, returns a (pseudo)random number from the uniform
     * distribution over the interval (0, 1), after advancing the state
     * by one step. The returned number has 32 bits of precision in the
     * sense that it is always a multiple of 1/(232 − 208). However, if
     * <code>IncreasedPrecis(true)</code> has been called for this stream,
     * the state is advanced by two steps and the returned number has 53
     * bits of precision.
     */
    double RandU01 ();

    /**
     * Returns a (pseudo)random number from the discrete uniform distribution
     * over the integers \f${i, i + 1, . . . , j}\f$.
     * Makes one call to RandU01().
     */
    int RandInt (int i, int j);

    /**
     * Returns a (pseudo)random number from the normal distribution with mean 0
     * and standard deviation 1 after advancing the state identically to a
     * single call to RandU01().  The value is using an inverse CDF technique.
     * The inverse CDF of the standard normal is approximated following <a
     * href="http://home.online.no/~pjacklam/notes/invnorm/"> Peter John
     * Acklam</a> and refined using a single step of Halley's rational method
     * as Acklam describes.
     */
    double RandN01 ();

private:

    /** Stores the current seed \f$C_g\f$. */
    double Cg[6];

    /** Stores the beginning of the current block (substream) \f$B_g\f$. */
    double Bg[6];

    /** Stores the beginning of the current stream \f$I_g\f$. */
    double Ig[6];

    /** Indicates whether to generate antithetic random numbers. */
    bool anti;

    /** Indicates whether to generate increased precision random numbers. */
    bool incPrec;

    /** String to store the optional name of the current RngStream object. */
    std::string name;

    /**
     * Static vector to store the beginning state of the next RngStream
     * to be created (instantiated).
     */
    static double nextSeed[6];

    /** The backbone uniform random number generator. */
    double U01 ();

    /**
     * The backbone uniform random number generator
     * with increased precision.
     */
    double U01d ();

};

}  // end namespace suzerain

#endif  /* SUZERAIN_RNGSTREAM_HPP */
