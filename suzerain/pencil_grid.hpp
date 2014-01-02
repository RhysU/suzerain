//--------------------------------------------------------------------------
//
// Copyright (C) 2008-2014 Rhys Ulerich
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

#ifndef SUZERAIN_PENCIL_GRID_HPP
#define SUZERAIN_PENCIL_GRID_HPP

/** @file
 * Class to manage distributed data concerns for P3DFFT or underling usage.
 */

#include <suzerain/common.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/timers.h>

#ifdef SUZERAIN_HAVE_UNDERLING
#include <underling/underling.hpp>
#include <underling/underling_fftw.hpp>
#endif

namespace suzerain {

// Forward declarations for different implementations
class pencil_grid_p3dfft;
class pencil_grid_underling;

// Choose default pencil grid implementation based on build-time availability
#if defined(SUZERAIN_HAVE_P3DFFT)
/** Use P3DFFT-based \c pencil_grid implementation as the default */
typedef pencil_grid_p3dfft pencil_grid_default;
#elif defined(SUZERAIN_HAVE_UNDERLING)
/** Use underling-based \c pencil_grid implementation as the default */
typedef pencil_grid_underling pencil_grid_default;
#else
# error "Found neither suzerain-p3dfft nor underling libraries"
#endif

/**
 * An abstract base class for %pencil grid implementations atop various
 * communications libraries.
 */
class pencil_grid : public boost::noncopyable
{
public:
#ifndef SUZERAIN_PARSED_BY_DOXYGEN
    // See http://eigen.tuxfamily.org/dox/TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#endif

    /** Default constructor for abstract base class */
    pencil_grid() : rank_zero_zero_modes(/* deliberately invalid */-1) {}

    /** Virtual destructor for an abstract base class */
    virtual ~pencil_grid() {};

    /** Human-readable name of concrete implementations */
    virtual const char * implementation() const = 0;

    /** Global grid extents using physical space sizes. */
    Array3i global_physical_extent;

    /** Global grid extents using wave space sizes. */
    Array3i global_wave_extent;

    /**
     * The processor grid extents.
     *
     * In physical space, \f$ P_0 \f$ is the grid extent in the Z direction
     * and \f$ P_1 \f$ is the grid extent in the Y direction.
     * In wave space, \f$ P_0 \f$ is the grid extent in the X direction
     * and \f$ P_1 \f$  is the grid extent in the Z direction.
     *
     * @return the processor grid extents in the \f$ P_0 \f$
     *         and \f$ P_1 \f$ directions as indices 0 and 1, respectively.
     */
    Array2i processor_grid;

    /**
     * Local pencil physical space starting indices (inclusive) within
     * the global extents.
     */
    Array3i local_physical_start;

    /**
     * Local pencil physical space ending indices (exclusive)
     * within the global extents.
     */
    Array3i local_physical_end;

    /**
     * Local pencil physical space extents.
     */
    Array3i local_physical_extent;

    /**
     * Retrieve the number of contiguous real scalars required to store a
     * pencil's worth of data.  This accounts for any padding required due to
     * differences in the local physical and wave space extents, any padding
     * required by the transformation routines, and the fact that wave storage
     * requires complex scalars.
     *
     * @return Number of real-valued scalars (i.e. <tt>suzerain::real_t</tt>s)
     *         required to store one pencil's contiguous data.
     */
    virtual std::size_t local_physical_storage() const = 0;

    /**
     * Local pencil wave space starting indices (inclusive) within
     * the global extents.
     */
    Array3i local_wave_start;

    /**
     * Local pencil wave space ending indices (exclusive)
     * within the global extents.
     */
    Array3i local_wave_end;

    /**
     * Local pencil wave space extents.
     */
    Array3i local_wave_extent;

    /**
     * Which single rank within MPI_COMM_WORLD reports has_zero_zero_modes()?
     *
     * \warning MPI_COMM_WORLD rank zero will not always report that it
     * has_zero_zero_modes()!
     */
    int rank_zero_zero_modes;

    /**
     * The inverse of the product of the global physical extents in X and Z.
     *
     * This factor shows up repeatedly as a normalization constant.
     * See the model document's "Spatial discretization" section for
     * why and it's "Combined space and time discretization" for how
     * this factor \f$\chi\f$ enters into, e.g., timestepping.
     */
    real_t chi() const
    {
        return static_cast<real_t>(1) / (   global_physical_extent.x()
                                          * global_physical_extent.z());
    }

    /**
     * Does this rank possess the "zero-zero" Fourier modes?
     *
     * More concretely, is <tt>local_wave_start.x() == 0 &&
     * local_wave_start.z() == 0 && local_wave_extent.prod() > 0</tt> true?
     * The last condition is required to detect degenerate conditions where
     * some ranks may not contain any wave space data.
     */
    bool has_zero_zero_modes() const
    {
        return    local_wave_start.x() == 0     // Zero mode in X?
               && local_wave_start.z() == 0     // Zero mode in Z?
               && local_wave_extent.prod() > 0; // Nontrivial wave data?
    }

    /**
     * Retrieve the number of contiguous complex scalars required to store a
     * pencil's worth of data.  This accounts for any padding required due to
     * differences in the local physical and wave space extents, any padding
     * required by the transformation routines, and the fact that physical
     * space storage requires real-valued scalars.
     *
     * @return Number of complex-valued scalars (i.e. <tt>suzerain::real_t[2]
     *         </tt>s) required to store one pencil's contiguous data.
     */
    virtual std::size_t local_wave_storage() const = 0;

    /**
     * Collectively transform a field from wave space to physical space.  The
     * input will be complex-valued and stored according to local_wave_start()
     * and local_wave_end().  The output will be real-valued and stored
     * according to local_physical_state() and local_physical_end().
     *
     * @param inout Field to transform in place.
     */
    void transform_wave_to_physical(real_t * inout) const
    {
        SUZERAIN_TIMER_SCOPED("transform_wave_to_physical");

        return  transform_wave_to_physical_(inout);
    }

    /**
     * Collectively transform a field from physical space to wave space.  The
     * input will be complex-valued and stored according to
     * local_physical_start() and local_physical_end().  The output will be
     * real-valued and stored according to local_wave_state() and
     * local_wave_end().
     *
     * @param inout Field to transform in place.
     */
    void transform_physical_to_wave(real_t * inout) const
    {
        SUZERAIN_TIMER_SCOPED("transform_physical_to_wave");

        return transform_physical_to_wave_(inout);
    }

protected:

    /**
     * Compute the value of \ref rank_zero_zero_modes.  Should be used by
     * subclasses after the pencil decomposition is established.
     */
    int compute_rank_zero_zero_modes_() const;

private:

    /** Subclasses override to implement transform_wave_to_physical(). */
    virtual void transform_wave_to_physical_(real_t * inout) const = 0;

    /** Subclasses override to implement transform_physical_to_wave(). */
    virtual void transform_physical_to_wave_(real_t * inout) const = 0;

};

#ifdef SUZERAIN_HAVE_P3DFFT
/**
 * Encapsulates P3DFFT %pencil grid details, including the global grid extent
 * and the processor grid decomposition parameters.  Appropriately handles
 * munging this information to obtain P3DFFT calls which are stride one in Y in
 * wave space.  Unless otherwise noted, all indices start from zero with X, Y,
 * and Z having indices 0, 1, and 2, respectively.  Under the covers, P3DFFT's
 * <tt>p3dfft_setup</tt> and <tt>p3dfft_clean</tt> are invoked on construction
 * and destruction of an instance.
 */
class pencil_grid_p3dfft : public pencil_grid
{
public:

    /**
     * Constructs an instance for the given global physical grid extents.
     *
     * @param global_physical_extent Global physical grid extents in the
     *        streamwise (X), wall-normal (Y), and spanwise (Z)
     *        directions.
     * @param processor_grid The processor grid decomposition to use in the
     *        \f$ P_0 \f$ and \f$ P_1 \f$ directions.  Providing a zero
     *        for either value causes the value to be determined automatically.
     */
    template < typename RandomAccessContainer1,
               typename RandomAccessContainer2 >
    pencil_grid_p3dfft(const RandomAccessContainer1 &global_physical_extent,
                       const RandomAccessContainer2 &processor_grid,
                       unsigned rigor_fft = 0,
                       unsigned rigor_mpi = 0)
    {
        construct_(boost::numeric_cast<int>(global_physical_extent[0]),
                   boost::numeric_cast<int>(global_physical_extent[1]),
                   boost::numeric_cast<int>(global_physical_extent[2]),
                   boost::numeric_cast<int>(processor_grid[0]),
                   boost::numeric_cast<int>(processor_grid[1]),
                   rigor_fft, rigor_mpi);
    }

    virtual ~pencil_grid_p3dfft();

    virtual const char * implementation() const;

    virtual std::size_t local_physical_storage() const;

    virtual std::size_t local_wave_storage() const;

private:

    virtual void transform_wave_to_physical_(real_t * inout) const;

    virtual void transform_physical_to_wave_(real_t * inout) const;

    /** Was p3dfft_setup successfully called during construction? */
    bool p3dfft_setup_called_;

    /** Internal routine performing construction-like tasks */
    void construct_(int Nx, int Ny, int Nz, int Pa, int Pb,
                    unsigned rigor_fft, unsigned rigor_mpi);
};
#endif /* SUZERAIN_HAVE_P3DFFT */

#ifdef SUZERAIN_HAVE_UNDERLING
/**
 * Encapsulates underling %pencil grid details, including the global grid
 * extent and the processor grid decomposition parameters.  Appropriately
 * handles munging this information to obtain underling calls which are stride
 * one in Y in wave space.  Unless otherwise noted, all indices start from zero
 * with X, Y, and Z having indices 0, 1, and 2, respectively.
 */
class pencil_grid_underling : public pencil_grid
{
public:

    /**
     * Constructs an instance for the given global physical grid extents.
     *
     * @param global_physical_extent Global physical grid extents in the
     *        streamwise (X), wall-normal (Y), and spanwise (Z)
     *        directions.
     * @param processor_grid The processor grid decomposition to use in the
     *        \f$ P_0 \f$ and \f$ P_1 \f$ directions.  Providing a zero
     *        for either value causes the value to be determined automatically.
     */
    template < typename RandomAccessContainer1,
               typename RandomAccessContainer2 >
    pencil_grid_underling(const RandomAccessContainer1 &global_physical_extent,
                          const RandomAccessContainer2 &processor_grid,
                          unsigned rigor_fft = 0,
                          unsigned rigor_mpi = 0)
    {
        construct_(boost::numeric_cast<int>(global_physical_extent[0]),
                   boost::numeric_cast<int>(global_physical_extent[1]),
                   boost::numeric_cast<int>(global_physical_extent[2]),
                   boost::numeric_cast<int>(processor_grid[0]),
                   boost::numeric_cast<int>(processor_grid[1]),
                   rigor_fft, rigor_mpi);
    }

    virtual ~pencil_grid_underling();

    virtual const char * implementation() const;

    virtual std::size_t local_physical_storage() const;

    virtual std::size_t local_wave_storage() const;

private:

    virtual void transform_wave_to_physical_(real_t * inout) const;

    virtual void transform_physical_to_wave_(real_t * inout) const;

    /** Internal routine performing construction-like tasks */
    void construct_(int Nx, int Ny, int Nz, int Pa, int Pb,
                    unsigned rigor_fft, unsigned rigor_mpi);

/** @{ */

    scoped_ptr<underling::grid>       grid;
    scoped_ptr<underling::problem>    problem;
    scoped_ptr<underling::plan>       transpose;
    scoped_ptr<underling::fftw::plan> n1_c2c_backward;
    scoped_ptr<underling::fftw::plan> n2_c2r_backward;
    scoped_ptr<underling::fftw::plan> n2_r2c_forward;
    scoped_ptr<underling::fftw::plan> n1_c2c_forward;
    shared_array<underling::real>     buf;

/** @} */

};
#endif /* SUZERAIN_HAVE_UNDERLING */

} // namespace suzerain

#endif // SUZERAIN_PENCIL_GRID_HPP
