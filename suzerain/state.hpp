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
 * state.hpp: Class to manage mutable state vectors
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_STATE_HPP
#define __SUZERAIN_STATE_HPP

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al.h>

namespace suzerain
{

template< typename FPT >
class IState
{
public:
    typedef boost::multi_array_types::index index;
    typedef boost::multi_array_types::size_type size_type;
    typedef boost::multi_array_types::difference_type difference_type;

    const size_type variable_count; /**< Number of state variables */
    const size_type vector_length;  /**< Length of each state vector */
    const size_type vector_count;   /**< Number of state vectors */

    IState(size_type variable_count,
           size_type vector_length,
           size_type vector_count)
        :   variable_count(variable_count),
            vector_length(vector_length),
            vector_count(vector_count) {};

    virtual ~IState() {};

    virtual bool isConformant(const IState * const other)
    {
        return    variable_count == other->variable_count
               && vector_length  == other->vector_length
               && vector_count   == other->vector_count;
    }

    virtual void scaleAddScaled(const FPT thisScale,
                                const FPT otherScale,
                                IState<FPT> * const other)
                                throw(std::bad_cast) = 0;
};

template< typename FPT >
class RealState : public IState<FPT>, private boost::noncopyable
{
public:
    typedef typename boost::multi_array_ref<FPT, 3> state_type;

    explicit RealState(typename IState<FPT>::size_type variable_count,
                       typename IState<FPT>::size_type vector_length,
                       typename IState<FPT>::size_type vector_count)
                       throw(std::bad_alloc);

    virtual void scaleAddScaled(const FPT thisScale,
                                const FPT otherScale,
                                IState<FPT> * const other)
                                throw(std::bad_cast);
private:
    boost::shared_array<double> raw;  /**< Raw, aligned storage */

public:
    state_type data;  /**< MultiArray of real-valued state */
};

template< typename FPT >
RealState<FPT>::RealState(typename IState<FPT>::size_type variable_count,
                          typename IState<FPT>::size_type vector_length,
                          typename IState<FPT>::size_type vector_count)
throw(std::bad_alloc)
    : IState<FPT>(variable_count, vector_length, vector_count),
      raw(reinterpret_cast<FPT *>(suzerain_blas_malloc(
                sizeof(FPT)*variable_count*vector_length*vector_count)),
          std::ptr_fun(free)),
      data(raw.get(),
           boost::extents[variable_count][vector_length][vector_count],
           boost::fortran_storage_order())
{
    if (!raw) throw std::bad_alloc();
}

namespace detail {

    template<typename FPT>
    void axpby(const int n,
               const FPT alpha, const FPT *x, const int incx,
               const FPT beta,        FPT *y, const int incy);

    template<>
    void axpby<float>(const int n,
                      const float alpha, const float *x, const int incx,
                      const float beta,        float *y, const int incy)
    {
        return suzerain_blas_saxpby(n, alpha, x, incx, beta, y, incy);
    }

    template<>
    void axpby<double>(const int n,
                       const double alpha, const double *x, const int incx,
                       const double beta,        double *y, const int incy)
    {
        return suzerain_blas_daxpby(n, alpha, x, incx, beta, y, incy);
    }
}

template< typename FPT >
void RealState<FPT>::scaleAddScaled(const FPT thisScale,
                                    const FPT otherScale,
                                    IState<FPT> * const other)
throw(std::bad_cast)
{
}

template< typename FPT >
class ComplexState : public IState<FPT>, private boost::noncopyable
{
public:
    typedef typename boost::multi_array_ref<std::complex<FPT>, 3> state_type;
    typedef typename boost::multi_array_ref<FPT, 3> component_type;

    explicit ComplexState(typename IState<FPT>::size_type variable_count,
                          typename IState<FPT>::size_type vector_length,
                          typename IState<FPT>::size_type vector_count)
    throw(std::bad_alloc);

    virtual void scaleAddScaled(const FPT thisScale,
                                const FPT otherScale,
                                IState<FPT> * const other);
private:
    boost::shared_array<FPT> raw;  /**< Raw, aligned storage */

public:
    state_type data;      /**< MultiArray of complex-valued state */
    component_type real;  /**< MultiArray of state's real part */
    component_type imag;  /**< MultiArray of state's imaginary part */
};

} // namespace suzerain

#endif // __SUZERAIN_STATE_HPP
