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

/** @file
 * @copydoc mpi_datatype.cpp
 */

#include <suzerain/mpi_datatype.hpp>

namespace suzerain {

namespace mpi {

// Boilerplate common to each template specialization
#define SPECIALIZE(T,C)                         \
    const MPI_Datatype datatype< T >::value = C

// Template specializations translating from C++ to MPI_Datatypes
SPECIALIZE(char,                      MPI_CHAR);
SPECIALIZE(signed char,               MPI_SIGNED_CHAR);
SPECIALIZE(unsigned char,             MPI_UNSIGNED_CHAR);
SPECIALIZE(wchar_t,                   MPI_WCHAR);
SPECIALIZE(signed short int,          MPI_SHORT);
SPECIALIZE(signed int,                MPI_INT);
SPECIALIZE(signed long int,           MPI_LONG);
SPECIALIZE(signed long long int,      MPI_LONG_LONG);
SPECIALIZE(unsigned short int,        MPI_UNSIGNED_SHORT);
SPECIALIZE(unsigned int,              MPI_UNSIGNED);
SPECIALIZE(unsigned long int,         MPI_UNSIGNED_LONG);
SPECIALIZE(unsigned long long int,    MPI_UNSIGNED_LONG_LONG);
SPECIALIZE(float,                     MPI_FLOAT);
SPECIALIZE(double,                    MPI_DOUBLE);
SPECIALIZE(long double,               MPI_LONG_DOUBLE);
SPECIALIZE(std::complex<float>,       MPI_C_FLOAT_COMPLEX);
SPECIALIZE(std::complex<double>,      MPI_C_DOUBLE_COMPLEX);
SPECIALIZE(std::complex<long double>, MPI_C_LONG_DOUBLE_COMPLEX);

#undef SPECIALIZE

} // namespace mpi

} // namespace suzerain
