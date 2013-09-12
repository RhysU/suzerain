/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2011, 2012, 2013 Rhys Ulerich
 * Copyright (C) 2012, 2013 The PECOS Development Team
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

#ifndef SUZERAIN_COUNTOF_H
#define SUZERAIN_COUNTOF_H

/** @file
 * Macros to count the number of array elements at compile time.  Provides both
 * a type-unsafe C version and a type-safe C++ version.  The C++ version is
 * from Ivan J. Johnson's article "Counting Array Elements at Compile Time"
 * published 6 March 2007 in Dr. Dobb's (http://drdobbs.com/cpp/197800525).
 */

#ifndef __cplusplus

/** Count the number of elements in an array at compile time */
#define SUZERAIN_COUNTOF(x) (sizeof(x)/sizeof((x)[0]))

#else /* __cplusplus */

/** Type safe count the number of elements in an array at compile time */
#define SUZERAIN_COUNTOF(x)  (                                                  \
    0*sizeof(reinterpret_cast<const ::suzerain::BAD_ARGUMENT_TO_COUNTOF*>(x)) + \
    0*sizeof(::suzerain::BAD_ARGUMENT_TO_COUNTOF::check_type((x), &(x))     ) + \
    sizeof(x)/sizeof((x)[0])   )

#ifndef SUZERAIN_PARSED_BY_DOXYGEN
namespace suzerain {

class BAD_ARGUMENT_TO_COUNTOF
{
public:
    class Is_pointer;
    class Is_array {};
    template<typename T>
    static Is_pointer check_type(const T*, const T* const*);
    static Is_array check_type(const void*, const void*);
};

}
#endif /* SUZERAIN_PARSED_BY_DOXYGEN */

#endif /* __cplusplus */

#endif /* SUZERAIN_COUNTOF_H */
