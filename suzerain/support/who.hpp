//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#ifndef SUZERAIN_SUPPORT_WHO_HPP
#define SUZERAIN_SUPPORT_WHO_HPP

/** @file
 * The <tt>Who<T></tt> template to track a current name for logging purposes.
 */

#include <string>

namespace suzerain {

namespace support {

/**
 * Who tracks a current name for logging purposes.  Subclasses may wish to
 * privately inherit from it to pervasively provide a mutable logger name for
 * use by utilities within \ref logger.hpp.  Non-virtual inheritance is useful
 * when each subclass should have a separately usable name.
 *
 * @tparam T A marker type to prevent ambiguous inheritance issues when
 *           multiple types within a hierarchy privately inherit from Who.
 */
template <class T>
struct Who
{
public:

    /**
     * Construct an instance.
     *
     * @param who Initial name to be used.
     */
    Who(const std::string& who) : who(who) {}

    /** The currently specified name. */
    std::string who;
};

} // end namespace support

} // end namespace suzerain

#endif // SUZERAIN_SUPPORT_WHO_HPP
