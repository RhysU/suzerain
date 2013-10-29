//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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
 * @copydoc largo_formulation.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/largo_formulation.hpp>

#include <boost/assign/list_of.hpp>

#include <suzerain/common.hpp>

namespace suzerain {

// Must be constructed before any static instances are constructed
std::map<std::string,const largo_formulation*> largo_formulation::by_name;

// BEGIN Add known Largo formulations here
const largo_formulation largo_formulation::disable(
        0, "disable", false,
        "No slow growth formulation is in use");

const largo_formulation largo_formulation::temporal(
        1, "bl_temporal", true,
        "Original temporal formulation by Topalian et al.",
        boost::assign::list_of("temporal")
            .convert_to_container<std::vector<std::string> >());

// Notice "Unimplemented placeholder" below.  Model does not exist.
const largo_formulation largo_formulation::spatial(
        2, "bl_spatial", false,
        "Unimplemented placeholder for spatial formulation by Topalian et al.",
        boost::assign::list_of("spatial")
            .convert_to_container<std::vector<std::string> >());

const largo_formulation largo_formulation::temporal_tensor_consistent(
        3, "bl_temporal_tensor-consistent", true,
        "Temporal tensor-consistent formulation by Topalian et al.",
        boost::assign::list_of("bl_temporal_tensor_consistent")
                              ("temporal_tensor_consistent")
                              ("temporal_tensor-consistent")
            .convert_to_container<std::vector<std::string> >());

const largo_formulation largo_formulation::spatiotemporal(
        4, "bl_spatiotemporal", false,
        "Spatiotemporal formulation by Topalian et al.",
        boost::assign::list_of("spatiotemporal")
            .convert_to_container<std::vector<std::string> >());

const largo_formulation largo_formulation::temporal_consistent(
        5, "bl_temporal_consistent", true,
        "Temporal consistent formulation with baseflow support "
        "by Topalian et al.",
        boost::assign::list_of("temporal_consistent")
            .convert_to_container<std::vector<std::string> >());

const largo_formulation largo_formulation::spatiotemporal_consistent(
        6, "bl_spatiotemporal_consistent", false,
        "Spatiotemporal consistent formulation with baseflow support "
        "by Topalian et al.",
        boost::assign::list_of("spatiotemporal_consistent")
            .convert_to_container<std::vector<std::string> >());

// END Add known Largo formulations here

largo_formulation::largo_formulation(
        const int   v,
        const char *n,
        const bool  t,
        const char *d)
    : v(v), n(n), t(t), d(d)
{
    register_name(this->n, this);
}

largo_formulation::largo_formulation(
        const int   v,
        const char *n,
        const bool  t,
        const char *d,
        const std::vector<std::string>& misspellings)
    : v(v), n(n), t(t), d(d)
{
    register_name(this->n, this);

    for (std::size_t i = 0; i < misspellings.size(); ++i) {
        register_name(misspellings[i], this);
    }
}

void
largo_formulation::register_name(const std::string& name,
                                 const largo_formulation* instance)
{
    const std::string& trimmed = boost::algorithm::trim_copy(name);
    if (by_name.find(trimmed) != by_name.end()) {
        throw std::logic_error(std::string("Name collision on '")
            + trimmed + "' when registering with largo_formulation::by_name");
    }
    by_name[trimmed] = instance;
}

const largo_formulation&
largo_formulation::lookup(const std::string& name)
{
    using namespace std;
    using namespace boost::algorithm;
    map<string,const largo_formulation*>::const_iterator i = by_name.find(trim_copy(name));
    if (i == by_name.end()) {
        ostringstream oss;
        oss << "Unknown largo_formulation '" << name << "'";
        throw invalid_argument(oss.str());
    } else {
        return *((*i).second);
    }
}

// Map by_name potentially contains multiple spellings of each formulation.
// Use name() member with set behavior to return only official spellings.
std::set<std::string>
largo_formulation::names()
{
    using namespace std;
    set<string> retval;
    map<string,const largo_formulation*>::const_iterator i   = by_name.begin();
    map<string,const largo_formulation*>::const_iterator end = by_name.end();
    while (i != end) {
        retval.insert((*i).second->name());
        ++i;
    }
    return retval;
}

} // namespace suzerain
