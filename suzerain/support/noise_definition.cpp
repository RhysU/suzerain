//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

/** @file
 * @copydoc noise_definition.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/support/noise_definition.hpp>

#include <suzerain/exprparse.hpp>
#include <suzerain/validation.hpp>

namespace suzerain {

namespace support {

noise_definition::noise_definition(
            real_t        fluct_percent,
            unsigned long fluct_seed,
            real_t        kxfrac_min,
            real_t        kxfrac_max,
            real_t        kzfrac_min,
            real_t        kzfrac_max)
    : noise_specification(fluct_percent,
                          fluct_seed,
                          kxfrac_min,
                          kxfrac_max,
                          kzfrac_min,
                          kzfrac_max)
{
}

boost::program_options::options_description
noise_definition::options_description()
{
    using boost::bind;
    using boost::program_options::options_description;
    using boost::program_options::value;
    using std::bind2nd;
    using std::pointer_to_binary_function;
    using std::string;
    using validation::ensure_nonnegative;
    using validation::ensure_positive;

    // For brevity below
    pointer_to_binary_function<unsigned long,const char*,void>
            ptr_fun_ensure_positive_ulint(ensure_positive<unsigned long>);
    pointer_to_binary_function<real_t,const char*,void>
            ptr_fun_ensure_nonnegative_real(ensure_nonnegative<real_t>);

    options_description retval(
            "Additive random velocity perturbations on startup");

    retval.add_options()
    ("fluct_percent", value(&percent)->default_value(percent)
        ->notifier(bind2nd(ptr_fun_ensure_nonnegative_real,
                            "fluct_percent")),
        "Maximum fluctuation magnitude to add as a percentage of"
        " centerline mean streamwise velocity")
    ("fluct_kxfrac", value<string>()->default_value("0:1")
        ->notifier(bind(&exprparse_range<const string&,real_t>, _1,
                        &kxfrac_min, &kxfrac_max,
                        0, 1, 0, 1, "fluct_kxfrac")),
        "Range of X wavenumbers in which to generate fluctuations")
    ("fluct_kzfrac", value<string>()->default_value("0:1")
        ->notifier(bind(&exprparse_range<const string&,real_t>, _1,
                        &kzfrac_min, &kzfrac_max,
                        0, 1, 0, 1, "fluct_kzfrac")),
        "Range of Z wavenumbers in which to generate fluctuations")
    ("fluct_seed", value(&seed)->default_value(seed)
        ->notifier(bind2nd(ptr_fun_ensure_positive_ulint,
                            "fluct_seed")),
        "rngstream generator seed (L'Ecuyer et al. 2002)")
    ;

    return retval;
}

} // end namespace support

} // end namespace suzerain
