/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * fftw.cpp: miscellaneous utilities for working with FFTW's C interface
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <suzerain/fftw.hpp>
#include <cstdlib>

namespace suzerain {

namespace fftw {

rigor rigor_from(const char *name)
{
    if (name) {
        switch (tolower(name[0])) {
            default:  break;
            case 'm': return measure;
            case 'p': return patient;
            case 'e':
                switch (tolower(name[1])) {
                    default:  // Assume exhaustive
                    case 's': return estimate;
                    case 'x': return exhaustive;
                }
        }
    }

    return measure;
}

rigor rigor_from(unsigned flags)
{
    if (flags & FFTW_EXHAUSTIVE) {
        return exhaustive;
    } else if (flags & FFTW_PATIENT) {
        return patient;
    } else if (flags & FFTW_ESTIMATE) {
        return estimate;
    } else {
        return measure; // Default
    }
}

const char * c_str(rigor r)
{
    switch (r) {
        case estimate:   return "estimate";
        case patient:    return "patient";
        case exhaustive: return "exhaustive";
        default:
        case measure:    return "measure";    // Default
    }
}

void FFTWDefinition::normalize_rigor_string(std::string input)
{
    this->rigor_string_ = c_str(rigor_from(input.c_str()));
}

FFTWDefinition::FFTWDefinition()
    : options_("FFTW definition")
{
    namespace po = ::boost::program_options;

    std::string rigor_description;
    rigor_description += "Planning rigor; one of {";
    rigor_description += c_str(estimate);
    rigor_description += ", ";
    rigor_description += c_str(measure);
    rigor_description += ", ";
    rigor_description += c_str(patient);
    rigor_description += ", ";
    rigor_description += c_str(exhaustive);
    rigor_description += "}";

    std::string rigor_default(c_str(measure));

    std::string nthreads_description(
            "Number of threads to use for FFTW planning");
#ifdef HAVE_FFTW3_THREADS
#if defined HAVE_OPENMP
    nthreads_description += " (OpenMP per OMP_NUM_THREADS)";
    const int nthreads_default =
        (getenv("OMP_NUM_THREADS") && atoi(getenv("OMP_NUM_THREADS")) > 0)
        ? atoi(getenv("OMP_NUM_THREADS")) : 1;
#elif defined HAVE_PTHREAD
    // TODO Provide sane nthreads default for FFTW pthread environment
    nthreads_description += " (pthread)";
    const int nthreads_default = 1;
#else
#error "Sanity check failed; unknown FFTW threading model in use."
#endif
#else  /* HAVE_FFTW3_THREADS not defined */
    nthreads_description += " (Disabled)";
    const int nthreads_default = 1;
#endif /* HAVE_FFTW3_THREADS */

    this->nthreads_ = nthreads_default;

    options_.add_options()
        ("rigor",
         po::value<std::string>(&rigor_string_)
                ->notifier(
                    std::bind1st(
                        std::mem_fun(&FFTWDefinition::normalize_rigor_string),
                        this))
                ->default_value(rigor_default),
         rigor_description.c_str())
        ("nthreads",
         po::value<int>(&nthreads_)
                ->default_value(nthreads_default),
         nthreads_description.c_str())
    ;
}

} // namespace fftw

} // namespace suzerain
