//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2010 The PECOS Development Team
//
// Please see http://pecos.ices.utexas.edu for more information.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------
//
// fftw.cpp: miscellaneous utilities for working with FFTW's C interface
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/fftw.hpp>
#include <suzerain/validation.hpp>

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
            case 'w': return wisdom_only;
        }
    }

    return measure;
}

rigor rigor_from(const unsigned flags)
{
    if (flags & FFTW_WISDOM_ONLY) {
        return wisdom_only;
    } else if (flags & FFTW_EXHAUSTIVE) {
        return exhaustive;
    } else if (flags & FFTW_PATIENT) {
        return patient;
    } else if (flags & FFTW_ESTIMATE) {
        return estimate;
    } else {
        return measure; // Default
    }
}

const char * c_str(const rigor r)
{
    switch (r) {
        case estimate:    return "estimate";
        case patient:     return "patient";
        case exhaustive:  return "exhaustive";
        case wisdom_only: return "wisdom_only";
        default:
        case measure:     return "measure";    // Default
    }
}

int default_nthreads()
{
    int retval = 1;

#ifdef HAVE_FFTW3_THREADS
#if defined HAVE_OPENMP || defined HAVE_PTHREAD
    // Use OMP_NUM_THREADS for both OpenMP and pthreads
    const char * const envstr = getenv("OMP_NUM_THREADS");
    if (envstr) {
        const int envnum = atoi(envstr);
        if (envnum > 0) retval = envnum;
    }
#else
# error "Sanity check failed; unknown FFTW threading model in use."
#endif
#endif

    return retval;
}

void FFTWDefinition::normalize_rigor_fft(std::string input)
{
    this->rigor_fft_ = c_str(rigor_from(input.c_str()));
}

void FFTWDefinition::normalize_rigor_mpi(std::string input)
{
    this->rigor_mpi_ = c_str(rigor_from(input.c_str()));
}

FFTWDefinition::FFTWDefinition(
        const rigor rigor_fft,
        const rigor rigor_mpi)
    : IDefinition("FFTW definition"),
      rigor_fft_(c_str(rigor_fft)),
      rigor_mpi_(c_str(rigor_mpi)),
      nthreads_(default_nthreads())
{
    namespace po = ::boost::program_options;
    using ::std::bind2nd;
    using ::std::ptr_fun;
    using ::suzerain::validation::ensure_nonnegative;

    std::string rigor_options;
    rigor_options += " {";
    rigor_options += c_str(estimate);
    rigor_options += ", ";
    rigor_options += c_str(measure);
    rigor_options += ", ";
    rigor_options += c_str(patient);
    rigor_options += ", ";
    rigor_options += c_str(exhaustive);
    rigor_options += ", or ";
    rigor_options += c_str(wisdom_only);
    rigor_options += "}";

    std::string rigor_fft_description;
    rigor_fft_description += "FFTW FFT planning rigor ";
    rigor_fft_description += rigor_options;

    std::string rigor_mpi_description;
    rigor_mpi_description += "FFTW MPI planning rigor ";
    rigor_mpi_description += rigor_options;

    std::string nthreads_description("Number of FFTW threads to use");
#ifdef HAVE_FFTW3_THREADS
#if defined HAVE_OPENMP
    nthreads_description += " (OpenMP per OMP_NUM_THREADS)";
#elif defined HAVE_PTHREAD
    nthreads_description += " (pthread per OMP_NUM_THREADS)";
#else
#error "Sanity check failed; unknown FFTW threading model in use."
#endif
#else  // HAVE_FFTW3_THREADS not defined
    nthreads_description += " (Disabled)";
#endif // HAVE_FFTW3_THREADS

    this->add_options()
        ("rigor_fft", po::value(&rigor_fft_)
            ->notifier(
                std::bind1st(
                    std::mem_fun(&FFTWDefinition::normalize_rigor_fft),
                    this)
                )
            ->default_value(rigor_fft_),
         rigor_fft_description.c_str())
        ("rigor_mpi", po::value(&rigor_mpi_)
            ->notifier(
                std::bind1st(
                    std::mem_fun(&FFTWDefinition::normalize_rigor_mpi),
                    this)
                )
            ->default_value(rigor_mpi_),
         rigor_mpi_description.c_str())
        ("nthreads", po::value(&nthreads_)
                ->default_value(nthreads_),
         nthreads_description.c_str())
        ("wisdom", po::value(&wisdom_),
         "File used for accumulating FFTW wisdom")
        ("timelimit", po::value(&timelimit_)
                ->default_value(FFTW_NO_TIMELIMIT, "unlimited")
                ->notifier(bind2nd(ptr_fun(ensure_nonnegative<double>),
                                           "timelimit")),
         "Maximum time allowed for preparing any plan")
    ;
}

} // namespace fftw

} // namespace suzerain
