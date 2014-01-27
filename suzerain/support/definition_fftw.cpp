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
 * @copydoc definition_fftw.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/support/definition_fftw.hpp>

namespace suzerain {

namespace support {

void definition_fftw::normalize_rigor_fft(std::string input)
{
    this->rigor_fft = fftw::rigor_from(input.c_str());
}

void definition_fftw::normalize_rigor_mpi(std::string input)
{
    this->rigor_mpi = fftw::rigor_from(input.c_str());
}

definition_fftw::definition_fftw(
        const fftw::rigor rigor_fft,
        const fftw::rigor rigor_mpi)
    : rigor_fft(rigor_fft)
    , rigor_mpi(rigor_mpi)
{
}

boost::program_options::options_description
definition_fftw::options_description()
{
    namespace po = boost::program_options;
    using std::bind1st;
    using std::bind2nd;
    using std::mem_fun;
    using std::ptr_fun;
    using std::string;
    using std::ostringstream;

    ostringstream oss;
    oss << " {"    << fftw::c_str(fftw::estimate)
        << ", "    << fftw::c_str(fftw::measure)
        << ", "    << fftw::c_str(fftw::patient)
        << ", "    << fftw::c_str(fftw::exhaustive)
        << ", or " << fftw::c_str(fftw::wisdom_only)
        << "}";
    const string rigor_options = oss.str();

    string rigor_fft_description = "FFTW FFT planning rigor ";
    rigor_fft_description += rigor_options;

    string rigor_mpi_description = "FFTW MPI planning rigor ";
    rigor_mpi_description += rigor_options;

    po::options_description retval("FFTW planning options:");

    retval.add_options()
    ("rigor_fft",
     po::value<std::string>(NULL)
     ->notifier(bind1st(mem_fun(&definition_fftw::normalize_rigor_fft), this))
     ->default_value(fftw::c_str(rigor_fft)),
     rigor_fft_description.c_str())
    ("rigor_mpi",
     po::value<std::string>(NULL)
     ->notifier(bind1st(mem_fun(&definition_fftw::normalize_rigor_mpi), this))
     ->default_value(fftw::c_str(rigor_mpi)),
     rigor_mpi_description.c_str())
    ("plan_wisdom",
     po::value(&plan_wisdom),
     "File for accumulating FFTW planning wisdom")
    ("plan_timelimit",
     po::value(&plan_timelimit)
     ->default_value(FFTW_NO_TIMELIMIT, "unlimited"),
     "Maximum seconds for creating any FFTW plan")
    ;

    return retval;
}

} // end namespace support

} // end namespace suzerain
