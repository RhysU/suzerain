//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
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
// fftw_definition.cpp: A definition_base to configure FFTW from the CLI
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/fftw_definition.hpp>

namespace suzerain {

void fftw_definition::normalize_rigor_fft(std::string input)
{
    this->rigor_fft = fftw::rigor_from(input.c_str());
}

void fftw_definition::normalize_rigor_mpi(std::string input)
{
    this->rigor_mpi = fftw::rigor_from(input.c_str());
}

fftw_definition::fftw_definition(
        const fftw::rigor rigor_fft,
        const fftw::rigor rigor_mpi)
    : definition_base("FFTW planning options:"),
      rigor_fft(rigor_fft),
      rigor_mpi(rigor_mpi)
{
    namespace po = ::boost::program_options;
    using ::std::bind1st;
    using ::std::bind2nd;
    using ::std::mem_fun;
    using ::std::ptr_fun;

    std::string rigor_options;
    rigor_options += " {";
    rigor_options += fftw::c_str(fftw::estimate);
    rigor_options += ", ";
    rigor_options += fftw::c_str(fftw::measure);
    rigor_options += ", ";
    rigor_options += fftw::c_str(fftw::patient);
    rigor_options += ", ";
    rigor_options += fftw::c_str(fftw::exhaustive);
    rigor_options += ", or ";
    rigor_options += fftw::c_str(fftw::wisdom_only);
    rigor_options += "}";

    std::string rigor_fft_description;
    rigor_fft_description += "FFTW FFT planning rigor ";
    rigor_fft_description += rigor_options;

    std::string rigor_mpi_description;
    rigor_mpi_description += "FFTW MPI planning rigor ";
    rigor_mpi_description += rigor_options;

    this->add_options()
        ("rigor_fft",
         po::value<std::string>(NULL)
                ->notifier(bind1st(
                        mem_fun(&fftw_definition::normalize_rigor_fft), this))
                ->default_value(fftw::c_str(rigor_fft)),
         rigor_fft_description.c_str())
        ("rigor_mpi",
          po::value<std::string>(NULL)
                ->notifier(bind1st(
                        mem_fun(&fftw_definition::normalize_rigor_mpi), this))
                ->default_value(fftw::c_str(rigor_mpi)),
         rigor_mpi_description.c_str())
        ("plan_wisdom",
          po::value(&plan_wisdom),
         "File used for accumulating FFTW planning wisdom")
        ("plan_timelimit",
          po::value(&plan_timelimit)
                ->default_value(FFTW_NO_TIMELIMIT, "unlimited"),
         "Maximum number of seconds allowed for creating any single FFTW plan")
    ;
}

} // namespace suzerain
