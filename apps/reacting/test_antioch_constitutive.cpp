//--------------------------------------------------------------------------
//
// Copyright (C) 2013 The PECOS Development Team
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

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#ifdef HAVE_ANTIOCH // only makes sense when antioch is available

#include <suzerain/support/support.hpp>
#include <suzerain/support/logging.hpp>

#include "antioch_constitutive.hpp"

/*
 * Initialize minimial test environment: MPI and logging
 */
void initialize_test_env(int argc, char **argv)
{
    MPI_Init(&argc, &argv);                     // Initialize MPI...
    atexit((void (*) ()) MPI_Finalize);         // ...finalize at exit
    suzerain::support::logging::initialize(     // Initialize logging
        MPI_COMM_WORLD,    
        suzerain::support::log4cxx_config_console); // correct???

}

/*
 * Check antioch_constitutive save/load for internal consistency
 */
int test_save_load(MPI_Comm comm)
{
    using suzerain::reacting::antioch_constitutive;
    using suzerain::real_t;

    // Prepare input data
    std::vector<std::string> species_names;
    const unsigned int Ns=2;
    species_names.reserve(Ns);
    species_names.push_back( "N2" );
    species_names.push_back( "O2" );

    // Not a real file... only checking name consistency
    std::string chem_input_file("fake_filename.xml");
        
    const real_t Le = 0.7;
    const real_t alpha = 0.5;

    antioch_constitutive acl1(species_names, chem_input_file, Le, alpha);
    antioch_constitutive acl2;

    
    // Create hdf5 file and save antioch data
    {
        esio_handle h = esio_handle_initialize(comm);
        esio_file_create(h, "test_save_load_antioch.h5", 1 /* overwrite */);
        acl1.save(h);
        esio_file_close(h);
        esio_handle_finalize(h);
    }
    

    // Open hdf5 file and read antioch data
    {
        esio_handle h = esio_handle_initialize(comm);
        esio_file_open(h, "test_save_load_antioch.h5", 0 /* read-only */);
        acl2.load(h);
        esio_file_close(h);
        esio_handle_finalize(h);
    }
    

    // check consistency
    int nerr=0;

    // ...between acl1 and expected (i.e., ctor behavior)
    {
        if (acl1.species_names.size() != species_names.size()) {
            std::cerr << "Error: acl1.species_names.size() = " << acl1.species_names.size()
                      << " where " << Ns << " was expected" << std::endl;
            nerr += 1;
        }
        
        for (unsigned int i=0; i<acl1.species_names.size(); ++i) {
            if( acl1.species_names[i] != species_names[i] ) {
                std::cerr << "Error: acl1.species_names[" << i << "] = " 
                          << acl1.species_names[i] << " but expected"
                          << species_names[i] << std::endl;
                nerr += 1;
            }
        }

        if (acl1.chem_input_file != "fake_filename.xml") {
            std::cerr << "Error: acl1.chem_input_file = " << acl1.chem_input_file
                      << " where fake_filename.xml was expected" << std::endl;
            nerr += 1;
        }
        
        if (acl1.Le != Le) {
            std::cerr << std::setprecision(20);
            std::cerr << "Error: acl1.Le = " << acl1.Le
                      << " where " << Le << " was expected" << std::endl;
            nerr += 1;
        }

        if (acl1.alpha != alpha) {
            std::cerr << std::setprecision(20);
            std::cerr << "Error: acl1.alpha = " << acl1.alpha
                      << " where " << alpha << " was expected" << std::endl;
            nerr += 1;
        }
    }


    // ...between acl2 and acl1 (i.e., save/load consistency)
    {
        if (acl2.species_names.size() != acl1.species_names.size()) {
            std::cerr << "Error: acl2.species_names.size() = " << acl2.species_names.size()
                      << " where " << acl1.species_names.size() << " was expected" << std::endl;
            nerr += 1;
        }
        
        for (unsigned int i=0; i<acl1.species_names.size(); ++i) {
            if( acl2.species_names[i] != acl1.species_names[i] ) {
                std::cerr << "Error: acl2.species_names[" << i << "] = " 
                          << acl2.species_names[i] << " but acl1.species_names["
                          << i << "] = " << acl1.species_names[i] << std::endl;
                nerr += 1;
            }
        }

        if (acl2.chem_input_file != acl1.chem_input_file) {
            std::cerr << "Error: acl2.chem_input_file = " << acl2.chem_input_file
                      << " where " << acl1.chem_input_file << " was expected" << std::endl;
            nerr += 1;
        }

        if (acl2.Le != acl1.Le) {
            std::cerr << std::setprecision(20);
            std::cerr << "Error: acl2.Le = " << acl2.Le
                      << " where " << acl1.Le << " was expected" << std::endl;
            nerr += 1;
        }

        if (acl2.alpha != acl1.alpha) {
            std::cerr << std::setprecision(20);
            std::cerr << "Error: acl2.alpha = " << acl2.alpha
                      << " where " << acl1.alpha << " was expected" << std::endl;
            nerr += 1;
        }
    }


    return nerr;
}


int main(int argc, char **argv)
{
    // MPI_Init(&argc, &argv);                     // Initialize MPI...
    // atexit((void (*) ()) MPI_Finalize);         // ...finalize at exit
    // suzerain::support::logging::initialize(     // Initialize logging
    //     MPI_COMM_WORLD,    
    //     suzerain::support::log4cxx_config_console); // correct???

    initialize_test_env(argc, argv);

    std::cout << "Running the antioch_constitutive unit test suite" << std::endl;

    std::cout << "Running test_save_load..." << std::endl;
    int ierr = test_save_load(MPI_COMM_WORLD);
    if (ierr!=0) std::cout << " FAILED." << std::endl;
    else std::cout << " passed." << std::endl;
                     

    return ierr;
}

#endif // HAVE_ANTIOCH