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

#include "antioch_constitutive.hpp"

#include <esio/esio.h>

#include <suzerain/common.hpp>
#include <suzerain/support/support.hpp>
#include <suzerain/support/logging.hpp>

#ifdef SUZERAIN_HAVE_ANTIOCH

/**
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

/**
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

/**
 * Check that call to init_antioch puts objects in desired state
 */
int test_init_antioch(const std::string& chem_xml_file)
{

    using suzerain::reacting::antioch_constitutive;
    using suzerain::real_t;

    // Prepare input data
    std::vector<std::string> species_names;
    const unsigned int Ns=5;
    species_names.reserve(Ns);
    species_names.push_back( "N2" );
    species_names.push_back( "O2" );
    species_names.push_back( "N" );
    species_names.push_back( "O" );
    species_names.push_back( "NO" );

    const real_t Le = 0.7;
    const real_t alpha = 0.5;

    antioch_constitutive acl1(species_names, chem_xml_file, Le, alpha);

    acl1.init_antioch();

    // Check that initialized objects appear to have valid information...
    // but not exhaustively... that's for antioch to test
    int nerr=0;

    // The mixture object
    // ... number of species
    if (acl1.mixture->n_species() != Ns) {
        std::cerr << "Error: acl1.mixture->n_species() = " << acl1.mixture->n_species()
                  << " but should be " << Ns << std::endl;
        nerr += 1;
    }

    // ... species name map
    const std::vector<Antioch::Species> species_list = acl1.mixture->species_list();
    for (unsigned int i=0; i<Ns; i++) {
        // convenience
        const std::map<std::string,Antioch::Species>& smap = acl1.mixture->species_name_map();

        if( smap.find( species_names[i] )->second != species_list[i] ){
            std::cerr << "Error: species name map and species list ordering mismatch" << std::endl
                      << "species_name_map = " << smap.find( species_names[i] )->second
                      << ", species_list = " << species_list[i] << std::endl;
            nerr += 1;
        }
    }

    // The reaction object
    // ... number of species
    if (acl1.reactions->n_species() != acl1.mixture->n_species()) {
        std::cerr << "Error: acl1.reactions->n_species() = " << acl1.mixture->n_species()
                  << " but should be " << Ns << std::endl;
        nerr += 1;
    }

    // ... number of reactions (NOTE: result specific to air_sp5.xml and species selected!)
    if (acl1.reactions->n_reactions() != 5) {
        std::cerr << "Error: acl1.reactions->n_reactions() = " << acl1.reactions->n_reactions()
                  << " but should be " << 5 << std::endl;
        nerr += 1;
    }

    // The cea_thermo object has required curve fits
    if (!acl1.cea_thermo->check()) {
        std::cerr << "Error: acl1.cea_thermo does not have curve fits for all specis (why not?)." << std::endl;
        nerr += 1;
    }

    // the kinetics evaluator
    if (acl1.kinetics->n_species()!=Ns) {
        std::cerr << "Error: acl1.kinetics->n_species() = " << acl1.kinetics->n_species()
                  << " but should be " << Ns << std::endl;
        nerr += 1;
    }

    if (acl1.kinetics->n_reactions()!=5) {
        std::cerr << "Error: acl1.kinetics->n_reactions() = " << acl1.kinetics->n_reactions()
                  << " but should be " << 5 << std::endl;
        nerr += 1;
    }


    return nerr;
}


/**
 * Check that call to init_antioch puts objects in desired state
 */
int test_init_antioch_CPAir()
{

    using suzerain::reacting::antioch_constitutive;
    using suzerain::real_t;

    // Prepare input data
    std::vector<std::string> species_names;
    const unsigned int Ns=1;
    species_names.reserve(Ns);
    species_names.push_back( "CPAir" );

    const real_t Le = 0.7;
    const real_t alpha = 0.5;

    //antioch_constitutive acl1(species_names, chem_xml_file, Le, alpha);
    antioch_constitutive acl1;

    acl1.species_names = species_names;
    acl1.Le = Le;
    acl1.alpha = alpha;

    acl1.init_antioch();

    // Check that initialized objects appear to have valid information...
    // but not exhaustively... that's for antioch to test
    int nerr=0;

    // The mixture object
    // ... number of species
    if (acl1.mixture->n_species() != Ns) {
        std::cerr << "Error: acl1.mixture->n_species() = " << acl1.mixture->n_species()
                  << " but should be " << Ns << std::endl;
        nerr += 1;
    }

    // ... species name map
    const std::vector<Antioch::Species> species_list = acl1.mixture->species_list();
    for (unsigned int i=0; i<Ns; i++) {
        // convenience
        const std::map<std::string,Antioch::Species>& smap = acl1.mixture->species_name_map();

        if( smap.find( species_names[i] )->second != species_list[i] ){
            std::cerr << "Error: species name map and species list ordering mismatch" << std::endl
                      << "species_name_map = " << smap.find( species_names[i] )->second
                      << ", species_list = " << species_list[i] << std::endl;
            nerr += 1;
        }
    }

    // The reaction object
    // ... number of species
    if (acl1.reactions->n_species() != acl1.mixture->n_species()) {
        std::cerr << "Error: acl1.reactions->n_species() = " << acl1.mixture->n_species()
                  << " but should be " << Ns << std::endl;
        nerr += 1;
    }

    // ... number of reactions (NOTE: result specific to air_sp5.xml and species selected!)
    if (acl1.reactions->n_reactions() != 0) {
        std::cerr << "Error: acl1.reactions->n_reactions() = " << acl1.reactions->n_reactions()
                  << " but should be " << 0 << std::endl;
        nerr += 1;
    }

    // The cea_thermo object has required curve fits
    if (!acl1.cea_thermo->check()) {
        std::cerr << "Error: acl1.cea_thermo does not have curve fits for all specis (why not?)." << std::endl;
        nerr += 1;
    }

    return nerr;
}


/**
 * Check that call to antioch_constitutive::evaluate "works"
 */
int test_evaluate(const std::string& chem_xml_file)
{

    using suzerain::reacting::antioch_constitutive;
    using suzerain::real_t;

    // Prepare input data
    std::vector<std::string> species_names;
    const unsigned int Ns=5;
    species_names.reserve(Ns);
    species_names.push_back( "N2" );
    species_names.push_back( "O2" );
    species_names.push_back( "N" );
    species_names.push_back( "O" );
    species_names.push_back( "NO" );

    const real_t Le = 0.7;
    const real_t alpha = 0.5;

    antioch_constitutive acl1(species_names, chem_xml_file, Le, alpha);

    acl1.init_antioch();

    // Set up state (not physically meaningful yet)
    real_t e          = 5717500; // a really big number s.t. T isn't really, really small
    real_t m[3]       = {0.0, 0.0, 0.0};
    real_t rho        = 1.0;
    real_t species[5] = {0.5, 0.2, 0.1, 0.1, 0.1};
    real_t cs[5]      = {0.5, 0.2, 0.1, 0.1, 0.1};

    real_t Tguess     = -1; 

    // Storage for computed quantities
    real_t T=-1, p=-1, Ds[5], mu, kap,
        hs[5], om[5]={1.0, 1.0, 1.0, 1.0, 1.0}, a=-1, Cp=0;

    // Eval rxn sources, trans, thermo
    acl1.evaluate(e, m, rho, species, cs, Tguess,  /* input  */
                  T, p, Ds, mu, kap, hs, om, a, Cp /* output */);


    // check that it did something potentially sane
    //
    // FIXME: make these checks stronger once final version of
    // evaluate function is complete

    int nerr=0;

    //... T and p are positive
    if (T<=0.0) {
        std::cerr << "Error: T = " << T << " <= 0 is not allowed!" << std::endl;
        nerr += 1;
    }

    if (p<=0.0) {
        std::cerr << "Error: p = " << p << " <= 0 is not allowed!" << std::endl;
        nerr += 1;
    }

    //... source terms sum to zero
    real_t som = 0.0;
    for (unsigned int i=0; i<5; ++i) som += om[i];
    if (std::abs(som)>5e-10) { // Tolerance empirical, max src term has magn 1.8e6
        std::cerr << "Error: reaction source terms do not sum to zero. "
                  << "Computed sum = " << som << std::endl;
        nerr += 1;
    }

    return nerr;
}

/**
 * Check that call to antioch_constitutive::evaluate "works"
 */
int test_evaluate_eigen(const std::string& chem_xml_file)
{

    using suzerain::reacting::antioch_constitutive;
    using suzerain::real_t;
    using suzerain::Vector3r;
    using suzerain::VectorXr;

    // Prepare input data
    std::vector<std::string> species_names;
    const unsigned int Ns=5;
    species_names.reserve(Ns);
    species_names.push_back( "N2" );
    species_names.push_back( "O2" );
    species_names.push_back( "N" );
    species_names.push_back( "O" );
    species_names.push_back( "NO" );

    const real_t Le = 0.7;
    const real_t alpha = 0.5;

    antioch_constitutive acl1(species_names, chem_xml_file, Le, alpha);

    acl1.init_antioch();

    // Set up state (not physically meaningful yet)
    real_t e   = 5717500; // a really big number s.t. T isn't really, really small
    Vector3r m (0.0, 0.0, 0.0);
    real_t rho = 1.0;
    VectorXr species(Ns), cs(Ns);

    species(0) = cs(0) = 0.5;
    species(1) = cs(1) = 0.2;
    species(2) = cs(2) = 0.1;
    species(3) = cs(3) = 0.1;
    species(4) = cs(4) = 0.1;

    real_t Tguess= -1; 

    // Storage for computed quantities
    real_t T=-1, p=-1, mu, kap, a=-1, Cp=0;
    VectorXr Ds(Ns), hs(Ns), om(Ns);
    om(0) = om(1) = om(2) = om(3) = om(4) = 1.0;

    // Eval rxn sources, trans, thermo
    acl1.evaluate(e, m, rho, species, cs, Tguess,  /* input  */
                  T, p, Ds, mu, kap, hs, om, a, Cp /* output */);


    // check that it did something potentially sane
    //
    // FIXME: make these checks stronger once final version of
    // evaluate function is complete

    int nerr=0;

    //... T and p are positive
    if (T<=0.0) {
        std::cerr << "Error: T = " << T << " <= 0 is not allowed!" << std::endl;
        nerr += 1;
    }

    if (p<=0.0) {
        std::cerr << "Error: p = " << p << " <= 0 is not allowed!" << std::endl;
        nerr += 1;
    }

    //... source terms sum to zero
    real_t som = 0.0;
    for (unsigned int i=0; i<Ns; ++i) som += om(i);
    if (std::abs(som)>5e-10) { // Tolerance empirical, max src term has magn 1.8e6
        std::cerr << "Error: reaction source terms do not sum to zero. "
                  << "Computed sum = " << som << std::endl;
        nerr += 1;
    }

    return nerr;
}

/**
 * Check that call to evaluate_pressure_derivs_and_trans "works"
 */
int test_evaluate_pressure_derivs_and_trans(const std::string& chem_xml_file)
{

    using suzerain::reacting::antioch_constitutive;
    using suzerain::real_t;
    using suzerain::Vector3r;
    using suzerain::VectorXr;

    // Prepare input data
    std::vector<std::string> species_names;
    const unsigned int Ns=5;
    species_names.reserve(Ns);
    species_names.push_back( "N2" );
    species_names.push_back( "O2" );
    species_names.push_back( "N" );
    species_names.push_back( "O" );
    species_names.push_back( "NO" );

    const real_t Le = 0.7;
    const real_t alpha = 0.5;

    antioch_constitutive acl1(species_names, chem_xml_file, Le, alpha);

    acl1.init_antioch();

    // Set up state (not physically meaningful yet)
    real_t e   = 5717500; // a really big number s.t. T isn't really, really small
    Vector3r m (1.0, 2.0, 3.0);
    real_t rho = 1.5;
    VectorXr species(Ns), cs(Ns);

    cs(0) = 0.5;
    cs(1) = 0.2;
    cs(2) = 0.1;
    cs(3) = 0.1;
    cs(4) = 0.1;

    species = rho*cs;

    real_t Tguess = -1;

    // Storage for computed quantities
    real_t T=-1, p=-1, p_rho=-1, p_rsum=-1, p_e=-1, mu=0, kap=0;
    real_t gamma=-1, a=-1;
    Vector3r p_m;
    VectorXr Ds(Ns);

    // Eval rxn sources, trans, thermo
    acl1.evaluate_pressure_derivs_and_trans(
        e, m, rho, species, cs, Tguess,  /* input */
        T, p, p_rho, p_rsum, p_m, p_e, mu, kap, Ds, gamma, a /* output */);


    // TODO: Add a finite difference check against acl1.evaluate

    int nerr=0;

    // FIXME: For now this passes as long as above call doesn't die somehow

    return nerr;
}


/**
 * Check that call to etots_from_T "works"
 */
int test_etots_from_T(const std::string& chem_xml_file)
{

    using suzerain::reacting::antioch_constitutive;
    using suzerain::real_t;
    using suzerain::Vector3r;
    using suzerain::VectorXr;

    // Prepare input data
    std::vector<std::string> species_names;
    const unsigned int Ns=5;
    species_names.reserve(Ns);
    species_names.push_back( "N2" );
    species_names.push_back( "O2" );
    species_names.push_back( "N" );
    species_names.push_back( "O" );
    species_names.push_back( "NO" );

    const real_t Le = 0.7;
    const real_t alpha = 0.5;

    antioch_constitutive acl1(species_names, chem_xml_file, Le, alpha);

    acl1.init_antioch();

    // Set up temperature
    real_t T=1000;

    // Storage for computed quantities
    VectorXr etots(Ns);

    // Eval rxn sources, trans, thermo
    acl1.etots_from_T(
        T,   /* input */
        etots /* output */);


    // TODO: Add check

    int nerr=0;

    // FIXME: For now this passes as long as above call doesn't die somehow

    return nerr;
}



int main(int argc, char **argv)
{

    // Check command line count.
    if( argc != 2 ) {
        std::cerr << "Error: Must specify reaction set XML input file." << std::endl;
        return 1;
    }

    std::string chem_xml_file(argv[1]);

    initialize_test_env(argc, argv);

    int etot=0;

    std::cout << "Running test_save_load..." << std::endl;
    int ierr = test_save_load(MPI_COMM_WORLD); etot += ierr;
    if (ierr!=0) std::cout << " FAILED." << std::endl;
    else std::cout << " passed." << std::endl;

    std::cout << "Running test_init_antioch..." << std::endl;
    ierr = test_init_antioch(chem_xml_file); etot += ierr;
    if (ierr!=0) std::cout << " FAILED." << std::endl;
    else std::cout << " passed." << std::endl;

    std::cout << "Running test_init_antioch_CPAir..." << std::endl;
    ierr = test_init_antioch_CPAir(); etot += ierr;
    if (ierr!=0) std::cout << " FAILED." << std::endl;
    else std::cout << " passed." << std::endl;

    std::cout << "Running test_evaluate..." << std::endl;
    ierr = test_evaluate(chem_xml_file); etot += ierr;
    if (ierr!=0) std::cout << " FAILED." << std::endl;
    else std::cout << " passed." << std::endl;

    std::cout << "Running test_evaluate..." << std::endl;
    ierr = test_evaluate_eigen(chem_xml_file); etot += ierr;
    if (ierr!=0) std::cout << " FAILED." << std::endl;
    else std::cout << " passed." << std::endl;

    std::cout << "Encountered " << etot << " total errors." << std::endl;

    return etot;
}

#endif // HAVE_ANTIOCH
