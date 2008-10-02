#include <stdio.h>
#include <string.h>
#include <math.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <fstream>
#include <iterator>
using namespace std;

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#include "legendrePoly.h"
#include "residual.h"
#include "solver.h"
#include "tools.h"
#include "burgers.h"



//-------------------------------------------------
// Evaluate BC function
//
int
boundaryCondition(const char *fname, const double time, double *UB)
{ 
  int ierr;
  const double PI = 3.141592653589793e+00;
  
  // pick requested bcs
  if( strcmp(fname, "nwave")==0 ){
    ierr = nWave(-1.0, time, 1e-2, 1.0, 16.0, 1.0, &(UB[0]));
    ierr = nWave( 1.0, time, 1e-2, 1.0, 16.0, 1.0, &(UB[1]));
  } else if( strcmp(fname, "zeroFunction")==0 ){
    UB[0] = UB[1] = 0.0;
  } else if( strcmp(fname, "steadyShock")==0 ){
    UB[0] =  1.0;
    UB[1] = -1.0;
  } else{
    printf("Unknown BC function.  Exiting.\n"); fflush(stdout);
    return -1;
  }

//    // 1 + sine
//    UB[0] =  1.0 + 0.5*sin(10.0*2.0*PI*time);
//    UB[1] = -1.0 - 0.5*cos(10.0*2.0*PI*time);

  return 0;
}

//-------------------------------------------------
// Set initial state
//
int
initialCondition(const char *fname, const double UB[2], const int N, double *U)
{ 

  // L2 projection of ic function onto solution space
  int ierr;
  const int solnOrder = N-1+2;
  const int quadOrder = 4*solnOrder + 1; // exact for poly of order 3*N 
                                         // (of course, ic may not be poly, but hopefully error is small enough)
  const int Nquad = quadOrder/2 + 1;

  double xq[Nquad], wq[Nquad];
  double phi[N];
  double rhs[N], u, u0;

  gsl_matrix *iMM;

  // Get quadrature points
  ierr = legendreGaussQuad(Nquad, xq, wq);
  if( ierr != 0 ) return ierr;

  // zero rhs
  for( int ii=0; ii<N; ii++ ) rhs[ii] = 0.0;

  // Loop over quadrature points to form right hand side vector
  for( int iquad=0; iquad<Nquad; iquad++ ){
  
    // Evaluate IC
    if( strcmp(fname, "nwave")==0 ){
      ierr = nWave(xq[iquad], 0.0, 1e-2, 1.0, 16.0, 1.0, &u);
      if( ierr != 0 ) return ierr;
    } else if( strcmp(fname, "zeroFunction")==0 ){
      u = 0.0;
    } else{
      printf("Unknown IC function.  Exiting.\n"); fflush(stdout);
      return -1;
    }
      

    // Write solution to screen         
    FILE *fp = fopen ("nWave.dat", "a");
    fprintf(fp, "%.15E, %.15E\n", xq[iquad], u);
    fclose(fp);

    u0 = u - (UB[0]*0.5*(1.0-xq[iquad]) + UB[1]*0.5*(1.0+xq[iquad]));

    // Evaluate basis
    ierr = legendrePolyZero(N, xq[iquad], phi, (double *)NULL);
    if( ierr != 0 ) return ierr;

    // Add to right hand side
    for( int ii=0; ii<N; ii++ ) rhs[ii] += wq[iquad]*phi[ii]*u0;
  }

  // Invert system to compute initial state
  iMM = gsl_matrix_calloc((size_t)N, (size_t)N);
  ierr = legendre0InverseMassMatrix(N, iMM->data);

  gsl_vector_view gslU = gsl_vector_view_array(U  , N);
  gsl_vector_view gslR = gsl_vector_view_array(rhs, N);

  gsl_blas_dgemv(CblasNoTrans, 1.0, iMM, &gslR.vector, 0.0, &gslU.vector);

  // Clean up
  gsl_matrix_free(iMM);

  return 0;
}

// A helper function to simplify the main part.
template<class T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(cout, " ")); 
    return os;
}

//-------------------------------------------------
// Reads ascii restart file
//
int
readRestartFile(const char *fname, double *time, double UB[2], int Nmode, double *U)
{
  int ierr, NmodeRead, minNmode; 

  FILE *fp = fopen(fname, "r");
  if( fp == NULL ){
    printf("Could not read file %s.  Exiting.\n", fname); fflush(stdout);
    return -1;
  }

  ierr = fscanf(fp, "%lf\n", time);
  if( ierr != 1 ) return -1;

  ierr = fscanf(fp, "%lf %lf\n", &(UB[0]), &(UB[1]));
  if( ierr != 2 ) return -1;

  ierr = fscanf(fp, "%d\n", &NmodeRead);
  if( ierr != 1 ) return -1;

  minNmode = NmodeRead;
  if( Nmode < minNmode ) minNmode = Nmode;

  for( int ii=0; ii<minNmode; ii++ ){
    ierr = fscanf(fp, "%lf\n", &(U[ii]));
    if( ierr != 1 ) return -1;
  }

  return 0;
}


//-------------------------------------------------
// Main routine.  Reads command line inputs
// and calls appropriate solver.
//
int
main( int argc, char * argv[] )
{
  int ierr, Nmode, nIter, nStep, nwrite;
  double nu, dt;
  burgersSteady Bsteady, *pBsteady;
  burgersUnsteady Bunsteady, *pBunsteady;
  bool steady, restart;

  pBsteady = &Bsteady;
  pBunsteady = &Bunsteady;

  // Declare a group of options that will be 
  // allowed only on command line
  po::options_description generic("Generic options");
  generic.add_options()
    ("help", "produce help message")    
    ;
    
  // Declare a group of options that will be 
  // allowed both on command line and in
  // config file
  po::options_description config("Configuration options");
  config.add_options()
    ("equation", po::value< string >()->default_value("Burgers"), 
     "Burgers (Burgers direct numerical simulation) or RAB (Reynolds-averaged Burgers)")

    ("viscosity", po::value< double >(&nu)->default_value(1e-1), 
     "fluid viscosity")

    ("steady", po::value< bool >(&steady)->default_value(false), 
     "true (for steady soln) or false (for unsteady)")

    ("nmode", po::value< int >(&Nmode)->default_value(2), 
     "number of spatial basis functions")

    ("nstep", po::value< int >(&nStep)->default_value(2), 
     "number of time steps (unsteady ONLY)")

    ("dt", po::value< double >(&dt)->default_value(1e-3), 
     "time step size (unsteady ONLY)")

    ("nwrite", po::value< int >(&nwrite)->default_value(nStep), 
     "write solution every nwrite steps (unsteady ONLY)")

    ("niter", po::value< int >(&nIter)->default_value(2), 
     "maximum number of Newton iterations (steady ONLY)")

    ("ic", po::value< string >()->default_value("zeroFunction"), 
     "name of function that sets initial condition (NOTE: not yet enabled)")

    ("bc", po::value< string >()->default_value("zeroFunction"), 
     "name of function that sets boundary condition (NOTE: not yet enabled)")

    ("output-file", po::value< string >()->default_value("burgers.dat"), 
     "name of output file for final solution dump")

    ("restart", po::value< bool >(&restart)->default_value(false), 
     "true for restarted case, false for fresh start")

    ("restart-file", po::value< string >()->default_value("burgers.dat"), 
     "name of file to restart from (restart=true ONLY)")

    ;
  
  // Hidden options, will be allowed both on command line and
  // in config file, but will not be shown to the user.
  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", po::value< string >(), "input file")
    ;
  
  po::options_description cmdline_options;
  cmdline_options.add(generic).add(config).add(hidden);
  
  po::options_description config_file_options;
  config_file_options.add(config).add(hidden);
  
  po::options_description visible("Allowed options");
  visible.add(generic).add(config);
  
  po::positional_options_description p;
  p.add("input-file", -1);
  
  // parse command line options
  po::variables_map vm;
  store(po::command_line_parser(argc, argv).
	options(cmdline_options).positional(p).run(), vm);
    
  // help if requested
  if (vm.count("help")) {
    cout << "\n";
    cout << "Usage: burgers [FILE] [options]" << "\n";
    cout << "\n";
    cout << "Description:\n";
    cout << "  Burgers (or Reynolds-average Burgers) equation solver.\n";
    cout << "  All configuration options described below may be set in\n";
    cout << "  FILE or on the command line.\n";
    cout << "\n";
    cout << "  Spatial discretization is Legendre-Galerkin with\n";
    cout << "  user-specified order.  For unsteady problems, temporal\n";
    cout << "  discretization is explicit 4th-order Runge-Kutta scheme.\n";
    cout << "  For steady problems, Newton-Raphson method is used.\n";
    cout << "\n";
    cout << "  For more information, contact <oliver@ices.utexas.edu>";
    cout << visible << "\n";
    return 0;
  }

  // dump out all option values being used
  if (vm.count("input-file")){
    cout << "\n";
    cout << "Reading input file " 
	 << vm["input-file"].as< string >() << "\n";
    cout << "\n";

    // parse configure file options
    ifstream ifs( (vm["input-file"].as< string >()).c_str() );
    store(parse_config_file(ifs, config_file_options), vm);
    notify(vm);
    
  } else{
    cout << "\n";
    cout << "No input file specified.\n";
    cout << "Using default settings/command line options only.\n";
    cout << "\n";
  }

  if (vm.count("equation")){
    cout << "equation is " 
	 << vm["equation"].as< string >() << "\n";
  }
  
  if (vm.count("viscosity")){
    cout << "viscosity is " 
	 << vm["viscosity"].as< double >() << "\n";
  }

  if (vm.count("steady")){
    if( steady ){
      cout << "solution is steady\n";
      cout << "using Newton solver with maximum of " << nIter << " iterations\n";
    } else{
      cout << "solution is unsteady\n";
      cout << "using RK4 time integration with dt = " << dt 
	   << " and nstep = " << nStep 
	   << " and nwrite = " << nwrite << "\n";
    }
  }

  
  if (vm.count("nmode")){
    cout << "nmode is " 
	 << vm["nmode"].as< int >() << "\n";
  }
   
  if (vm.count("output-file")){
    cout << "Output file is " 
	 << vm["output-file"].as< string >() << "\n";
  }

  if (vm.count("restart")){
    if( restart ){
      cout << "Restarting from " << vm["restart-file"].as< string >() << "\n";
    } else{
      cout << "This is a fresh start.\n";

      if (vm.count("bc")){
	cout << "using function " 
	     << vm["bc"].as< string >() << " to set boundary condition\n";
      }
      
      if (vm.count("ic")){
	cout << "using function " 
	     << vm["ic"].as< string >() << " to set initial condition\n";
      }

    }
  }

  cout << "\n";

  if( steady ){
    pBsteady->Nmode = Nmode;
    pBsteady->nu = nu;

    pBsteady->U = gsl_vector_calloc( (size_t)Nmode );
    pBsteady->R = gsl_vector_calloc( (size_t)Nmode );

    // Set initial condition
    if( restart ){
      printf("Reading restart file!\n"); fflush(stdout);
      double tmptime;
      readRestartFile((vm["restart-file"].as< string >()).c_str(), 
		      &tmptime, pBsteady->UB, Nmode, pBsteady->U->data);
      if( ierr != 0 ) return ierr;

    } else {
      ierr = boundaryCondition((vm["bc"].as< string >()).c_str(), 0.0, pBsteady->UB);
      if( ierr != 0 ) return ierr;

      ierr = initialCondition((vm["ic"].as< string >()).c_str(), pBsteady->UB, Nmode, pBsteady->U->data);
      if( ierr != 0 ) return ierr;
    }

    ierr = steadyNewton(nIter, pBsteady);
    if( ierr != 0 ) return ierr;

    // Dump final solution to file         
    FILE *fp = fopen ((vm["output-file"].as< string >()).c_str(), "w");
    
    printf("\n");
    printf("Writing %s... ", (vm["output-file"].as< string >()).c_str() );

    fprintf(fp, "%.15E\n", 0.0); // time
    fprintf(fp, "%.15E %.15E\n", pBsteady->UB[0], pBsteady->UB[1]); // BCs
    fprintf(fp, "%d\n", Nmode); // number of modes
    gsl_vector_fprintf(fp, pBsteady->U, "%.15E"); // solution coefficients
    fclose(fp);
    printf("done.\n"); fflush(stdout);

    // Clean up memory
    gsl_vector_free(pBsteady->U);
    gsl_vector_free(pBsteady->R);

    
  } else{
    pBunsteady->Nmode = Nmode;

    pBunsteady->Nstep = nStep;
    pBunsteady->Nwrite = nwrite;
    pBunsteady->Nstat = nStep/2; // convert to user set soon

    pBunsteady->nu = nu;
    pBunsteady->time = 0.0; // may be overwritten by restart file read
    pBunsteady->dt = dt;

    pBunsteady->U    = gsl_vector_calloc( (size_t)Nmode );
    pBunsteady->Uavg = gsl_vector_calloc( (size_t)Nmode );

    if( (vm["bc"].as< string >()).length() > 79 ){
      printf("BC function name is too long!\n"); fflush(stdout);
      return -1;
    }
    strcpy(pBunsteady->BCname, (vm["bc"].as< string >()).c_str());

    // Set initial condition
    if( restart ){
      printf("Reading restart file!\n"); fflush(stdout);
      readRestartFile((vm["restart-file"].as< string >()).c_str(), 
		      &(pBunsteady->time), pBunsteady->UB, Nmode, pBunsteady->U->data);
      if( ierr != 0 ) return ierr;
    } else {
      ierr = boundaryCondition((vm["bc"].as< string >()).c_str(), pBunsteady->time, pBunsteady->UB);
      if( ierr != 0 ) return ierr;

      ierr = initialCondition((vm["ic"].as< string >()).c_str(), pBunsteady->UB, Nmode, pBunsteady->U->data);
      if( ierr != 0 ) return ierr;
    }

    pBunsteady->UBavg[0] = pBunsteady->UBavg[1] = 0.0;

    ierr = unsteadyRK4(pBunsteady);
    if( ierr != 0 ) return ierr;

    // Dump final solution to file         
    FILE *fp = fopen ((vm["output-file"].as< string >()).c_str(), "w");
    
    printf("\n");
    printf("Writing %s... ", (vm["output-file"].as< string >()).c_str() );

    fprintf(fp, "%.15E\n", pBunsteady->time); // time
    fprintf(fp, "%.15E %.15E\n", pBunsteady->UB[0], pBunsteady->UB[1]); // BCs
    fprintf(fp, "%d\n", Nmode); // number of modes
    gsl_vector_fprintf(fp, pBunsteady->U, "%.15E"); // solution coefficients
    fclose(fp);
    printf("done.\n"); fflush(stdout);

    // Clean up memory
    gsl_vector_free(pBunsteady->U   );
    gsl_vector_free(pBunsteady->Uavg);
  }

  printf("\n");
  return 0;
}
