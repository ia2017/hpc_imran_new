#include <iostream>
#include "LidDrivenCavity.h"
#include "PoissonSolver.h"
#include <chrono>
#include <mpi.h>
#include <boost/program_options.hpp>
using namespace std;

namespace po=boost::program_options;    //To reduce typing
 

int main(int argc, char **argv){
    
    // ------------- Boost options part ------------- //
    po::options_description opts(
    "To run lid-driven cavity problem");
    opts.add_options()
    ("Lx", po::value<double>()->default_value(1),"Length of the domain in the x-direction.")
    ("Ly", po::value<double>()->default_value(1),"Length of the domain in the y-direction.")
    ("Nx", po::value<int>()->default_value(20),"Number of grid points in the x-direction.")
    ("Ny", po::value<int>()->default_value(20),"Number of grid points in the y-direction.")
    ("Px", po::value<int>()->default_value(1),"Number of partitions in the x-direction.")
    ("Py", po::value<int>()->default_value(1),"Number of partitions in the y-direction.")
    ("dt", po::value<double>()->default_value(0.001),"Time step size.")
    ("T", po::value<double>()->default_value(5),"Final time.")
    ("Re", po::value<double>()->default_value(100),"Reynolds number.")    
    ("help","Print help message.");
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc,argv,opts),vm);
    po::notify(vm);
    
    if(vm.count("help")){
        cout<<"This program solves the lid-driven cavity problem."<<endl;
        cout<<opts<<endl;
        return 0;
    }
    
    auto start = chrono::steady_clock::now(); // Start timer
    
    // Getting variables from command line
    const double Lx = vm["Lx"].as<double>();
    const double Ly = vm["Ly"].as<double>();
    const int Nx = vm["Nx"].as<int>();
    const int Ny = vm["Ny"].as<int>();
    const int Px = vm["Px"].as<int>();
    const int Py = vm["Py"].as<int>();
    const double dt = vm["dt"].as<double>();
    const double T = vm["T"].as<double>();
    const double Re = vm["Re"].as<double>();
    
    
    // ---------- Initialising parallel program (MPI) -----------//
    //Initialising parallel program
    int rank, size, retval_rank, retval_size;
    MPI_Init(&argc, &argv);
    retval_rank = MPI_Comm_rank(MPI_COMM_WORLD, &rank); // zero-based
    retval_size = MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (retval_rank == MPI_ERR_COMM || retval_size == MPI_ERR_COMM) {
        cout << "Invalid communicator" << endl;
    }
     
    // Checking if CFL condition is met
    try{
        
        double dx=Lx/(Nx-1);
        double dy=Ly/(Ny-1);
    
        if (dt>=Re*dx*dy/4) throw logic_error("dt too large");

    }
    catch(const exception& e){
        
        if(rank==0){
            cout<< "Error in input: "<<e.what()<<endl;  // Catching errors
        }
        return 1;
    }
    
    // Checking partitioning and no. of processes
    if(Px*Py!=size){
        if(rank==0){
            cout<<"Error: Partitions do not correspond with number of processes chosen."<<endl;
        }
        
        return 1;
    }
    
    // Checking if partitioning is less than no. of points
    if(Px>=(Nx-2) || Py>=(Ny-2)){
        if(rank==0){
            cout<<"Error: Not enough points for partitions."<<endl;
        }
        
        return 1;
    }
    
     
    //----------------------------START OF SOLVER------------------------------//
    
    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();

    // Configure the solver here...
    solver->SetDomainSize(Lx, Ly);
    solver->SetGridSize(Nx, Ny);
    solver->SetReynoldsNumber(Re); 
    solver->SetTimeStep(dt);
    solver->SetFinalTime(T);
    solver->SetPartitions(Px,Py);
     
    // Setting initial and boundary conditions
    solver->Initialise(rank,size); // Includes MPI and Cart Create
    solver->BoundaryVorticity();
    solver->InteriorVorticity();

    // Run the solver
    solver->Integrate();
    
    if (rank==0){
        solver->OutputFile();   // Output data to .txt files
    }  
    
    //----------------------------END OF SOLVER------------------------------//  
  
    auto end = chrono::steady_clock::now(); //End timer
    
    // Output code runtime
    if(rank==0){
        cout<<endl;
        cout<<"Elapsed time in milliseconds : " 
            <<chrono::duration_cast<chrono::milliseconds>(end - start).count()
            <<" ms" << endl;

        cout<<"Elapsed time in seconds : " 
            <<chrono::duration_cast<chrono::seconds>(end - start).count()
            <<" sec"<<endl;
    }
    
    MPI_Finalize(); // Finalizing MPI

	return 0;
}