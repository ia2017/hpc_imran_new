#include "LidDrivenCavity.h"
#include "PoissonSolver.h"
#include <iostream>
#include <cstdlib>
#include <cblas.h>
#include <mpi.h>
#include <fstream>

using namespace std;


//---Constructors---
LidDrivenCavity::LidDrivenCavity(){
}   // nothing to do here

//Destructor
LidDrivenCavity::~LidDrivenCavity(){
    
    delete[] s;
    delete[] v;
    delete[] s_nplus1;
    delete[] v_nplus1;
}

//------Functions------

// Setting domain size (Lx * Ly)
void LidDrivenCavity::SetDomainSize(double xlen, double ylen){
    
    Lx=xlen;
    Ly=ylen;
}

// Setting no. of grid points (Nx * Ny)
void LidDrivenCavity::SetGridSize(int nx, int ny){
    
    Nx=nx;
    Ny=ny;
}

// Setting delta t
void LidDrivenCavity::SetTimeStep(double deltat){
    
    dt=deltat;
    
}

// Setting final time T
void LidDrivenCavity::SetFinalTime(double finalt){
    
    T=finalt;
}

// Setting Reynolds number Re
void LidDrivenCavity::SetReynoldsNumber(double re){
    
    Re=re;
}

// Setting Reynolds number Re
void LidDrivenCavity::SetPartitions(int px, int py){
    
    Px=px;
    Py=py;
}

// Setting initial conditions
void LidDrivenCavity::Initialise(int rk, int sz){    
    
    rank=rk;
    size=sz;

    // Calculating dx and dy
    dx=Lx/(Nx-1);
    dy=Ly/(Ny-1);
    
    /* Allocating memory for matrices and
     * setting vorticity and steamfunction to zero*/
     
    v = new double[Nx*Ny]();
    v_nplus1 = new double[Nx*Ny]();
    s = new double[Nx*Ny]();
    s_nplus1 = new double[Nx*Ny]();
    
    if(rank==0){
        cout<<"Initialised lid driven cavity."<<endl;   // Output message
    }
    
    //----Organising processes as a Cartesian grid-----
        
    MPI_Comm mygrid;
    
    const int dims = 2;         // 2-D grid
    int sizes[dims] = {Px, Py};     // Size of grid based on partitioning
    int periods[dims] = {0, 0};     // Non-cyclic
    int reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, dims, sizes, periods,
                    reorder, &mygrid);  //Creating grid
    
    MPI_Cart_shift(mygrid,0,1,&north,&south);   // Getting north and south values
    MPI_Cart_shift(mygrid,1,1,&west,&east);     // Getting east and west values
    
}

void LidDrivenCavity::BoundaryVorticity(){
    
    //------Updating boundary vorticity-------
        
    for(int i=0;i<Nx;i++){
        
        v[i*Ny+Ny-1]=(s[i*Ny+Ny-1]-s[i*Ny+Ny-2])*2/(dy*dy)-2*U/dy; //Top
        v[i*Ny]=(s[i*Ny]-s[i*Ny+1])*2/(dy*dy);                     //Bottom
    }
    
    for(int j=0;j<Ny;j++){
        
        v[j]=(s[j]-s[Ny+j])*2/(dx*dx);                             //Left
        v[(Nx-1)*Ny+j]=(s[(Nx-1)*Ny+j]-s[(Nx-2)*Ny+j])*2/(dx*dx);  //Right
        
    }
    
    
}

void LidDrivenCavity::InteriorVorticity(){

    
    //------Updating interior vorticity (at time t)-------
    
    for(int i=1;i<Nx-1;i++){
        for(int j=1;j<Ny-1;j++){
        
        v[i*Ny+j] = -((s[(i+1)*Ny+j]-2.0*s[i*Ny+j]+s[(i-1)*Ny+j])/(dx*dx)
                   +(s[i*Ny+j+1]-2.0*s[i*Ny+j]+s[i*Ny+j-1])/(dy*dy));
        
        }
    }    
    
}

void LidDrivenCavity::Integrate(){
    
    double t=0;     // Setting t=0
    
    // Setting up for poisson
    PoissonSolver* pos = new PoissonSolver();
    pos->SetGrid(Nx,Ny);
    pos->SetDeltas(dx,dy);
    pos->SetPartitions(Px,Py);
    pos->GetRankSize(rank,size);             // Getting rank and size
    pos->GetCart(north,south,east,west);     // Getting north, east, south and west values
    
    if(rank==0){
        cout<<"Integration start runtime: "<<endl;      // Output message
    }
    
    
    //-----------------------------START OF LOOP------------------------------------------//
    do{
    
        //------Computing interior vorticity (at time t+dt)-------
    
        for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny;j++){
        
            v_nplus1[i*Ny+j]=v[i*Ny+j]+dt*
                        (1.0/Re*((v[(i+1)*Ny+j]-2.0*v[i*Ny+j]+v[(i-1)*Ny+j])/(dx*dx)+
                        (v[i*Ny+j+1]-2.0*v[i*Ny+j]+v[i*Ny+j-1])/(dy*dy))
                         -((s[i*Ny+j+1]-s[i*Ny+j-1])/(2*dy))*((v[(i+1)*Ny+j]-v[(i-1)*Ny+j])/(2*dx))
                         +((s[(i+1)*Ny+j]-s[(i-1)*Ny+j])/(2*dx))*((v[i*Ny+j+1]-v[i*Ny+j-1])/(2*dy)));
        
            }
        }
    
    
        /*------Using Poisson Solver to get --------
        * -----new stream function (at time t+dt)----*/
        pos->UpdateVorticity(v_nplus1);    // Update vorticity into Poisson Solver
        pos->UpdateStream(s);          // Update stream function
     
        // Solving for Serial Case
        if(size==1){
            pos->SolveSerial();         // includes setting BC
        }
        // Solving for Pararllel Case
        else{
            pos->SetBoundaryConditions();  // boundaries streams = 0
            pos->LocalVorticity();         // partition vorticity
            pos->LocalStream();            // partition stream function
            pos->SolveParallel();
        }
     
        pos->GetStream(s_nplus1);   // Extracting solved stream function
        
        //-----Updating for next time iterations-------
    
        cblas_dcopy(Nx*Ny,s_nplus1,1,s,1);      // s=s_nplus1;
        cblas_dcopy(Nx*Ny,v_nplus1,1,v,1);      // v=v_nplus1
        t=t+dt;                                 // updating time

        LidDrivenCavity::BoundaryVorticity();   // Updating Boundary Voricity
    
        if(rank==0){
            cout<<"t: "<<t<<endl;               // Output runtime
        }
    
    
    }while(t<T);
    //-----------------------------END OF LOOP------------------------------------------//

    if(rank==0){
        cout<<"Integration Complete!"<<endl;
    }
        
}

// Print vorticity
void LidDrivenCavity::printv(){
    
    cout<<endl<<"Vorticity: "<<endl<<endl;
        
        for (int j=Ny-1;j>=0;j--){
            cout<<"[";
            
            for (int i=0;i<Nx;i++){
                
            cout.precision(3);
            cout << v[i*Ny+j];
            if (i==Nx-1){
                break;
            }
            cout<<" ";
        }
        cout<<"]"<<endl;
    
    }
    cout<<endl<<endl;
}

// Print stream function
void LidDrivenCavity::prints(){
            
    cout<<endl<<"Stream function: "<<endl<<endl;
        
        for (int j=Ny-1;j>=0;j--){
            cout<<"[";
            
            for (int i=0;i<Nx;i++){
            
            cout.precision(3);
            cout << s[i*Ny+j];
            if (i==Nx-1){
                break;
            }
            cout<<" ";
        }
        cout<<"]"<<endl;
    
    }
    cout<<endl<<endl;
}

// Output data into .txt files
void LidDrivenCavity::OutputFile(){
    
    int Re1=int(Re);    // For file labelling
    
    // Creating new file for vorticity
    ofstream fout("lidvortex_re" + std::to_string(Re1) + ".txt",ios::out | ios::trunc);
    fout.precision(5);
    
    
    for (int j=Ny-1;j>=0;j--){        
        for (int i=0;i<Nx;i++){
        
            fout.precision(5);
            fout.width(10);
            fout << v[i*Ny+j];        
        
        }
        fout<<endl;
    }
    
    
    fout.close();
    fout.clear();
    
    // Data for stream function
    ofstream fout1("lidstream_re" + std::to_string(Re1) + ".txt",ios::out | ios::trunc);
    fout1.precision(5);
    
    for (int j=Ny-1;j>=0;j--){        
        for (int i=0;i<Nx;i++){
        
            fout1.precision(5);
            fout1.width(10);
            fout1 << s[i*Ny+j];        
        
        }
        fout1<<endl;
    }
    
    
    fout1.close();
    fout1.clear();
    
    // Data for parameters
    ofstream fout2("lidparameters_re" + std::to_string(Re1) + ".txt",ios::out | ios::trunc);
    fout2.precision(3);
    
    fout2<<"Nx: "<<Nx<<endl;
    fout2<<"Ny: "<<Ny<<endl;
    fout2<<"Lx: "<<Lx<<endl;
    fout2<<"Ly: "<<Ly<<endl;
    fout2<<"Px: "<<Px<<endl;
    fout2<<"Py: "<<Py<<endl;
    fout2<<"dt: "<<dt<<endl;
    fout2<<"T: "<<T<<endl;
    fout2<<"Re: "<<Re<<endl;
    
    fout2.close();
    fout2.clear();
    
}