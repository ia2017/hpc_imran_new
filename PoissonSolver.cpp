/*---- HPC Coursework ---------
----Class to solve Poisson Equations ------
-----Created by: Imran bin Ahmad Azhar -----
-----------------------------------------*/

#include "PoissonSolver.h"
#include <cmath>
#include <cblas.h>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <mpi.h>

using namespace std;

// Constructor
PoissonSolver::PoissonSolver(){
}   // Nothing to do here

// Destructor
PoissonSolver::~PoissonSolver(){
    delete[] s;
    delete[] v;
    delete[] s_local;
    delete[] v_local;
}

// Setting Nx and Ny
void PoissonSolver::SetGrid(int nx, int ny){
    
    Nx=nx;
    Ny=ny;
}

// Setting dx and dy
void PoissonSolver::SetDeltas(double deltax,double deltay){
    
    dx=deltax;
    dy=deltay;
    
}   

// Setting partitioning
void PoissonSolver::SetPartitions(int px, int py){
    
    Px=px;
    Py=py;
}

// Getting rank and size from MPI
void PoissonSolver::GetRankSize(int rk, int sz){
    
    rank=rk;
    size=sz;

}

// Getting parameters for MPI communication
void PoissonSolver::GetCart(int n, int s, int e, int w){
    
    north=n;
    south=s;
    east=e;
    west=w;
    
}

// Updating vorticity
void PoissonSolver::UpdateVorticity(double* vort){
    
    v=new double[Nx*Ny]();
    
    for (int i=0;i<Nx*Ny;i++){
        v[i]=vort[i];
    }
    
} 

// Updating stream function
void PoissonSolver::UpdateStream(double* stream){
    
    s=new double[Nx*Ny]();
    
    for (int i=0;i<Nx*Ny;i++){
        s[i]=stream[i];
    }
    
}     

void PoissonSolver::SetBoundaryConditions(){
    
    //------Setting BCs at the wall where s=0-------
    
    for(int i=0;i<Nx;i++){
        s[i*Ny+Ny-1]=0.0;       //Top
        s[i*Ny]=0.0;           //Bottom
    }
    
    for(int j=0;j<Ny;j++){
        
        s[j]=0.0;              //Left
        s[(Nx-1)*Ny+j]=0.0;    //Right
        
    }
    
}

void PoissonSolver::GetStream(double* stream){
    
    cblas_dcopy(Nx*Ny,s,1,stream,1); // Copy v into s
    
}

// Partition vorticity into sub domains
void PoissonSolver::LocalVorticity(){
    
    // Checking if no. of points equates to partitioning 
    
    double Nxlc,Nylc;   //double to check decimals
    padx=0,pady=0;

    Nxlc=(double(Nx)-2.0)/double(Px);
    Nylc=(double(Ny)-2.0)/double(Py);
    
    if((Nx-2)%Px==0) Nxl=(Nx-2)/Px;
    else{ 
        Nxl=int(ceil(Nxlc)); //rounding up and converting to int
        padx=Px-(Nx-2)%Px;
    }
    if((Ny-2)%Py==0) Nyl=(Ny-2)/Py;
    else{
        Nyl=int(ceil(Nylc)); //rounding up and converting to int
        pady=Py-(Ny-2)%Py;
    }
    
    //Determining cases for padded rows
    if(padx==0 && pady==0){
        cases=0;
    }
    else{
        cases=1;
    }
    
    
    //-----Reduce to unknown points GLOBAL points-------
    
    int nsv=Nxl*Px*Nyl*Py;  //No. of unknown points
    Nsd=Nxl*Nyl;           //No. of sub domain points
    
    double* v_sd=new double[nsv]();
    
    for(int i=1;i<Nx-1;i++){
        for(int j=1;j<Ny-1;j++){
            v_sd[(i-1)*(Nyl*Py)+j-1]=v[i*Ny+j];
        }
    }
    
    //-------Partitioning to subdomains and sending to processes------
    
    double* temp=new double[Nsd]();
    v_local=new double[Nsd]();
    
    //--------------COMM STARTS---------------------
    MPI_Barrier(MPI_COMM_WORLD);    //Synchronising processes
    
    if(rank==0){
    for(int i=0;i<Nxl;i++){
        for(int j=0;j<Nyl;j++){
            v_local[i*Nyl+j]=v_sd[i*(Nyl*Py)+j];
        }
    }
    
    int l=0,m=1;

    for(int i=1;i<size;i++){
        if(i%Py==0){
            l++;   //reset
            m=0;
        }
        for(int j=0;j<Nxl;j++){
            for (int k=0;k<Nyl;k++){
                temp[j*Nyl+k]=v_sd[(j+Nxl*l)*(Nyl*Py)+k+Nyl*m];
            }
        }
        m++;
        MPI_Send(temp, Nsd, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    }
    
    }
    else{
        MPI_Recv(v_local, Nsd, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    //--------------COMM ENDS---------------------
    
    /*
    if(rank==1){
        
        printVec1(v_local,Nsd);
    }*/
    
    delete[] temp;
    delete[] v_sd;
}

// Partition stream functions into sub domains
void PoissonSolver::LocalStream(){
    
    
    //-----Reduce to unknown points GLOBAL points-------
    
    int nsv=Nxl*Px*Nyl*Py;  //No. of unknown points
    
    double* s_sd=new double[nsv]();
    
    for(int i=1;i<Nx-1;i++){
        for(int j=1;j<Ny-1;j++){
            s_sd[(i-1)*(Nyl*Py)+j-1]=s[i*Ny+j];
        }
    }
    
    //-------Partitioning to subdomains and sending to processes------
    
    double* temp=new double[Nsd]();
    s_local=new double[Nsd]();
    
    //--------------COMM STARTS---------------------
    MPI_Barrier(MPI_COMM_WORLD);    //Synchronising processes
    
    if(rank==0){
    for(int i=0;i<Nxl;i++){
        for(int j=0;j<Nyl;j++){
            s_local[i*Nyl+j]=s_sd[i*(Nyl*Py)+j];
        }
    }
    
    int l=0,m=1;

    for(int i=1;i<size;i++){
        if(i%Py==0){
            l++;   //reset
            m=0;
        }
        for(int j=0;j<Nxl;j++){
            for (int k=0;k<Nyl;k++){
                temp[j*Nyl+k]=s_sd[(j+Nxl*l)*(Nyl*Py)+k+Nyl*m];
            }
        }
        m++;
        MPI_Send(temp, Nsd, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    }
    
    }
    else{
        MPI_Recv(s_local, Nsd, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    //--------------COMM ENDS---------------------
    
    /*
    if(rank==1){
        
        printVec1(v_local,Nsd);
    }*/
    
    delete[] temp;
    delete[] s_sd;
}

// Function to print out vector vertically (for checking only)
void PoissonSolver::printVec1(double v[],int size){
        
    cout<<endl; //leaving spacing
    
    for (int i=0;i<size;i++){
        
        cout.precision(5);          // Setting precision to 5
        cout<<"["<<v[i]<<"]"<<endl; // Output vector v
        
    }
    cout<<endl; //leaving spacing
}

// Function to solve linear system of equations using CG in parallel
void PoissonSolver::ConjugateGradient(double A[],double x[],double b[]){
    
    //-------Declaring CG parameters---------
    
    int off=0, k=0, m=0;   // Initialise variables
    int nsv=Nxl*Px*Nyl*Py;  // Unknown variables for global
    const double eps=0.00000001;  // Defining error constant
    int spc=Nyl+2;                  // Spacing for calculations
    double alpha_top,alpha_bottom,alpha,beta_top,beta_bottom,beta,
            alpha_top_r,alpha_bottom_r,beta_top_r,beta_bottom_r,
            sum, res, res_local;  // Declaring variables
    
    // Temporary memory for west and east transfers
    double temp_send[Nxl]; 
    double temp_recv[Nxl];
    
    // Memory allocations
    double* r_k=new double[Nsd]();
    double* r_kplus1=new double[Nsd]();
    double* p_k=new double[Nsd]();
    double* p_bound=new double[(Nxl+2)*(Nyl+2)]();
    double* s_bound=new double[(Nxl+2)*(Nyl+2)]();    // local stream + its boundaries    
    double* s_k=new double[Nsd]();        // stream at k iteration
    double* Apk=new double[Nsd]();
    double* z_k=new double[Nsd]();    // Preconditions
    double* z_kplus1=new double[Nsd]();
    //double* r_kplus1_global=new double[nsv]();    // r_k+1 global for error calculation   

    //----------------START of Conjugate Gradient Method--------------------//
    
    cblas_dcopy(Nsd,x,1,s_k,1);   //Copying initial s into s_k
    
    //Inserting LOCAL s into s_k
    
    for(int i=0;i<Nxl;i++){
        for(int j=0;j<Nyl;j++){
            s_bound[Nyl+2+i*(Nyl+2)+j+1]=s_k[i*Nyl+j];
        }
    }
    
    /*--------Communicate to find Initial Boundary points--------
     * --------once other processes are solved-----------
     * ---------simultaneously----------------------------*/
        
    //MPI_Barrier(MPI_COMM_WORLD);    // Synchronising processes

    // ---------------North Transfers---------------
    if(north>=0){
        
        MPI_Sendrecv(s_bound+(Nyl+2)+1, Nyl, MPI_DOUBLE, north, 1, 
                    s_bound+1, Nyl, MPI_DOUBLE, north, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    else{
        for(int j=0;j<Nyl;j++){
            s_bound[j]=0;   // if nothing is north
        }
    }
        
    // ---------------West Transfers---------------
    if(west>=0){
        
        memset(temp_send,0,sizeof(double)*Nxl);
        memset(temp_recv,0,sizeof(double)*Nxl);
            
        for(int i=0;i<Nxl;i++){
            temp_send[i]=s_bound[Nyl+2+i*(Nyl+2)+1];     // Preparing west transfers
        }
        
        MPI_Sendrecv(temp_send, Nxl, MPI_DOUBLE, west, 2, 
                    temp_recv, Nxl, MPI_DOUBLE, west, 3, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            
        for(int i=0;i<Nxl;i++){
                s_bound[Nyl+2+i*(Nyl+2)]=temp_recv[i];     // receive from west
        }
            
    }
    else{
        for(int i=0;i<Nxl;i++){
                s_bound[Nyl+2+i*(Nyl+2)]=0;     // if nothing is west
        }
    }
    // ---------------East Transfers---------------
    if(east>=0){
        
        memset(temp_send,0,sizeof(double)*Nxl);
        memset(temp_recv,0,sizeof(double)*Nxl);
            
        for(int i=0;i<Nxl;i++){
            temp_send[i]=s_bound[Nyl+2+i*(Nyl+2)+Nyl];   // Preparing east transfers
        }
        
        
    MPI_Sendrecv(temp_send, Nxl, MPI_DOUBLE, east, 3, 
                temp_recv, Nxl, MPI_DOUBLE, east, 2, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    for(int i=0;i<Nxl;i++){
        s_bound[Nyl+2+i*(Nyl+2)+Nyl+1]=temp_recv[i];   // receive from east
    }
    }
    else{
        for(int i=0;i<Nxl;i++){
            s_bound[Nyl+2+i*(Nyl+2)+Nyl+1]=0;   // if nothing is east
        }
    }
         
    // ---------------South Transfers---------------
    if(south>=0){
        
        MPI_Sendrecv(s_bound+(Nyl+2)*Nxl+1, Nyl, MPI_DOUBLE, south, 0, 
                    s_bound+(Nyl+2)*(Nxl+1)+1, Nyl, MPI_DOUBLE, south, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        
    }
    else{
        for(int j=0;j<Nyl;j++){
            s_bound[(Nyl+2)*(Nxl+1)+j+1]=0;     // if nothing is south
        }
    }
        
    //----------------END OF COMMUNICATION------------------//

    //------Calculating residual----------
    
    // Parallel multiply A and s0
    // For evenly distributed points
    if(cases==0){
        
        m=0, off=0; // Initialising local iterations
        
        // Finding sum
        for(int i=0;i<Nsd;i++){
            
            if(i%Nyl==0 && i!=0){
                off++;   // array offset
                m=0;
            }
            // Performing matrix vector multiplication and summing
            sum=0;
            sum=sum+A[0]*s_bound[m+1+off*spc];
            sum=sum+A[1]*s_bound[m+1+spc-1+off*spc];
            sum=sum+A[2]*s_bound[m+1+spc+off*spc];
            sum=sum+A[3]*s_bound[m+1+spc+1+off*spc];
            sum=sum+A[4]*s_bound[m+1+2*spc+off*spc]; 
                
            r_k[i]=b[i]-sum;   // r0=b-A*s0
            m++;
        }
    }
    // For padded rows/columns case (due to uneven no. of points)
    else{
        
        //---------Top side padded rows---------------
        if(east<0 && south>0){
        
        // Finding sum
        for(int i=0;i<Nxl;i++){
            for(int j=0;j<Nyl-pady;j++){
                    
                // Performing matrix vector multiplication and summing
                sum=0;
                sum=sum+A[0]*s_bound[j+1+i*spc];
                sum=sum+A[1]*s_bound[j+1+spc-1+i*spc];
                sum=sum+A[2]*s_bound[j+1+spc+i*spc];
                sum=sum+A[3]*s_bound[j+1+spc+1+i*spc];
                sum=sum+A[4]*s_bound[j+1+2*spc+i*spc]; 
                
                r_k[i*Nyl+j]=b[i*Nyl+j]-sum;   // r0=b-A*s0
                    
            }
        }
        
        }
        //----------Top right corner padded rows-------------
        else if(east<0 && south<0){
            
        // Finding sum
        for(int i=0;i<Nxl-padx;i++){
            for(int j=0;j<Nyl-pady;j++){
                    
                // Performing matrix vector multiplication and summing
                sum=0;
                sum=sum+A[0]*s_bound[j+1+i*spc];
                sum=sum+A[1]*s_bound[j+1+spc-1+i*spc];
                sum=sum+A[2]*s_bound[j+1+spc+i*spc];
                sum=sum+A[3]*s_bound[j+1+spc+1+i*spc];
                sum=sum+A[4]*s_bound[j+1+2*spc+i*spc]; 
            
                r_k[i*Nyl+j]=b[i*Nyl+j]-sum;   // r0=b-A*s0
                    
            }
        }
            
        }
        //----------Right side padded rows-------------
        else if(east>0 && south<0){
            
        // Finding sum
        for(int i=0;i<Nxl-padx;i++){
            for(int j=0;j<Nyl;j++){
                    
                // Performing matrix vector multiplication and summing
                sum=0;
                sum=sum+A[0]*s_bound[j+1+i*spc];
                sum=sum+A[1]*s_bound[j+1+spc-1+i*spc];
                sum=sum+A[2]*s_bound[j+1+spc+i*spc];
                sum=sum+A[3]*s_bound[j+1+spc+1+i*spc];
                sum=sum+A[4]*s_bound[j+1+2*spc+i*spc]; 
            
                r_k[i*Nyl+j]=b[i*Nyl+j]-sum;   // r0=b-A*s0
                    
            }
        }
            
        }
        else{
        
        m=0, off=0; // Initialising local iterations
        
        // Finding sum
        for(int i=0;i<Nsd;i++){
            
            if(i%Nyl==0 && i!=0){
                off++;   // array offset
                m=0;
            }
            
            // Performing matrix vector multiplication and summing
            sum=0;
            sum=sum+A[0]*s_bound[m+1+off*spc];
            sum=sum+A[1]*s_bound[m+1+spc-1+off*spc];
            sum=sum+A[2]*s_bound[m+1+spc+off*spc];
            sum=sum+A[3]*s_bound[m+1+spc+1+off*spc];
            sum=sum+A[4]*s_bound[m+1+2*spc+off*spc]; 
                
            r_k[i]=b[i]-sum;   // r0=b-A*s0
            m++;
        }
        
        }
        
    }
    
    //----------Initial Pre-conditioning Step---------------
    
    // Using Jacobi conditioner
    for(int i=0;i<Nsd;i++){
        z_k[i]=r_k[i]/A[2];
    }
       
    cblas_dcopy(Nsd,z_k,1,p_k,1);     //Inserting LOCAL r into LOCAL p
    
    //--------------------------START OF LOOP-----------------------------//
    do{
        // Reset for next iteration
        memset(p_bound,0,sizeof(double)*(Nxl+2)*(Nyl+2));
        
        // Getting local p + its boundaries for A*p multiplication
            
        for(int i=0;i<Nxl;i++){
            for(int j=0;j<Nyl;j++){
                p_bound[Nyl+2+i*(Nyl+2)+j+1]=p_k[i*Nyl+j];
            }
        }
        
        /*---------Communicate to find Boundary points--------
         * --------once other processes are solved-----------
         * ---------simultaneously----------------------------*/
        
        //MPI_Barrier(MPI_COMM_WORLD);    // Synchronising processes

        // ---------------North Transfers---------------
        if(north>=0){
        
            MPI_Sendrecv(p_bound+(Nyl+2)+1, Nyl, MPI_DOUBLE, north, 1, 
                        p_bound+1, Nyl, MPI_DOUBLE, north, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        else{
            for(int j=0;j<Nyl;j++){
                p_bound[j]=0;   // if nothing is north
            }
        }
        
        
        // ---------------West Transfers---------------
        if(west>=0){
            
            // Reset values to 0
            memset(temp_send,0,sizeof(double)*Nxl);
            memset(temp_recv,0,sizeof(double)*Nxl);
            
            for(int i=0;i<Nxl;i++){
                temp_send[i]=p_bound[Nyl+2+i*(Nyl+2)+1];     // Preparing west transfers
            }
        
            MPI_Sendrecv(temp_send, Nxl, MPI_DOUBLE, west, 2, 
                        temp_recv, Nxl, MPI_DOUBLE, west, 3, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            
            for(int i=0;i<Nxl;i++){
                p_bound[Nyl+2+i*(Nyl+2)]=temp_recv[i];     // receive from west
            }
            
        }
        else{
            for(int i=0;i<Nxl;i++){
                p_bound[Nyl+2+i*(Nyl+2)]=0;     // if nothing is west
            }
        }
        // ---------------East Transfers---------------
        if(east>=0){
            
            // Reset values to 0
            memset(temp_send,0,sizeof(double)*Nxl);
            memset(temp_recv,0,sizeof(double)*Nxl);
            
            for(int i=0;i<Nxl;i++){
                temp_send[i]=p_bound[Nyl+2+i*(Nyl+2)+Nyl];   // Preparing east transfers
            }
                    
            MPI_Sendrecv(temp_send, Nxl, MPI_DOUBLE, east, 3, 
                        temp_recv, Nxl, MPI_DOUBLE, east, 2, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

            for(int i=0;i<Nxl;i++){
                p_bound[Nyl+2+i*(Nyl+2)+Nyl+1]=temp_recv[i];   // receive from east
            }
        }
        else{
            for(int i=0;i<Nxl;i++){
                p_bound[Nyl+2+i*(Nyl+2)+Nyl+1]=0;       // if nothing is east
            }
        }
         
        // ---------------South Transfers---------------
        if(south>=0){
        
            MPI_Sendrecv(p_bound+(Nyl+2)*Nxl+1, Nyl, MPI_DOUBLE, south, 0, 
                        p_bound+(Nyl+2)*(Nxl+1)+1, Nyl, MPI_DOUBLE, south, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        
        }
        else{
            for(int j=0;j<Nyl;j++){
                p_bound[(Nyl+2)*(Nxl+1)+j+1]=0;     // if nothing is south
            }
        }
        
        //----------------END OF COMMUNICATION------------------//
        
        
        
        //-----------------Finding Alpha_k-----------------//
        
        // Finding alpha bottom
        // Parallel multiply A and s0
        
        // For non-padded cases
        if(cases==0){
        m=0, off=0; // Initialising local iterations
        
        // Finding sum
        for(int i=0;i<Nsd;i++){
            
            if(i%Nyl==0 && i!=0){
                off++;   // array offset
                m=0;
            }
            
            // Performing matrix vector multiplication and summing
            sum=0;
            sum=sum+A[0]*p_bound[m+1+off*spc];
            sum=sum+A[1]*p_bound[m+1+spc-1+off*spc];
            sum=sum+A[2]*p_bound[m+1+spc+off*spc];
            sum=sum+A[3]*p_bound[m+1+spc+1+off*spc];
            sum=sum+A[4]*p_bound[m+1+2*spc+off*spc]; 
            
            Apk[i]=sum;   //need to check if it works
            m++;
        }
        }
        // For padded cases
        else{
        
            //---------Top side padded rows---------------
            if(east<0 && south>0){
                
            // Finding sum
            for(int i=0;i<Nxl;i++){
                for(int j=0;j<Nyl-pady;j++){
                    
                    // Performing matrix vector multiplication and summing
                    sum=0;
                    sum=sum+A[0]*p_bound[j+1+i*spc];
                    sum=sum+A[1]*p_bound[j+1+spc-1+i*spc];
                    sum=sum+A[2]*p_bound[j+1+spc+i*spc];
                    sum=sum+A[3]*p_bound[j+1+spc+1+i*spc];
                    sum=sum+A[4]*p_bound[j+1+2*spc+i*spc]; 
                
                    Apk[i*Nyl+j]=sum; 
    
                }
            }
        
            }
            //----------Top right corner padded rows-------------
            else if(east<0 && south<0){
            
            // Finding sum
            for(int i=0;i<Nxl-padx;i++){
                for(int j=0;j<Nyl-pady;j++){
                    
                    // Performing matrix vector multiplication and summing
                    sum=0;
                    sum=sum+A[0]*p_bound[j+1+i*spc];
                    sum=sum+A[1]*p_bound[j+1+spc-1+i*spc];
                    sum=sum+A[2]*p_bound[j+1+spc+i*spc];
                    sum=sum+A[3]*p_bound[j+1+spc+1+i*spc];
                    sum=sum+A[4]*p_bound[j+1+2*spc+i*spc]; 
            
                    Apk[i*Nyl+j]=sum; 
                    
                }
            }
            
            }
            //----------Right side padded rows-------------
            else if(east>0 && south<0){
            
            // Finding sum
            for(int i=0;i<Nxl-padx;i++){
                for(int j=0;j<Nyl;j++){
                    
                    // Performing matrix vector multiplication and summing
                    sum=0;
                    sum=sum+A[0]*p_bound[j+1+i*spc];
                    sum=sum+A[1]*p_bound[j+1+spc-1+i*spc];
                    sum=sum+A[2]*p_bound[j+1+spc+i*spc];
                    sum=sum+A[3]*p_bound[j+1+spc+1+i*spc];
                    sum=sum+A[4]*p_bound[j+1+2*spc+i*spc]; 
            
                    Apk[i*Nyl+j]=sum; 
                    
                }
            }
            
            }
            else{
        
            m=0, off=0; // Initialising local iterations
        
            // Finding sum
            for(int i=0;i<Nsd;i++){
            
                if(i%Nyl==0 && i!=0){
                    off++;   // array offset
                    m=0;
                }
            
                // Performing matrix vector multiplication and summing
                sum=0;
                sum=sum+A[0]*p_bound[m+1+off*spc];
                sum=sum+A[1]*p_bound[m+1+spc-1+off*spc];
                sum=sum+A[2]*p_bound[m+1+spc+off*spc];
                sum=sum+A[3]*p_bound[m+1+spc+1+off*spc];
                sum=sum+A[4]*p_bound[m+1+2*spc+off*spc]; 
                
                Apk[i]=sum; 
                m++;
            }
            }
        }
        
        alpha_bottom=cblas_ddot(Nsd,p_k,1,Apk,1);   //Finding alpha bottom
        alpha_top=cblas_ddot(Nsd,r_k,1,z_k,1);  // Finding alpha top
        
        // Reduce to get global alpha
        
        MPI_Allreduce(&alpha_top,&alpha_top_r,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);  //alpha top
        MPI_Allreduce(&alpha_bottom,&alpha_bottom_r,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);  //alpha top
        
        alpha=alpha_top_r/alpha_bottom_r; // Dividing for global case
        
        //-------Finding LOCAL s_k+1----------
        
        cblas_daxpy(Nsd,alpha,p_k,1,s_k,1);  //Calculating s_k+1
    
        //-------Finding LOCAL r_k+1-----------
        
        cblas_dcopy(Nsd,r_k,1,r_kplus1,1);  //Copying values to k+1
        cblas_daxpy(Nsd,-alpha,Apk,1,r_kplus1,1);  //Calculating r_k+1
        
        //---------Calculating residual for error-----------
        
        //memset(r_kplus1_global,0,sizeof(double)*nsv);     // Reset values to 0
        //MPI_Allgather(r_kplus1,Nsd,MPI_DOUBLE,r_kplus1_global,Nsd,MPI_DOUBLE,MPI_COMM_WORLD);   // Getting global r_k+1
        res_local=cblas_dnrm2(Nsd,r_kplus1,1);  // Calculating local norms
        res_local=res_local*res_local;          //squaring
        MPI_Allreduce(&res_local,&res,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);  // Reducing to get global squares
        res=sqrt(res);  // sqrt to get global norm
       
        if(abs(res)<eps) break;  // BREAKING loop if error is below eps
        
        //----------Pre-conditioning for next iteration---------
        
        for(int i=0;i<Nsd;i++){
            z_kplus1[i]=r_kplus1[i]/A[2];
        }
        
        //----------Calculating beta-------------------
        
        beta_top=cblas_ddot(Nsd,z_kplus1,1,r_kplus1,1);   //Local calculations
        beta_bottom=cblas_ddot(Nsd,z_k,1,r_k,1);
        MPI_Allreduce(&beta_top,&beta_top_r,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);   //Reducing to get global beta
        MPI_Allreduce(&beta_bottom,&beta_bottom_r,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        beta=beta_top_r/beta_bottom_r;  //Dividing for global case
        
        cblas_dcopy(Nsd,z_kplus1,1,z_k,1);          // Storing z_k+1 into z_k for next iteration
        cblas_daxpy(Nsd,beta,p_k,1,z_kplus1,1);     // Calculating p_k+1 and storing in z_k+1
        cblas_dcopy(Nsd,z_kplus1,1,p_k,1);          // Storing p_k+1 into p_k for next iteration
        cblas_dcopy(Nsd,r_kplus1,1,r_k,1);          // Storing r_k+1 into r_k for next iteration
        
        k++;
        
    }while(k<5000);
    
    //--------------------END OF LOOP------------------------//
    
    // Wrapping up
    cblas_dcopy(Nsd,s_k,1,x,1);
    
    //----------------END of Conjugate Gradient Method--------------------//
    
    if (rank==0){
        cout<<"k: "<<k<<endl;       // Output no. of iterations
    }

    // Clearing memory
    delete[] r_k;    
    delete[] r_kplus1;
    delete[] p_k;
    delete[] p_bound;
    delete[] s_bound;
    delete[] s_k;
    delete[] Apk;
    delete[] z_k;
    delete[] z_kplus1;
    //delete[] r_kplus1_global;


}

// *Jacobi method (unused, only for reference)
void PoissonSolver::Jacobi(double* A,double* s, double* v){
    
    //----Firstly, organising processes as a Cartesian grid-----
    
    MPI_Comm mygrid;
    
    const int dims = 2;
    int sizes[dims] = {Px, Py};
    int periods[dims] = {0, 0};
    int reorder = 1;
    int rank_grid[dims];
    MPI_Cart_create(MPI_COMM_WORLD, dims, sizes, periods,
                    reorder, &mygrid);  //Creating grid
    MPI_Cart_coords(mygrid,rank,dims,rank_grid);    //Getting coordinates
    
    int north,south,east,west;
    MPI_Cart_shift(mygrid,0,1,&north,&south);
    MPI_Cart_shift(mygrid,1,1,&west,&east);
    
    //-------Declaring Jacobi parameters---------
    
    int mod, off=0, k=0, m=0;
    int nsv=Nxl*Px*Nyl*Py;  // Unknown variables for global
    const double eps=0.00001;  // Defining error constant
    double nrm_kplus1, nrm_k, err,sum;
    double* s_k=new double[(Nxl+2)*(Nyl+2)];      // stream at k iteration
    double* s_kplus1=new double[Nsd]; // stream at k plus 1 iteration
    double* s_kplus1_global;    
    double* s_k_global;
    if(rank==0){
        s_kplus1_global=new double[nsv];    //allocating memory for rank 0 only
        s_k_global=new double[nsv];
    }

    int kl=Nyl;                     // Lower diagonal
    int spc=Nyl+2;                  // Spacing for calculations
    double* temp_send=new double[Nxl](); //for west and east transfer
    double* temp_recv=new double[Nxl]();
    
    //Gather inital global s for error calculation later
    //MPI_Gather(&s_local,Nsd,MPI_DOUBLE,s_k_global,Nsd,MPI_DOUBLE,0,MPI_COMM_WORLD); 
    
    // Jacobi iteration
    
    //------------------------START OF LOOP---------------------------//
    do{
        
        /*---------Communicate to find Boundary points--------
         * --------once other processes are solved-----------
         * ---------simultaneously----------------------------*/
        
        // Looking at all directions for neighbouring points
        //**********CHECK FOR DEADLOCKKKKKK**************
        
        
        MPI_Barrier(MPI_COMM_WORLD);    // Synchronising processes

        
        // ---------------North Transfers---------------
        if(north>=0){
        
            MPI_Sendrecv(s_k+(Nyl+2)+1, Nyl, MPI_DOUBLE, north, 1, 
                        s_k+1, Nyl, MPI_DOUBLE, north, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        else{
            for(int j=0;j<Nyl;j++){
                s_k[j]=0;
            }
        }
        
        
        // ---------------West Transfers---------------
        if(west>=0){
            
            memset(temp_send,0,sizeof(double)*Nxl);
            memset(temp_recv,0,sizeof(double)*Nxl);
            
            for(int i=0;i<Nxl;i++){
                temp_send[i]=s_k[Nyl+2+i*(Nyl+2)+1];     // Preparing west transfers
            }
        
            MPI_Sendrecv(temp_send, Nxl, MPI_DOUBLE, west, 2, 
                        temp_recv, Nxl, MPI_DOUBLE, west, 3, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            
            for(int i=0;i<Nxl;i++){
                s_k[Nyl+2+i*(Nyl+2)]=temp_recv[i];     // receive from west
            }
            
        }
        else{
            for(int i=0;i<Nxl;i++){
                s_k[Nyl+2+i*(Nyl+2)]=0;
            }
        }
        // ---------------East Transfers---------------
        if(east>=0){
            
            memset(temp_send,0,sizeof(double)*Nxl);
            memset(temp_recv,0,sizeof(double)*Nxl);
            
            for(int i=0;i<Nxl;i++){
                temp_send[i]=s_k[Nyl+2+i*(Nyl+2)+Nyl];   // Preparing east transfers
            }
        
        
            MPI_Sendrecv(temp_send, Nxl, MPI_DOUBLE, east, 3, 
                        temp_recv, Nxl, MPI_DOUBLE, east, 2, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

            for(int i=0;i<Nxl;i++){
                s_k[Nyl+2+i*(Nyl+2)+Nyl+1]=temp_recv[i];   // receive from east
            }
        }
        else{
            for(int i=0;i<Nxl;i++){
                s_k[Nyl+2+i*(Nyl+2)+Nyl+1]=0;
            }
        }
         
        // ---------------South Transfers---------------
        if(south>=0){
        
            MPI_Sendrecv(s_k+(Nyl+2)*Nxl+1, Nyl, MPI_DOUBLE, south, 0, 
                        s_k+(Nyl+2)*(Nxl+1)+1, Nyl, MPI_DOUBLE, south, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        
        }
        else{
            for(int j=0;j<Nxl;j++){
                s_k[(Nyl+2)*(Nxl+1)+j+1]=0;
            }
        }
        
        //------------Non Padded rows case---------------
        if(cases==0){
        
        //Initialising local iterations
        m=0;
        off=0;
        
        // Finding sum
        for(int i=0;i<Nsd;i++){
            
            if(i%Nyl==0 && i!=0){
                off++;   // array offset
                m=0;
            }
                // Performing matrix vector multiplication and summing
                sum=0;
                sum=sum+A[0]*s_k[m+1+off*spc];
                sum=sum+A[1]*s_k[m+1+spc-1+off*spc];
                sum=sum+A[3]*s_k[m+1+spc+1+off*spc];
                sum=sum+A[4]*s_k[m+1+2*spc+off*spc]; 
                
                s_kplus1[i]=1.0/(A[2])*(v_local[i]-sum);   // Finding s at k plus 1
            m++;
       }
       
       }
       //------------Padded rows case---------------
       else{
           
            //---------Top side padded rows---------------
            if(east<0 && south>0){
        
            // Finding sum
            for(int i=0;i<Nxl;i++){
                for(int j=0;j<Nyl-pady;j++){
                    
                    // Performing matrix vector multiplication and summing
                    sum=0;
                    sum=sum+A[0]*s_k[j+1+i*spc];
                    sum=sum+A[1]*s_k[j+1+spc-1+i*spc];
                    sum=sum+A[3]*s_k[j+1+spc+1+i*spc];
                    sum=sum+A[4]*s_k[j+1+2*spc+i*spc]; 
                
                    s_kplus1[i*Nyl+j]=1.0/(A[2])*(v_local[i*Nyl+j]-sum);   // Finding s at k plus 1
                    
                }
            }
        
            }
            //---------Top right corner side padded rows---------------
            else if(east<0 && south<0){
        
            // Finding sum
            for(int i=0;i<Nxl-padx;i++){
                for(int j=0;j<Nyl-pady;j++){
                    
                    // Performing matrix vector multiplication and summing
                    sum=0;
                    sum=sum+A[0]*s_k[j+1+i*spc];
                    sum=sum+A[1]*s_k[j+1+spc-1+i*spc];
                    sum=sum+A[3]*s_k[j+1+spc+1+i*spc];
                    sum=sum+A[4]*s_k[j+1+2*spc+i*spc]; 
                
                    s_kplus1[i*Nyl+j]=1.0/(A[2])*(v_local[i*Nyl+j]-sum);   // Finding s at k plus 1
                    
                }
            }
        
            }
            //---------Right corner side padded rows---------------
            else if(east>0 && south<0){
        
            // Finding sum
            for(int i=0;i<Nxl-padx;i++){
                for(int j=0;j<Nyl;j++){
                    
                    // Performing matrix vector multiplication and summing
                    sum=0;
                    sum=sum+A[0]*s_k[j+1+i*spc];
                    sum=sum+A[1]*s_k[j+1+spc-1+i*spc];
                    sum=sum+A[3]*s_k[j+1+spc+1+i*spc];
                    sum=sum+A[4]*s_k[j+1+2*spc+i*spc]; 
                
                    s_kplus1[i*Nyl+j]=1.0/(A[2])*(v_local[i*Nyl+j]-sum);   // Finding s at k plus 1
                    
                }
            }
        
            }
            else{
                
            //Initialising local iterations
            m=0;
            off=0;
        
            // Finding sum
            for(int i=0;i<Nsd;i++){
            
                if(i%Nyl==0 && i!=0){
                    off++;   // array offset
                    m=0;
                }
                    // Performing matrix vector multiplication and summing
                    sum=0;
                    sum=sum+A[0]*s_k[m+1+off*spc];
                    sum=sum+A[1]*s_k[m+1+spc-1+off*spc];
                    sum=sum+A[3]*s_k[m+1+spc+1+off*spc];
                    sum=sum+A[4]*s_k[m+1+2*spc+off*spc]; 
                
                    s_kplus1[i]=1.0/(A[2])*(v_local[i]-sum);   // Finding s at k plus 1
                m++;
            }
       
            }
           
           
       }
    
        
        //Gather all points into global s in rank 0
        MPI_Gather(s_kplus1,Nsd,MPI_DOUBLE,s_kplus1_global,Nsd,MPI_DOUBLE,0,MPI_COMM_WORLD);   // for k plus 1
        //MPI_Allgather(s_kplus1,Nsd,MPI_DOUBLE,s_kplus1_global,Nsd,MPI_DOUBLE,MPI_COMM_WORLD);
        memset(s_k,0,sizeof(double)*(Nxl+2)*(Nyl+2));    //clearing values for next iteration


        if(rank==0){
            //Calculating Norm of both solutions
            nrm_kplus1=cblas_dnrm2(nsv,s_kplus1_global,1);     // Norm of current iteration.
            nrm_k=cblas_dnrm2(nsv,s_k_global,1);     // Norm of previous.
            err=nrm_kplus1-nrm_k;           // Calculating Error
        
            cblas_dcopy(nsv,s_kplus1_global,1,s_k_global,1);  // Copying global for next iteration
        }
        
        MPI_Bcast(&err,1,MPI_DOUBLE,0,MPI_COMM_WORLD);    //Broadcasting error term for other ranks
        
        // Placing back k+1 into k for next iteration
        for(int i=0;i<Nxl;i++){
            for(int j=0;j<Nyl;j++){
                s_k[Nyl+2+i*(Nyl+2)+j+1]=s_kplus1[i*Nyl+j];
            }
        }
        
        k++;    // for checking
        
    }while(abs(err)>=eps);
    //-------------------------END OF LOOP------------------------
    
    //Extracting interior points only
    for(int i=0;i<Nxl;i++){
        for(int j=0;j<Nyl;j++){
            s_local[i*Nyl+j]=s_k[Nyl+2+i*(Nyl+2)+j+1];
        }
    }

    
    if(rank==0){
        //printVec1(s_kplus1,Nsd);
        //printVec1(s_k,(Nxl+2)*(Nyl+2));
        cout<<"k: "<<k<<endl;
        //cout<<"Error: "<<err<<endl;

    }
    
    delete[] temp_recv;
    delete[] temp_send;
   
    delete[] s_k;
    delete[] s_kplus1;
    
    if(rank==0){
        delete[] s_kplus1_global;
        delete[] s_k_global;
    }
    
    
}
 
// Solving in serial using LAPACK functions
void PoissonSolver::SolveSerial(){
     
    // -------- Defining parameters ---------- //
    
    int nsv=(Nx-2)*(Ny-2); // No. of unknown rows
    int kl=Ny-2;    // Lower diagonal
    int ku=kl;      // Upper diagonal
    int ldh=1+2*kl+ku;  // Leading dimension of banded matrix
    double* A=new double[ldh*nsv]();  // Initialising A matrix
    int info;   // Info for lapack
    int* ipiv=new int[nsv](); // Pivot array for lapack
    double* s_nsv=new double[nsv]();  // Unknown stream functions
    
    // Getting coefficients for A matrix
    double a=-1.0/(dx*dx);
    double b=-1.0/(dy*dy);
    double c=(2.0/(dx*dx)+2.0/(dy*dy));
    
    // Assembling A matrix
    
    for (int i=0;i<nsv;i++){
        
        A[i*ldh+kl]=a;
        A[i*ldh+2*kl]=c;
        A[i*ldh+3*kl]=a; 
        
        int j=i%kl;
        int k=(i+1)%kl;
        
        if(j==0){
            
        A[i*ldh+2*kl-1]=0;
        A[i*ldh+2*kl+1]=b;

        }
        else if(k==0){
        A[i*ldh+2*kl-1]=b;
        A[i*ldh+2*kl+1]=0;           
            
        }
        else{
        A[i*ldh+2*kl-1]=b;
        A[i*ldh+2*kl+1]=b;           
            
        }
    }  

    cblas_dcopy(Nx*Ny,v,1,s,1); // Copying values from vorticity into stream function
    SetBoundaryConditions();    // Setting BCs
    
    //Eliminating known variables from stream function vector
        
    for (int i=1;i<Nx-1;i++){
        for(int j=1;j<Ny-1;j++){
            s_nsv[(i-1)*(Ny-2)+(j-1)]=s[i*Ny+j];
    
        }
    }
    
    //Factorising A matrix
    F77NAME(dgbtrf)(nsv,nsv,kl,ku,A,ldh,ipiv,&info);
    
    //Solving
    F77NAME(dgbtrs)('N',nsv,kl,ku,1,A,ldh,ipiv,s_nsv,nsv,&info);

    //Placing back solved stream functions into main vector s
    for (int i=1;i<Nx-1;i++){
        for(int j=1;j<Ny-1;j++){
            s[i*Ny+j]=s_nsv[(i-1)*(Ny-2)+(j-1)];
        }
    }
    
    // Clearing memory
    delete[] A;
    delete[] ipiv;
    delete[] s_nsv;

}

// Solving in parallel using CG method
void PoissonSolver::SolveParallel(){
    
    // ------Defining parameters------ //
    
    int nsv=Nxl*Px*Nyl*Py;  // No. of unknown rows including padded
    double A[5];  // Initialising A matrix
    double a=-1.0/(dx*dx);  // Initialising coefficients within matrix
    double b=-1.0/(dy*dy);
    double c=(2.0/(dx*dx)+2.0/(dy*dy));
    
    // -------------Assembling LOCAL A------------- //
    
    A[0] = a;
    A[1] = b;
    A[2] = c;
    A[3] = b;  
    A[4] = a;         
            
    // ------Solving unknown streams using Parallel Jacobi iteration------//

    //Jacobi(A,s_local,v_local);  //Using Parallel Jacobi
    ConjugateGradient(A,s_local,v_local);
    
    // ------Gather all results and placing into s in rank 0-------//
    
    double* temp=new double[nsv]();   // temporary storage
    
    //--------------COMM STARTS---------------------
    
    //MPI_Barrier(MPI_COMM_WORLD);    //Synchronising processes
    
    // Gathering all local points into rank 0
    MPI_Allgather(s_local,Nsd,MPI_DOUBLE,temp,Nsd,MPI_DOUBLE,MPI_COMM_WORLD);
    
    //--------------COMM ENDS---------------------
    
    // Rearranging to correct order in memory
    
        
    double* s_nsv=new double[(Nx-2)*(Ny-2)]();    // To store interior stream functions
    int k=0,l=0,m=0,p=0;

    for(int i=0;i<Nx-2;i++){
            
        l=0, m=0;
            
        if(i%Nxl==0 && i!=0){
            k++;
            p=0;
        }
        for(int j=0;j<Ny-2;j++){
                
            if(j%Nyl==0 && j!=0){
                l++;
                m=0;
                    
            }
            s_nsv[i*(Ny-2)+j]=temp[l*(Nsd)+m+k*(Nyl*Py*Nxl)+p*Nyl];
            m++;
        }
        p++;
    }
    
    // Inserting into global s (i.e. + boundaries)
    
    for(int i=1;i<Nx-1;i++){
        for(int j=1;j<Ny-1;j++){
            s[i*Ny+j]=s_nsv[(i-1)*(Ny-2)+j-1];
        }
    }
    
    // Clearing memory
    delete[] s_nsv;
    delete[] temp;
}
