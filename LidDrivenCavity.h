#pragma once

#include <string>
using namespace std;

/**
 * @class LidDrivenCavity
 * @author Imran Ahmad Azhar
 * @date 25/03/20
 * @file LidDrivenCavity.h
 * @brief A class to solve lid driven cavity problem using navier-stokes
 */

class LidDrivenCavity
{
public:
    
    LidDrivenCavity();
    ~LidDrivenCavity();

    /**
     * @brief Functions to set up Lid Driven Cavity
     * @param rk - Rank of process
     * @param sz - Total number of processes
     */
    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
    void SetPartitions(int px, int py);
    void Initialise(int rk, int sz);
    
    //MPI Functions
    void LocalVorticity();
    void LocalStream();

    
    void BoundaryVorticity();
    void InteriorVorticity();
    
    /**
     * @brief Integrating lid driven cavity problem
     * if size=1: Solve in serial
     * if size>1: Solve in parallel
     */
    void Integrate();
    void printv();
    void prints();
    void OutputFile();

private:
    
    double* v = nullptr;    /**< Vorticity */
    double* s = nullptr;    /**< Stream function */

    double dt;
    double T;
    int    Nx;
    int    Ny;
    double Lx;
    double Ly;
    double Re;
    
    // New variables
    double dx; 
    double dy;
    double* v_nplus1 = nullptr;
    double* s_nplus1 = nullptr;
    double U = 1.0;
    
    // For MPI
    int Px;
    int Py;
    int Nxl;
    int Nyl;
    int rank;
    int size;
    int north;      /**< rank of process to the north */
    int south;      /**< rank of process to the south */
    int east;       /**< rank of process to the east */
    int west;       /**< rank of process to the west */
    
};