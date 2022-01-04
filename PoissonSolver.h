#ifndef POISSONSOLVER_H
#define POISSONSOLVER_H

#include <string>
#include <iostream>

using namespace std;

// Declaring LAPACK functions to solve linear system of equations in serial
    
#define F77NAME(x) x##_
extern "C" {
        
        
    // LAPACK routine for solving systems of linear equations
    void F77NAME(dgbtrf)(const int& m,const int& n, const int& kl, const int& ku, double * A,
                         const int& lda, int * ipiv, int* info);
       
    void F77NAME(dgbtrs)(const int& trans,const int& n, const int& kl, const int& ku, const int& nrhs,
                         double * A, const int& lda, const int * ipiv, double* b, const int&ldb, int* info);
        
}

/**
 * @class PoissonSolver
 * @author Imran Ahmad Azhar
 * @date 25/03/20
 * @file PoissonSolver.h
 * @brief A class to solve the Poisson problem using LU Decomposition in serial
 * and Conjugate Gradient in Parallel
 * 
 */

class PoissonSolver{
    
public:

    // Constructors and destructors
    PoissonSolver();
    ~PoissonSolver();
    
    /**
     * @brief Functions to set up poisson problem
     * @param stream - public stream to extract from private
     */
    void SetGrid(int nx, int ny);
    void SetDeltas(double deltax,double deltay);
    void UpdateVorticity(double* vort);
    void UpdateStream(double* stream);
    void SetPartitions(int px, int py);
    void SetBoundaryConditions();
    void GetStream(double* stream);
    void GetRankSize(int rk, int sz);
    void GetCart(int n, int s, int e, int w);
    void LocalVorticity();
    void LocalStream();
    
    /**
     * @brief Function to print out vector (only for checking
     * @param v - vector array
     * @param size - size of vector
     */
    void printVec1(double v[], int size);   // -------for checking
    
    /**
     * @brief Function to solve linear system of equations in parallel
     * using Conjugate Gradient method
     * @param A - A matrix
     * @param s - stream function vector
     * @param v - vorticity vector
     */
    void ConjugateGradient(double A[],double s[],double v[]);
    
    /**
     * @brief Function to solve linear system of equation in parallel
     * using Gauss Jacobi method (not used, only for reference)
     * @param A - A matrix
     * @param s - stream function vector
     * @param v - vorticity vector
     */
    void Jacobi(double* A, double* s, double* v);

    /**
     * @brief Solving poisson problem in serial
     */
    void SolveSerial();
    
    /**
     * @brief Solving poisson problem in parallel
     */
    void SolveParallel();
    
    
private:
    
    double* v = nullptr;        /**< Vorticity */
    double* s = nullptr;        /**< Stream function */
    double* v_local = nullptr;  /**< Local Vorticity function */
    double* s_local = nullptr;  /**< Local Stream function */

    int    Nx;      /**< No. of x grid points */
    int    Ny;      /**< No. of y grid points */
    int    Nxl;     /**< No. of local x grid points */
    int    Nyl;     /**< No. of local y grid points */
    int    Nsd;     /**< Size of local domain */
    int    Px;
    int    Py;
    double dx;
    double dy;
    int    cases=0; // Cases: 0 if not padded, 1 if padded
    int    padx;    // Padded rows in the x direction
    int    pady;    // Padded rows in the y direction
    int    rank;    // Rank of process
    int    size;    // Total number of processes
    int    north;
    int    south;
    int    east;
    int    west;
    

};

#endif // POISSONSOLVER_H