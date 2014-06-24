#ifndef SIMHPP
#define SIMHPP

#ifndef DEBUG
#define DEBUG 0
#endif //DEBUG
// Stephen Inglis, 2014.06.23
// sim.hpp, header file for the main simulation part of the Random 2D Ising model
// Will eventually have two replicas so we can calculate the entanglement entropy

#include <string>
#include <fstream>
#include <iostream>
#include <math.h>
#include "MersenneTwister.h"
#include "Eigen/Core"
#include "stat.hpp"

class Sim{
    private:
        // Placeholder
        Eigen::Matrix<double, Eigen::Dynamic, 1> spins; // All spins +- 1
        Eigen::Matrix<double, Eigen::Dynamic, 1> Jmat; // All elements +- 1
        void loadParams(); // Loads parameters from param.dat
        MTRand* rand; // MersenneTwister pseudorandom number generator
        int seed; // Seed for the random number generator
        int L; // Size of the lattice
        double P; // Probability of 
        double beta; // Inverse temperature
        int sweeps; // Number of MC sweeps per bin
        int bins; // Number of total bins
        int regionA; // Number of rows in regionA (for simplicity, we will only add one entire row at a time)
    public:
        Sim(); //Default constructor
};

Sim::Sim(){
    loadParams();
}

void Sim::loadParams(){
    std::string filename = "param.dat";
    std::string g; // Garbage string for going through param file
    std::fstream inFile(filename.c_str());
    inFile >> g >> L;
    inFile >> g >> P;
    inFile >> g >> beta;
    inFile >> g >> seed;
    inFile >> g >> sweeps;
    inFile >> g >> bins;
    inFile >> g >> regionA;
    if(DEBUG){
        std::cout << "L = " << L << std::endl;
        std::cout << "P = " << P << std::endl;
        std::cout << "beta = " << beta << std::endl;
        std::cout << "seed = " << seed << std::endl;
        std::cout << "sweeps = " << sweeps << std::endl;
        std::cout << "bins = " << bins << std::endl;
        std::cout << "regionA = " << regionA << std::endl;
    }
}

#endif //SIMHPP
