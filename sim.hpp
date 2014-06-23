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
    public:
        Sim(); //Default constructor
};

Sim::Sim(){
}

void Sim::loadParams(){
}

#endif //SIMHPP
