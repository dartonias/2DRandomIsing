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
#include <iomanip>
#include <math.h>
#include "MersenneTwister.h"
#include "Eigen/Core"
#include "stat.hpp"
#include <assert.h>

class Sim{
    private:
        static const int binSize = 1000; // Averages per line (larger --> smaller files)
        static const int bufferSize = 10; // Lines held per write (larger --> less often writes to file)
        Eigen::Matrix<int, Eigen::Dynamic, 1> cluster; // For the cluster update
        Eigen::Matrix<int, Eigen::Dynamic, 1> Jmat; // All elements +- 1
        Eigen::Matrix<int, Eigen::Dynamic, 1> spins; // All spins +- 1
        void loadParams(); // Loads parameters from param.dat
        MTRand* rand; // MersenneTwister pseudorandom number generator
        int seed; // Seed for the random number generator
        int L; // Size of the lattice
        int Nspins; // Derived from L above
        int Nbonds; // Derived from L above
        double P; // Probability of disorder bond, 0 = pure ferromagnet
        double beta; // Inverse temperature
        int regionA; // Number of rows in regionA (for simplicity, we will only add one entire row at a time)
        double Padd; // Probability of adding to a cluster
        double Energy; // Total energy
        observable<double>* obs_mag2; // For binder cumulant
        observable<double>* obs_mag4;
    public:
        Sim(); //Default constructor
        Sim(MTRand* _rand, Eigen::Matrix<int, Eigen::Dynamic, 1> _Jmat, double _beta); //Parallel temperting constructor
        Eigen::Matrix<int, Eigen::Dynamic, 1> getSpins();
        void setSpins(Eigen::Matrix<int, Eigen::Dynamic, 1> _spins);
        int sweeps; // Number of MC sweeps per bin
        int bins; // Number of total bins
        void singleUpdate(); // Single spin update attempt
        int adjS(int x, int y); // Find the y'th neighbor to spin x
        int adjJ(int x, int y); // Find the y'th bond of spin x (it connect's with the y'th spin above)
        void addNeighbours(int z, int& s);
        void clusterUpdate(); // Swendsen-Wang update
        int getNspins();
        int getNbonds();
        void updateE(); // Updates the Energy variable above
        void updateBinder(); // Updates magnetization for Binder Ratio
        void saveJ();
        void loadJ();
        Eigen::Matrix<int, Eigen::Dynamic, 1> getJ();
        MTRand* getRand();
        double getE();
        double getB();
};

int Sim::getNspins(){
    return Nspins;
}

int Sim::getNbonds(){
    return Nbonds;
}

Sim::Sim(){
    loadParams();
    obs_mag2 = new observable<double>("mag2_"+std::to_string(beta),bufferSize,binSize,0);
    obs_mag4 = new observable<double>("mag4_"+std::to_string(beta),bufferSize,binSize,0);
    Padd = 1. - exp(-2*beta); // Probability of adding a spin when they are satisfied for the cluster move
    rand = new MTRand(seed);
    Nspins = L*L;
    Nbonds = 2*Nspins;
    spins.resize(Nspins);
    cluster.resize(Nspins);
    for(int i=0;i<Nspins;i++){
        spins(i) = rand->randInt(1)*2 - 1; // +- 1 random initial state
    }
    Jmat.resize(Nbonds);
    for(int i=0;i<Nbonds;i++){
        if(rand->randExc() < P){
            Jmat(i) = 1; // If P is small, this usually won't happen
        }
        else{
            Jmat(i) = -1; // Usually we are ferromagnetic
        }
    }
    if(DEBUG){
        std::cout << Jmat.transpose() << std::endl;
    }
}

Sim::Sim(MTRand* _rand, Eigen::Matrix<int, Eigen::Dynamic, 1> _Jmat, double _beta){
    loadParams();
    beta = _beta;
    obs_mag2 = new observable<double>("mag2_"+std::to_string(beta),bufferSize,binSize,0);
    obs_mag4 = new observable<double>("mag4_"+std::to_string(beta),bufferSize,binSize,0);
    Padd = 1. - exp(-2*beta); // Probability of adding a spin when they are satisfied for the cluster move
    rand = _rand;
    Nspins = L*L;
    Nbonds = 2*Nspins;
    spins.resize(Nspins);
    cluster.resize(Nspins);
    for(int i=0;i<Nspins;i++){
        spins(i) = rand->randInt(1)*2 - 1; // +- 1 random initial state
    }
    Jmat.resize(Nbonds);
    Jmat = _Jmat;
    if(DEBUG){
        std::cout << Jmat.transpose() << std::endl;
    }
}

void Sim::loadParams(){
    std::string filename = "param.dat";
    std::string g; // Garbage string for going through param file
    std::fstream inFile(filename.c_str());
    inFile >> g >> L;
    inFile >> g >> P;
    inFile >> g >> beta; //beta_min
    inFile >> g >> g; //beta_max
    inFile >> g >> g; //num_temps
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

int Sim::adjS(int x, int y){
    // We're finding the y'th neighbor to spin x
    assert((y>=0)&&(y<4));
    assert((x>=0)&&(x<Nspins));
    if (y==0){ // Spin to the right
        if ((x%L)==(L-1)){ // We are a boundary spin on the right
            return x+1-L;
        }
        return x+1;
    }
    else if(y==1){ // Spin above
        if ((x/L)==(L-1)){ // We are a boundary spin at the top
            return x+L-Nspins;
        }
        return x+L;
    }
    else if(y==2){ // Spin to the left
        if ((x%L)==(0)){ // We are a boundary spin on the left
            return x-1+L;
        }
        return x-1;
    }
    else if(y==3){ // Spin below
        if ((x/L)==(0)){ // We are a boundary spin on the bottom
            return x-L+Nspins;
        }
        return x-L;
    }
}

int Sim::adjJ(int x, int y){
    // We're finding the y'th bond to spin x
    assert((y>=0)&&(y<4));
    assert((x>=0)&&(x<Nspins));
    if (y==0){ // Spin to the right
        return 2*x+0;
    }
    else if(y==1){ // Spin above
        return 2*x+1;
    }
    else if(y==2){ // Spin to the left
        if ((x%L)==(0)){ // We are a boundary spin on the left
            return 2*(x-1+L)+0;
        }
        return 2*(x-1)+0;
    }
    else if(y==3){ // Spin below
        if ((x/L)==(0)){ // We are a boundary spin on the bottom
            return 2*(x-L+Nspins)+1;
        }
        return 2*(x-L)+1;
    }
}

void Sim::singleUpdate(){
    // Choose a random spin
    int z = rand->randInt(Nspins-1);
    double field = 0;
    for(int i=0;i<4;i++){ // Loop over neighbours
        field += spins(z) * spins(adjS(z,i)) * Jmat(adjJ(z,i));
    }
    if (field>=0) spins(z) = spins(z) * -1;
    else{ // Probability to flip
        if(rand->randExc() < exp(2*field*beta)){ // field is negative, so the right hand term is less than 1
            spins(z) = spins(z) * -1;
        }
    }
}

void Sim::addNeighbours(int z,int& s){ // Recursive part
    s = s+1;
    cluster(z) = -1;
    for(int i=0;i<4;i++){
        if(cluster(adjS(z,i))==1){ // Only try if it's not in the cluster
            if((spins(z) * spins(adjS(z,i)) * Jmat(adjJ(z,i)))==-1){ // If the spins match ...
                if(rand->randExc() < Padd){ // Probability to add
                    addNeighbours(adjS(z,i),s);
                }
            }
        }
    }
}

void Sim::clusterUpdate(){
    // Choose a random spin
    int z = rand->randInt(Nspins-1);
    cluster.fill(1); // No spins are in the cluster
    if(DEBUG){
        int s = 0;
        addNeighbours(z,s);
        std::cout << s << std::endl;
        for(int i=0;i<Nspins;i++){
            if((i%L)==0) std::cout << std::endl;
            std::cout << std::setw(2);
            std::cout << cluster(i) << " ";
        }
        std::cout << std::endl;
    }
    else{
        int s = 0;
        addNeighbours(z,s);
    }
    // Now that it's finished, flip all spins in the cluster
    // All elements in "cluster" are -1, and will flip spins in "spins"
    if(DEBUG){
        updateE();
        double tE1 = Energy;
        spins = spins.array() * cluster.array();
        updateE();
        double tE2 = Energy;
        std::cout << "dE = " << (tE2 - tE1) << std::endl << std::endl;
    }
    else spins = spins.array() * cluster.array();
}

void Sim::updateE(){
    Energy = 0;
    // Loop over spins
    for(int i=0;i<Nspins;i++){
        // Loop over unique bonds, 2 per spin
        for(int b=0;b<2;b++){
            Energy += spins(i) * spins(adjS(i,b)) * Jmat(adjJ(i,b));
        }
    }
}

void Sim::updateBinder(){
    double mag = spins.sum()*1.0/Nspins;
    obs_mag2->pe(pow(mag,2));
    obs_mag4->pe(pow(mag,4));
}

void Sim::saveJ(){
    ofstream jfile;
    jfile.open("Jmat.dat");
    for(int i=0;i<Nbonds;i++){
        jfile << Jmat(i);
    }
    jfile.close();
}

void Sim::loadJ(){
    ifstream jfile;
    jfile.open("Jmat.dat");
    for(int i=0;i<Nbonds;i++){
        jfile >> Jmat(i);
    }
    jfile.close();
}

Eigen::Matrix<int, Eigen::Dynamic, 1> Sim::getJ(){
    return Jmat;
}

MTRand* Sim::getRand(){
    return rand;
}

double Sim::getE(){
    return Energy;
}

double Sim::getB(){
    return beta;
}

Eigen::Matrix<int, Eigen::Dynamic, 1> Sim::getSpins(){
    return spins;
}

void Sim::setSpins(Eigen::Matrix<int, Eigen::Dynamic, 1> _spins){
    spins = _spins;
}

// Some helper functions, for flow control

bool fexists(const char* filename){
    ifstream ifile(filename);
    return ifile;
}

void runParams(double& beta_min, double& beta_max, int& num_temps){
    std::string filename = "param.dat";
    std::string g; // Garbage string for going through param file
    std::fstream inFile(filename.c_str());
    inFile >> g >> g;
    inFile >> g >> g;
    inFile >> g >> beta_min; //beta_min
    inFile >> g >> beta_max; //beta_max
    inFile >> g >> num_temps; //num_temps
}

#endif //SIMHPP
