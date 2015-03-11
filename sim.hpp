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

#define PI 3.1415926535897

// Some helper functions, for flow control

bool fexists(const char* filename){
    ifstream ifile(filename);
    return ifile;
}

void getTemps(vector<double>& temps){
    std::string filename = "temps.dat";
    std::fstream inFile(filename.c_str());
    temps.resize(0);
    double cur_temp;
    while (inFile >> cur_temp){
        temps.push_back(cur_temp);
    }
}

class Sim{
    private:
        static const int binSize = 1000; // Averages per line (larger --> smaller files)
        static const int bufferSize = 10; // Lines held per write (larger --> less often writes to file)
        Eigen::Matrix<int, Eigen::Dynamic, 2> cluster; // For the cluster update
        Eigen::Matrix<double, Eigen::Dynamic, 1> Jmat; // All elements +- 1
        Eigen::Matrix<int, Eigen::Dynamic, 2> spins; // All spins +- 1
        Eigen::Matrix<double, Eigen::Dynamic, 1> chiwave; // For calculating the correlation length
        Eigen::Matrix<double, Eigen::Dynamic, 1> chiwave2; // For calculating the correlation length
        Eigen::Matrix<double, 2, 2> trans_mat; // For the transfer matrix calculation
        Eigen::Matrix<double, 2, 2> trans_prod_top; // For the transfer matrix calculation
        Eigen::Matrix<double, 2, 2> trans_prod_bot; // For the transfer matrix calculation
        Eigen::Matrix<double, 2, 2> trans_prod_con; // For the transfer matrix calculation
        void loadParams(); // Loads parameters from param.dat
        MTRand* rand; // MersenneTwister pseudorandom number generator
        int seed; // Seed for the random number generator
        int L; // Size of the lattice
        int Nspins; // Derived from L above
        int Nbonds; // Derived from L above
        double P; // Probability of disorder bond, 0 = pure ferromagnet
        double beta; // Inverse temperature
        int regionA; // Number of rows in regionA (for simplicity, we will only add one entire row at a time)
        double Energy; // Total energy
        observable<double>* obs_mag2; // For binder cumulant
        observable<double>* obs_mag4;
        observable<double>* obs_F;
        observable<double>* obs_ratio;
        int reverse; // Variable for calculating the regions in reverse size, since we lack symmetry
    public:
        Sim(double _beta=0.1); //Default constructor
        Sim(MTRand* _rand, Eigen::Matrix<double, Eigen::Dynamic, 1> _Jmat, double _beta); //Parallel temperting constructor
        Sim& operator=(const Sim& rh);
        Eigen::Matrix<int, Eigen::Dynamic, 2> getSpins();
        void setSpins(Eigen::Matrix<int, Eigen::Dynamic, 2> _spins);
        int sweeps; // Number of MC sweeps per bin
        int bins; // Number of total bins
        void singleUpdate(); // Single spin update attempt
        int adjS(int x, int y); // Find the y'th neighbor to spin x
        int adjJ(int x, int y); // Find the y'th bond of spin x (it connect's with the y'th spin above)
        void addNeighbours(int z, int zz, int& s);
        void clusterUpdate(); // Swendsen-Wang update
        int getNspins();
        int getNbonds();
        void updateE(); // Updates the Energy variable above
        void updateBinder(); // Updates magnetization for Binder Ratio
        void updateRatio(); // Updates ratio measurement from the current regionA choice to the next one (adding a row)
        void saveJ();
        void loadJ();
        Eigen::Matrix<double, Eigen::Dynamic, 1> getJ();
        MTRand* getRand();
        double getE();
        void setE(double newE);
        double getB();
        void printSpins();
};

Sim& Sim::operator=(const Sim& rh){
    Energy = rh.Energy;
    spins = rh.spins;
    return *this;
}

int Sim::getNspins(){
    return Nspins;
}

int Sim::getNbonds(){
    return Nbonds;
}

Sim::Sim(double _beta){
    loadParams();
    beta = _beta;
    obs_mag2 = new observable<double>("mag2_"+std::to_string(beta),bufferSize,binSize,0);
    obs_mag4 = new observable<double>("mag4_"+std::to_string(beta),bufferSize,binSize,0);
    obs_F = new observable<double>("F_"+std::to_string(beta),bufferSize,binSize,0);
    obs_ratio = new observable<double>("ratio_"+std::to_string(beta),bufferSize,binSize,0);
    //Padd = 1. - exp(-2*beta); // Probability of adding a spin when they are satisfied for the cluster move
    rand = new MTRand(seed);
    Nspins = L*L;
    Nbonds = 2*Nspins;
    spins.resize(Nspins,2);
    cluster.resize(Nspins,2);
    chiwave.resize(Nspins,1);
    chiwave2.resize(Nspins,1);
    for(int i=0;i<Nspins;i++){
        spins(i,0) = rand->randInt(1)*2 - 1; // +- 1 random initial state
        if(reverse){
            if ((Nspins-1-i) < regionA*L){
                spins(i,1) = spins(i,0); // spins must match in region A
            }
            else{
                spins(i,1) = rand->randInt(1)*2 - 1; // +- 1 random initial state
            }
        }
        else{
            if (i < regionA*L){
                spins(i,1) = spins(i,0); // spins must match in region A
            }
            else{
                spins(i,1) = rand->randInt(1)*2 - 1; // +- 1 random initial state
            }
        }
        chiwave(i) = cos(2*PI/L*(i%L));
        chiwave2(i) = sin(2*PI/L*(i%L));
    }
    // Even though we allow real random variables, the default constructor will still do +-1 for the variables
    Jmat.resize(Nbonds);
    for(int i=0;i<Nbonds;i++){
        if(rand->randExc() < P){
            Jmat(i) = 1*0.75; // If P is small, this usually won't happen
        }
        else{
            Jmat(i) = -1*0.65; // Usually we are ferromagnetic
        }
    }
    if(DEBUG){
        std::cout << Jmat.transpose() << std::endl;
    }
}

Sim::Sim(MTRand* _rand, Eigen::Matrix<double, Eigen::Dynamic, 1> _Jmat, double _beta){
    loadParams();
    beta = _beta;
    obs_mag2 = new observable<double>("mag2_"+std::to_string(beta),bufferSize,binSize,0);
    obs_mag4 = new observable<double>("mag4_"+std::to_string(beta),bufferSize,binSize,0);
    obs_F = new observable<double>("F_"+std::to_string(beta),bufferSize,binSize,0);
    obs_ratio = new observable<double>("ratio_"+std::to_string(beta),bufferSize,binSize,0);
    //Padd = 1. - exp(-2*beta); // Probability of adding a spin when they are satisfied for the cluster move
    rand = _rand;
    Nspins = L*L;
    Nbonds = 2*Nspins;
    spins.resize(Nspins,2);
    cluster.resize(Nspins,2);
    chiwave.resize(Nspins,1);
    chiwave2.resize(Nspins,1);
    for(int i=0;i<Nspins;i++){
        spins(i,0) = rand->randInt(1)*2 - 1; // +- 1 random initial state
        if(reverse){
            if ((Nspins-1-i) < regionA*L){
                spins(i,1) = spins(i,0); // spins must match in region A
            }
            else{
                spins(i,1) = rand->randInt(1)*2 - 1; // +- 1 random initial state
            }
        }
        else{
            if (i < regionA*L){
                spins(i,1) = spins(i,0); // spins must match in region A
            }
            else{
                spins(i,1) = rand->randInt(1)*2 - 1; // +- 1 random initial state
            }
        }
        chiwave(i) = cos(2*PI/L*(i%L));
        chiwave2(i) = sin(2*PI/L*(i%L));
    }
    // Here we read in Jmat
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
    inFile >> g >> seed;
    inFile >> g >> sweeps;
    inFile >> g >> bins;
    inFile >> g >> regionA;
    inFile >> g >> reverse;
    if(DEBUG){
        std::cout << "L = " << L << std::endl;
        std::cout << "P = " << P << std::endl;
        std::cout << "beta = " << beta << std::endl;
        std::cout << "seed = " << seed << std::endl;
        std::cout << "sweeps = " << sweeps << std::endl;
        std::cout << "bins = " << bins << std::endl;
        std::cout << "regionA = " << regionA << std::endl;
        std::cout << "reverse = " << reverse << std::endl;
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
    if(reverse){
        if ((Nspins-1-z) < regionA*L){ // We are in region A
            for(int i=0;i<4;i++){ // Loop over neighbours
                field += spins(z,0) * spins(adjS(z,i),0) * Jmat(adjJ(z,i));
                field += spins(z,1) * spins(adjS(z,i),1) * Jmat(adjJ(z,i));
            }
            if (field>=0){
                spins(z,0) = spins(z,0) * -1;
                spins(z,1) = spins(z,1) * -1;
            }
            else{ // Probability to flip
                if(rand->randExc() < exp(2*field*beta)){ // field is negative, so the right hand term is less than 1
                    spins(z,0) = spins(z,0) * -1;
                    spins(z,1) = spins(z,1) * -1;
                }
            }
        }
        else{ // We are not in region A
            int zz = rand->randInt(1); // Choose a layer
            for(int i=0;i<4;i++){ // Loop over neighbours
                field += spins(z,zz) * spins(adjS(z,i),zz) * Jmat(adjJ(z,i));
            }
            if (field>=0){
                spins(z,zz) = spins(z,zz) * -1;
            }
            else{ // Probability to flip
                if(rand->randExc() < exp(2*field*beta)){ // field is negative, so the right hand term is less than 1
                    spins(z,zz) = spins(z,zz) * -1;
                }
            }
        }
    }
    else{
        if (z < regionA*L){ // We are in region A
            for(int i=0;i<4;i++){ // Loop over neighbours
                field += spins(z,0) * spins(adjS(z,i),0) * Jmat(adjJ(z,i));
                field += spins(z,1) * spins(adjS(z,i),1) * Jmat(adjJ(z,i));
            }
            if (field>=0){
                spins(z,0) = spins(z,0) * -1;
                spins(z,1) = spins(z,1) * -1;
            }
            else{ // Probability to flip
                if(rand->randExc() < exp(2*field*beta)){ // field is negative, so the right hand term is less than 1
                    spins(z,0) = spins(z,0) * -1;
                    spins(z,1) = spins(z,1) * -1;
                }
            }
        }
        else{ // We are not in region A
            int zz = rand->randInt(1); // Choose a layer
            for(int i=0;i<4;i++){ // Loop over neighbours
                field += spins(z,zz) * spins(adjS(z,i),zz) * Jmat(adjJ(z,i));
            }
            if (field>=0){
                spins(z,zz) = spins(z,zz) * -1;
            }
            else{ // Probability to flip
                if(rand->randExc() < exp(2*field*beta)){ // field is negative, so the right hand term is less than 1
                    spins(z,zz) = spins(z,zz) * -1;
                }
            }
        }
    }
}

void Sim::addNeighbours(int z,int zz,int& s){ // Recursive part
    s = s+1;
    cluster(z,zz) = -1;
    if(reverse){
        if ((Nspins-1-z) < regionA*L){
            if(cluster(z,(zz+1)%2)==1) addNeighbours(z,(zz+1)%2,s);
        }
    }
    else{
        if (z < regionA*L){
            if(cluster(z,(zz+1)%2)==1) addNeighbours(z,(zz+1)%2,s);
        }
    }
    for(int i=0;i<4;i++){
        if(cluster(adjS(z,i),zz)==1){ // Only try if it's not in the cluster
            double bondE = spins(z,zz) * spins(adjS(z,i),zz) * Jmat(adjJ(z,i));
            if(bondE <= 0){ // If the spins match ...
                if(rand->randExc() < (1-exp(2*bondE*beta))){ // Probability to add, bondE is negative so this prob is greater than zero
                    addNeighbours(adjS(z,i),zz,s);
                }
            }
        }
    }
}

void Sim::clusterUpdate(){
    // Choose a random spin
    int z = rand->randInt(Nspins-1);
    int zz = rand->randInt(1);
    cluster.fill(1); // No spins are in the cluster
    if(DEBUG){
        int s = 0;
        addNeighbours(z,zz,s);
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
        addNeighbours(z,zz,s);
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
            // For each bond, calculate contribution from each layer
            Energy += spins(i,0) * spins(adjS(i,b),0) * Jmat(adjJ(i,b));
            Energy += spins(i,1) * spins(adjS(i,b),1) * Jmat(adjJ(i,b));
        }
    }
}

void Sim::updateBinder(){
    double mag = spins.col(0).sum()*1.0/(Nspins);
    obs_mag2->pe(pow(mag,2));
    obs_mag4->pe(pow(mag,4));
    double F =  (chiwave.array() * (spins.col(0).array()).cast<double>()).sum();
    double F2 =  (chiwave2.array() * (spins.col(0).array()).cast<double>()).sum();
    obs_F->pe(1.0/Nspins * (pow(F,2) + pow(F2,2)));
}

void Sim::updateRatio(){
    double field_1_top = 0; // 1 and 2 are the first and second spin of the transfer matrix
    double field_2_top = 0; // While top and bot represent the two replicas
    double field_1_bot = 0;
    double field_2_bot = 0;
    trans_prod_top.setIdentity(2,2);
    trans_prod_bot.setIdentity(2,2);
    trans_prod_con.setIdentity(2,2);
    int s1,s2;
    if(reverse){
        // Spins we are integrating over are spin N-1-(regionA*L) to N-1-(regionA*L + (L-1))
        for(int i=0;i<L;i++){
            s1 = Nspins-1-regionA*L-i;
            if (i==(L-1)){
                s2 = Nspins-1-regionA*L;
            }
            else{
                s2 = s1-1;
            }
            if(i==0){
                field_1_top = Jmat(adjJ(s1,1)) * spins(adjS(s1,1),0) + Jmat(adjJ(s1,3)) * spins(adjS(s1,3),0);
                field_1_bot = Jmat(adjJ(s1,1)) * spins(adjS(s1,1),1) + Jmat(adjJ(s1,3)) * spins(adjS(s1,3),1);
            }
            else{
                field_1_top = field_2_top;
                field_1_bot = field_2_bot;
            }
            field_2_top = Jmat(adjJ(s2,1)) * spins(adjS(s2,1),0) + Jmat(adjJ(s2,3)) * spins(adjS(s2,3),0);
            field_2_bot = Jmat(adjJ(s2,1)) * spins(adjS(s2,1),1) + Jmat(adjJ(s2,3)) * spins(adjS(s2,3),1);

            trans_mat(0,0) = exp(-1.*beta*(Jmat(adjJ(s1,0)) + field_1_top/2. + field_2_top/2.));
            trans_mat(0,1) = exp(-1.*beta*(-1.*Jmat(adjJ(s1,0)) + field_1_top/2. - field_2_top/2.));
            trans_mat(1,0) = exp(-1.*beta*(-1.*Jmat(adjJ(s1,0)) - field_1_top/2. + field_2_top/2.));
            trans_mat(1,1) = exp(-1.*beta*(Jmat(adjJ(s1,0)) - field_1_top/2. - field_2_top/2.));
            trans_prod_top *= trans_mat;

            trans_mat(0,0) = exp(-1.*beta*(Jmat(adjJ(s1,0)) + field_1_bot/2. + field_2_bot/2.));
            trans_mat(0,1) = exp(-1.*beta*(-1.*Jmat(adjJ(s1,0)) + field_1_bot/2. - field_2_bot/2.));
            trans_mat(1,0) = exp(-1.*beta*(-1.*Jmat(adjJ(s1,0)) - field_1_bot/2. + field_2_bot/2.));
            trans_mat(1,1) = exp(-1.*beta*(Jmat(adjJ(s1,0)) - field_1_bot/2. - field_2_bot/2.));
            trans_prod_bot *= trans_mat;

            trans_mat(0,0) = exp(-1.*beta*(2.*Jmat(adjJ(s1,0)) + field_1_top/2. + field_1_bot/2. + field_2_top/2. + field_2_bot/2.));
            trans_mat(0,1) = exp(-1.*beta*(-2.*Jmat(adjJ(s1,0)) + field_1_top/2. + field_1_bot/2. - field_2_top/2. - field_2_bot/2.));
            trans_mat(1,0) = exp(-1.*beta*(-2.*Jmat(adjJ(s1,0)) - field_1_top/2. - field_1_bot/2. + field_2_top/2. + field_2_bot/2.));
            trans_mat(1,1) = exp(-1.*beta*(2.*Jmat(adjJ(s1,0)) - field_1_top/2. - field_1_bot/2. - field_2_top/2. - field_2_bot/2.));
            trans_prod_con *= trans_mat;
        }
    }
    else{
        // Spins we are integrating over are spin regionA*L to regionA*L + (L-1)
        for(int i=0;i<L;i++){
            s1 = regionA*L+i;
            if (i==(L-1)){
                s2 = regionA*L;
            }
            else{
                s2 = s1+1;
            }
            if(i==0){
                field_1_top = Jmat(adjJ(s1,1)) * spins(adjS(s1,1),0) + Jmat(adjJ(s1,3)) * spins(adjS(s1,3),0);
                field_1_bot = Jmat(adjJ(s1,1)) * spins(adjS(s1,1),1) + Jmat(adjJ(s1,3)) * spins(adjS(s1,3),1);
            }
            else{
                field_1_top = field_2_top;
                field_1_bot = field_2_bot;
            }
            field_2_top = Jmat(adjJ(s2,1)) * spins(adjS(s2,1),0) + Jmat(adjJ(s2,3)) * spins(adjS(s2,3),0);
            field_2_bot = Jmat(adjJ(s2,1)) * spins(adjS(s2,1),1) + Jmat(adjJ(s2,3)) * spins(adjS(s2,3),1);

            trans_mat(0,0) = exp(-1.*beta*(Jmat(adjJ(s1,0)) + field_1_top/2. + field_2_top/2.));
            trans_mat(0,1) = exp(-1.*beta*(-1.*Jmat(adjJ(s1,0)) + field_1_top/2. - field_2_top/2.));
            trans_mat(1,0) = exp(-1.*beta*(-1.*Jmat(adjJ(s1,0)) - field_1_top/2. + field_2_top/2.));
            trans_mat(1,1) = exp(-1.*beta*(Jmat(adjJ(s1,0)) - field_1_top/2. - field_2_top/2.));
            trans_prod_top *= trans_mat;

            trans_mat(0,0) = exp(-1.*beta*(Jmat(adjJ(s1,0)) + field_1_bot/2. + field_2_bot/2.));
            trans_mat(0,1) = exp(-1.*beta*(-1.*Jmat(adjJ(s1,0)) + field_1_bot/2. - field_2_bot/2.));
            trans_mat(1,0) = exp(-1.*beta*(-1.*Jmat(adjJ(s1,0)) - field_1_bot/2. + field_2_bot/2.));
            trans_mat(1,1) = exp(-1.*beta*(Jmat(adjJ(s1,0)) - field_1_bot/2. - field_2_bot/2.));
            trans_prod_bot *= trans_mat;

            trans_mat(0,0) = exp(-1.*beta*(2.*Jmat(adjJ(s1,0)) + field_1_top/2. + field_1_bot/2. + field_2_top/2. + field_2_bot/2.));
            trans_mat(0,1) = exp(-1.*beta*(-2.*Jmat(adjJ(s1,0)) + field_1_top/2. + field_1_bot/2. - field_2_top/2. - field_2_bot/2.));
            trans_mat(1,0) = exp(-1.*beta*(-2.*Jmat(adjJ(s1,0)) - field_1_top/2. - field_1_bot/2. + field_2_top/2. + field_2_bot/2.));
            trans_mat(1,1) = exp(-1.*beta*(2.*Jmat(adjJ(s1,0)) - field_1_top/2. - field_1_bot/2. - field_2_top/2. - field_2_bot/2.));
            trans_prod_con *= trans_mat;
        }
    }
    obs_ratio->pe(trans_prod_con.trace() / (trans_prod_top.trace() * trans_prod_bot.trace()));
}

void Sim::saveJ(){
    ofstream jfile;
    jfile.open("Jmat.dat");
    for(int i=0;i<Nbonds;i++){
        jfile << Jmat(i) << " ";
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

Eigen::Matrix<double, Eigen::Dynamic, 1> Sim::getJ(){
    return Jmat;
}

MTRand* Sim::getRand(){
    return rand;
}

double Sim::getE(){
    return Energy;
}

void Sim::setE(double newE){
    Energy = newE;
}

double Sim::getB(){
    return beta;
}

Eigen::Matrix<int, Eigen::Dynamic, 2> Sim::getSpins(){
    return spins;
}

void Sim::setSpins(Eigen::Matrix<int, Eigen::Dynamic, 2> _spins){
    spins = _spins;
}

void Sim::printSpins(){
    std::cout << spins.transpose() << std::endl;
}

#endif //SIMHPP
