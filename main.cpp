#include "sim.hpp"
#include <vector>
#include <cmath>

int main(int argc, char** argv){
    Eigen::Matrix<int, Eigen::Dynamic, 1> temp_spins;
    double beta_min;
    double beta_max;
    int num_temps;
    runParams(beta_min,beta_max,num_temps);
    vector<Sim> simVec;
    simVec.resize(0);
    simVec.push_back(Sim());
    if(fexists("LOADJ")){
        simVec[0].loadJ();
    }
    if(fexists("SAVEJ")){
        simVec[0].saveJ();
        return 0;
    }
    for(int i=0;i<(num_temps-1);i++){
        simVec.push_back(Sim(simVec[0].getRand(),simVec[0].getJ(),(beta_max-beta_min)/(num_temps-1)*(i+1)+beta_min));
    }
    for(int i_bins=0;i_bins<simVec[0].bins;i_bins++){
        for(int i_sweeps=0;i_sweeps<simVec[0].sweeps;i_sweeps++){
            for(int z=0;z<num_temps;z++){
                for(int i_onesweep=0;i_onesweep<simVec[0].getNspins();i_onesweep++){
                    simVec[z].singleUpdate();
                }
                simVec[z].clusterUpdate();
                simVec[z].updateBinder();
                simVec[z].updateE();
            }
            // Parallel tempering step
            for(int z=0;z<num_temps;z++){
                int flipme = simVec[0].getRand()->randInt(num_temps-2);
                double dE = simVec[flipme+1].getE() - simVec[flipme].getE();
                double dB = simVec[flipme+1].getB() - simVec[flipme].getB();
                if(simVec[0].getRand()->randExc() < exp(dE*dB)){
                    temp_spins = simVec[flipme].getSpins();
                    simVec[flipme].setSpins(simVec[flipme+1].getSpins());
                    simVec[flipme+1].setSpins(temp_spins);
                }
            }
        }
    }
    return 0;
}
