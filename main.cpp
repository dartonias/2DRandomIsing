#include "sim.hpp"
#include <vector>
#include <cmath>

int main(int argc, char** argv){
    Sim* mySim = new Sim();
    for(int i_bins=0;i_bins<mySim->bins;i_bins++){
        for(int i_sweeps=0;i_sweeps<mySim->sweeps;i_sweeps++){
            for(int i_onesweep=0;i_onesweep<mySim->getNspins();i_onesweep++){
                mySim->singleUpdate();
            }
            mySim->clusterUpdate();
            mySim->updateE();
            mySim->updateBinder();
        }
    }
    return 0;
}
