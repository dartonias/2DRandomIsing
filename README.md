# 2DRandomIsing
Code for calculating the Renyi entanglement entropy in the 2D Ising model with arbitrary couplings.

## Parameter file
The way the code is currently structured is most clearly seen by looking at the parameter file, `param.dat`
```
Lattice_size_L 16
Disorder_prob_P 0.02
Random_seed 12345
Sweeps_per_bin 1000
Num_bins 100
Strips_in_regionA 0
reverse 0
```
The `Lattice_size_L` determines the size of the lattice, so it is `LxL` in the end.
The `Disorder_prob_P` determines what fraction of the bonds are antiferromagnetic.
This setting is overridden if the file `LOADJ` is present, in which case the program loads the couplings from `Jmat.dat` in the same directory as the executable.
The `Random_seed` seeds the Mersenne Twister random number generator.
A *sweep* in the program corresponds to

1. `Nspins` single spin updates
1. One cluster update
1. One measurement of the energy `updateE()`
1. One measurement of the ratio of partition functions `updateRatio()`
1. `num_temps` swap attempts for randomly chosen adjacent temperatures for parallel tempering

`Num_bins` is the number is the number of total bins, each one a single measurement in the file, and each bin is compried of `Sweeps_per_bin` sweeps and measurements.
`Strips_in_regionA` describes how many `Lx1` strips of spins are in region A.
If this is set to zero, the two system are disconnected, and `updateRatio()` is measuring `Z[A,2,T]/Z^2[T]` where `A` is the first strip of spins, i.e. the spins from index `0` to index `L-1`.
To measure the next strip (spins from index `L` to index `2L-1`) we only change the param file such that `Strips_in_regionA` is 1.
The largest sensible value for `Strips_in_regionA` is `L-1`, since in this case `updateRatio()` will measure the ratio of adding the last strip of spins.
If we tried setting `Strips_in_regionA` to `L`, this would cause the code to error unless `updateRatio()` was turned off.
Lastly, `reverse` determines whether we are building up region A from index `0` of the lattice or from index `L*L-1`.
Both are needed to calculate the mutual information.

## Basic run cycle of the code
The basic procedure for running the code (as I have been so far) is the following:

1. Create all the `Jmat.dat` files needed for simulation. If only +-1 values are needed, this can be done by running the code with the file `SAVEJ` present in the same location as the executable.
1. Choose the temperatures to simulate, and put these in `temps.dat`. `Temps.py` has an example script for creating evenly spaced beta points.
1. For each realization of disorder (or `Jmat.dat` file):
  1. Run a simulation for each valid value of `Strips_in_regionA`. These run from `0` to `L-1`
  1. For each value of `Strips_in_regionA`, we must run once with `reverse` set to `0`, and once with `reverse` set to `1`

This means, for each realization of disorder, we need to run the code 2L times, with each run internally containing `num_temps` simulations.
Each simulation will create data files according to the beta used for each simulation.
The data one is left with is samples of `Z[A',2,T]/Z[A,2,T]` where `A` is specific in the parameter file, and `A'` is region A plus one more strip.
> This value is stored directly, so in larger systems there may be problems with numerical overflow that require some slight modifications

To calculate the mutual information, one must calculate as usual

`I_2(A:B) = S_2(A) + S_2(B) - S_2(A \union B)`

`S_2(A) = \sum_{A_i=0}^{A_i+1=A} log[Z[A_i+1,2,T]/Z[A_i,2,T]]`

i.e., if we want to calculate `S_2(A)` for the half strip, we have to take the sum of the logs of `Z[A_i+1,2,T]/Z[A_i,2,T]` for the first `L/2` simulation.
A concrete example -- with `L=8`, this will be the product of the simulations with `Strips_in_regionA` set from 0..3 (so 4 total simulations).
To calculate `S_2(B)`, we use the complimentary number of simulations with `reverse` set to 1.
Again, for `L=8`, if we wanted to calculate `S_2(A)` for the quarter region, we would use the firt 2 simulations with `reverse` set to 0 (this would be 0 and 1) and the first 6 simulations with `reverse` set to 1 (0..5)
The final value `S_2(A \union B)` can be calculated from the simulations with `reverse` set to either 0 or 1, it makes no difference, and in my processing I usually symmetrize over both methods to reduce error.

I usually organize runs by: `REVERSE/DISORDER/REGION_A`
but any reordering of that should work just as well.

## Getting started

Using the following scripts to test the code.
Please note, the scripts are organized for use on a cluster, that means jobs are submitted through `sqsub` and we tar the final results to keep things organized (plus, some clusters have file number restrictions, and this circumvents that problem).
It is also set up in such a way that the number of jobs submitted to the cluster simultaneously is never too many.
After crashing a cluster once with 80,000 jobs, such measures are necessary.

* Compile the code (with make) -- you should get an executable called `Ising.out`.
* Move into to the directory `sharcnet_test`.
* Copy `Ising,out` into the this folder.
* Type `python makedirs.py` -- this will create two folders based on the size inside the file `SIZE`.
* Go into the two folders (called `L#/` and `L#r/`) and run the command `python submit_sharcnet.py` -- this will sbumit the actual jobs, and id currently hardcoded to only submit 200 at a time. Running the command again will submit the next 200 jobs, and so on and so on.
* When the jobs are finished, run the command `. ./tar_stuff` in each of the `L` folders -- this will create `data.tar` containing all the simulation data for that set.
* In the folder `sharcnet_test` run the command `python analyze_MI.py` -- this will create `MI_L#.dat` and `MIE_L#.dat` containing the mutual information and error. Each line is a temperature (it matches up with `temps.dat`) and each column is a different size of region A. The middle column is for the half cut.
* Finally, running `python plot_MI.py` will produce a crude plot.

Currently the files are setup to do a quick test of L=4 with only 16 realizations of disorder.
The system size can be changed by modifing the number in the file `SIZE` (don't change the first part in the file) and the number of disorders must be chnged in the following files

* `analyze_MI.py`
* `genJs.py`
* `submit_sharcnet.py`
