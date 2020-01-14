# TowerBuilding

## Overview
This repository provides the code used to study fire ant tower building, and its possible robotic implementation, in [Nave et al. (2019)](https://www.biorxiv.org/content/10.1101/864306v1). The primary function provided in this repository, `TowerSimulation.m`, provides a MATLAB simulation of agents climbing on top of one another to build towers.

### Packages required
 - Image Processing Toolbox
 - Statistics and Machine Learning Toolbox 

## Function usage
**`TowerSimulation(c, p_u, k_nl)`** - Runs a simulation and generates an output file based on input parameter names
### Required parameters
 - `c`: Attraction parameter, represents the ratio of the attraction force to randomness
 - `p_u`: Unlock probability, probability that an locked agent chooses to unlock
 - `k_nl`: Neighbor-influenced locking factor, increase in locking probability per adjacent neighbor
These are all discussed in detail in [Nave et al. (2019)](https://www.biorxiv.org/content/10.1101/864306v1) section 2.2.

### Optional parameters
 - `'NumAgents'`: Number of agents in the simulation, default: `1,000`
 - `'NumSteps'`: Number of simulation time steps, default: `20,000`
 - `'SideLength'`: Arena size is SideLength x SideLength, default: `100`
 - `'ClimbHeight'`: Maximum climb height of each agent, default: `1`
 - `'Psl'`: Probability that an agent spontaneously locks, default: `0.00005`
 - `'Periodic'`: Whether or not boundary conditions are periodic
 - `'OutputFolder'`: Location of output data, default: `'Data/Output'`
 - `'SaveFreq'`: Frequency of data output to save file, default: `250`
 - `'RNGseed'`: Seed value for random number generator, default: `'shuffle'`

### Output format
The output file names will be in the form `'cx.xx.pux.xxx.knlx.xxx.vx'` where `x`s are numbers corresponding to the file inputs.

Two output files will be generated. The first, ending in `_parameters.mat`, will contain all of the required and optional parameters used for the simulation. The second, ending in `_output.txt`, will be a comma-separated variable with the raw output data. It will have `NumAgents` columns and 5\*(1+`NumSteps`/`SaveFreq`) rows. On each saved frame, including the initialized data, 5 rows are saved: `X`, `Y`, `Z`, `locked`, and `covered`. Each entry in the `locked` and `covered` rows is either 1 or 0, signifying True or False.

**`TowerProperties(c, p_u, k_nl, v)`** - 
### Required parameters
 - `c`: Attraction parameter, represents the ratio of the attraction force to randomness
 - `p_u`: Unlock probability, probability that an locked agent chooses to unlock
 - `k_nl`: Neighbor-influenced locking factor, increase in locking probability per adjacent neighbor
These are all discussed in detail in [Nave et al. (2019)](https://www.biorxiv.org/content/10.1101/864306v1) section 2.2.
 - `v`: 

### Optional parameters
 - `'NumAgents'`: Number of agents in the simulation, default: `1,000`
 - `'NumSteps'`: Number of simulation time steps, default: `20,000`
 - `'SideLength'`: Arena size is SideLength x SideLength, default: `100`
 - `'ClimbHeight'`: Maximum climb height of each agent, default: `1`
 - `'Psl'`: Probability that an agent spontaneously locks, default: `0.00005`
 - `'Periodic'`: Whether or not boundary conditions are periodic
 - `'OutputFolder'`: Location of output data, default: `'Data/Output'`
 - `'SaveFreq'`: Frequency of data output to save file, default: `250`
 - `'RNGseed'`: Seed value for random number generator, default: `'shuffle'`

### Output format
The output file names will be in the form `'cx.xx.pux.xxx.knlx.xxx.vx'` where `x`s are numbers corresponding to the file inputs.

Two output files will be generated. The first, ending in `_parameters.mat`, will contain all of the required and optional parameters used for the simulation. The second, ending in `_output.txt`, will be a comma-separated variable with the raw output data. It will have `NumAgents` columns and 5\*(1+`NumSteps`/`SaveFreq`) rows. On each saved frame, including the initialized data, 5 rows are saved: `X`, `Y`, `Z`, `locked`, and `covered`. Each entry in the `locked` and `covered` rows is either 1 or 0, signifying True or False.

## References:
1. Gary K. Nave Jr., Nelson T. Mitchell, Jordan A. Chan Dick, Tyler Schuessler, Joaquin A. Lagarrigue, Orit Peleg. (2019) bioRxiv 864306; doi: https://doi.org/10.1101/864306.
2. Phonekeo, S., Mlot, N., Monaenkova, D., Hu, D. L., & Tovey, C. (2017). Fire ants perpetually rebuild sinking towers. Royal Society open science, 4(7), 170475. 
