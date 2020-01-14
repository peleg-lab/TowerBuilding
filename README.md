# TowerBuilding

## Overview
This repository provides the code used to study fire ant tower building, and its possible robotic implementation, in [Nave et al. (2019)](https://www.biorxiv.org/content/10.1101/864306v1). The primary function provided in this repository, `TowerSimulation.m`, provides a MATLAB simulation of agents climbing on top of one another to build towers.

### Packages required
 - Image Processing Toolbox
 - Statistics and Machine Learning Toolbox 

## Function usage
### `TowerSimulation(c, p_u, k_nl)` 
Runs a simulation and generates an output file based on the chosen input parameters.
#### Required parameters
 - `c`: Attraction parameter, represents the ratio of the attraction force to randomness
 - `p_u`: Unlock probability, probability that an locked agent chooses to unlock
 - `k_nl`: Neighbor-influenced locking factor, increase in locking probability per adjacent neighbor

These are all discussed in detail in [Nave et al. (2019)](https://www.biorxiv.org/content/10.1101/864306v1) section 2.2.

#### Optional parameters
 - `'NumAgents'`: Number of agents in the simulation, default: `1,000`
 - `'NumSteps'`: Number of simulation time steps, default: `20,000`
 - `'SideLength'`: Arena size is SideLength x SideLength, default: `100`
 - `'ClimbHeight'`: Maximum climb height of each agent, default: `1`
 - `'Psl'`: Probability that an agent spontaneously locks, default: `0.00005`
 - `'Periodic'`: Whether or not boundary conditions are periodic
 - `'OutputFolder'`: Location of output data, default: `'Data/Output'`
 - `'SaveFreq'`: Frequency of data output to save file, default: `250`
 - `'RNGseed'`: Seed value for random number generator, default: `'shuffle'`

#### Output format
The output file names will be in the form `'cX.XX.puX.XXX.knlX.XXX.vX'` where `X`s are numbers corresponding to the file inputs. `TowerSimulation` automatically increments trial number `vX` as more simulations are run.

Two output files will be generated. The first, ending in `_parameters.mat`, will contain all of the required and optional parameters used for the simulation. The second, ending in `_output.txt`, will be a comma-separated variable with the raw output data. It will have `NumAgents` columns and 5\*(1+`NumSteps`/`SaveFreq`) rows. On each saved frame, including the initialized data, 5 rows are saved: `X`, `Y`, `Z`, `locked`, and `covered`. Each entry in the `locked` and `covered` rows is either 1 or 0, signifying True or False.

### `TowerProperties(c, p_u, k_nl, v)`
Assesses the results of `TowerSimulation` to measure several properties of the towers in the simulation.
#### Required parameters
 - `c`: Attraction parameter, represents the ratio of the attraction force to randomness
 - `p_u`: Unlock probability, probability that an locked agent chooses to unlock
 - `k_nl`: Neighbor-influenced locking factor, increase in locking probability per adjacent neighbor

The above are all discussed in detail in [Nave et al. (2019)](https://www.biorxiv.org/content/10.1101/864306v1) section 2.2.

 - `v`: Simulation trial number. This is generated as a part of the output of `TowerSimulation`, and is used in the naming process.

#### Optional parameters
 - `'InputFolder'`: Location of input data, default: `'Data/Output'`
 - `'OutputFolder'`: Location of output data, default: `'Data/Output'`

#### Output format
One output file will be generated, `'cX.XX.puX.XXX.knlX.XXX.vX_output.mat'`, containing a MATLAB struct variable `out`, which is also returned via `out = TowerProperties(c, p_u, k_nl, v)`. 

`out` contains the following fields, each of which is a 1-dimensional array of that property at each saved time step:
 - `out.NumTowers`: the number of towers in the simulation
 - `out.AverageArea`: the average base area of all towers
 - `out.AverageDiameter`: the average effective diameter of all towers
 - `out.AverageHeight`: the average height of all towers
 - `out.AverageRatio`: the average height-diameter ratio of all towers
 - `out.AverageNumAnts`: the average number of agents of all towers
 - `out.MaxArea`: the base area of the largest tower
 - `out.MaxHeight`: the height of the largest tower
 - `out.MaxRatio`: the height-diameter ratio of the largest tower
 - `out.MaxNumAnts`: the number of agents in the largest tower
 
Note that in all of these properties, a tower is defined as any group of agents with a base area greater than 1, to remove free individuals from analysis. The largest tower is defined as the tower containing the most agents.

### `TowerVideo(c, p_u, k_nl, v)`
Visualizes the results of `TowerSimulation` in video format.
#### Required parameters
 - `c`: Attraction parameter, represents the ratio of the attraction force to randomness
 - `p_u`: Unlock probability, probability that an locked agent chooses to unlock
 - `k_nl`: Neighbor-influenced locking factor, increase in locking probability per adjacent neighbor

The above are all discussed in detail in [Nave et al. (2019)](https://www.biorxiv.org/content/10.1101/864306v1) section 2.2.

 - `v`: Simulation trial number. This is generated as a part of the output of `TowerSimulation`, and is used in the naming process.

#### Optional parameters
 - `'InputFolder'`: Location of input data, default: `'Data/Output'`
 - `'OutputFolder'`: Location of output data, default: `'Data/Output'`
 - `'FrameRate'`: The frame rate of the output video, default: `30`
 - `'Quality'`: The quality of the output video on a scale of 1-100, default: `75`
 - `'ImageFrequency'`: The frequency of time steps loaded into the video in terms of saved time steps, default: `1`
 
The `'ImageFrequency'` is in terms of saved time steps, and is therefore related to the `'SaveFreq'` input of `TowerSimulation`. If `TowerSimulation` saves once every 250 time steps and `'ImageFrequency'` is chosen to be 2, then each frame of the output video will be once every 500 simulation time steps. These will be displayed at the chosen `'FrameRate'`.

#### Output format
One output file will be generated, `'cX.XX.puX.XXX.knlX.XXX.vX.avi'`. This video will contain each plotted time step, displayed at the chosen frame rate.

## References:
1. Gary K. Nave Jr., Nelson T. Mitchell, Jordan A. Chan Dick, Tyler Schuessler, Joaquin A. Lagarrigue, Orit Peleg. (2019) bioRxiv 864306; doi: https://doi.org/10.1101/864306.
