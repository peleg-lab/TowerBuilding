% Choose relevant parameters:
c = 2.0;     % c defines the ratio of attractive force to random ness
p_u = 0.200; % p_u is the probability that a locked agent unlocks at each time step 
k_nl = 1/26; % k_nl is the probability of locking per neighbor

% Run the tower simulation code, choosing a few input parameters:
TowerSimulation(c, p_u, k_nl, 'NumAgents', 400, 'NumSteps', 10000, 'SideLength', 50, 'RNGseed', 1);

v = 1;       % Only one simulation has been run, so we are interested in version 1

% Analyze the tower properties in the resulting simulation:
TowerProperties(c, p_u, k_nl, v, 'InputFolder', 'Data/Output')

% Generate a video of the simulation results
TowerVideo(c, p_u, k_nl, v, 'InputFolder', 'Data/Output', 'OutputFolder', 'Data/Videos')
