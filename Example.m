
c = 2.0;
p_u = 0.200;
k_nl = 1/26;

TowerSimulation(c, p_u, k_nl, 'NumAgents', 400, 'NumSteps', 10000, 'SideLength', 50, 'RNGseed', 1);
TowerProperties('Data/Output', c, p_u, k_nl, 1)
TowerVideo('Data/Output', c, p_u, k_nl, 1)
