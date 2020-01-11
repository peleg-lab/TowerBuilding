function Value = fitnessfunc(x, mins, maxs)
    y = x;
    y(x<mins)=mins(x<mins);
    y(x>maxs)=maxs(x>maxs);

    N = 1000;    

    TowerSimulation(y(1),y(2),y(3), 'NumAgents', N, 'TimeSteps', 50000);
    out = TowerProperties(folder, y(1), y(2), y(3));
    Value = (N-out.MaxNumAnts(end))/N
		+max([(14-out.MaxHeight(end))/14,0])
		+100*sum((x-y).^2);
end
