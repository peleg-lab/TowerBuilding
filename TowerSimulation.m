%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Agent-based modeling code for tower-building            %
% inspired by fire ants by:                               %
% Gary Nave, Postdoc, University of Colorado Boulder      %
% Updated Jan. 9, 2020                                    %

% Initial agent-based code template by:                   %
% Kirstin Petersen, Asst. Prof. ECE, Cornell University   %
% Homework for SYSEN 6000 May 5th 2017                    %
% Updated May 17th 2017                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TowerSimulation(c,p_u,k_nl,varargin)
global PERIODIC L


%%% Required input parameters
% c     %Ratio of attraction to randomness
% p_u   %probability agent spontaneously unlocks position
% k_nl  %probability agent locks position 
        %  depending on num. of locked neighbors

%%% Optional parameters, parsed using inputParser
% NumAgents    % Number of agents in the simulation
               % default: 1,000
% NumSteps     % Number of simulation time steps
               % default: 20,000
% SideLength   % Arena size is SideLength x SideLength
               % default: 100
% ClimbHeight  % Max climb height of each agent
               % default: 1
% Psl          % Probability of spontaneous locking
               % default: 0.00005
% Periodic     % Whether or not to use periodic boundary conditions
               % Use 1 for periodic or 0 for free
               % default: 1
% OutputFolder % Folder to save output data
               % default: 'Output'
% SaveFreq     % Frequency of saving to output file, once every SaveFreq frames
               % default: 250
% RNGseed      % Seed value for random number generator.
               % Can be replaced with a scalar when testing
               % default: 'shuffle'

parser = inputParser;

addParameter(parser, 'NumAgents', 1000); 
addParameter(parser, 'NumSteps', 20000); 
addParameter(parser, 'SideLength', 100); 
addParameter(parser, 'ClimbHeight', 1); 
addParameter(parser, 'Psl', 0.00005); 
addParameter(parser, 'Periodic', 1); 
addParameter(parser, 'OutputFolder', 'Data/Output'); 
addParameter(parser, 'SaveFreq', 250);
addParameter(parser, 'RNGseed', 'shuffle'); 

parse(parser, varargin{:});

%%% Internal parameters
N = parser.Results.NumAgents;           % No. of agents, default 1,000
frames = parser.Results.NumSteps;       % No. of time steps, default 500,000
L = parser.Results.SideLength;          % Axis limits, default 100
max_climb = parser.Results.ClimbHeight; % Max height of climbing agents, default 1
p_s_l = parser.Results.Psl;             % Probability agent spontaneously locks
% default, 0.00005

% 1: periodic boundary condition, 0: unlimited
% default, 1
PERIODIC = parser.Results.Periodic; 

%% Seed the random number generator
s = rng(parser.Results.RNGseed); % Initialize the random number generator
% default, 'shuffle'

% List of distances to 8 neighboring pixels for any given pixel.
spacelist = [-1, -1; 
             -1,  0; 
             -1,  1; 
              0, -1; 
              0,  1; 
              1, -1; 
              1,  0; 
              1,  1];

%%% Output file parameters
% Subfolder to output data to:
folder = parser.Results.OutputFolder;

% Make sure that the folder exists 
if ~exist(folder)
    mkdir(folder))
end

% Generate output filename
savename = sprintf('c%4.2f.pu%5.3f.knl%5.3f',c,p_u,k_nl);

% Frequency of file output, once every save_freq time steps
save_freq = parser.Results.SaveFreq;

% Generate version number by checking existing files
for i=1:25
    temp = strcat(savename,sprintf('.v%1d', i)); % Append version number
    % Check for existing filename
    fnamecheck = strcat(folder,'/',temp,'_parameters.mat');
    % Stop at the first version number that does not exist
    if exist(fnamecheck,'file')==0 
        savename = temp; 
        break
    end
end

% Save .mat file with existing parameters
save(strcat(folder,'/',savename,'_parameters.mat'))

%%% Initialize variables
p = L*rand(2,N)-L/2; % List of all individual (x,y) positions, randomized under uniform distribution from -L/2 to L/2
indlist = floor(p+L/2+1); % Positions as indices 1 through L in the lattice, rounding down to the nearest whole number
p = indlist-0.5-L/2; % Convert back to (x, y) positions, with each now at the center of a pixel

h_level = zeros(1,N); % List of all vertical positions, initialized to 0
locked=zeros(1,N); % List of locked states (1=locked, 0=unlocked), all start unlocked
covered=zeros(1,N); % List of covered state (1=covered, 0=uncovered), all start uncovered
h_map = zeros(L, L); % Initial height map

% Save data: for each frame, 5 rows are written: x, y, z, locked, covered
p_frame = [p;h_level;locked;covered];
dlmwrite(strcat(folder,'/',savename,'_output.txt'),p_frame,'-append')

% Calculate inter-particle distances
d2 = pdist([p(1,:)' p(2,:)' h_level'],@distfun3d);
D2 = squareform(d2); %Matrix representation for the distance between particles
% Generate lists of adjacent neighbors
[l1,l2]=find(0<D2 & D2<=1);


%%% Main loop
for k=1:frames
    %% 1. Initialize data vectors
    v = zeros(2,N); % List of all velocities, initialized to 0
    locked_neigh=zeros(1,N); % Number of locked neighbors for each individual
    free_list = find(locked==0); % List of free agents
    
    %% 2. Calculate motion
    if ~isempty(free_list) % Only execute if there are free agents

        %% First, we generate a map of locations for all free individuals, prevents to prevent position conflicts, step_map
        step_map = zeros(L,L); 
        indlist = floor(p+L/2+1); % Indexed locations of all positions, p
        for fi = 1:length(free_list) % Iterate through free agents
            i = free_list(fi); % Index of current agent
            ind = indlist(:,i); % Indexed position of current agent
            step_map(ind(1),ind(2)) = step_map(ind(1),ind(2))+1; % Add one to step_map at current agent's position
        end
        
        %% Next, calculate attraction
        if ~isempty(free_list)
        for i = free_list % Iterate through free agents
            ind = indlist(:,i); % Indexed position of current agent
            list = l1(l2==i); % List of neighbors to current agent
            v_attraction=zeros(2,1); % Initialize attraction vector
            
            if ~isempty(list)
            for neigh_i = list % Iterate through neighbors
                dx = p(:,neigh_i)-p(:,i); %Vector from one agent to the next

                % Correction of distances for periodic boundary
                if PERIODIC==1
                     if dx(1) > L/2
                         dx(1) = dx(1) - L;
                     elseif dx(1) < -L/2
                        dx(1) = dx(1) + L;
                     end
                     if dx(2) > L/2
                        dx(2) = dx(2) - L;
                     elseif dx(2) < -L/2
                        dx(2) = dx(2) + L;
                    end
                end
                    
                v_attraction = v_attraction + dx; %Add distance vector dx to attraction vector
            end
            v_attraction = v_attraction/length(list); % Divide by number of neighbors
            end
           
            % Generate random vector
            randang = 2*pi*rand(); % First, generate random angle
            v_rand = [cos(randang); sin(randang)]; % Then, take sin and cos to generate XY vector
            
            % Generate velocity vector to neighboring square
            tempv = c*v_attraction+v_rand; % Sum attraction vector (scaled by c) and random vector
            ang = mod(atan2(tempv(2),tempv(1)),2*pi); % Take angle of tempv as between 0 and 2*pi
            DiscAng = floor((ang+pi/8)/(pi/4))*pi/4; % Discretize that angle to the nearest of 8 direction vectors
            v(:,i) = round([cos(DiscAng), sin(DiscAng)]);  %Update velocity of current agent with each component magnitude 0 or 1
        end
        end

        %% Finally, check that no two free agents occupy the same square
        for fi=randperm(length(free_list)) % Iterate through free_list in random order
            i = free_list(fi); % Index of current agent, i
            ind = indlist(:,i); % Indexed position of current agent

            p_temp=p(:, i)+v(:, i); % Intended new position
            % Adjust new position for periodic boundary
            if PERIODIC==1
                p_temp(p_temp>L/2) = p_temp(p_temp>L/2) - L;
                p_temp(p_temp<-L/2)  = p_temp(p_temp<-L/2) + L;
            end
            
            newind = floor((p_temp+L/2)+1); % New indexed position of current agent
            step_map(ind(1),ind(2)) = step_map(ind(1), ind(2))-1; % Leave the old position, set step_map to 0
            success = 0; % Marker for successful movement
            
            if step_map(newind(1),newind(2)) == 0 % Check if new position is unoccupied
                if h_map(newind(1),newind(2)) - h_level(i) <= max_climb % Check if the height difference is too tall
                    p(:, i) = p_temp; % Successfully moved to new position
                    step_map(newind(1),newind(2)) = step_map(newind(1),newind(2))+1; % Set step_map to 1 at new position
                else % If climb to new position fails
                    step_map(ind(1),ind(2)) = step_map(ind(1),ind(2))+1; % Current agent stays at old position, set step_map to 1
                end
                success = 1; % Movement successfully determined
            else % If new position is occupied
                for r=randperm(8) % Check each of 8 neighbor positions, in random order
                    % Generate random neighbor position index, from index r
                    temp_ind = newind+spacelist(r,:)';
                    % Ensure neighbor position obeys periodic boundary condition
                    if PERIODIC == 1
                        temp_ind = mod(temp_ind, L); 
                        temp_ind(temp_ind==0)=L; 
                    end

                    if step_map(temp_ind(1), temp_ind(2)) == 0 % Check if new position is unoccupied
                        if h_map(temp_ind(1), temp_ind(2)) - h_level(i) <= max_climb % Check if the height difference is too tall
                            p(:,i) = [temp_ind(1); temp_ind(2)]-L/2-0.5; % Successfully moved to new position
                            step_map(temp_ind(1), temp_ind(2)) = step_map(temp_ind(1), temp_ind(2))+1; % New position is occupied
                        else % If climb to new position fails
                            step_map(ind(1),ind(2)) = step_map(ind(1),ind(2))+1; % Current agent stays at old position, set step_map to 1
                        end
                        success = 1; % Movement successfully determined
                        break
                    end
                end
            end
            if success == 0 % Check if all attempted positions are occupied
                step_map(ind(1),ind(2)) = step_map(ind(1),ind(2))+1; % Agent remains in place
            end
        end
    end
    
    % Check periodic boundary condition
    tmp_p = p;
    if PERIODIC==1
        tmp_p(1,p(1,:)>L/2)  = tmp_p(1,p(1,:)>L/2) - L;
        tmp_p(2,p(2,:)>L/2)  = tmp_p(2,p(2,:)>L/2) - L;
        tmp_p(1,p(1,:)<-L/2) = tmp_p(1,p(1,:)<-L/2) + L;
        tmp_p(2,p(2,:)<-L/2) = tmp_p(2,p(2,:)<-L/2) + L;
    end
    p = tmp_p;
    
    %% Recalculate interagent distances 
    d2 = pdist([p(1, :)' p(2, :)' h_level'], @distfun3d);
    D2 = squareform(d2);
    [l1,l2]=find(0<D2 & D2<1.5);
    
    %%% 4. Determine status changes (free, locked, covered) of each agent
    all_map = zeros(L, L); % Initialize all_map, a height map of all agents
    indlist = floor(p+L/2+1); % Determine indexed locations of all agents

    % Build all_map by iterating through each agent and adding its height to its position
    for i=1:N 
        ind = indlist(:, i);
        all_map(ind(1), ind(2)) = all_map(ind(1), ind(2))+1;

        % While iterating, also count the number of locked neighbors for each agent
        list2 = l1(l2==i); 
        locked_neigh(i) = sum(locked(list2));
    end
    for i=randperm(N) % For each agent in random order
        ind = indlist(:, i);
        % Determine z position for free agents
        if locked(i) == 0
            h_level(i) = h_map(ind(1), ind(2)); % z position is equal to number of locked agents in its square
        end
        % Check if each locked agent is covered, based on comparison of its z position and all_map
        if locked(i) == 1 && all_map(ind(1), ind(2)) > h_level(i)+1
            covered(i) = 1;
        else
            covered(i) = 0;

            % If uncovered and unlocked
            if locked(i)==0
                % Randomly lock position, based on spontaneous lock probability and neighbor lock probability
                if rand()< p_s_l || rand() < locked_neigh(i)*k_nl
                    locked(i) = 1;
                    h_map(ind(1), ind(2)) = h_map(ind(1), ind(2))+1; % Increase h_map at locked agent location
                end
            elseif locked(i)==1 % If uncovered and locked
                % Spontaneously unlock position:
                if rand()<p_u
                    locked(i) = 0;
                    h_map(ind(1), ind(2)) = h_map(ind(1), ind(2))-1; % Decrease h_map at unlocked agent location
                end
            end
        end
    end

    % Save positions of all agents once every 250 frames
    p_frame = [p;h_level;locked;covered];
    if mod(k, save_freq) == 0 
    dlmwrite(strcat(folder,'/',savename,'_output.txt'),p_frame,'-append')
    end
end
end



% Function for calculation of distance
function D3 = distfun3d(XI,XJ)
global PERIODIC L
dx = (XI-XJ); % Calculate interparticle distance vectors

% Adjust distances for periodic boundary conditions
tmp_p = dx;
if PERIODIC == 1
    tmp_p(find(dx(:,1)>L/2),1) = tmp_p(find(dx(:,1)>L/2),1) - L;
    tmp_p(find(dx(:,2)>L/2),2) = tmp_p(find(dx(:,2)>L/2),2) - L;
    tmp_p(find(dx(:,1)<-L/2),1)  = tmp_p(find(dx(:,1)<-L/2),1) + L;
    tmp_p(find(dx(:,2)<-L/2),2)  = tmp_p(find(dx(:,2)<-L/2),2) + L;
end

% Take the max of the absolute value of each distance
% Two particles with a distance vector of [1,1,1] are considered adjacent (distance 1)
D3 = max([abs(tmp_p)]');
end
