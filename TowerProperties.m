function out = TowerProperties(c, p_u, k_nl, v, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB code for analyzing properties of towers %
% based on the results of TowerSimulation.m      %
%                                                %
% Code by Gary K. Nave, Jr.                      %
% University of Colorado Boulder                 %
% gary.k.nave@gmail.com                          %
% Last updated: January, 2020                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Input parameters
    % c, p_u, and k_nl are simulation parameters
    % that are used in the filename
    % For details, see TowerSimulation.m
    % v is the version (or trial) number, which
    % inreases automatically with successive trials

    % All of these define the filename:
    filename = sprintf('c%4.2f.pu%5.3f.knl%5.3f.v%d', c, p_u, k_nl, v);

%%% Optional parameters, parsed using inputParser
    % InputFolder    % Folder to read data from
                     % default: 'Data/Output'
    % InputFolder    % Folder to save data to
                     % default: 'Data/Output'
    % FrameRate      % Frame rate of output video
                     % default: 30
    % Quality        % Quality of output video, 1-100
                     % default: 75
    % ImageFrequency % Frequency of shown images
                     % In terms of *saved* frames
                     % default: 1

%%% Parse input parameters %%%%%
    parser = inputParser;

    addParameter(parser, 'InputFolder', 'Data/Output');  
    addParameter(parser, 'OutputFolder', 'Data/Output'); 

    parse(parser, varargin{:});

    folder = parser.Results.InputFolder;
    outfolder = parser.Results.OutputFolder;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Load necessary parameters from parameters.mat file
    load(sprintf('%s/%s_parameters.mat', folder, filename), 'L', 'N', 'frames', 'save_freq')
    % L: arena side length
    % N: number of agents
    % frames: number of time steps
    % save_freq: frequency of saving

    % Load data from output.txt file
    p_history = dlmread(sprintf('%s/%s_output.txt', folder, filename));

    % Number of recorded frames (including initial points) is:
    num_frames = frames/save_freq + 1;
    
    % out is the structure containing all tower property data
    out = struct;

    %%% Initialize all tracked tower properties
    out.NumTowers = zeros(1, num_frames);       % Number of towers in the simulation
    
    % Average tower properties show the mean of all towers in the simulation 
    out.AverageArea = zeros(1, num_frames);     % Base area
    out.AverageDiameter = zeros(1, num_frames); % Base effective diameter
    out.AverageHeight = zeros(1, num_frames);   % Tower height
    out.AverageRatio = zeros(1, num_frames);    % Tower height-effective diameter ratio
    out.AverageNumAnts = zeros(1, num_frames);  % Number of agents per tower
    
    % Max tower properties show the properties of the largest tower in the simulation
    % based on number of agents
    out.MaxArea = zeros(1, num_frames);         % Base area
    out.MaxHeight = zeros(1, num_frames);       % Tower height
    out.MaxRatio = zeros(1, num_frames);        % Tower height-effective diameter ratio
    out.MaxNumAnts = zeros(1, num_frames);      % Number of agents per tower
    
    % Iterate through all simulation frames
    for f = 0:num_frames-1
        % Load x, y positions from p_history
        p = p_history((5*f+1):(5*f+2), :);
        
        % Find index-based positions
        indlist = floor((p+L/2)+1);

        % Initialize height map
        h_map = zeros(L, L);
        % Add one to height map for every agent
        for i=1:N
            ind = indlist(:, i);
            h_map(ind(1), ind(2)) = h_map(ind(1), ind(2))+1; 
        end

        % Binarize height map to occupied/unoccupied pixels
        binary = imbinarize(h_map, 0.5);
        % Label connected components (which correspond to tower bases)
        labeled = bwlabel(binary, 8);
        % Identify equivalent diameter and area of each connected component/tower base
        stats = regionprops(labeled, 'EquivDiameter', 'Area');
        
        % Identify max height in pixel and tower number for each agent
        towers = zeros(1, N);
        h_list = zeros(1, N);
        for i=1:N
            ind = indlist(:, i);
            h_list(i) = h_map(ind(1), ind(2));
            towers(i) = labeled(ind(1), ind(2));
        end
        
        % Number of towers corresponds to the number of connected components
        out.NumTowers(f+1) = length(stats);
        
        area = zeros(1,out.NumTowers(f+1));
        diameter = zeros(1,out.NumTowers(f+1));
        height = zeros(1,out.NumTowers(f+1));
        numAnts = zeros(1,out.NumTowers(f+1));
        
        for k = 1:out.NumTowers(f+1)
            area(k) = stats(k).Area;
            diameter(k) = stats(k).EquivDiameter;
            height(k) = max(h_list(find(towers==k)));
            numAnts(k) =  length(h_list(find(towers==k)));
            
            if numAnts(k) > out.MaxNumAnts(f+1)
                out.MaxHeight(f+1) = height(k);
                out.MaxArea(f+1) = stats(k).Area;
                out.MaxNumAnts(f+1) = numAnts(k);
                out.MaxRatio(f+1) = height(k)/diameter(k);
            end
        end

        ratio = height/diameter;
        out.AverageArea(f+1) = mean(area);
        out.AverageDiameter(f+1) = mean(diameter);
        out.AverageHeight(f+1) = mean(height);
        out.AverageRatio(f+1) = mean(ratio);
        out.AverageNumAnts(f+1) = mean(numAnts);
    end
    save(strcat(outfolder,'/',filename,'_out.mat'),'out')
end
