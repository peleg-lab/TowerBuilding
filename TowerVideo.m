function TowerVideo(c, p_u, k_nl, v, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB code for making movies based %
% on the results of TowerSimulation.m %
%                                     %
% Code by Gary K. Nave, Jr.           %
% University of Colorado Boulder      %
% gary.k.nave@gmail.com               %
% Last updated: January, 2020         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    % InputFolder    % Folder to save video to
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
    addParameter(parser, 'FrameRate', 30); 
    addParameter(parser, 'Quality', 75); 
    addParameter(parser, 'ImageFrequency', 1); 

    parse(parser, varargin{:});

    folder = parser.Results.InputFolder;
    outfolder = parser.Results.OutputFolder;
    framerate = parser.Results.FrameRate;
    quality = parser.Results.Quality;
    framefreq = parser.Results.ImageFrequency;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Load necessary parameters from parameters.mat file
    load(sprintf('%s/%s_parameters.mat', folder, filename), 'L', 'N', 'frames', 'save_freq')
    % L: arena side length
    % N: number of agents
    % frames: number of time steps
    % save_freq: frequency of saving

    % Load data from output.txt file
    p_history = dlmread(sprintf('%s/%s_output.txt', folder, filename));
    
    % Make sure that the output folder exists 
    if ~exist(outfolder)
        mkdir(outfolder)
    end

    % Video savename
    myVideo = VideoWriter(sprintf('%s/%s.avi', outfolder, filename));
    % Video Frame Rate
    myVideo.FrameRate = framerate;
    % Video Quality, scale 1-100
    myVideo.Quality = quality;

    % Frame list, in terms of number of saved frames
    framelist = 0:framefreq:frames/save_freq;

    % Initialize number of frames
    M=moviein(length(framelist));
    
    % Initialize figure
    figure('Units', 'normalized', 'Position', [0.1 0.1 0.5 0.5])
    for k=framelist

        % Load x, y, z position from frame of interest
        p = p_history((5*k+1):(5*k+2), :);
        z = p_history(5*k+3, :);
        
        % Clear figure
        clf
        % Run the matlab_viz function (see matlab_viz.m)
        matlab_viz([p(1,:)',p(2,:)',z'], 10);
        % Generate title, including frame number
        title(sprintf('Time Step: %1d',k*save_freq))
        % Set axis limits
        xlim([-L/2,L/2]); ylim([-L/2,L/2])
        % Save frame
        M(:,(k-1)/1+2) = getframe(gcf);
    end
    % Write data to video
    open(myVideo);
    writeVideo(myVideo, M);
end
