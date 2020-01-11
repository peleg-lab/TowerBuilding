function TowerVideo(folder, c, p_u, k_nl, v)
    filename = sprintf('c%4.2f.pu%5.3f.knl%5.3f.v%d', c, p_u, k_nl, v);
    load(sprintf('%s/%s_parameters.mat', folder, filename), 'L', 'N', 'frames', 'save_freq')
    p_history = dlmread(sprintf('%s/%s_output.txt', folder, filename));
    
    myVideo = VideoWriter(sprintf('%s/%s.avi', folder, filename));
    myVideo.FrameRate = 10;
    myVideo.Quality = 75;
    M=moviein(length(10:10:frames));
    
    figure('Units', 'normalized', 'Position', [0.1 0.1 0.5 0.5])
    for k=0:1:frames
    
        p = p_history((5*k+1):(5*k+2), :);
        z = p_history(5*k+3, :);
        
        clf
        matlab_viz([p(1,:)',p(2,:)',z'], 10);
        title(sprintf('Time Step: %1d',k*save_freq))
        xlim([-50,50]); ylim([-50,50])
        M(:,(k-10)/10+1) = getframe(gcf);
    end
    open(myVideo);
    writeVideo(myVideo, M);
end
