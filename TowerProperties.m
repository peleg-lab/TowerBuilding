function out = TowerProperties(folder, c, p_u, k_nl, v)
    filename = sprintf('c%4.2f.pu%5.3f.knl%5.3f.v%d', c, p_u, k_nl, v);
    
    load(sprintf('%s/%s_parameters.mat', folder, filename), 'L', 'frames', 'N')
    p_history = dlmread(sprintf('%s/%s_output.txt', folder, filename));
    
    out = struct;

    out.NumTowers = zeros(1, frames);
    out.AverageArea = zeros(1, frames);
    out.AverageDiameter = zeros(1, frames);
    out.AverageHeight = zeros(1, frames);
    out.AverageRatio = zeros(1, frames);
    out.AverageNumAnts = zeros(1, frames);
    
    out.MaxArea = zeros(1, frames);
    out.MaxHeight = zeros(1, frames);
    out.MaxRatio = zeros(1, frames);
    out.MaxNumAnts = zeros(1, frames);
    
    for f = 1:frames
        p = p_history((5*f+1):(5*f+2), :);
        
        indlist = floor((p+L/2)+1);
        h_map = zeros(L, L);
        for i=1:N
            ind = indlist(:, i);
            h_map(ind(1), ind(2)) = h_map(ind(1), ind(2))+1;
        end
        binary = imbinarize(h_map, 0.5);
        labeled = bwlabel(binary, 8);
        stats = regionprops(labeled, 'EquivDiameter', 'Area');
        
        towers = zeros(1, N);
        h_list = zeros(1, N);
        for i=1:N
            ind = indlist(:, i);
            h_list(i) = h_map(ind(1), ind(2));
            towers(i) = labeled(ind(1), ind(2));
        end
        
        out.NumTowers(f) = length(stats);
        
        area = zeros(1,out.NumTowers(f));
        diameter = zeros(1,out.NumTowers(f));
        height = zeros(1,out.NumTowers(f));
        numAnts = zeros(1,out.NumTowers(f));
        
        for k = 1:out.NumTowers(f)
            area(k) = stats(k).Area;
            diameter(k) = stats(k).EquivDiameter;
            height(k) = max(h_list(find(towers==k)));
            numAnts(k) =  length(h_list(find(towers==k)));
            
            if numAnts(k) > out.MaxNumAnts(f)
                out.MaxHeight(f) = height(k);
                out.MaxArea(f) = stats(k).Area;
                out.MaxNumAnts(f) = numAnts(k);
                out.MaxRatio(f) = height(k)/diameter(k);
            end
        end

        ratio = height/diameter;
        out.AverageArea(f) = mean(area);
        out.AverageDiameter(f) = mean(diameter);
        out.AverageHeight(f) = mean(height);
        out.AverageRatio(f) = mean(ratio);
        out.AverageNumAnts(f) = mean(numAnts);
    end
    save(strcat(folder,'/',filename,'_out.mat'),'out')
end