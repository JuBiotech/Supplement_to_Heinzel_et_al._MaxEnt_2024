function plotSingleMetaData( meta, xInd, dataInd )
%PLOTSINGLEMETADATA Summary of this function goes here
%   Detailed explanation goes here

    n = length(meta);
    colors = distinguishable_colors(n, 'w');
    lStr = cell(n, 1);
    
    h = plot(meta{1}(:, xInd), meta{1}(:, dataInd));
    set(h, 'Color', colors(1, :));
    lStr{1} = 'Gruppe 1';
    
    hold on;
    for i = 2:n
        h = plot(meta{i}(:, xInd), meta{i}(:, dataInd));
        set(h, 'Color', colors(i, :));
        lStr{i} = ['Gruppe ' num2str(i)];
    end
    hold off;
    
    legend(lStr, 'Location', 'NorthWest');
end

