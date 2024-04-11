function compArea = calc2DPolyCoverage( run )
%CALC2DPOLYCOVERAGE Summary of this function goes here
%   Detailed explanation goes here

    compArea = cell(length(run.polys2D)+1, 6);
    compArea(1, :) = [{'Group 1'}, {'Group 2'}, {'Area'}, {'Clamped Area'}, {'Ratio'}, {'Ratio 100%'}];
    for i = 1:length(run.polys2D)
    
        pos = cell2mat(run.polys2D{i}(2:end, 1:3));
        metaNames = run.polys2D{i}(1, 2:3);

        bounds = [pos(pos(1:end-1,1) == 0, 2), pos(pos(:,1) == 180, 2); ...
            pos(pos(:,1) == 90, 3), pos(pos(:,1) == 270, 3)];
        
        maxs = [max(bounds(1, :)); max(bounds(2, :))];
        mins = [min(bounds(1, :)); min(bounds(2, :))];
        
        % Clamp the polygon
        pos = pos(:, 2:3);
        pos(:, 1) = min(max(pos(:, 1), mins(1)), maxs(1));
        pos(:, 2) = min(max(pos(:, 2), mins(2)), maxs(2));
        
        clampedArea = polyarea(pos(:,1), pos(:,2));
        
        % Full polygon size
        fullArea = run.areaTable2D{i, 3};
        
        compArea(i+1, 1:2) = metaNames;
        compArea(i+1, 3:6) = num2cell([fullArea, clampedArea, clampedArea / fullArea, clampedArea / fullArea * 100]);
    end

end

