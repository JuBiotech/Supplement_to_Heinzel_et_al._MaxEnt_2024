function compArea = calc2DPolyStats( run, origObjVal )
%CALC2DPOLYSTATS Calculates some key values for each polygon in the given
%run. Calculates the 1D bounding box area, the polygon area, the clamped
%area, the coverage, objective value (min, max, range, relative range to
%original objective value) and the length of the principal axis.
%
%The clamped polygon is created by taking the intersection of the polygon
%and the 1D bounding box.
%
%The coverage is the percentage of the area the clamped polygon covers with
%respect to the original polygon.
%
%The principal axis is the axis of the polygon with the longest length.
%
% Parameters:
%   - run: The run field of the generateRobustRanges() output.
%   - origObjVal: Original objective value. Optional, is determined by
%       run.resSphere, if possible. This field only affects the relative
%       objective range.
%
% Returns a cell table with headings and the values.

    if (nargin < 2) || isempty(origObjVal)
        if isfield(run, 'resSphere')
            origObjVal = run.resSphere.origVal;
        end
    end

    compArea = cell(length(run.polys2D)+1, 13);
    compArea(1, :) = [{'Group 1'}, {'Group 2'}, {'1D BB'}, {'Area'}, {'Clamped Area'}, {'Coverage'}, {'ObjMin'}, {'ObjMax'}, {'ObjRange'}, {'ObjRelRange'}, {'MaxAxis'}, {'MinAxis'}, {'AxisRatio'}];
    for i = 1:length(run.polys2D)
    
        pos = cell2mat(run.polys2D{i}(2:end, 2:3));
        angles = cell2mat(run.polys2D{i}(2:end, 1));
        objVals = cell2mat(run.polys2D{i}(2:end, 4));
        metaNames = run.polys2D{i}(1, 2:3);

        objVals = [min(objVals), max(objVals)];
        
        bounds = [pos(angles(1:end-1) == 0, 1), pos(angles == 180, 1); ...
            pos(angles == 90, 2), pos(angles == 270, 2)];
        
        % Full polygon size
        fullArea = run.areaTable2D{i, 3};
        
        % Principal axis
        axisLen = zeros(ceil(size(pos, 1) / 2), 1);
        for j = 1:size(pos, 1)
            if (angles(j) >= angles(1) + 180)
                break;
            end
            
            if angles(j)+180 > 360
                secInd = abs(angles - angles(j)+180) < 1e-6;
            else
                secInd = abs(angles - angles(j)-180) < 1e-6;
            end
            
            bla = sum(secInd);
            if bla == 0
                continue;
            end
            
            axisLen(j) = norm(pos(j, :) - pos(secInd, :));
        end
        axisLen = axisLen(~isnan(axisLen));

        % 1D bounding box
        maxs = [max(bounds(1, :)); max(bounds(2, :))];
        mins = [min(bounds(1, :)); min(bounds(2, :))];
        
        % Clamp the polygon
        pos(:, 1) = min(max(pos(:, 1), mins(1)), maxs(1));
        pos(:, 2) = min(max(pos(:, 2), mins(2)), maxs(2));
        
        clampedArea = polyarea(pos(:,1), pos(:,2));
        
        compArea(i+1, 1:2) = metaNames;
        compArea(i+1, 3:13) = num2cell([prod(maxs - mins), fullArea, clampedArea, clampedArea / fullArea, ...
            objVals(1), objVals(2), objVals(2)-objVals(1), (objVals(2)-objVals(1)) / abs(origObjVal), max(axisLen), ...
            min(axisLen), max(axisLen) / min(axisLen)]);
    end

end

