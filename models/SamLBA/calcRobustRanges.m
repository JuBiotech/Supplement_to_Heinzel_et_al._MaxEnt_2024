function [ranges, table] = calcRobustRanges( model, groups, normPC, singlePC, maxIter )
%CALCROBUSTRANGES Calculates robust ranges for each coefficient group.
% A robust range is the range in which a coefficent group can vary without
% the solution of the FBA with the uncertain data violating the conditions
%	norm(v_nom - v_sol) <= norm(v_nom) * normPC / 100
%	abs(v_nom - v_sol) <= v_nom * singlePC / 100
% So there are boundaries on each flux vector component and on the flux vector
% as a whole.
%
% Algorithm: This method uses a simple bisection search for the min / max range.
%
% Parameters:
%	- model: Model
%	- groups: Cell array with group structures. See MCGroups.m for hints.
%	- normPC: Allowed maximum relative deviation of the flux vector norm
%	- singlePC: Allowd maximum relative deviation of each flux vector component.
%	- maxIter: Maximum iteration count used for bisection search.
%
% Returns:
%	- ranges: A matrix with the ranges for each group. Each line corresponds to
%		a metabolite of a group (in order of groups cell array and group.row vector).
%		The columns give min, max, interval width and nominal coefficient value (in this order)
%	- table: A cell array holding the ranges matrix and additional captions.

    % Check constraints
    if ~isfield(model, 'cM')
        model.cM = [];
    end
    if ~isfield(model, 'cB')
        model.cB = [];
    end

    if nargin <= 4
        maxIter = 100;
    end
    
	[origSol] = solveLPProblem(-1, model.c, model.cM, model.cB, model.S, model.b, model.lb, model.ub);
    [groups, sampleSize, origData] = checkGroupsAndSaveOriginalData(model, groups, 'none');    
    
    diffSol = abs(origSol);
    diffSol(diffSol == 0) = 1;
    
    ranges = zeros(sampleSize, 4);
    table = cell(sampleSize, 1);
    
    ranges(:, 4) = origData;
    
    tHandle = tic;
    index = 1;
    for i = 1:length(groups)
		if lower(groups{i}.type) == 's'
        	mets = model.mets(groups{i}.row);
		else
			mets = cell(length(groups{i}.row), 1);
			for j = 1:length(groups{i}.row)
    			mets{j} = ['Group ' groups{i}.type ' ' num2str(groups{i}.row(j))];
            end
		end
		if size(mets,2) > 1
			mets = mets';
		end
        table(index:index+length(mets)-1) = mets;

        [ranges(index:index+length(mets)-1, 1), ranges(index:index+length(mets)-1, 2)] = ...
            robustRangeSingle(model, groups(i), norm(origSol) * normPC / 100, diffSol * singlePC / 100, maxIter);
        
        ranges(index:index+length(mets)-1, 3) = abs(ranges(index:index+length(mets)-1,2) - ranges(index:index+length(mets)-1, 1));
        index = index + length(mets);
    end
    table = [[{'Metabolite', 'Min', 'Max', 'Interval', 'Nominal'}]; [table num2cell(ranges)]];

    tElapsed = toc(tHandle);
    fprintf('Elapsed time: %g sec = %g h %g min\n', tElapsed, floor(tElapsed / 3600), floor((tElapsed - floor(tElapsed / 3600)*3600) / 60));
end
