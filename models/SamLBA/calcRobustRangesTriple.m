function calcRobustRangesTriple( model, groups,  angPerHalfCircle, normPC, singlePC, maxIter, baseName, nAnglesRun, nWorkers )
%CALCROBUSTRANGESDOUBLE Calculates robust ranges or polygon for each pair of coefficient groups.
% A robust range is the range in which a coefficent group can vary without
% the solution of the FBA with the uncertain data violating the conditions
%	norm(v_nom - v_sol) <= norm(v_nom) * normPC / 100
%	abs(v_nom - v_sol) <= v_nom * singlePC / 100
% So there are boundaries on each flux vector component and on the flux vector
% as a whole.
%
% In 2D a pair of groups is taken and the search direction is given by an angle (polar coordinates).
%
% Algorithm: This method uses a simple bisection search for the min / max in the direction of
%	[cos(angle); sin(angle)] thus weighting the groups. Instead of looking for absolute values
%	for the coefficients of the groups it uses a multiplicator. Thus maximizing / minimizing the factor
%	r in coefficients(group) * (1+r).
%
% Parameters:
%	- model: Model
%	- groups: Cell array with group structures. See MCGroups.m for hints.
%	- angles: A vector with angles (in degree) which are used for search directions.
%	- normPC: Allowed maximum relative deviation of the flux vector norm
%	- singlePC: Allowd maximum relative deviation of each flux vector component.
%	- maxIter: Maximum iteration count used for bisection search.
%	- nWorkers: Optional. Count of workers the Parallel Computing Toolbox uses.
%
% Returns:
%	- areaTable: A matrix with the areas for each pair of groups.
%		Note that the area here is computed by taking the relative group range. The
%		areas do not have to be weighted for comparison.
%	- polys: A cell array holding all computed relative robust range polygons. Each cell
%		is a table with headings. The columns are (in this order): 
%		Angle, first metabolite of group 1, first metabolite of group 2, objective function value

    % Check constraints
    if ~isfield(model, 'cM')
        model.cM = [];
    end
    if ~isfield(model, 'cB')
        model.cB = [];
    end

    if nargin <= 5
        maxIter = 5000;
    end
    
    groupPairs = nchoosek(1:length(groups), 3);

    tHandle = tic;
    
    for i = 1:size(groupPairs, 1)
        gp = groups(groupPairs(i, :));

        mets = getMetaboliteHeading(model, gp);
        
        tSub = tic;
        [vals, table, relVals] = robustRangeTriple( model, gp, angPerHalfCircle, normPC, singlePC, maxIter, true, nAnglesRun, nWorkers );
        tElapsed = toc(tSub);
        
        save([baseName, '_', mets{1}, '_', mets{2}, '_', mets{3}, '.mat'], 'mets', 'vals', 'table', 'relVals', 'tElapsed', 'nWorkers', 'angPerHalfCircle');
        fprintf('Status: %g \n', floor(i / size(groupPairs, 1) * 100));
    end

    tElapsed = toc(tHandle);
    fprintf('Elapsed time: %g sec = %g h %g min\n', tElapsed, floor(tElapsed / 3600), floor((tElapsed - floor(tElapsed / 3600)*3600) / 60));
end

function mets = getMetaboliteHeading(model, gp)
    mets = [];
	for i = 1:length(gp)
		if lower(gp{i}.type) == 's'
		   	mets = [mets, model.mets(gp{i}.row)];
		else
			for j = 1:length(gp{i}.row)
    			mets = [mets, {['Group ' gp{i}.type ' ' num2str(gp{i}.row(j))]}];
            end
		end
	end
	
	if size(mets, 1) > 1
		mets = mets';
	end
end
