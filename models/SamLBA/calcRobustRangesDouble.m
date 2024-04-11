function [areaTable, polys,absData] = calcRobustRangesDouble( model, groups, angles, normPC, singlePC, maxIter, nWorkers )
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
%		Angle, first metabolite of group 1, first metabolite of group 2,
%		objective function value, multiplicator

    % Check constraints
    if ~isfield(model, 'cM')
        model.cM = [];
    end
    if ~isfield(model, 'cB')
        model.cB = [];
    end

    if nargin <= 5
        maxIter = 100;
    end
    
    if mod(angles(1) - angles(end), 360) == 0
        % Remove last angle, as it is the same as the first (mod 360)
        angles = angles(1:end-1);
    end
    
	[origSol] = solveLPProblem(-1, model.c, model.cM, model.cB, model.S, model.b, model.lb, model.ub);
    
    diffSol = abs(origSol);
    diffSol(diffSol == 0) = 1;
    diffSol = diffSol * singlePC / 100;    
    
    normPC = norm(origSol) * normPC / 100;
    
    groupPairs = nchoosek(1:length(groups), 2);
    
    areaTable = cell(size(groupPairs, 1), 2);
    areas = zeros(size(groupPairs, 1), 1);
%    areasWeighted = zeros(size(groupPairs, 1), 1);
    polys = cell(size(groupPairs, 1), 1);
    
    nPairs = size(groupPairs, 1);
    nParRun = 4;

    tHandle = tic;
    
    for j = 0:floor(nPairs / nParRun)
%        parfor i = (j*nParRun+1):min(nPairs, ((j+1)*nParRun))
        for i = (j*nParRun+1):min(nPairs, ((j+1)*nParRun))
    %    parfor i = 1:size(groupPairs, 1)
    %    for i = 1:size(groupPairs, 1)
            gp = groups(groupPairs(i, :));

            areaTable(i,:) = getMetaboliteHeading(model, gp);

            if (nargin <= 6) || isempty(nWorkers)
%                matlabpool open;
            else
%                matlabpool('open', nWorkers);
            end
            
            [absData, polyCoords] = robustRangeDouble(model, gp, angles, ...
                normPC, diffSol, maxIter);

%             matlabpool close;
            
            polyCoords = [polyCoords; polyCoords(1, :)];
            polys{i} = [{'angle'}, getMetaboliteHeading(model, gp), {'ObjVal'}, {'Multiplicator'}; num2cell(polyCoords)];

            % Calculate area
            areas(i) = polyarea(polyCoords(:, 2), polyCoords(:, 3));
%            areasWeighted(i) = areas(i) / abs(model.S(gp{1}.row, gp{1}.column) * model.S(gp{2}.row, gp{2}.column));
        end
        fprintf('Status: %g \n', min(nPairs, ((j+1)*nParRun)) * 100 / nPairs);
    end
    areaTable = [areaTable num2cell(areas)]; % num2cell(areasWeighted)];

    tElapsed = toc(tHandle);
    fprintf('Elapsed time: %g sec = %g h %g min\n', tElapsed, floor(tElapsed / 3600), floor((tElapsed - floor(tElapsed / 3600)*3600) / 60));
end

function mets = getMetaboliteHeading(model, gp)
    mets = [];
	for i = 1:length(gp)
		if lower(gp{i}.type) == 's'
            if length(gp{i}.row) > 1
                temp = '';
    			for j = 1:length(gp{i}.row)
        			temp = [temp, model.mets{gp{i}.row(j)} ' & '];
                end
                temp = temp(1:end-3);
                mets = [mets, {temp}];
            else
    		   	mets = [mets, model.mets(gp{i}.row)];
            end
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
