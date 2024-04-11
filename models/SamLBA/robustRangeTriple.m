function [vals, table, relVals] = robustRangeTriple( model, groups, angPerHalfCircle, normBoundary, singleBoundary, maxIter, relative, nAnglesRun, nWorkers )
%ROBUSTRANGETRIPLE Calculates the robust coefficient ranges of the three given coefficient groups.
% Computes the robust polyhedron of the given coefficients. A robust extreme, which forms the
% vertices of a robust polyhedron, is the maximum / minimum coefficient value for which the optimal solution
% of the FBA problem subject to the modified coefficients suffices the following conditions:
%	norm(v_nom - v_sol) <= normBoundary  and  abs(v_nom - v_sol) <= singleBoundary
%
% See getRobustExtremeConstrainedFlux for a description of the method used.
%
% This method uses spherical coordinates two create a grid on the unit sphere. The robust extreme 
% is computed for each direction in the three groups specified by the grid points.
% The parameter angPerHalfCircle gives the amount of steps used to partition a 180 degree angle.
% So there are angPerHalfCircle * (angPerHalfCircle - 1) * 2 grid points (leaving one out to avoid doubling).
%
% Note: Using spherical coordinates to create the grid (that is partitioning the theta and phi angles)
%	does not yield a uniform grid. This method leads to a concentration of the grid points at the poles.
%
% Parameters:
%	- model: Model.
%	- groups: Cell array with two groups. See MCGroups for hints.
%	- angPerHalfCircle: Count of partitions of a 180 degree angle.
%	- normBoundary: Upper bound on the solution deviation in the norm.
%	- singleBoundary: Upper bound on the absolute flux vector component deviation.
%	- maxIter: Maximum number of iterations. See getRobustExtremeConstrainedFlux for hints.
%	- relative: Set to true to treat the parameters normBoundary and singleBoundary as percentage of
%		the nominal values. Defaults to false.
%
% Returns:
%	- vals: A matrix in which each line is a ray cast in the direction of the i-th angle in the
%		angles vector. The columns are (in this order): 
%		theta in rad, phi in rad, absolute coefficient values, objective function value
%	- table: A table (cell array) of the matrix vals with headings.
%	- relVals: The same as vals except that the middle columns contain the relative robust multiplicator.
%		So there are always three columns (one for each group) in the middle of the matrix.
%       The last column contains the relative multiplicator of the current ray.

    % Check constraints
    if ~isfield(model, 'cM')
        model.cM = [];
    end
    if ~isfield(model, 'cB')
        model.cB = [];
    end

	[origSol] = solveLPProblem(-1, model.c, model.cM, model.cB, model.S, model.b, model.lb, model.ub);

    if (nargin > 6) && relative
        normBoundary = norm(origSol) * normBoundary / 100;
        diffSol = abs(origSol);
        diffSol(diffSol == 0) = 1;
        singleBoundary = diffSol * singleBoundary / 100;
    end
    
    % Set single boundaries
    if length(singleBoundary) ~= length(origSol)
        singleBoundary = ones(length(origSol), 1) .* singleBoundary;
    end
   
    if any(model.ub < model.lb)
        error('ROBUSTRANGETRIPLE:BoundaryCheck', 'Problem is infeasible. Please check your boundary conditions.');
    end
    
	nomVal = getModelCoefficients(model, groups);
    vals = zeros(angPerHalfCircle * (angPerHalfCircle - 1) * 2, length(nomVal) + 3);
    relVals = zeros(angPerHalfCircle * (angPerHalfCircle - 1) * 2, 7);
    phi = linspace(-pi / 2, pi / 2, angPerHalfCircle);
    theta = linspace(0, 2 * pi, 2 * angPerHalfCircle);
    theta = theta(1:end-1);

    curInd = 1;
    nCurRun = nAnglesRun;
    
    for i = 1:length(theta)

        tempVals = zeros(length(phi), length(nomVal) + 3);
        tempRelVals = zeros(length(phi), 7);

        if (nCurRun+1 == i) || (i == 1)
            nCurRun = nCurRun + nAnglesRun;
            
            if (nargin <= 8) || isempty(nWorkers)
                matlabpool open;
            else
                matlabpool('open', nWorkers);
            end
        end
        
        parfor k = 1:length(phi)
            dir = [cos(phi(k)) * cos(theta(i)); cos(phi(k)) * sin(theta(i)); sin(phi(k))];

            [resLine, objVal, relLine, multi] = getRobustExtremeConstrainedFlux(dir, nomVal, origSol, model, groups, normBoundary, singleBoundary, maxIter);
            tempVals(k, :) = [theta(i), phi(k), resLine', objVal];
            tempRelVals(k, :) = [theta(i), phi(k), relLine', objVal, multi];
        end

        if (i == nCurRun) && (matlabpool('size') > 0)
            matlabpool close;
        end

        vals(curInd:curInd+length(phi)-1, :) = tempVals;
        relVals(curInd:curInd+length(phi)-1, :) = tempRelVals;
        curInd = curInd + length(phi);
    end
    
    if matlabpool('size') > 0
        matlabpool close;
    end
    
    if nargout > 1
        metNames = {};
        for i = 1:length(groups)
            metNames = [metNames model.mets(groups{i}.row)'];
        end
        table = [[{'theta'}, {'phi'}, metNames , {'objVal'}] ; num2cell(vals)];
    end
end
