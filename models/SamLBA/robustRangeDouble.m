function [vals, relVals] = robustRangeDouble( model, groups, angles, normBoundary, singleBoundary, maxIter )
%ROBUSTRANGEDOUBLE Calculates the robust coefficient ranges of the two given coefficient groups.
% Computes the robust polygon of the given coefficients. A robust extreme, which forms the
% vertices of a robust polygon, is the maximum / minimum coefficient value for which the optimal solution
% of the FBA problem subject to the modified coefficients suffices the following conditions:
%	norm(v_nom - v_sol) <= normBoundary  and  abs(v_nom - v_sol) <= singleBoundary
%
% See getRobustExtremeConstrainedFlux for a description of the method used.
%
% This method determines the robust extreme for each direction in the two groups specified
% by the angles (in degree) given in the angles vector.
%
% Note: By default this function runs in parallel using the Parallel Computing Toolbox.
%	However, a worker pool will not be started or closed by this function, so the caller has to
%	open and close the pool himself.
%
% Parameters:
%	- model: Model.
%	- groups: Cell array with two groups. See MCGroups for hints.
%	- angles: Vector with angles in degree.
%	- normBoundary: Upper bound on the solution deviation in the norm.
%	- singleBoundary: Upper bound on the absolute flux vector component deviation.
%	- maxIter: Maximum number of iterations. See getRobustExtremeConstrainedFlux for hints.
%
% Returns:
%	- vals: A matrix in which each line is a ray cast in the direction of the i-th angle in the
%		angles vector. The columns are (in this order): angle, absolute coefficient values, objective function value
%	- relVals: The same as vals except that the middle columns contain the relative robust multiplicator.
%		So there are always two columns (one for each group) in the middle of the matrix.
%       The last column contains the relative multiplicator of the current ray.

    % Check constraints
    if ~isfield(model, 'cM')
        model.cM = [];
    end
    if ~isfield(model, 'cB')
        model.cB = [];
    end

	[origSol] = solveLPProblem(-1, model.c, model.cM, model.cB, model.S, model.b, model.lb, model.ub);

    % Set single boundaries
    if length(singleBoundary) ~= length(origSol)
        singleBoundary = ones(length(origSol), 1) .* singleBoundary;
    end
    
    if any(model.ub < model.lb)
        error('ROBUSTRANGEDOUBLE:BoundaryCheck', 'Problem is infeasible. Please check your boundary conditions.');
    end
    
	nomVal = getModelCoefficients(model, groups);
    vals = zeros(length(angles), length(nomVal) + 2);
    relVals = zeros(length(angles), 5);

% Uncomment the following line (and comment the parfor) to run in serial mode
    for i = 1:length(angles)

% Uncomment the following line (and comment the for in the previous line) to run in parallel mode
%    parfor i = 1:length(angles)
		phi = (pi/180) * angles(i);
		phi = [cos(phi); sin(phi)];
		[resLine, objVal, relLine, multi] = getRobustExtremeConstrainedFlux(phi, nomVal, origSol, model, groups, normBoundary, singleBoundary, maxIter);
        vals(i, :) = [angles(i), resLine', objVal];
        relVals(i, :) = [angles(i), relLine', objVal, multi];
    end
end
