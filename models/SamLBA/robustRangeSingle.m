function [ minVal, maxVal, relVals ] = robustRangeSingle( model, group, normBoundary, singleBoundary, maxIter )
%ROBUSTRANGESINGLE Calculates the robust coefficient range of the given coefficient group.
% Computes the robust interval (min, max) of the given coefficient. A robust extreme, which forms the
% endpoint of a robust interval, is the maximum / minimum coefficient value for which the optimal solution
% of the FBA problem subject to the modified coefficient suffices the following conditions:
%	norm(v_nom - v_sol) <= normBoundary  and  abs(v_nom - v_sol) <= singleBoundary
%
% See getRobustExtremeConstrainedFlux for a description of the method used.
%
% Parameters:
%	- model: Model.
%	- group: Group structure. See MCGroups for hints.
%	- normBoundary: Upper bound on the solution deviation in the norm.
%	- singleBoundary: Upper bound on the absolute flux vector component deviation.
%	- maxIter: Maximum number of iterations. See getRobustExtremeConstrainedFlux for hints.
%
% Returns:
%	- minVal: Minimum robust flux coefficient values.
%	- maxVal: Maximum robust flux coefficient values.
%	- relVals: Vector with min and max relative robust multiplicator. See getRobustExtremeConstrainedFlux. 

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
        error('ROBUSTRANGESINGLE:BoundaryCheck', 'Problem is infeasible. Please check your boundary conditions.');
    end
    
    relVals = zeros(2,1);
	nomVal = getModelCoefficients(model, group);
	[minVal, temp, relVals(1)] = getRobustExtremeConstrainedFlux(-1, nomVal, origSol, model, group, normBoundary, singleBoundary, maxIter);
	[maxVal, temp, relVals(2)] = getRobustExtremeConstrainedFlux(1, nomVal, origSol, model, group, normBoundary, singleBoundary, maxIter);
end
