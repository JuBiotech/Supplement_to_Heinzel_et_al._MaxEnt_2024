function [minFlux, maxFlux] = fluxVariabilityAnalysis(model, minNorm, dpcDevFromObjVal)
%FLUXVARIABILITYANALYSIS Performs a flux variability analysis yielding an axis aligned bounding box of the optimal flux configurations.
% FVA constrains the model to the optimal objective function value (or a certain amount of it) and tries
% to maximize / minimize the other flux vector components (one at a time).
% This yields an axis aligned bounding box of the hyperplane with the optimal objective function value.
%
% Note: If the optimum is unique the FVA coefficient intervals are zero.
%
% Parameters:
%	- model: Model.
%   - minNorm: Optional. Set to false to disable 1-Norm minimization (cycle
%       avoiding). Set to true to enable norm minimization. Default is
%       false.
%	- dpcDevFromObjVal: Optional. If > 0 the original objective function value has to be attained up to
%		a deviation given by dpcDevFromObjVal in percent. For example, if dpcDevFromObjVal = 5 the
%		objective function value may be in the range of 95 - 105% of the nominal objective function value.
%
% Returns:
%	- minFlux: Vector with minimum flux components.
%	- maxFlux: Vector with maximum flux components.

    % Check constraints
    if ~isfield(model, 'cM')
        model.cM = [];
    end
    if ~isfield(model, 'cB')
        model.cB = [];
    end
    
    if size(model.c, 1) == 1
        model.c = model.c';
    end
    
    if (nargin <= 1) || isempty(minNorm)
        minNorm = false;
    end
    params.minNorm = minNorm;
    
    [origSol, origObjVal] = solveLPProblem(-1, model.c, model.cM, model.cB, model.S, model.b, model.lb, model.ub, params);
    
    if (nargin <= 2) || (dpcDevFromObjVal <= 0)
        S = [model.S; model.c'];
        b = [model.b; origObjVal];
        cM = model.cM;
        cB = model.cB;
    else
        dpcDevFromObjVal = dpcDevFromObjVal / 100;
        S = model.S;
        b = model.b;
        cM = [model.cM; model.c'; -model.c'];
        cB = [model.cB; (1+dpcDevFromObjVal) * origObjVal; -(1-dpcDevFromObjVal) * origObjVal];
    end
    
    inds = abs(model.ub - model.lb) > 0;
    maxFlux = zeros(length(model.c), 1);
    
    maxFlux(~inds) = model.lb(~inds);
    minFlux = maxFlux;
    
    inds = find(inds);
    
    for i = 1:length(inds)
        c = zeros(length(model.c), 1);
        c(inds(i)) = 1;
        [sol, maxFlux(inds(i))] = solveLPProblem(-1, c, cM, cB, S, b, model.lb, model.ub, params);
        [sol, minFlux(inds(i))] = solveLPProblem(1, c, cM, cB, S, b, model.lb, model.ub, params);
    end
end
