function [ sol, val, status ] = quickSolveFBA( model, minNorm )
%QUICKSOLVEFBA Solves the FBA problem associated with the given model.
%
% Parameters:
%   - model: Model.
%   - minNorm: Set to true to minimize the 1-norm as a second stage to
%       remove cycles. Optional, defaults to true.
%
% Returns:
%   - sol: Vector with flux distribution
%   - val: Objective function value
%   - status: Status (1 if optimum has been found)

    if (nargin < 2) || isempty(minNorm)
        minNorm = true;
    end
    
    params.minNorm = minNorm;

    % Check constraints
    if ~isfield(model, 'cM')
        model.cM = [];
    end
    if ~isfield(model, 'cB')
        model.cB = [];
    end

    [sol, val, status, meta] = solveLPProblem(-1, model.c, model.cM, model.cB, model.S, model.b, model.lb, model.ub, params);
end

