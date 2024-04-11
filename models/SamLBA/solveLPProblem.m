function [sol, val, status, meta] = solveLPProblem(osense, target, Ain, bin, Aeq, beq, lb, ub, params)
%SOLVELPPROBLEM Solves the supplied LP: max / min target' * x
%subject to Ain * x <= bin, Aeq * x = beq, lb <= x <= ub
%
% The solver is specified by the global variable "LPsolver".
%   LPsolver = 1: GLPK (default)
%              2: MOSEK
%
% Parameters:
%   - osense: Objective sense (-1 max, +1 min)
%   - target: Objective function
%   - Ain: Inequality constraint matrix
%   - bin: Inequality right hand side
%   - Aeq: Equality constraint matrix
%   - beq: Equality right hand side
%   - lb: Lower bound of the solution vector
%   - ub: Upper bound of the solution vector
%   - params: Optional. See glpk.m for further information.
%       Set iterCount to true to count iterations.
%       Set minNorm to true to minimize the 1-norm of the flux solution.
%       Set objTol for objective value tolerance in norm minimization:
%           c' * v >= objVal - objTol
%       Defaults to 1e-8.
%
%       Note: minNorm defaults to true. Params may be used as scalar 
%           boolean indicating params.minNorm.
%
% Returns
%   - sol: Solution vector
%   - val: Objective value
%   - status: Status (1 if optimum has been found)
%   - meta: Meta data like reduced costs, shadow prices etc.

    if nargin <= 8
        params.minNorm = true;
    elseif ~isstruct(params)
        minNorm = params;
        params = [];
        params.minNorm = minNorm;
    elseif ~isfield(params, 'minNorm')
        params.minNorm = true;
    end
    
    global LPsolver;
    
    if isempty(LPsolver)
        LPsolver = 1;
    end
    
    % Sanity check
    if any(any(~isfinite(Ain))) || any(~isfinite(bin)) || any(any(~isfinite(Aeq))) || any(~isfinite(beq))
        status = 2;
        val = nan;
        sol = nan(size(target));
        return;
    end
    
    if LPsolver == 1
        [sol, val, status, meta] = solveLPProblemGLPK(osense, target, Ain, bin, Aeq, beq, lb, ub, params);
    else
        [sol, val, status, meta] = solveLPProblemMOSEK(osense, target, Ain, bin, Aeq, beq, lb, ub, params);
    end
    
end
