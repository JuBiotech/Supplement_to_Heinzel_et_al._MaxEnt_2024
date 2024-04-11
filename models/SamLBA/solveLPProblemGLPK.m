function [sol, val, status, meta] = solveLPProblemGLPK(osense, target, Ain, bin, Aeq, beq, lb, ub, params)
%SOLVELPPROBLEM uses GLPK to solve the supplied LP: max / min target' * x
%subject to Ain * x <= bin, Aeq * x = beq, lb <= x <= ub
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
    
    % Solve standard problem
    [sol, val, status, meta] = glpkInterface(osense, target, Ain, bin, Aeq, beq, lb, ub, params);
    
    % Minimize 1-norm
    if params.minNorm
        
        if isfield(params, 'objTol')
            objTol = params.objTol;
        else
            objTol = 1e-8;
        end
%        params = [];
%        params.dual = 2;
        params.msglev = 0;
        
        % Minimize |v|_1
        % Subject to
        %   c^T * v >= val
        %   S * v = 0
        %   N * v <= b
        %   lb <= v <= ub
        %
        % Now transform this into an LP-Problem:
        % Layout: [v ; mu]
        %
        % Minimize sum( mu )
        % Subject to
        %   c^T * v >= val
        %   S * v = 0
        %   N * v <= b
        %   lb <= v <= ub
        %   v <= mu
        %   -v <= mu
        %   mu >= 0
        
        nVars = length(target);
        
        if size(target, 1) == 1
            target = target';
        end
        
        if ~isempty(Ain)
            Ain = [Ain, zeros(size(Ain)) ; ...
                [osense * target', zeros(1, nVars)] ; ...
                eye(nVars), -eye(nVars); ...
                -eye(nVars), -eye(nVars)];
            bin = [bin ; osense * val + objTol; zeros(2 * nVars, 1)];
        else
            Ain = [[osense * target', zeros(1, nVars)] ; ...
                eye(nVars), -eye(nVars); ...
                -eye(nVars), -eye(nVars)];
            bin = [osense * val + objTol ; zeros(2 * nVars, 1)];
        end
        
        Aeq = [Aeq, zeros(size(Aeq))];

        target = [zeros(nVars, 1) ; ones(nVars, 1)];
        
        if ~isempty(ub)
            ub = [ub ; max(abs(ub), abs(lb))];
        end
        lb = [lb ; zeros(length(lb), 1)];
        
        sol = glpkInterface(1, target, Ain, bin, Aeq, beq, lb, ub, params);
        sol = sol(1:nVars);
    end
end
