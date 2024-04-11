function [sol, val, status, meta] = solveLPProblemMOSEK(osense, target, Ain, bin, Aeq, beq, lb, ub, params)
%SOLVELPPROBLEMMOSEK uses MOSEK to solve the supplied LP: max / min target' * x
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

    meta = [];
    [sol, val, status] = solveMosekLP(osense, target, Ain, bin, Aeq, beq, lb, ub, params);
    
    % Minimize 1-norm
    if params.minNorm
        
        if isfield(params, 'objTol')
            objTol = params.objTol;
            params = rmfield(params, 'objTol');
        else
            objTol = 0;% 1e-8;
        end
        
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
        
        sol = solveMosekLP(1, target, Ain, bin, Aeq, beq, lb, ub, params);
        sol = sol(1:nVars);
    end
end

function [sol, val, stat] = solveMosekLP(osense, c, Ain, bin, S, b, lb, ub, params)
%MSK_OPTIMIZER_FREE_SIMPLEX    

    prob.c        = c;
    prob.blx      = lb;
    prob.bux      = ub;
    prob.a        = [];
    prob.blc      = [];
    prob.buc      = [];
    if ~isempty(S) && ~isempty(b)
        prob.a        = [prob.a ; sparse(S)];
        prob.blc      = [prob.blc; b];
        prob.buc      = [prob.buc; b];
    end
    if ~isempty(Ain) && ~isempty(bin)
        prob.a        = [prob.a ; sparse(Ain)];
        prob.blc      = [prob.blc; -inf(length(bin), 1)];
        prob.buc      = [prob.buc; bin];
    end

    % Specify parameters
    params.MSK_IPAR_OPTIMIZER = 7; %'MSK_OPTIMIZER_FREE_SIMPLEX';
    params.MSK_IPAR_LOG = 0;
    
    if isfield(params, 'minNorm')
        params = rmfield(params, 'minNorm');
    end
    
    % Optimize the problem.
    if osense == -1
        [r,res] = mosekopt('maximize', prob, params);
    else
        [r,res] = mosekopt('minimize', prob, params);
    end

    try 
      % Display the optimal solution.
      val = res.sol.bas.pobjval;
      sol = res.sol.bas.xx;
      stat = (r == 0);
    catch
      fprintf('MSKERROR: Could not get solution')
      stat = false;
    end
end