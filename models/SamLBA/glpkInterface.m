function [sol, val, status, meta] = glpkInterface(osense, target, Ain, bin, Aeq, beq, lb, ub, params)
%GLPKINTERFACE uses GLPK to solve the supplied LP: max / min target' * x
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
%
% Returns
%   - sol: Solution vector
%   - val: Objective value
%   - status: Status (1 if optimum has been found)
%   - meta: Meta data like reduced costs, shadow prices etc.

    % Build matrix
    A = [Aeq; Ain];
    b = [beq; bin];
    b = full(b);
    
    % Check values
    if any(any(~isfinite(A))) || any(~isfinite(b))
        error('GLPK:NOTFINITE', 'A matrix or righthandside is invalid. It contains NaN or Inf.');
    end
    
    csense = [repmat('S', size(Aeq, 1), 1); repmat('U', size(Ain, 1), 1)];

    if (nargin <= 8) || ~isfield(params, 'msglev')
        params.msglev = 0;
    end
    
    if (nargin <= 8) || ~isfield(params, 'lpsolver')
        params.lpsolver = 1;
    end

    % Old way of calling glpk
    %[x,f,stat,extra] = glpkmex(osense,c,A,b,csense,lb,ub,[],params);

    if (nargin <= 8) || ~isfield(params, 'iterCount')
        [sol,val,origStat,extra] = glpk(target,A,b,lb,ub,csense,[],osense,params);
        meta.iterCount = NaN;
        timeMatlab = 0;
    else
        params.msglev = 3;
        
        tStart = tic;
        % Extract iteration count
        t = evalc('[sol,val,origStat,extra] = glpk(target,A,b,lb,ub,csense,[],osense,params);');
        timeMatlab = toc(tStart);
        
        indA = strfind(t, 'obj =');
        if ~isempty(indA)
            indA = indA(end);
            indNL = strfind(t, 10);
            if ~isempty(indNL)
                indNL = max(indNL(indNL < indA));
                t = t((indNL+2):(indA-3));
                meta.iterCount = str2double(t);
            else
                meta.iterCount = NaN;
            end
        else
            meta.iterCount = NaN;
        end
    end

    % Note that status handling may change (see glplpx.h)
    if (origStat == 180 || origStat == 5)
        stat = 1; % Optimal solution found
    elseif (origStat == 182 || origStat == 183 || origStat == 3 || origStat == 110)
        stat = 0; % Infeasible
    elseif (origStat == 184 || origStat == 6)
        stat = 2; % Unbounded
    else
        stat = -1; % Solution not optimal or solver problem
    end

    
    status = stat;
	if isfield(extra, 'lambda')
    	meta.dual = extra.lambda;
	else
		meta.dual = NaN;
	end

	if isfield(extra, 'redcosts')
	    meta.reduced = extra.redcosts;
	else
		meta.reduced = NaN;
	end

	if isfield(extra, 'time')
	    meta.time = max(timeMatlab, extra.time);
	else
		meta.time = timeMatlab;
	end

	if isfield(extra, 'memory')
	    meta.memory = extra.memory;
	else
		meta.memory = NaN;
	end

end
