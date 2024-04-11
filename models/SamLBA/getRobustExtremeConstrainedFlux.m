function [val, objVal, relVal, multi] = getRobustExtremeConstrainedFlux(dir, nomVal, vNom, model, groups, normBound, singleBound, maxIter, tolerance)
%GETROBUSTEXTREMECONSTRAINEDFLUX Does a line search to find the coefficient values for which the optimal solution to the modified problem lies in the neighbourhood of the nominal solution.
%
% This method employs the bisection method to perform a line search on the model 
% coefficients given in the groups cell array. The direction of the search is given by the unit vector dir.
% The method looks for the point where the optimal solution of the FBA problem subject 
% to the modified group coefficients violates the constraints:
%	norm(v_sol - v_nom) <= normBound  and  all(abs(v_sol - v_nom) <= singleBound)
% The constraints are checked subject to a given tolerance (to account for numerics).
%
% The returned coefficient values, objective value function and relative group multiplicators are
% taken at the last feasible position. If the method was not able to reach the robust extreme, it
% returns inf or -inf (or inf in direction of dir). Maybe increasing the maxIter parameter can
% help to attain a robust extreme.
%
% The function works with a relative multiplicator. Instead of maximizing / minimizing the 
% absolute coefficient values, the method max's / min's the factor r in coefficients(group) * (1 + r * dir).
%
% Parameters:
%	- dir: Direction vector (norm = 1). Direction in which the method searches for a robust extreme.
%	- vNom: Nominal flux solution.
%	- model: Model.
%	- groups: Groups to use. See MCGroups.m for hints.
%	- normBound: Constraint on the norm of the flux: norm(vNom - vSample) <= normBound.
%	- singleBound: Constraint on the fluxes: abs(vNom - vSample) <= normBound.
%	- maxIter: Maximum iteration count.
%   - tolerance: Tolerance for feasibility check. Defaults to 1e-12.
%
% Returns:
%	- val: Vector with the absolute coefficient values (in order of the given groups) of the robust extreme.
%	- objVal: Objective function value at the robust extreme.
%	- relVal: Vector with the relative group multiplicators (in order of the given groups) of the robust extreme.
%   - multi: Multiplicator of the current ray.

    % ------------------
%    glpkParams.dual = 1;
%    glpkParams.lpsolver = 3;
    glpkParams = [];
    % ------------------

    if (nargin <= 8) || isempty(tolerance)
        tolerance = 1e-12;
    end
    
    precisionTol = 1e-14;
    
    if any(~isfinite(dir))
        val = nomVal;
        return;
    end

	% State: 
    %   0 - push
    %   1 - bisection to get the precise boundary
    
    state = 0;
    iter = 0;

    % Increment-increase: Paper version = 1
    pIncMult = 10;
    
    lastMult = 0;
    pMult = 0;
    pInc = 0.05;
    a = 0;
    insideObjVal = 0;
    b = precisionTol + 1;
    
    [vNom, ~]  = solveLPProblem(-1, model.c, model.cM, model.cB, model.S, model.b, model.lb, model.ub, glpkParams);
    
    while (iter <= maxIter) && (abs(a - b) >= precisionTol)
        
        % Do a step
		curModel = modifyModelCoefficientsRelative(model, groups, ones(size(dir)) + pMult * dir);
        
    	[sol, objVal] = solveLPProblem(-1, curModel.c, curModel.cM, curModel.cB, curModel.S, curModel.b, curModel.lb, curModel.ub, glpkParams);
        
        if isnan(objVal) && all(isnan(sol))
            break;
        end
        
%        fprintf('Mult: %g ; Norm %g ; Single %g\n', pMult, norm(vNom - sol), any(abs(vNom - sol) > singleBound + tolerance));
        
        % Check conditions
        if isnan(objVal) || (norm(vNom - sol) > normBound) || any(abs(vNom - sol) > singleBound + tolerance)
            % Infeasible: We've gone too far
            if (state == 0)
                % Pull mode
                state = 1;
                a = lastMult;
                b = pMult;
                lastMult = pMult;
                pMult = (a+b) / 2;
            elseif state == 1
				% Bisection mode
                lastMult = pMult;
                b = pMult;
                pMult = (a + pMult) / 2;
            end
        else
            % Still valid: Push to the boundary
            if (state == 0)
				% Push further
                lastMult = pMult;
                pMult = pMult + pInc;
                pInc = pInc * pIncMult;
            elseif state == 1
                a = pMult;
                pMult = (b + pMult) / 2;
            end
            insideObjVal = objVal;
        end
        
%        fprintf('Iter %g - Phase: %g - curVal: %g\n', iter, state, curVal);
        iter = iter + 1;
    end

    if state == 0
		% Did not reach boundary 
        val = getModelCoefficients(model, groups, inf(length(dir), 1));
        relVal = inf(length(dir), 1);
		objVal = inf;
        multi = inf;
    else
        objVal = insideObjVal;
        val = getModelCoefficients(model, groups, ones(size(dir)) + dir * a);
        relVal = dir * a;
        multi = a;
    end
    
%    fprintf('Total iter %g - In Phase 1: %g  - in Phase2: %g \n', iter, temp, iter - temp);
end
