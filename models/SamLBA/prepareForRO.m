function [ roProblem ] = prepareForRO( mode, Sfile, target, C, lower, upper, U, params )
%PREPAREFORRO translates the given problem to an Robust Optimization
%Problem.
%
% Background: Problem max f'*u subject to Su = 0 (Steady state), Cu <= b,
% lower <= u <= upper is prepared for robust optimization.
%
% Parameters:
%   - mode: String. Modes are "box" for box-uncertainty, "budget" for
%       budget-uncertainty.
%   - Sfile: Filename of HDF5 file with stoichiometric matrix. Or structure 
%       with stoichiometry (S), reactions (rxns, cell array with names) and 
%       pools (mets, cell array with names). Accepts SBML models.
%   - target: Target function given as string (linear combination of
%       fluxes) or given as vector.
%   - C: Constraints. String format: '12.7v - 2e-1 u_3 >= 3; 2x >= 3' using
%       the variable names (v, u_3, x) given by Sfile.
%       Structure format: C.C * x <= C.b
%   - lower: Lower bounds on variables
%   - upper: Upper bounds on variables
%   - U: Uncertainties structure (see parseUncertainty.m for details)
%       |- Target function (Index) = U.target
%       |- Stoichiometric matrix (Reaction, Pool) = U.stoichiometry
%       |- Constraints (Index, Index) = U.constreq, U.constrrhs
%       |- Lower and upper bounds (Index each) = U.upper, U.lower
%   - params: Only needed for budget mode. Epsilon represents the accepted 
%       maximum probability of constraint violation.
%
% Each uncertainty with bounds.
%
% Returns a structure containing all necessary information to solve the
% problem using a LP-solver:
%   max objective^T * x
%   subject to  eqA * x = eqB 
%               iqA * x <= iqB
%               lower <= x <= upper
%   The solution is given by x(realSolutionIndices), which discards slack
%   variables.

    % Step 1: Fetch matrix, pools and reactions
    if ~ischar(Sfile)
		
		% Check for SBML file
		if isfield(Sfile, 'rxns') && isfield(Sfile, 'S') && isfield(Sfile, 'mets')
			S = Sfile.S;
			pools = Sfile.mets;
			reactions = Sfile.rxns;
		else
			S = Sfile.stoichiometry;
			pools = Sfile.pools;
			reactions = Sfile.reactions;
		end
    else
        S = hdf5read(Sfile, '/stoichiometry/matrix');
        [pools, reactions] = getPoolsAndReactions(Sfile);
    end

    % Step 2: Build target function
    if ischar(target)
        f = parseTargetFunction(target, reactions);
    else
        f = target;
    end

    % Step 3: Parse Constraints
    [cons, consB] = parseConstraints(C, reactions);

	roProblem.rawUpper = upper;
	roProblem.rawLower = lower;
	roProblem.S = S;

    % Harden boundaries
    if isfield(U, 'upper') && ~isempty(U.upper)
        upper = applyStrictBounds(upper, U.upper, @min);
    end
    if isfield(U, 'lower') && ~isempty(U.lower)
        lower = applyStrictBounds(lower, U.lower, @max);
    end
    
    if strcmp(mode, 'box')
        [roProblem.objective, roProblem.iqA, roProblem.iqB, roProblem.eqA, ...
            roProblem.eqB, roProblem.lower, roProblem.upper, roProblem.realSolutionIndices] = ...
            prepareForBoxRO(U, f, S, reactions, pools, cons, consB, lower, upper);
    elseif strcmp(mode, 'budget')
        [roProblem.objective, roProblem.iqA, roProblem.iqB, roProblem.eqA, ...
            roProblem.eqB, roProblem.lower, roProblem.upper, roProblem.realSolutionIndices] = ...
            prepareForBudgetRO(U, f, S, reactions, pools, cons, consB, lower, upper, params);
    end

    % Sanitize matrices
    if isempty(roProblem.iqA)
        roProblem.iqA = [];
        roProblem.iqB = [];
    end
    if isempty(roProblem.eqA)
        roProblem.eqA = [];
        roProblem.eqB = [];
    end
    
end

function bounds = applyStrictBounds(bounds, uncertain, evaluator)
%APPLYSTRICTBOUNDS applies the strictest boundaries to the given bounds
%vector using the information supplied by the uncertainty cell array.
%The evaluator is used to select the strictest boundary.

    for i = 1:length(uncertain)
        bounds(uncertain{i}.index) = evaluator(uncertain{i}.bounds);
    end
end
