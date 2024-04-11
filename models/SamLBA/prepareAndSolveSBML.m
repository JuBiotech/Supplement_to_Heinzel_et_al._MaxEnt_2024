function roProb = prepareAndSolveSBML( mode, model, dpc, tol, epsilon, bmInd )
%PREPAREANDSOLVESBML prepares the given model (output of COBRA Toolbox) for
%robust optimization and solves the robust problem.
%The biomass equation (column 13 of S in e.coli) is robustified against 
%perturbation of each coefficient by dpc percent.
%
% Note: The steady-state equation S*v = b is replaced by the system
%   S * v <= b  and  -S * v <= -b
%   Instead of robustifying the equation, the system is robustified.
%
% Parameters:
%   - mode: 'box' or 'budget'. Specifies the uncertainty set.
%   - model: Model structure. Used fields of the structure:
%       o c: Objective function vector as column-vector
%       o S: Stoichiometric matrix
%       o b: Right hand side of the system S*v = b
%       o lb,ub: Lower and upper bound of flux vector v as column vector
%       o pools or mets: Cell array of pool names
%       o rxns or reactions: Cell array of reaction names
%   - dpc: Disturbance in percent
%   - tol: Tolerance for equality constraints. An uncertain equality
%       constraint is replaced by two inequalities with a tolerance:
%           s^T * v = 0  =>  s^T * v <= tol and s^T * v >= -tol
%   - epsilon: Used for budget mode
%   - bmInd: Index or name of biomass flux. (Optional)
%
% Returns: Structure with the fields
%   - sol: Robust solution with slack variables
%   - realSol: Robust solution without slack variables. Typically the
%       desired output of the robust optimization
%   - realSolutionIndices: Vector with indices of the real data of the
%       sol-vector. Other entries of the sol-vector are slack variables
%   - objVal: Robust objective value
%   - certSol: Solution using unperturbed data
%   - certObjVal: Objective value using unperturbed data
%   - nearestSol: Nearest (euclidean norm) unperturbed solution

	dpc = dpc / 100;

    % Check vectors
    if size(model.c, 2) > 1
        model.c = model.c';
    end
    if size(model.lb, 2) > 1
        model.lb = model.lb';
    end
    if size(model.ub, 2) > 1
        model.ub = model.ub';
    end

    % Transfer biomass equation from column to objective function
    if nargin < 6
        bmInd = findBiomassFlux(model);
    else
        bmInd = findBiomassFlux(model, bmInd);
    end
    
    % Biomass equation in column bmInd
	inds = find(model.S(:, bmInd) ~= 0);
	
    % We need to construct a structure U carrying all information of
    % uncertain coefficients
    % In this case, the structure has only one field: stoichiometry. The
    % stoichiometry-field is a cell-array of strucutres. One structure 
    % defines an uncertain coefficient by giving its reaction name (or 
    % index), its pool name (or index) and its boundaries (bounds: lower, 
    % upper).

    U.stoichiometry = cell(length(inds), 1);
	for i = 1:length(inds)
		element.reaction = model.rxns(bmInd);
		element.pool = model.mets(inds(i));
		if model.S(inds(i), bmInd) > 0
			element.bounds = model.S(inds(i), bmInd) * [1-dpc, 1+dpc];
		else
			element.bounds = model.S(inds(i), bmInd) * [1+dpc, 1-dpc];
        end
        element.rhsTolerance = [-tol; tol];
		U.stoichiometry{i} = element;
	end

	roProb = prepareForRO(mode, model, model.c, '', model.lb, model.ub, U, epsilon);

	% Solve
    params.iterCount = true;
    [roProb.sol, roProb.objVal, status, roProb.solverMeta] = glpkInterface(-1, roProb.objective, ...
        roProb.iqA, roProb.iqB, roProb.eqA, roProb.eqB, roProb.lower, roProb.upper, params);
    roProb.realSol = roProb.sol(roProb.realSolutionIndices);

    if status ~= 1
        roProb.objVal = NaN;
    end

    [roProb.certSol, roProb.certObjVal] = glpkInterface(-1, model.c, ...
		[], [], roProb.S, zeros(size(roProb.S, 1), 1), roProb.rawLower, roProb.rawUpper);
end

