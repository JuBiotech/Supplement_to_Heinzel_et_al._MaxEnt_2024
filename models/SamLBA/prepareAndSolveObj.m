function roProb = prepareAndSolveObj( mode, addMode, model, dpc, epsilon, bmInd )
%PREPAREANDSOLVEOBJ prepares the given model (output of COBRA Toolbox) for
%robust optimization and solves the robust problem.
%The objective function is robustified against perturbation of each
%coefficient by dpc percent.
%
% Parameters:
%   - mode: 'box' or 'budget'. Specifies the uncertainty set.
%   - addMode: 'all' to add all equations, 'bm' to add equations with
%       non-zero biomass coefficient only. Defaults to all.
%   - model: Model structure. Used fields of the structure:
%       o c: Objective function vector as column-vector
%       o S: Stoichiometric matrix
%       o b: Right hand side of the system S*v = b
%       o lb,ub: Lower and upper bound of flux vector v as column vector
%       o pools or mets: Cell array of pool names
%       o rxns or reactions: Cell array of reaction names
%   - dpc: Disturbance in percent
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
    
    if isempty(addMode)
        addMode = 'all';
    end
    
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

    if nargin < 6
        bmInd = findBiomassFlux(model);
    else
        bmInd = findBiomassFlux(model, bmInd);
    end
    
    % Solve the system for biomass flux
    if strcmpi(addMode, 'all')
        bmObj = sum(model.S);
    elseif strcmpi(addMode, 'bm')
        bmObj = sum(model.S(model.S(:, bmInd) ~= 0, :));
    end
    
    bmObj = full(bmObj)';
    bmObj = bmObj ./ bmObj(bmInd);
    bmObj(bmInd) = 0;
    model.c = -bmObj;
    
    % Disable biomass
    model.lb(bmInd) = 0;
    model.ub(bmInd) = 0;
    
    %-------------------------------------------------
    % Experimental: Add outflow for biomass components
    %-------------------------------------------------

%     bmInds = find(model.S(:, bmInd) ~= 0);
%     for i = 1:length(bmInds)
%         if isfield(model, 'rxns')
%             model.rxns{end+1} = ['vBMOut' num2str(i)];
%         else
%             model.reactions{end+1} = ['vBMOut' num2str(i)];
%         end
%         coeff = model.S(bmInds(i), bmInd);
%         
%         model.S(bmInds(i), end+1) = sign(coeff);
%         model.c(end+1) = 0;
% 
%         model.lb(end+1) = 0;
%         model.ub(end+1) = 0.837 * abs(coeff);
%     end
    
    % Biomass equation in objective function!
	inds = find(model.c ~= 0);
	
    % We need to construct a structure U carrying all information of
    % uncertain coefficients
    % In this case, the structure has only one field: target. The
    % target-field is a cell-array of strucutres. One structure defines an
    % uncertain coefficient by giving its index in the objective-vector and
    % its boundaries (bounds: lower, upper).
    
	U.target = cell(length(inds), 1);
	for i = 1:length(inds)
		element.index = inds(i);
		if model.c(inds(i)) > 0
			element.bounds = model.c(inds(i)) * [1-dpc, 1+dpc];
		else
			element.bounds = model.c(inds(i)) * [1+dpc, 1-dpc];
		end
		U.target{i} = element;
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
	
    %---------------------------------------------
    % The following part is optional.
    % Computation of the nearest nominal solution.
    %---------------------------------------------
    
	% Remove linear dependent rows for QP solver libQuadProgPP
% 	[R, lirows]=rref(model.S');
% 	S = model.S(lirows, :);
% 
% 	[ obj, roProb.nearestSol, status2 ] = libQuadProgPP( 2.*eye(length(roProb.realSol)), -2.*roProb.realSol, ...
% 		S, zeros(size(S, 1), 1), [], [], roProb.rawLower, roProb.rawUpper);
end

