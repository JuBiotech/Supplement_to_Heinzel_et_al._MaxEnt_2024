function model = reconstructProblemFromSample( model, groups, sample )
%RECONSTRUCTPROBLEMFROMSAMPLE Reconstructs an optimization problem from a given sample.
%
%	prob = reconstructProblemFromSample( model, groups, sample )
%
% Note: See mcGroups for hints on the parameters.
%
% This function essentially takes a sample and reconstructs the FBA problem
% with perturbed data solved by the mcGroups or mcGroupsParallel function.
%
% The groups have to be the same as in the call of mcGroups.
%
% Parameters:
%	- model: SBML model of a network. A structure with the following
%		fields:
%		o S: Stoichiometric matrix
%		o c: Target function vector
%		o b: Right hand side of the equation S * v = b. 
%			(Equals zero for FBA => Steady State)
%		o lb: Lower bound on the fluxes
%		o ub: Upper bound on the fluxes
%	- groups: Cell array of group description. A group description is a
%		structure with the following fields:
%		o row: Row of the affected coefficient(s) in the matrix (can be
%			vector).
%		o column: Column of the affected coefficient(s) in the matrix (can
%			be vector).
%		o relative: True / False. If true, the elements of the matrix get
%			multiplied by a drawn random value. Otherwise they will be
%			overwritten.
%		o mean: Mean of the normally distributed random variable or middle
%			of the range of the uniformly distributed variable.
%		o sigma: Standard deviation of the normally distributed variable or
%			half interval width of the uniformly distributed one.
%		o bounds: [bottom, top] Boundaries of the range of the uniform
%			distribution used when relative = false.
%		o isInteger: True / False. If true, the uniform random variable
%			will be integer. Effective when relative = false only.
%		o nonZero: True / False. If true, the uniform random variable will
%			not be zero. Effective when isInteger = true only.
%	- sample: Sample vector created by mcGroups or mcGroupsParallel.
%
% Returns a structure with the following fields:
%	- c: Objective function.
%	- S: Stoichiometric matrix.
%   - b: Right hand side of steady-state S * v = b
%   - cM: Inequality matrix.
%   - cB: Right hand side of inequality constraint cM * v <= cB
%	- lb: Lower flux bound.
%	- ub: Upper flux bound.
%	- sol: Solution of the problem.
%   - objVal: Objective value of the solution.

    % Check constraints
    if ~isfield(model, 'cM')
        model.cM = [];
    end
    if ~isfield(model, 'cB')
        model.cB = [];
    end
		
    offset = 1;

    % ------------------
    glpkParams.dual = 1;
    % ------------------

    % Obtain perturbed problem
    for k = 1:length(groups)
        g = groups{k};

        [model, data] = modifyGroupCoefficients(model, g, sample(k));

        data = reshape(data, numel(data), 1);
        offset = offset + length(data);
    end

    % Solve
    [model.solution, model.objVal, model.valid, model.solMetaInfo] = solveLPProblem(-1, model.c, ...
        model.cM, model.cB, model.S, model.b, model.lb, model.ub, glpkParams);
end

