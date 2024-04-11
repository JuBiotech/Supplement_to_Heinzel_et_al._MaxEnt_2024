function res = mcGroups( model, groups, nSamples, mode, plugins )
%MCGROUPS Performs a Monte-Carlo simulation on the given network.
%
%	res = mcGroups( model, groups, nSamples, mode )
%
% Some elements of the stoichiometric matrix (or boundary vectors,
% constraint system or objective vector) of the network are perturbed
% and the FBA problem is solved. The elements are specified by the groups
% allowing for disturbing groups of elements by the same factor.
%
% There are three modes available which determine the distribution of the
% disturbance: 'uniform', 'normal', 'truncnormal'.
%
% If the 'normal' mode is chosen, the elements will be multiplied by the
% factor epsilon where epsilon is normally distributed with mean and sigma
% given by the group.
%
% The 'truncnormal' behaves just like the 'normal' mode, but restricts the  
% samples to the interval mean +/- 2*sigma.
%
% The 'uniform' mode has several configuration options which are set by the 
% group:
%	- relative: If true, the elements will be multiplied by a random value 
%		of the interval (mean - 2*sigma; mean + 2*sigma).
%	- If relative is false, the elements will be set to a value determined
%		by the following options.
%	- bounds: The elements are set to a random value drawn from the
%		interval (bounds(1), bounds(2)). The former values of the elemnts
%		are discarded.
%	- isInteger: If true, the drawn random value will be integer.
%	- nonZero: If true, the drawn random value will not be zero. This
%		option is only activated when isInteger = true.
%
% If the mode parameter is left out, the function defaults to 'normal'.
%
% Parameters:
%	- model: SBML model of a network. A structure with the following
%		fields:
%		o S: Stoichiometric matrix
%		o c: Target function vector
%		o b: Right hand side of the equation S * v = b. 
%			(Equals zero for FBA => Steady State)
%       o cM: Constraint matrix (optional)
%       o cB: Right hand side of the inequality cM * v <= cB (optional)
%		o lb: Lower bound on the fluxes
%		o ub: Upper bound on the fluxes
%	- groups: Cell array of group description. A group description is a
%		structure with the following fields:
%       o type: Type string of coefficient group. Use name of field in
%           model. Example: Perturbation of constraint matrix cM needs type
%           string 'cM'. Default is 'S' for stoichiometric matrix.
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
%       o truncated: True / False. If true, the samples will be drawn from
%           a truncated normal distribution. The valid interval is set to 
%           mean +/- 2*sigma.
%           If false, the "normal" normal distribution will be used.
%       o truncateHalf: True / False. If true the standard deviation will
%           be halved and the valid interval is set to mean +/- sigma.
%           This implies trucated = true.
%	- nSamples: Number of samples of the MC simulation.
%	- mode: 'normal' / 'uniform' / 'truncnormal'. Sets the mode of the 
%       perturbance. Specifying 'truncnormal' overrides the 'truncated'
%       field in the groups structure.
%   - plugins: A cell array of plugins. A plugin is a function handle with
%       the following signature: data = plugin(S, solution, objValue, diff)
%       The result of the plugin will be saved in a cell array.
%       There are no requirements for the return values of a plugin as
%       everything is stored in a cell array.
%       This field is optional.
%
% Returns a structure with the following fields:
%	- origVal: Original objective function value of the unperturbed FBA problem.
%	- origSol: Original solution to unperturbed FBA problem.
%	- origData: Coefficients of the stoichiometric matrix before perturbance.
%	- vals: Vector of objective function values of the sampled problems.
%       NaN indicates an invalid solution (problem infeasible, unbounded
%       etc.).
%	- sols: Matrix in which each column contains a solution of a sampled 
%		problem.
%	- diffs: Vector with the Euclidian distances of the original solution
%		vector to a sampled solution. ||origSol, sols(:, i)||_2
%	- weightedDiffs: Matrix with the weighted Euclidian distances of the 
%       original solution vector to a sampled solution. First column
%       contains the diffs with nominal weights (nominal flux solution),
%       second column contains diffs with std.-weights (std of all fluxes
%       is taken).
%	- maxVal: The maximum absolute deviation of the objective function
%		values from the unperturbed objective function value. 
%		max( abs(vals(i) - objVal) )
%	- maxDiff: The maximum absolute deviation of the Euclidian distances.
%		max( diffs(i) )
%	- maxWeightedDiff: The maximum weighted Euclidean distance. Vector with
%       entries corresponding to columns of weightedDiffs.
%	- samples: Matrix in which each column contains the changed
%		coefficients of the stoichiometric matrix after perturbance.
%		The values are set in the same order as the groups. If a group 
%		perturbs more than one coefficient, the values are saved in 
%		column-major form.
%   - meta: Matrix in which each row contains the time needed to solve the
%       problem (first column) and the number of iterations (second
%       column).
%   - pluginData: Cell array of plugin data. Each cell holds the data of
%       the plugins returned for the current sample. Each of the nSamples
%       cells contains a cell array with length(plugins) cells which hold
%       the data of the plugin.
%       Note: This field only exists if there were plugins.

	if nargin == 3
		mode = 'normal';
    end
    
    if nargin <= 4
        plugins = [];
    end
    
    % Check constraints
    if ~isfield(model, 'cM')
        model.cM = [];
    end
    if ~isfield(model, 'cB')
        model.cB = [];
    end
    
	[res.origSol, res.origVal] = solveLPProblem(-1, model.c, model.cM, model.cB, model.S, model.b, model.lb, model.ub);
	
    % Calculate nominal weights
    wghtOrig = abs(res.origSol);
    wghtOrig(wghtOrig <= 1e-8) = 1;
    
	res.vals = zeros(nSamples, 1);
	res.sols = zeros(length(model.lb), nSamples);
	res.diffs = zeros(nSamples, 1);
    res.weightedDiffs = zeros(nSamples, 2);
    res.meta = zeros(nSamples, 2);
    
    if ~isempty(plugins)
        res.pluginData = cell(nSamples, 1);
        for i = 1:length(nSamples)
            res.pluginData{i} = cell(length(plugins), 1);
        end
    end
    
	% Calculate sample size and assign generators
	[groups, sampleSize, res.origData] = checkGroupsAndSaveOriginalData(model, groups, mode);
	res.samples = zeros(sampleSize, nSamples);
	
	% Sample main loop
    curPrg = -1;
    glpkParams.iterCount = true;
	for i = 1:nSamples
        
        modelMod = model;
		
		offset = 1;
		
		% Obtain sample
		for k = 1:length(groups)
			g = groups{k};
			
			data = g.generator(g);
            [modelMod, data] = modifyGroupCoefficients(modelMod, g, data); 
            			
			% Save sample
			data = reshape(data, numel(data), 1);
			res.samples(offset:offset + length(data) - 1, i) = data;
			offset = offset + length(data);
		end
		
		% Solve
		[res.sols(:, i), res.vals(i), valid, metaInfo] = solveLPProblem(-1, modelMod.c, ...
            modelMod.cM, modelMod.cB, modelMod.S, modelMod.b, modelMod.lb, modelMod.ub, glpkParams);
		res.diffs(i) = norm(res.origSol - res.sols(:, i));
        res.weightedDiffs(i, 1) = sqrt(sum( ((res.origSol - res.sols(:, i)) ./ wghtOrig).^2 ));
        res.meta(i, :) = [metaInfo.time, metaInfo.iterCount];
        
        if valid ~= 1
            res.vals(i) = NaN;
        end
        
        % Plugins
        if ~isempty(plugins)
            for pI = 1:length(plugins)
                plug = plugins{pI};
                res.pluginData{i}{pI} = plug(S, res.sols(:, i), res.vals(i), res.diffs(i));
            end
        end
        
        % Progress
        prgVal = floor(i / nSamples * 100);
        if (mod(prgVal, 10) == 0) && (prgVal > curPrg)
            display(['Progress (%): ' num2str(prgVal)]);
            curPrg = prgVal;
        end
    end

    % Calculate sigma-weighted difference
    wghtSigma = ones(length(res.origSol), 1);
    for i = 1:length(res.origSol)
        compSigma = std(res.sols(i, :));
        if compSigma > 1e-8
            wghtSigma(i) = compSigma;
        end
    end
    
    for i = 1:nSamples
        cur = res.sols(:, i);
        res.weightedDiffs(i, 2) = sqrt(sum( ((res.origSol - cur) ./ wghtSigma).^2 ));
    end
    
	res.maxVal = max(abs(res.vals - res.origVal));
	res.maxDiff = max(res.diffs);
	res.maxWeightedDiff = [max(res.weightedDiffs(:, 1)), max(res.weightedDiffs(:, 2))];
end
