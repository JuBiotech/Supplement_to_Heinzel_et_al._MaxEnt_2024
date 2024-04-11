function res = inverseMCGroupsCoeffEstimation( model, groups, nSamples, mode, normPC, singlePC, nSlice, closePool, nWorkers )
%INVERSEMCGROUPSCOEFFESTIMATION Estimates an axis aligned bounding box of the robust coefficient space.
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
%   - closePool: True if the matlabpool shall be closed after computations
%       are finished. False to keep it as is. Optional, default is false.
%
% Returns a structure with the following fields:
%	- origVal: Original objective function value of the unperturbed FBA problem.
%	- origSol: Original solution to unperturbed FBA problem.
%	- origData: Coefficients of the stoichiometric matrix before perturbance.
%	- samples: Matrix in which each column contains the changed
%		coefficients of the stoichiometric matrix after perturbance.
%		The values are set in the same order as the groups. If a group 
%		perturbs more than one coefficient, the values are saved in 
%		column-major form.

	if isempty(mode)
		mode = 'normal';
	end

    if nargin <= 7
        closePool = true;
    end
    
    if nargin <= 6
        nSlice = nSamples / 100;
    end    
    
    % Check constraints
    if ~isfield(model, 'cM')
        model.cM = [];
    end
    if ~isfield(model, 'cB')
        model.cB = [];
    end
    
	[origSol, origVal] = solveLPProblem(-1, model.c, model.cM, model.cB, model.S, model.b, model.lb, model.ub);

    diffSol = abs(origSol);
    diffSol(diffSol == 0) = 1;
    diffSol = diffSol * singlePC / 100;    
    
    normPC = norm(origSol) * normPC / 100;
    
	% Calculate sample size and assign generators
	[groups, sampleSize, res.origData] = checkGroupsAndSaveOriginalData(model, groups, mode);
    
    % ------------------
    glpkParams.dual = 1;
    % ------------------
    
    tHandle = tic;
    maxVal = -inf(sampleSize, 1);
    minVal = inf(sampleSize, 1);
    nValid = 0;
    
    % Slice
    for j = 1:floor(nSamples / nSlice)

        if matlabpool('size') == 0
            if (nargin <= 8) || isempty(nWorkers)
                matlabpool open;
            else
                matlabpool('open', nWorkers);
            end
        end
    
        valid = false(nSlice, 1);
        samples = zeros(sampleSize, nSlice);
        
        % Sample main loop
        parfor i = 1:nSlice

            modelMod = model;
            offset = 1;
            sample = zeros(sampleSize, 1);

            % Obtain sample
            for k = 1:length(groups)
                g = groups{k};

                data = g.generator(g);
                [modelMod, data] = modifyGroupCoefficients(modelMod, g, data); 

                % Save sample
                data = reshape(data, numel(data), 1);
                sample(offset:offset + length(data) - 1) = data;
                offset = offset + length(data);
            end

            % Solve
            [sol, ~, valid(i)] = solveLPProblem(-1, modelMod.c, ...
                modelMod.cM, modelMod.cB, modelMod.S, modelMod.b, modelMod.lb, modelMod.ub, glpkParams);

            % Check solution
            valid(i) = valid(i) && all(abs(sol - origSol) <= diffSol);

            if isfinite(normPC)
                valid(i) = valid(i) && (norm(sols - origSol) <= cmpNorm);
            end
            
            samples(:,i) = sample;
        end
        matlabpool close;
        
        for i = 1:nSlice
            if valid(i)
                maxVal = max(samples(:,i), maxVal);
                minVal = min(samples(:,i), minVal);
            end
        end

        nValid = nValid + sum(valid);
        
        fprintf('Status: %f  - Valid up to now: %g\n', j * nSlice / nSamples * 100, nValid / nSamples * 100);
        
    end
    
    table = cell(sampleSize, 1);
    for i = 1:sampleSize
        table(i) = model.mets(groups{i}.row);
    end
    
    tElapsed = toc(tHandle);
    fprintf('Elapsed time: %g sec = %g h %g min\n', tElapsed, floor(tElapsed / 3600), floor((tElapsed - floor(tElapsed / 3600)*3600) / 60));
    
    if closePool && (matlabpool('size') > 0)
        matlabpool close;
    end
    
    % Assign res structure
    res.maxVal = maxVal;
    res.minVal = minVal;
    res.validPC = nValid / nSamples * 100;
    res.origSol = origSol;
    res.origVal = origVal;
    res.table = [{'Group'}, {'Nominal'}, {'Min'}, {'Max'}, {'Width'}, {'RelWidthPC'}; ...
        table, num2cell(res.origData), num2cell([res.minVal, res.maxVal, abs(res.maxVal - res.minVal), abs(res.maxVal - res.minVal) ./ abs(res.origData) * 100])];
end
