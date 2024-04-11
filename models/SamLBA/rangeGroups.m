function res = rangeGroups( model, groups, parallel, closePool, nWorkers )
%RANGEGROUPS Scans a parameter range solving an FBA problem for each
%parameter set.
%
%	res = rangeGroups( model, groups )
%
% Some elements of the problem (stoichiometry, objective function, boundary
% vectors or constraint system) of the network are scanned in a given range
% and for each parameter set the corresponding FBA problem is solved.
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
%			multiplied by the current value of the range. Otherwise they 
%			will be overwritten. F.e. use relative = true and [0.8, 1.2, 3]
%           to scan the coefficient group from 80% to 120% in 3 steps. Set
%           relative = false to scan the values 0.8, 1, 1.2 (absolute).
%		o range: [bottom, top, steps] Range which is scanned given by low
%			value, high value and steps to scan. linspace() is used to
%			generate the scanned values.
%		o rawRange: Optional. Vector with values to be scanned. Will ignore
%           range if this field is set.
%   - parallel: Optional. If true parallelization is used. Default is true.
%   - closePool: True if the matlabpool shall be closed after computations
%       are finished. False to keep it as is. Optional, default is false.
%   - nWorkers: Number of tasks the Parallel Computation Toolbox uses. This
%       is optional. Leave out to use the default number of tasks.
%
% Returns a structure with the following fields:
%	- origVal: Original objective function value of the unperturbed FBA problem.
%	- origSol: Original solution to unperturbed FBA problem.
%	- origData: Coefficients of the stoichiometric matrix before perturbance.
%	- vals: nD array of objective function values of the ranged problems.
%       NaN indicates an invalid solution (problem infeasible, unbounded
%       etc.). Dimensions are given by the steps (or length of rawRange) of
%       the groups. The range is scanned in ascending order.
%	- diffs: nD array of the Euclidian distances of the original solution
%		vector to a ranged solution. ||origSol, sols(:, i)||_2
%       See vals for explanation of nD array.
%	- sols: Matrix in which each column contains a solution of a ranged 
%		problem. Use ind2sub( lengthOfRanges ) to get the corresponding
%		parameter values.
%	- samples: Matrix in which each column contains the changed
%		coefficients of the stoichiometric matrix after ranging.
%		The values are set in the same order as the groups. If a group 
%		perturbs more than one coefficient, the values are saved in 
%		column-major form. See sols on how to convert linear index to
%		parameters.

    if nargin <= 2
        parallel = true;
    end
    if nargin <= 3
        closePool = false;
    end

    % Check constraints
    if ~isfield(model, 'cM')
        model.cM = [];
    end
    if ~isfield(model, 'cB')
        model.cB = [];
    end
    
    % ------------------
    glpkParams.dual = 1;
    % ------------------
	[origSol, origVal] = solveLPProblem(-1, model.c, model.cM, model.cB, model.S, model.b, model.lb, model.ub, glpkParams);

    % Calculate nominal weights
    wghtOrig = abs(origSol);
    wghtOrig(wghtOrig <= 1e-8) = 1;
    
	% Assign generators
	[groups, sampleSize, res.origData] = checkGroupsAndSaveOriginalData(model, groups, 'uniform');
	
    % Extract ranges
    ranges = cell(length(groups), 1);
    dimSizes = zeros(1, length(groups));
    linearizedLength = 1;
    for i = 1:length(groups)
        g = groups{i};
        
        if ~isfield(g, 'rawRange')
            ranges{i} = linspace(g.range(1), g.range(2), g.range(3));
        else
            ranges{i} = g.rawRange;
        end
        dimSizes(i) = length(ranges{i});
        linearizedLength = linearizedLength * length(ranges{i});
    end

    % Create linear arrays to store function values, samples, solutions,
    % diffs
    vals = zeros(linearizedLength, 1);
    diffs = zeros(linearizedLength, 1);
    samples = zeros(sampleSize, linearizedLength);
    sols = zeros(length(model.lb), linearizedLength);
    weightedDiffs = zeros(linearizedLength, 1);

    if parallel
        if matlabpool('size') == 0
            if (nargin <= 4) || isempty(nWorkers)
                matlabpool open;
            else
                matlabpool('open', nWorkers);
            end
        end
        
        parfor i = 1:linearizedLength
            % Map linear index to multidim array index
            indices = cell(length(ranges), 1);
            [indices{:}] = ind2sub(dimSizes, i);
            indices = cell2mat(indices);

            modelMod = model;
            sample = zeros(sampleSize, 1);

            offset = 1;
            % Obtain grid point
            for k = 1:length(groups)
                g = groups{k};

                data = ranges{k}(indices(k));
                [modelMod, data] = modifyGroupCoefficients(modelMod, g, data); 

                % Save point
                data = reshape(data, numel(data), 1);
                sample(offset:offset + length(data) - 1) = data;
                offset = offset + length(data);
            end

            samples(:, i) = sample;
            
            % Solve
            [sols(:, i), vals(i), valid, metaInfo] = solveLPProblem(-1, modelMod.c, ...
                modelMod.cM, modelMod.cB, modelMod.S, modelMod.b, modelMod.lb, modelMod.ub, glpkParams);
            diffs(i) = norm(origSol - sols(:, i));
            weightedDiffs(i) = sqrt(sum( ((origSol - sols(:, i)) ./ wghtOrig).^2 ));
%            res.meta(i, :) = [metaInfo.time, metaInfo.iterCount];

            if valid ~= 1
                vals(i) = NaN;
            end
        end
        
        if (matlabpool('size') > 0) && closePool
            matlabpool close;
        end
    else
        for i = 1:linearizedLength
            % Map linear index to multidim array index
            indices = cell(length(ranges), 1);
            [indices{:}] = ind2sub(dimSizes, i);
            indices = cell2mat(indices);

            modelMod = model;

            offset = 1;
            % Obtain grid point
            for k = 1:length(groups)
                g = groups{k};

                data = ranges{k}(indices(k));
                [modelMod, data] = modifyGroupCoefficients(modelMod, g, data); 

                % Save point
                data = reshape(data, numel(data), 1);
                samples(offset:offset + length(data) - 1, i) = data;
                offset = offset + length(data);
            end

            % Solve
            [sols(:, i), vals(i), valid, metaInfo] = solveLPProblem(-1, modelMod.c, ...
                modelMod.cM, modelMod.cB, modelMod.S, modelMod.b, modelMod.lb, modelMod.ub, glpkParams);
            diffs(i) = norm(origSol - sols(:, i));
            weightedDiffs(i) = sqrt(sum( ((origSol - sols(:, i)) ./ wghtOrig).^2 ));
%            res.meta(i, :) = [metaInfo.time, metaInfo.iterCount];

            if valid ~= 1
                vals(i) = NaN;
            end
        end
    end
    
    if length(dimSizes) > 1
        res.diffs = reshape(diffs, dimSizes);
        res.weightedDiffs = reshape(weightedDiffs, dimSizes);
        res.vals = reshape(vals, dimSizes);
    else
        res.diffs = diffs;
        res.weightedDiffs = weightedDiffs;
        res.vals = vals;
    end
    res.samples = samples;
    res.origSol = origSol;
    res.origVal = origVal;
    res.sols = sols;
end
