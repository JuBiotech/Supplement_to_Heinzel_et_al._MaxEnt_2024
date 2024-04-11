function res = inverseMCGroupsCoeffEstimationSphere( model, groups, nSamples, normPC, singlePC, maxIter, nSlice, closePool, nWorkers )
	%INVERSEMCGROUPSCOEFFESTIMATIONSPHERE Estimates the robust coefficient
	%space by sampling directions from the unit sphere and using the ray
	%bisection method of getRobustExtremeConstrainedFlux().
	%
	% To account for memory leaks in some optimizers (GLPK) the samples are
	% partitioned into slices. Each slice reopens the MATLAB pool for parallel
	% computation.
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
	%	- nSamples: Number of samples.
	%	- normPC: Allowed 2-norm deviation of the flux vector from the nominal 
	%       one.
	%   - singlePC: Allowed deviation of each flux from the nominal one.
	%   - maxIter: Maximum iteration count for the bisection.
	%   - nSlice: Number of samples in a slice.
	%   - closePool: True if the matlabpool shall be closed after computations
	%       are finished. False to keep it as is.
	%   - nWorkers: Optional. Amount of workers used for the computation.
	%
	% Returns a structure with the following fields:
	%	- origVal: Original objective function value of the unperturbed FBA problem.
	%	- origSol: Original solution to unperturbed FBA problem.
	%	- samples: Matrix in which each column contains the changed
	%		coefficients of the stoichiometric matrix after perturbance.
	%		The values are set in the same order as the groups. If a group 
	%		perturbs more than one coefficient, the values are saved in 
	%		column-major form.
	%   - maxVal: Maximum coefficient value (see res.samples for ordering).
	%   - minVal: Minimum coefficient value (see res.samples for ordering).
	%   - validPC: Percentage of valid samples.
	%   - table: Cell array with summarized, readable data. Contains absolute 
	%       values.

    if nargin <= 5
        maxIter = 200;
    end
    
    if nargin <= 7
        closePool = true;
    end
    if nargin <= 6
        nSlice = nSamples / 100;
    else
        nSlice = min(nSlice, nSamples);
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
	[groups, sampleSize, res.origData] = checkGroupsAndSaveOriginalData(model, groups, 'uniform');
    
    tHandle = tic;
    res.maxVal = -inf(sampleSize, 1);
    res.minVal = inf(sampleSize, 1);
    res.samples = zeros(sampleSize, nSamples);
    nValid = 0;

    nomVal = full(res.origData);

    % Slice
    for j = 1:floor(nSamples / nSlice)

%         if matlabpool('size') == 0
%             if (nargin <= 8) || isempty(nWorkers)
%                 matlabpool open;
%             else
%                 matlabpool('open', nWorkers);
%             end
%         end
        
        sampSlice = zeros(sampleSize, nSlice);
        
        % Sample main loop
%         parfor i = 1:nSlice
       for i = 1:nSlice

            % Sample uniform on sphere
            sphereDir = randn(length(groups), 1);
            sphereDir = sphereDir ./ norm(sphereDir);
            
			[sampSlice(:,i)] = getRobustExtremeConstrainedFlux(sphereDir, nomVal, origSol, model, groups, normPC, diffSol, maxIter);
        end
%         matlabpool close;
        
        for i = 1:nSlice
            if any(~isfinite(sampSlice(:,i)))
                continue;
            end
            
            res.maxVal = max(sampSlice(:,i), res.maxVal);
            res.minVal = min(sampSlice(:,i), res.minVal);
            nValid = nValid + 1;
        end
        
        res.samples(:, (j-1)*nSlice+1:j*nSlice) = sampSlice;
        
        fprintf('Status: %f  - Valid up to now: %g\n', j * nSlice / nSamples * 100, nValid / nSamples * 100);
    end
    
    table = cell(sampleSize, 1);
    index = 1;
    for i = 1:length(groups)
		if lower(groups{i}.type) == 's'
        	mets = model.mets(groups{i}.row);
		else
			mets = cell(length(groups{i}.row), 1);
			for j = 1:length(groups{i}.row)
    			mets{j} = ['Group ' groups{i}.type ' ' num2str(groups{i}.row(j))];
            end
		end
		if size(mets,2) > 1
			mets = mets';
		end
        table(index:index+length(mets)-1) = mets;
        index = index + length(mets);
    end
    
    tElapsed = toc(tHandle);
    fprintf('Elapsed time: %g sec = %g h %g min\n', tElapsed, floor(tElapsed / 3600), floor((tElapsed - floor(tElapsed / 3600)*3600) / 60));
    
    if closePool && (matlabpool('size') > 0)
        matlabpool close;
    end
    
    % Assign res structure
    res.validPC = nValid / nSamples * 100;
    res.origSol = origSol;
    res.origVal = origVal;
    res.table = [{'Group'}, {'Nominal'}, {'Min'}, {'Max'}, {'Width'}, {'RelWidthPC'}; ...
        table, num2cell(res.origData), num2cell([res.minVal, res.maxVal, abs(res.maxVal - res.minVal), ... 
        abs(res.maxVal - res.minVal) ./ abs(res.origData) * 100])];
end
