function [ minFluxes, maxFluxes ] = mcFVAWithSamples( model, samples, vals, groups, nWorkers )
%MCFVAWITHSAMPLES Summary of this function goes here
%   Detailed explanation goes here

    % Check constraints
    if ~isfield(model, 'cM')
        model.cM = [];
    end
    if ~isfield(model, 'cB')
        model.cB = [];
    end
    
    if size(model.c, 1) == 1
        model.c = model.c';
    end

    nonZeroRangeInds = abs(model.ub - model.lb) > 0;
    inds = find(nonZeroRangeInds);

    minFluxes = zeros(length(model.c), size(samples, 2));
    maxFluxes = zeros(length(model.c), size(samples, 2));
    
    if matlabpool('size') == 0
        if (nargin <= 4) || isempty(nWorkers)
            matlabpool open;
        else
            matlabpool('open', nWorkers);
        end
    end
    parfor i = 1:size(samples,2)
        modelMod = model;
        curSamp = samples(:,i);
        
        % Obtain sample
        for k = 1:length(groups)
            g = groups{k};
            g.relative = false;
            modelMod = modifyGroupCoefficients(modelMod, g, curSamp(k));
        end
        
        [minFluxes(:, i), maxFluxes(:,i)] = localFVA(modelMod, nonZeroRangeInds, inds, vals(i));
    end
    
    matlabpool close;
end

function [minFlux, maxFlux] = localFVA(model, indsBin, inds, origObjVal)
%FLUXVARIABILITYANALYSIS Performs a flux variability analysis yielding an axis aligned bounding box of the optimal flux configurations.
% FVA constrains the model to the optimal objective function value (or a certain amount of it) and tries
% to maximize / minimize the other flux vector components (one at a time).
% This yields an axis aligned bounding box of the hyperplane with the optimal objective function value.
%
% Note: If the optimum is unique the FVA coefficient intervals are zero.
%
% Parameters:
%	- model: Model.
%   - minNorm: Optional. Set to false to disable 1-Norm minimization (cycle
%       avoiding). Set to true to enable norm minimization. Default is
%       false.
%	- dpcDevFromObjVal: Optional. If > 0 the original objective function value has to be attained up to
%		a deviation given by dpcDevFromObjVal in percent. For example, if dpcDevFromObjVal = 5 the
%		objective function value may be in the range of 95 - 105% of the nominal objective function value.
%
% Returns:
%	- minFlux: Vector with minimum flux components.
%	- maxFlux: Vector with maximum flux components.

    params.minNorm = false;
    
    S = [model.S; model.c'];
    b = [model.b; origObjVal];
    cM = model.cM;
    cB = model.cB;
    
    maxFlux = zeros(length(model.c), 1);
    maxFlux(~indsBin) = model.lb(~indsBin);
    minFlux = maxFlux;
    
    for i = 1:length(inds)
        c = zeros(length(model.c), 1);
        c(inds(i)) = 1;
        [sol, maxFlux(inds(i))] = solveLPProblem(-1, c, cM, cB, S, b, model.lb, model.ub, params);
        [sol, minFlux(inds(i))] = solveLPProblem(1, c, cM, cB, S, b, model.lb, model.ub, params);
    end
end
