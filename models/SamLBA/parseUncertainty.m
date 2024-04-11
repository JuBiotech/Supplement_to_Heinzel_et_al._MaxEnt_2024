function [dimSample, convMatrix, constMat, constB, mu, sigma] = parseUncertainty(U, C, S, reactions, pools, dimMap)
%PARSEUNCERTAINTY parses the given uncertainty set and calculates the
%some conversion matrices.
%
% Parameters:
%   - U: Uncertainties structure consisting of cell arrays.
%       |- Target function (Index) = U.target
%       |- Stoichiometric matrix (Reaction, Pool) = U.stoichiometry
%       |- Constraints (Index, Index) = U.constreq, U.constrrhs
%       |- Lower and upper bounds (Index each) = U.upper, U.lower
%       Each element of a cell array may have a bounds vector [lower, upper]
%   - C: Constraint system (C.C * x <= C.b)
%   - S: Stoichiometric matrix
%   - reactions: Cell array with reaction names
%   - dimMap: Data dimension = 3 * length(objective vector) + numel(S) +
%       numel(C.C) + length(C.b)
%
% Returns:
%   - dimSample: Sampling dimension
%   - convMatrix: Conversion matrix. Sample goes in and data vector comes
%       out. Data vector is composed of: 
%       * Target-fuction coefficients
%       * Stoichiometric coefficients (row-wise)
%       * Constraint-matrix coefficients (row-wise)
%       * Constraint-right-hand-side
%   - constMat: Constraint matrix for sampling.
%   - constB: Right-hand-side of constraints for sampling.
%   - mu: Vector of means for sampling.
%   - sigma: Covariance matrix for sampling.

    dimSample = getSampleDimension(U);
    
    convMatrix = zeros(dimMap, dimSample);
    constMat = zeros(dimSample * 2, dimSample);
    constB = zeros(dimSample * 2, 1);
    mu = zeros(dimSample, 1);
    sigma = zeros(dimSample, dimSample);

    curIndex = 1;
    base = 0;
    
    if isfield(U, 'target')
        % Add target function constraints
        for i = 1:length(U.target)
            tc = U.target{i};
            convMatrix(base + tc.index, curIndex) = 1;
            
            [constMat, constB, mu, sigma] = ...
                parseUncertaintyInfo(curIndex, constMat, constB, mu, sigma, tc);
            
            curIndex = curIndex + 1;
        end
        base = base + length(reactions);
    end
    
    if isfield(U, 'stoichiometry')
        % Add stoichiometric constraints
        for i = 1:length(U.stoichiometry)
            sto = U.stoichiometry{i};
            r = getIndexOfElement(reactions, sto.reaction);
            p = getIndexOfElement(pools, sto.pool) - 1;

            if strcmp(sto.reaction, '#row') && strcmp(sto.pool, '#column')
                % Perturb complete matrix
                inds = reshape(S' ~= 0, numel(S), 1)' .* [1:length(reactions) * length(pools)];
                inds(inds == 0) = [];
                convMatrix(base + inds, curIndex) = 1;
            else
                if strcmp(sto.reaction, '#row')
                    % Perturb one entire row
                    inds = [1:lenght(reactions)] .* (S(p+1, :) ~= 0);
                    inds(inds == 0) = [];
                    convMatrix(base + length(reactions) * p + inds, curIndex) = 1;
                elseif strcmp(sto.pool, '#column')
                    % Perturb one entire column
                    inds = [0:length(pools)-1] .* (S(:, r)' ~= 0);
                    inds(inds == 0) = [];
                    convMatrix(base + length(reactions) * inds + r, curIndex) = 1;
                else
                    % Perturb jut one coefficient
                    convMatrix(base + length(reactions) * p + r, curIndex) = 1;
                end
            end

            [constMat, constB, mu, sigma] = ...
                parseUncertaintyInfo(curIndex, constMat, constB, mu, sigma, sto);
            
            curIndex = curIndex + 1;
        end
        base = base + length(reactions) * length(pools);
    end
    
    if isfield(U, 'constreq')
        % Add constraints constraints
        for i = 1:length(U.constreq)
            ce = U.constreq{i};
            r = getIndexOfElement(reactions, ce.reaction);
            
            convMatrix(base + length(reactions) * (ce.index-1) + r, curIndex) = 1;

            [constMat, constB, mu, sigma] = ...
                parseUncertaintyInfo(curIndex, constMat, constB, mu, sigma, ce);
            
            curIndex = curIndex + 1;
        end
        base = base + numel(C);
    end

    if isfield(U, 'constrrhs')
        % Add constraint right-hand-side constraints
        for i = 1:length(U.constrrhs)
            cr = U.constrrhs{i};
                        
            convMatrix(base + cr.index, curIndex) = 1;
            [constMat, constB, mu, sigma] = ...
                parseUncertaintyInfo(curIndex, constMat, constB, mu, sigma, cr);

            curIndex = curIndex + 1;
        end
        base = base + size(C, 1);
    end

    if isfield(U, 'lower')
        % Add variable lower bound constraints
        for i = 1:length(U.lower)
            lb = U.lower{i};
                        
            convMatrix(base + lb.index, curIndex) = 1;
            [constMat, constB, mu, sigma] = ...
                parseUncertaintyInfo(curIndex, constMat, constB, mu, sigma, lb);

            curIndex = curIndex + 1;
        end
        base = base + length(reactions);
    end
        
    if isfield(U, 'upper')
        % Add variable upper bound constraints
        for i = 1:length(U.upper)
            ub = U.upper{i};
                        
            convMatrix(base + ub.index, curIndex) = 1;
            [constMat, constB, mu, sigma] = ...
                parseUncertaintyInfo(curIndex, constMat, constB, mu, sigma, ub);

            curIndex = curIndex + 1;
        end
    end
end

function [constMat, constB, mu, sigma] = ...
    parseUncertaintyInfo(curIndex, constMat, constB, mu, sigma, obj)
%PARSEUNCERTAINTYINFO parses the uncertainty-info given by the object.
%Mean and standard-deviation are saved into vectors and the
%constraint-matrix and -right-hand-side for sampling are built.

    constMat(2 * curIndex - 1 : 2 * curIndex, curIndex) = [1; -1];
    constB(2 * curIndex - 1 : 2 * curIndex) = [obj.bounds(2); -obj.bounds(1)];
    if isfield(obj, 'mean') && ~isempty(obj.mean)
        mu(curIndex) = obj.mean;
    end
    if isfield(obj, 'std') && ~isempty(obj.std)
        sigma(curIndex, curIndex) = obj.std;
    end
end

function index = getIndexOfElement(vals, e)
%GETINDEXOFELEMENT determines the index of the value e in the string cell
%array vals.
    inds = 1:length(vals);
    index = inds(strcmp(vals, e));
end

function dimSample = getSampleDimension(U)
%GETSAMPLEDIMENSION calculates the sampling dimension.
    dimSample = 0;

    if isfield(U, 'target')
        % Add target function constraints
        dimSample = dimSample + length(U.target);
    end
    
    if isfield(U, 'stoichiometry')
        % Add stoichiometric constraints
        dimSample = dimSample + length(U.stoichiometry);
    end
    
    if isfield(U, 'constreq')
        % Add constraints constraints
        dimSample = dimSample + length(U.constreq);
    end

    if isfield(U, 'constrrhs')
        % Add constraint right-hand-side constraints
        dimSample = dimSample + length(U.constrrhs);
    end
    
    if isfield(U, 'upper')
        % Add constraint right-hand-side constraints
        dimSample = dimSample + length(U.upper);
    end
    
    if isfield(U, 'lower')
        % Add constraint right-hand-side constraints
        dimSample = dimSample + length(U.lower);
    end
end