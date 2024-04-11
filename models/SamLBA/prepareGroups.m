function groups = prepareGroups( mode, model, sigmaPC, fluxes, independent, topFlux, topFluxExceptions )
%PREPAREGROUPS Prepares the coefficient groups for the given network.
% The groups created by this function preserve the correlation of the
% coefficients of a reaction (e.g. biomass reaction).
%
% Correlated coefficients will be put together in one group.
% The function assumes that coefficients are correlated if they have the same absolute value.
% For description of the group structure see mcGroups.m .
%
% This function can prepare the groups for MC and range methods (robust range, range scans).
%
% Parameters:
%   - mode: One of 'range', 'mc'. Creates groups for the given purpose.
%	- model: Model.
%	- sigmaPC: In mc-mode: Standard deviation of the biomass coefficients 
%       in percent(!). In range-mode: Vector [sigma, steps] with range and 
%       steps. The range is given by coefficient * [(1-sigma), (1+sigma)]
%	- fluxes: Vector with indices of perturbed fluxes.
%   - independent: Set to true if every coefficient shall form its own
%       group. Defaults to false.
%   - topFlux: Optional. If set ignores all metabolites in the top flux,
%       which also occur in other fluxes of the fluxes vector.
%   - topFluxExceptions: Optional. Vector of metbolite indices which will
%       be excluded from removal.
%
% Returns a cell array of group structures.

	groups = {};
    
    if nargin < 5
        independent = false;
    end
    
    if nargin < 6
        topFlux = [];
    else
        % Remove top flux
        fluxes(fluxes == topFlux) = [];
    end
    
    if nargin < 7
        topFluxExceptions = [];
    end
    
    % Base setting
    g.relative = true;
    if strcmpi('mc', mode)
        sigmaPC = sigmaPC / 100;
        g.mean = 1;
        g.sigma = sigmaPC;
    elseif strcmpi('range', mode)
        sigma = sigmaPC(1) / 100;
        g.range = [1-sigma, 1+sigma, sigmaPC(2)];
    end
    
    if independent
        % Create a single group for each coefficient neglecting
        % dependencies
        for i = 1:length(fluxes)
            col = fluxes(i);
            coeffInds = find(model.S(:, col) ~= 0);

            g.column = col;
            g.type = 'S';

            for j = 1:length(coeffInds)
                g.row = coeffInds(j);
                groups{end+1} = g;
            end
        end

	   if isfield(model, 'cM') && ~isempty(model.cM)
	       for i = 1:length(fluxes)
	            col = fluxes(i);
	            coeffInds = find(model.cM(:, col) ~= 0);

	            g.column = col;
	            g.type = 'cm';

	            for j = 1:length(coeffInds)
	                g.row = coeffInds(j);
	                groups{end+1} = g;
	            end
	        end
		end

    else
        % Create groups
        
        for i = 1:length(fluxes)
            col = fluxes(i);

            [corrs, uncorr] = determineGroupsInColumn(model.S(:, col));

            g.column = col;
            g.type = 'S';

            for j = 1:length(corrs)
                g.row = corrs{j};
                groups{end+1} = g;
            end

            % Save uncorrelated
            for j = 1:length(uncorr)
                g.row = uncorr(j);
                groups{end+1} = g;
            end
        end
    end
    
    % Handle top flux
    if ~isempty(topFlux)
        [corrs, uncorr] = determineGroupsInColumn(model.S(:, topFlux));
        
        g.column = topFlux;

        for j = 1:length(corrs)
            
            mets = corrs{j};
            
            % Check occurrence of metabolites in another flux
            delInd = [];
            for k = 1:length(mets)
                metInd = mets(k);
                if all(topFluxExceptions ~= metInd) && any(model.S(metInd, fluxes) ~= 0)
                    delInd = [delInd ; k];
                end
            end
            mets(delInd) = [];
            
            if isempty(mets)
                continue;
            end
            
            g.row = mets;
            groups{end+1} = g;
        end

        % Save uncorrelated
        for j = 1:length(uncorr)
            
            % Check occurrence of metabolite in another flux
            if all(topFluxExceptions ~= uncorr(j)) && any(model.S(uncorr(j), fluxes) ~= 0)
                continue;
            end
            
            g.row = uncorr(j);
            groups{end+1} = g;
        end
    end
end

function [corr, unCorr] = determineGroupsInColumn(col)
    
    corr = {};
    unCorr = [];
    vals = unique(abs(col));

    % Remove 0
    vals(vals == 0) = [];
    
    for i = 1:length(vals)
        inds = find(abs(col) == vals(i));
        if length(inds) > 1
            corr{end+1} = inds;
        else
            unCorr = [unCorr ; inds];
        end
    end
end
