function groups = prepareGroupsCoryneAll( model, sigmaPC, fluxes, topFlux, topFluxExceptions )
%PREPAREGROUPSCORYNE Prepares the groups for the Coryne Glutamicum network.
% The groups created by this function preserve the correlation of the
% coefficients of the biomass reaction.
%
% Correlated coefficients will be put together in one group.
% For description of the group structure see mcGroups.m .
%
% Parameters:
%	- model: Model.
%	- sigmaPC: Standard deviation of the biomass coefficients in percent(!)
%	- fluxes: Vector with indices of perturbed fluxes.
%   - topFlux: Optional. If set ignores all metabolites in the top flux,
%       which also occur in other fluxes of the fluxes vector.
%   - topFluxExceptions: Optional. Vector of metbolite indices which will
%       be excluded from removal.
%
% Returns a cell array of group structures.

	sigmaPC = sigmaPC / 100;
	groups = {};
    
    if nargin < 4
        topFlux = [];
    else
        % Remove top flux
        fluxes(fluxes == topFlux) = [];
    end
    
    if nargin < 5
        topFluxExceptions = [];
    end
    
    % Base setting
    g.relative = true;
    g.mean = 1;
    g.sigma = sigmaPC;

    for i = 1:length(fluxes)
        col = fluxes(i);
        
        [corrs, uncorr] = determineGroupsInColumn(model.S(:, col));
        
        g.column = col;

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
