function bmInd = findBiomassFlux( model, bmInd )
%FINDBIOMASSFLUX Looks for the biomass flux in the model and returns its
%index.
%
% The function looks for reactions with names like "biomass" or "cdw".
%
% Parameters:
%   - model: Model in which the biomass flux is searched
%   - bmInd: Optional. Name of the biomass flux or index.
%
% Returns: Index of the biomass flux

    if nargin < 2
        % Look for biomass flux in reaction names
        bmInd = lookForFluxName(model);

        % Biomass equation in column bmInd
        if isempty(bmInd)
            error('No biomass flux found. Please rename flux (name should be starting with "biomass" or "cdw").');
        end
        fprintf('Auto determined biomass flux: %s (Index: %g)\n', model.rxns{bmInd}, bmInd);
    else
        if ischar(bmInd)
            bmInd = lookForFluxName(model, bmInd);

            % Biomass equation in column bmInd
            if isempty(bmInd)
                error('No biomass flux with the given name found.');
            end
        else
            fprintf('Using %s (Index: %g) as biomass flux.\n', model.rxns{bmInd}, bmInd);
        end
    end

end

function bmInd = lookForFluxName(model, bmInd)
    if (nargin < 2) || isempty(bmInd) || ~ischar(bmInd)
        bmInd = [{'biomass'}, {'cdw'}];        
    else
        bmInd = [{bmInd}, {'biomass'}, {'cdw'}];
    end

    for i = 1:length(bmInd)
        bmName = bmInd{i};
        
        if isfield(model, 'rxns')
            bmIndicator = strncmpi(bmName, model.rxns, length(bmName));
        else
            bmIndicator = strncmpi(bmName, model.reactions, length(bmName));
        end

        if any(bmIndicator)
            bmInd = find(bmIndicator, 1, 'first');
            return;
        end
    end
    bmInd = [];
end