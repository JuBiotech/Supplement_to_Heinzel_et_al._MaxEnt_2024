function [origData] = getModelCoefficients(model, groups, relData)
%GETMODELCOEFFICIENTS Returns the model coefficients subject to the given groups.
% May also be used to apply a scaling factor to each group.
%
% Parameters:
%	- model: Model.
%	- groups: Groups. Have to be correct (check by checkGroupsAndSaveOriginalData() ).
%	- relData: Optional. Vector with a coefficient for each group which is applied to the returned data.
%
% Returns: Model coefficients given by the groups

    origData = [];

	for k = 1:length(groups)
		g = groups{k};
        
        if ~strcmpi(g.type, 'cM') && ~strcmpi(g.type, 'S')
            groupSize = length(g.index);
        else
            groupSize = length(g.row) * length(g.column);
        end
        
		if (nargin <= 2) || isempty(relData)
			origData = [origData; reshape(getCoefficientsOfGroup(model, g), groupSize, 1)];
		else
			origData = [origData; reshape(getCoefficientsOfGroup(model, g), groupSize, 1) .* relData(k)];
		end
    end
    
    origData = full(origData);
end
