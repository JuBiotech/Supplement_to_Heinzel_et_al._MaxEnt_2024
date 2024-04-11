function [groups, sampleSize, origData] = checkGroupsAndSaveOriginalData(model, groups, mode)
%CHECKGROUPSANDSAVEORIGINALDATA Checks the given groups against the model (for consistency) and saves the original coefficients.
% The given groups are checked for consistency. Assigns the correct random number generators according
% to the given mode. Collects the original coefficients and returns them in a vector.
%
% Parameters:
%	- model: Model.
%	- groups: Cell array with coefficient groups.
%	- mode: String. Choose from 'normal', 'truncnormal' and 'uniform'.
%
% Returns:
%	- groups: Cell array with checked groups.
%	- sampleSize: Size of the uncertain coefficient vector.
%	- origData: Original data of the given model compressed to a single vector (of length sampleSize).

    % Calculate sample size and check if groups are sane
	sampleSize = 0;
    origData = [];

	for k = 1:length(groups)
		g = groups{k};
        
        % Default type is 'S'
        if ~isfield(g, 'type')
            g.type = 'S';
        end
        
        if ~strcmpi(g.type, 'cM') && ~strcmpi(g.type, 'S')
            % Perturbation of vector (c, b, cB, lb, ub)
            if ~isfield(g, 'row') && ~isfield(g, 'column') && ~isfield(g, 'index')
                % Error
                error('MCGroups:checkGroupsAndSaveOriginalData', ... 
                    'Error: Group %d has to have a "row", "column" or "index" field!', k);
            end
            if isfield(g, 'row') && ~isfield(g, 'index')
                g.index = g.row;
            end
            if isfield(g, 'column') && ~isfield(g, 'index')
                g.index = g.column;
            end
            
            g = rmfield(g, 'row');
            g = rmfield(g, 'column');
            
            groupSize = length(g.index);
        else
            groupSize = length(g.row) * length(g.column);
        end

        % Assign generator
		if strcmp(mode, 'normal')
            if isfield(g, 'truncated') && g.truncated
                g.generator = @truncNormalGenerator;
            else
                g.generator = @normalGenerator;
            end
        elseif strcmp(mode, 'truncnormal')
            if isfield(g, 'truncatedHalf') && g.truncatedHalf
                g.generator = @truncNormalHalfGenerator;
            else
                g.generator = @truncNormalGenerator;
            end
		elseif strcmp(mode, 'uniform')
			if g.relative
                if isfield(g, 'uniformHalf') && g.uniformHalf
    				g.generator = @(x) relativeUniformGenerator(x, 1);
                else
    				g.generator = @(x) relativeUniformGenerator(x, 2);
                end
			else
				if isfield(g, 'isInteger') && g.isInteger
					if isfield(g, 'nonZero') && g.nonZero
						g.generator = @absoluteUniformNonNullIntegerGenerator;
					else
						g.generator = @absoluteUniformIntegerGenerator;
					end
				else
					g.generator = @absoluteUniformGenerator;
				end
			end
		end
        
        groups{k} = g;
        
		sampleSize = sampleSize + groupSize;
		origData = [origData; reshape(getCoefficientsOfGroup(model, g), groupSize, 1)];
	end
end
