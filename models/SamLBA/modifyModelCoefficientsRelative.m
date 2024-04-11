function [model] = modifyModelCoefficientsRelative(model, groups, data)
%MODIFYMODELCOEFFICIENTSRELATIVE Modifies the coefficients of the given model subject to the given coefficient groups.
% Multiplies each coefficient group by a factor given in the data vector and returns the modified model.
%
% Parameters:
%	- model: Model.
%	- groups: Cell array of groups. See MCGroups.m or rangeGroups.m for hints
%	- data: Vector of multiplicators. Each component corresponds to a group.
%
% Returns the modified model.

	for k = 1:length(groups)
		g = groups{k};

		switch lower(g.type)
		    case 's'
				model.S(g.row, g.column) = model.S(g.row, g.column) .* data(k);
		    case 'b'
				model.b(g.index) = model.b(g.index) .* data(k);
		    case 'cm'
				model.cM(g.row, g.column) = model.cM(g.row, g.column) .* data(k);
		    case 'cb'
				model.cB(g.index) = model.cB(g.index) .* data(k);
		    case 'c'
				model.c(g.index) = model.c(g.index) .* data(k);
		    case 'lb'
				model.lb(g.index) = model.lb(g.index) .* data(k);
		    case 'ub'
				model.ub(g.index) = model.ub(g.index) .* data(k);
		end
	end
end
