function [model, data] = modifyGroupCoefficients(model, g, data)
%MODIFYGROUPCOEFFICIENTS Modifies the model coefficients of a single group in relative or absolute fashion.
% Multiplies the coefficients of the given group by the given factor or assigns the given absolute data to the
% coefficients. The modified model is returned.
%
% Parameters:
%	- model: Model.
%	- g: Group to modify.
%	- data: Data (relative multiplicator or absolute values)
%
% Returns the modified model.

    switch lower(g.type)
        case 's'
			if g.relative                
				data = model.S(g.row, g.column) .* data;
				model.S(g.row, g.column) = data;
			else
				model.S(g.row, g.column) = ones(length(g.row), length(g.column)) .* data;
			end
        case 'b'
			if g.relative                
				data = model.b(g.index) .* data;
				model.b(g.index) = data;
			else
				model.b(g.index) = ones(length(g.index), 1) .* data;
			end
        case 'cm'
			if g.relative                
				data = model.cM(g.row, g.column) .* data;
				model.cM(g.row, g.column) = data;
			else
				model.cM(g.row, g.column) = ones(length(g.row), length(g.column)) .* data;
			end
        case 'cb'
			if g.relative                
				data = model.cB(g.index) .* data;
				model.cB(g.index) = data;
			else
				model.cB(g.index) = ones(length(g.index), 1) .* data;
			end
        case 'c'
			if g.relative                
				data = model.c(g.index) .* data;
				model.c(g.index) = data;
			else
				model.c(g.index) = ones(length(g.index), 1) .* data;
			end
        case 'lb'
			if g.relative                
				data = model.lb(g.index) .* data;
				model.lb(g.index) = data;
			else
				model.lb(g.index) = ones(length(g.index), 1) .* data;
			end
        case 'ub'
			if g.relative                
				data = model.ub(g.index) .* data;
				model.ub(g.index) = data;
			else
				model.ub(g.index) = ones(length(g.index), 1) .* data;
			end
    end
end
