function data = getCoefficientsOfGroup(model, g)
%GETCOEFFICIENTSOFGROUP Returns the model coefficients of the given group.
%
% Parameters:
%	- model: Model.
%	- g: Group. See MCGroups.m or rangeGroups.m for hints.
%
% Returns: Raw coefficients

    data = [];
    switch lower(g.type)
        case 's'
            data = model.S(g.row, g.column);
        case 'b'
            data = model.b(g.index);
        case 'cm'
            data = model.cM(g.row, g.column);
        case 'cb'
            data = model.cB(g.index);
        case 'c'
            data = model.c(g.index);
        case 'lb'
            data = model.lb(g.index);
        case 'ub'
            data = model.ub(g.index);
    end
end
