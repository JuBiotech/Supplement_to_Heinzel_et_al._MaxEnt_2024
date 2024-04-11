function M = visualizeNDArray( groups, data, makeMovie )
%VISUALIZENDARRAY Visualizes data obtained by rangeGroups().
%
% Parameters:
%   - groups: Cell array of groups used by rangeGroups. See rangeGroups.m
%       for information.
%   - data: Data to visualize. An nD-array where n is the number of groups.
%       For example use res.diffs or res.vals returned by rangeGroups().
%   - makeMovie: Optional. Makes a movie of 3D Plots in case of 3 groups.
%
% Returns: Movie (in case of 3D with makeMovie or 4D). Otherwise nothing.

    M = [];

    ranges = cell(length(groups), 1);
    dimSizes = zeros(1, length(groups));
    linearizedLength = 1;
    for i = 1:length(groups)
        g = groups{i};
        
        if ~isfield(g, 'rawRange')
            ranges{i} = linspace(g.range(1), g.range(2), g.range(3));
        else
            ranges{i} = g.rawRange;
        end
        dimSizes(i) = length(ranges{i});
        linearizedLength = linearizedLength * length(ranges{i});
    end

    switch length(dimSizes)
        case 1
            plot(ranges{1}, data);
        case 2
            [X,Y] = meshgrid(ranges{1}, ranges{2});
            surf(X,Y, data);
            shading interp;
        case 3
            if (nargin == 2) || ~makeMovie
                [X,Y,Z] = meshgrid(ranges{1}, ranges{2}, ranges{3});
                slice(X,Y,Z, data, mean(ranges{1}), mean(ranges{1}), mean(ranges{1}));
            else
                for i = 1:dimSizes(3)
                    [X,Y] = meshgrid(ranges{1}, ranges{2});
                    surf(X,Y, data(:, :, i));
                    M(i) = getframe;
                end
                movie(M);
            end
        case 4
            [X,Y,Z] = meshgrid(ranges{1}, ranges{2}, ranges{3});
            
            for i = 1:dimSizes(4)
                slice(X,Y,Z, data(:, :, :, ranges{4}(i)), mean(ranges{1}), mean(ranges{1}), mean(ranges{1}));
                M(i) = getframe;
            end
            movie(M);
    end

end
