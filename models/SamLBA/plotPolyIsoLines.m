function plotPolyIsoLines(runs, polInd)

    if (nargin < 2) || isempty(polInd)
        polInd = 1;
    end

    colOrder = get(gca,'ColorOrder');
    legendOrder = cell(length(runs), 1);
    
    for k =1:length(runs)

        curPoly = runs{k}.polys2D{polInd};
        polyData = cell2mat(curPoly(2:end, :));     
        
        polyData(:,2:3) = 100.*polyData(:, 2:3);
        
        plot(polyData(:, 2), polyData(:, 3), '.', 'Color', colOrder(mod(k, size(colOrder, 1))+1,:));
        
        if isfield(runs{k}, 'bound')
            legendOrder{k} = num2str(runs{k}.bound);
        else
            legendOrder{k} = num2str(runs{k}.singlePC);
        end
        
        grid on;
        hold on;

%        fill(polyData(:, 2), polyData(:, 3), polyData(:, 4));

%        angles = polyData(:, 1);
%        indA = find(angles == 0,1, 'first');
%        if isempty(indA)
%            indA = find(angles == 360,1, 'first');
%        end

%        indB = find(angles == 180,1, 'first');
%        plot(polyData([indA, indB], 2),polyData([indA, indB], 3), 'r');
%        indA = find(angles == 90,1, 'first')+1;
%        indB = find(angles == 270,1, 'first')+1;
%        plot(polyData([indA, indB], 2),polyData([indA, indB], 3), 'r');

        % Hier Save Befehl
        %saveas(gcf, ['polygone' num2str(k) '.png']);
    end

    poly = runs{1}.polys2D{polInd};
    legend(legendOrder);
    xlabel(poly{polInd,2});
    ylabel(poly{polInd,3});
    hold off;

end