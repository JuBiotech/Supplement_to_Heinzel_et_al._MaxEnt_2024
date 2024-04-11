function plotAndSavePolys( polys10pc)

    for k =1:length(polys10pc)

        curPoly = polys10pc{k};
        polyData = cell2mat(curPoly(2:end, :));     
        
        polyData(:,2:3) = 100.*polyData(:, 2:3);
        
        clf;
        plot(polyData(:, 2), polyData(:, 3), 'o-');

        grid on;
        hold on;

        fill(polyData(:, 2), polyData(:, 3), polyData(:, 4));

        xlabel(curPoly{1,2});
        ylabel(curPoly{1,3});

        angles = polyData(:, 1);
        indA = find(angles == 0,1, 'first');
        if isempty(indA)
            indA = find(angles == 360,1, 'first');
        end

        indB = find(angles == 180,1, 'first');
        plot(polyData([indA, indB], 2),polyData([indA, indB], 3), 'r');
        indA = find(angles == 90,1, 'first')+1;
        indB = find(angles == 270,1, 'first')+1;
        plot(polyData([indA, indB], 2),polyData([indA, indB], 3), 'r');

        colorbar;

        hold off;

        % Hier Save Befehl
        saveas(gcf, ['polygone' num2str(k) '.png']);
    end

end