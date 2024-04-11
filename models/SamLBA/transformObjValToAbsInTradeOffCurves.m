function transformObjValToAbsInTradeOffCurves(solution, dpc, folder, baseNameStr, figureFormats)
    
    if (nargin < 5) || isempty(figureFormats)
        figureFormats = [{'.fig'}, {'.tikz'}];
    else
        % Make sure each format has a leading dot
        for i = 1:length(figureFormats)
            if figureFormats{i}(1) ~= '.'
                figureFormats{i} = ['.' figureFormats{i}];
            end
        end
    end
    
    if (nargin < 4) || isempty(baseNameStr)
        baseNameStr = 'ro_';
    else
        % Take care of trailing underscore
        if baseNameStr(end) ~= '_'
            baseNameStr = [baseNameStr, '_'];
        end
    end
    
    if (nargin < 3) || isempty(folder)
        folder = '';
    else
        % Make sure the folder has a trailing / (filesep)
        if (~isempty(folder)) && (folder(end) ~= filesep)
            folder = [folder filesep];
        end
                
        % Create directory if it does not exist
        if ~isdir(folder)
            mkdir(folder);
        end
    end
    
    prefix = [folder, baseNameStr];
    
    rhos = unique(solution(1, :));
    eigenLbs = unique(solution(2, :));
    epsVals = sort(unique(solution(3, :)));
    
    for l = 1:length(rhos)

        for j = 1:length(eigenLbs)

            plotData = zeros(2, length(epsVals));
            for i = 1:length(epsVals)
                ind = (l-1) * length(eigenLbs) * length(epsVals) + (j-1) * length(epsVals) + i;
                plotData(:, i) = solution(4:5, ind);
            end

            semilogy(epsVals, abs(plotData(1,:)));
            hold on;
            grid on;
            semilogy(epsVals, plotData(2,:), 'r');
            hold off;

            title([num2str(dpc) '% Uncertainty - ' num2str(rhos(l)) ' Omega - ' num2str(eigenLbs(j)) ' Eigen LB']);
            legend('ObjFun', 'Lambda');
            saveGraph(gcf, [prefix, num2str(dpc) 'dpc_' num2str(rhos(l)) 'omega_' num2str(eigenLbs(j)) '_elb'], figureFormats);
        end
    end
end

function saveGraph(handle, fileName, formats)
    for i = 1:length(formats)
        if strcmpi('.tikz', formats{i})
            matlab2tikz([fileName, formats{i}], 'silent', true, 'height', '\figHeight', 'width', '\figWidth');
        else
            saveas(handle, [fileName, formats{i}]);
        end
    end
end
