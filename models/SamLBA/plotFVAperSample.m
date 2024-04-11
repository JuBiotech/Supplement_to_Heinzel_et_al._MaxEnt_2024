function plotFVAperSample( model, minFluxes, maxFluxes, folder, prefix, formats )
%PLOTFVAPERSAMPLE Summary of this function goes here
%   Detailed explanation goes here

    % Check constraints
    if ~isfield(model, 'cM')
        model.cM = [];
    end
    if ~isfield(model, 'cB')
        model.cB = [];
    end
    
    if size(model.c, 1) == 1
        model.c = model.c';
    end

    prefName = '';
    if (nargin >= 4) && (~isempty(folder)) 
        if folder(end) ~= filesep
            folder = [folder filesep];
        end

        % Create directory if it does not exist
        if ~isdir(folder)
            mkdir(folder);
        end
        
        prefName = folder;
    end
    
    if (nargin >= 5) && (~isempty(prefix))
        if prefix(end) ~= '_'
            prefix = [prefix '_'];
        end
        
        prefName = [prefName prefix];
    end
 
    if (nargin >= 6) && (~isempty(formats))
        % Make sure each format has a leading dot
        for i = 1:length(formats)
            if formats{i}(1) ~= '.'
                formats{i} = ['.' formats{i}];
            end
        end
    else
        formats = [{'.fig'}, {'.png'}];
    end
    
    nBuckets = 30;
    
    [origSol, origObjVal] = solveLPProblem(-1, model.c, model.cM, model.cB, model.S, model.b, model.lb, model.ub);
    [minOrig, maxOrig] = fluxVariabilityAnalysis(model);

    inds = find(abs(model.ub - model.lb) > 0);
    set(0, 'DefaultFigureVisible', 'off')

%     for i=1:length(inds)
%         
%         k = inds(i);
%         
%         plotFVAHist(nBuckets, maxFluxes(k,:), origSol(k), maxOrig(k), 'Maximum', model.rxnNames{k}, k);
%         saveGraph(gcf, [prefName 'max_' num2str(k) '_' model.rxnNames{k}], formats);
%         
%         plotFVAHist(nBuckets, minFluxes(k,:), origSol(k), minOrig(k), 'Minimum', model.rxnNames{k}, k);
%         saveGraph(gcf, [prefName 'min_' num2str(k) '_' model.rxnNames{k}], formats);
% 
%         hist(maxFluxes(k,:) - minFluxes(k,:), nBuckets);
%         hold on;
%         grid on;
% 
%         limits = get(gca, 'YLim');
%         line([maxOrig(k)-minOrig(k), maxOrig(k)-minOrig(k)],[0, limits(2)], 'Color', 'g', 'LineWidth', 1.5);
% 
%         hold off;
%         title(['FVA interval width - ' model.rxnNames{k} ' (' num2str(k) ')']);
%         legend('hist', 'origFVA');
%         saveGraph(gcf, [prefName 'intwidth_' num2str(k) '_' model.rxnNames{k}], formats);
%     end
    for i=1:length(inds)
        
        k = inds(i);
        figure;
        hold on;
        x=min(maxFluxes(k,:))-0.01:0.005:max(maxFluxes(k,:))+0.01;
        stairs(x,histc(maxFluxes(k,:),x),'-','Color','black','LineWidth',3);%,'Color',[0.8 0.8 0.8],
        
        y=min(minFluxes(k,:))-0.01:0.005:max(minFluxes(k,:))+0.01;
        stairs(y,histc(minFluxes(k,:),y),'--','Color','green','LineWidth',3);%,'Color',[0.8 0.8 0.8],
        
        limits = get(gca, 'YLim');
        line([maxOrig(k), maxOrig(k)],[0, limits(2)], 'Color', 'c', 'LineWidth', 3);
        line([minOrig(k), minOrig(k)],[0, limits(2)], 'Color', 'b', 'LineWidth', 3);
        line([origSol(k), origSol(k)],[0, limits(2)], 'Color', 'r', 'LineWidth', 1.5, 'LineStyle','--');        
        hold off;
        title(['FVA MinMax - ' model.rxnNames{k} ' (' num2str(k) ')']);
        saveGraph(gcf, [prefName 'MinMax_' num2str(k) '_' model.rxnNames{k}], formats);

    end
    set(0, 'DefaultFigureVisible', 'on')
end

function plotFVAHist(nBuckets, data, orig, fvaOrig, type, rxnName, k)
    hist(data, nBuckets);
    hold on;
    grid on;

    limits = get(gca, 'YLim');
    line([orig, orig],[0, limits(2)], 'Color', 'r', 'LineWidth', 1.5);
    line([fvaOrig, fvaOrig],[0, limits(2)], 'Color', 'g', 'LineWidth', 1.5);

    hold off;
    title(['FVA ' type ' - ' rxnName ' (' num2str(k) ')']);
    legend('hist', 'origSol', 'origFVA');
end

function saveGraph(handle, fileName, formats)
    for i = 1:length(formats)
        if strcmpi('.tikz', formats{i})
            matlab2tikz([fileName, formats{i}], 'silent', true, 'height', '\mcHistHeight', 'width', '\mcHistWidth');
        else
            saveas(handle, [fileName, formats{i}]);
        end
    end
end
