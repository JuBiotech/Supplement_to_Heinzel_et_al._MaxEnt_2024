function [metaUn, metaNo, metaTN] = generateDataAndPlots(model, dpcRange, opts)
%GENERATEDATAANDPLOTS Does MC simulations in which the uncertain
%coefficients are randomly disturbed by a given factor (dpcRange) and 
%conducts a FBA for each realisation of the uncertain coefficients.
%
% Results of each run are saved to mat-files. Histograms are generated and
% saved as figure and tikz-file.
%
% Parameters:
%   - model: Coryne model.
%   - dpcRange: Vector with disturbance in percent. F.e. [1 2 3] means, a
%       MC run is performed with disturbance of 1 percent, 2 percent, 3 
%       percent. Yielding a total of 3 MC runs.
%   - opts: Optional. Structure with options. Fields are:
%       o folders: Optional. Cell array with names of directories for data  
%           of uniform and normal MC runs. If not set or empty current 
%           directory will be used to dump the data. 
%               Order: normal, truncnormal, uniform.
%       o groupGenerator: Optional. Function handle to a function which
%           generates groups for a given perturbation level. Interface:
%               groupCellArray = groupGenerator(perturbanceLevel);
%           Defaults to prepareGroupsCoryne().
%       o nSamples: Number of samples for each run.
%       o typeStr: Type string. Used for file name composition.
%       o sampStr: Sample string. Used for file name composition.
%           File names: coryne_optBM_typeStr_DPCpc_sampStr_mode
%       o nBuckets: Number of buckets for histograms.
%       o enabledModes: Vector of booleans. If an entry is set to false,
%           the MC run with the corresponding perturbation mode will not be
%           performed. Order: [normal, truncnormal, uniform]
%
% Returns: Two matrices with meta data about the MC runs.
%   - metaUn: Contains information about the MC runs with uniform
%      distribution.
%   - metaNo: Contains information about the MC runs with Gaussian (normal)
%      distribution.
%   - metaTN: Contains information about the MC runs with truncated normal
%      distribution.
% Each matrix hast the columns dpc, objective value-(min, max, mean, std.),
%   Euclidean distance of flux vector to nominal flux-(min, max, mean,
%   std.)

    if (nargin <= 2)
        % Set defaults
        opts.groupGenerator = @prepareGroupsCglutSoxPhos;%@prepareGroupsCoryneTopBM;
        opts.nSamples = 1000;%10000;
        opts.sampStr = '10kSamp';
%         opts.groupGenerator = @prepareGroupsCglutSoxPhos;%@prepareGroupsCoryneTopBM;
        opts.nBuckets = 30;
        opts.enabledModes = [true, true, true];
        opts.typeStr = 'mcBM';
    end
    
    if ~isfield(opts, 'groupGenerator') || isempty(opts.groupGenerator)
        opts.groupGenerator = @prepareGroupsCoryneTopBM;
    end

    if ~isfield(opts, 'nBuckets') || (opts.nBuckets <= 0)
        opts.nBuckets = 30;
    end

    if ~isfield(opts, 'typeStr') || isempty(opts.typeStr)
        opts.typeStr = 'mcBM';
    end

    if ~isfield(opts, 'nSamples') || (opts.nSamples <= 0)
        opts.nSamples = 10000;
    end
    
    if ~isfield(opts, 'sampStr') || isempty(opts.sampStr)
        opts.sampStr = [num2str(floor(opts.nSamples / 1000)) 'kSamp'];
    end

    if ~isfield(opts, 'enabledModes') || isempty(opts.enabledModes)
        opts.enabledModes = [true, true, true];
    end
    
    if ~isfield(opts, 'folders') || isempty(opts.folders)
        opts.folders = {'', '', ''};
    else
        % Make sure each folder has a trailing / (filesep)
        for i = 1:length(opts.folders)
            if opts.folders{i}(end) ~= filesep
                opts.folders{i} = [opts.folders{i} filesep];
            end
                
            % Create directory if it does not exist
            if ~isdir(opts.folders{i}) && opts.enabledModes(i)
                mkdir(opts.folders{i});
            end
        end
    end
    
    %               -------- vals ---------   -------- diffs ---------
    % Columns: dpc, min, max, mean, std.dev., min, max, mean, std.dev.
    metaUn = zeros(length(dpcRange), 9);
    metaNo = zeros(length(dpcRange), 9);
    metaTN = zeros(length(dpcRange), 9);
    
    modes = {'normal', 'truncnormal', 'uniform'};

    tic;
    
%     % Open pool for parallel computation
%     % Comment if parallel execution is disabled!
%     if matlabpool('size') == 0
%         matlabpool;
%     end
    
    for i = 1:length(dpcRange)
        dpc = dpcRange(i);
        
        for j = 1:length(modes)
            
            if ~opts.enabledModes(j)
                continue;
            end
            
            mode = modes{j};
            groups = opts.groupGenerator(dpc);
%            groups = prepareGroupsCoryne(dpc);
            
            % Uncomment to select a mode (serial, parallel execution).
            % Parallel execution yields a major speedup.
            
            % =============
            % Parallel mode
%             res = mcGroupsParallel(model, groups, opts.nSamples, mode);
            
            % =============
            % Serial mode
            res = mcGroups(model, groups, opts.nSamples, mode);
            
            % speichert matfile mit Infos wie bei MCgroups
%             prefix = [opts.folders{j} 'coryne_optBM_' opts.typeStr '_' num2str(dpc) 'pc_' opts.sampStr '_' mode];
%             save([prefix '.mat'], 'res');

            % ===============
            % Plot generation
            %   Comment lines starting with "matlab2tikz" to disable
            %   tikz-file generation.

%             plotHistPercent(res.vals, opts.nBuckets);
% %             saveas(gcf, [prefix '_objValHist_pc.fig']);
% %             matlab2tikz([prefix '_objValHist_pc.tikz'], 'silent', true, 'height', '\mcHistHeight', 'width', '\mcHistWidth');
% 
%             plotHistPercent(res.diffs, opts.nBuckets);
% %             saveas(gcf, [prefix '_2normHist_pc.fig']);
% %             matlab2tikz([prefix '_2normHist_pc.tikz'], 'silent', true, 'height', '\mcHistHeight', 'width', '\mcHistWidth');
% 
%             plotHistPercent((res.vals - res.origSol(300)) / res.origSol(300)*100, opts.nBuckets);
% %             saveas(gcf, [prefix '_objValDeviation_pc.fig']);
% %             matlab2tikz([prefix '_objValDeviation_pc.tikz'], 'silent', true, 'height', '\mcHistHeight', 'width', '\mcHistWidth');
% 
%             plotHistPercent((res.diffs) / norm(res.origSol)*100, opts.nBuckets);
% %             saveas(gcf, [prefix '_2normDeviation_pc.fig']);
% %             matlab2tikz([prefix '_2normDeviation_pc.tikz'], 'silent', true, 'height', '\mcHistHeight', 'width', '\mcHistWidth');
            
            % ================
            % Meta information
            
            if strcmp(mode, 'uniform')
                metaUn(i, :) = [dpc, ...
                    min(res.vals), max(res.vals), mean(res.vals), std(res.vals, 1), ...
                    min(res.diffs), max(res.diffs), mean(res.diffs), std(res.diffs, 1)];
            elseif strcmp(mode, 'normal')
                metaNo(i, :) = [dpc, ...
                    min(res.vals), max(res.vals), mean(res.vals), std(res.vals, 1), ...
                    min(res.diffs), max(res.diffs), mean(res.diffs), std(res.diffs, 1)];
            else
                metaTN(i, :) = [dpc, ...
                    min(res.vals), max(res.vals), mean(res.vals), std(res.vals, 1), ...
                    min(res.diffs), max(res.diffs), mean(res.diffs), std(res.diffs, 1)];
            end
        end
    end

%     if opts.enabledModes(1)
%         writeMeta(['coryne_optBM_' opts.typeStr '_' opts.sampStr '_Normal.csv'], metaNo);
%     end
%     if opts.enabledModes(2)
%         writeMeta(['coryne_optBM_' opts.typeStr '_' opts.sampStr '_TruncNormal.csv'], metaTN);
%     end
%     if opts.enabledModes(3)
%         writeMeta(['coryne_optBM_' opts.typeStr '_' opts.sampStr '_Uniform.csv'], metaUn);
%     end
%     save(['coryne_optBM_' opts.typeStr '_' opts.sampStr '_Meta.mat'], 'metaUn', 'metaNo', 'metaTN');
    
%     % Comment if parallel execution is disabled!
%     matlabpool close;
        
    fprintf('Needed time: %g sec \n', toc());
end

function writeMeta(fileName, meta)
    fid = fopen(fileName, 'w');
    fprintf(fid, 'dpc;objMin;objMax;objMean;objStd;diffMin;diffMax;diffMean;diffStd\n');
    fclose(fid);
    dlmwrite(fileName, meta, '-append', 'delimiter', ';');
end
