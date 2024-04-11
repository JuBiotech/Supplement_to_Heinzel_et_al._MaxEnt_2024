function metaRes = generateDataAndPlotsSingle(model, dpcRange, opts)
%GENERATEDATAANDPLOTSSINGLE Does MC simulations in which the uncertain
%coefficients are randomly disturbed by a given factor (dpcRange) and 
%conducts a FBA for each realisation of the uncertain coefficients.
%
% Instead of perturbing every uncertain coefficient simultaneously, only
% only one coefficient at a time is perturbed and the others are kept
% nominal.
%
% Results of each run are saved to mat- and csv-files. Histograms are 
% generated and saved as figure and tikz-file if noPlots option is 
% disabled. So by default no plots are generated.
%
% Call: metaRes = generateDataAndPlotsSingle(model, dpcRange, opts)
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
%       o baseNameStr: Prefix of saved files. Could be network name.
%           Defaults to 'net'.
%       o sampStr: Sample string. Used for file name composition. Defaults
%           to '[nSamples / 1000]kSamp'.
%           File names: baseNameStr_mcGx_DPCpc_sampStr_mode
%       o nBuckets: Number of buckets for histograms.
%       o enabledModes: Vector of booleans. If an entry is set to false,
%           the MC run with the corresponding perturbation mode will not be
%           performed. Order: [normal, truncnormal, uniform]
%       o noPlots: If true no plots will be saved. Data will still be saved
%           to mat-Files. Default is true.
%       o dryRun: If true data will not be saved to files. Default is
%           false (save data to files and directories).
%       o parallel: If true the computation will be parallelized. The 
%           Parallel Computation Toolbox is required for this feature to
%           work. Default is true (use parallelization).
%       o nWorkers: Number of workers used by the Parallel Computation
%           Toolbox. Optional, leave out to use default configuration.
%       o figureFormats: Cell array with figure formats used by saveas().
%           A figure will be saved in each given format. Use 'tikz' to save
%           a tikz-file. Defaults to [{'fig'}, {'tikz'}]
%
% Returns: A structure with cell array of matrices with meta data about 
%   the MC runs. The structure has the fields:
%       - uniform: Cell array with information about the uniform
%          distributed samples
%       - normal: Cell array with information about the normal (Gaussian)
%          distributed samples
%       - truncNormal: Cell array with information about the truncated
%          normal distributed samples
%
%   Each cell array{i} contains information about the run in which the 
%   uncertain coefficient (group) groups{i} is varyied, where groups is
%   returned by the prepareGroupsCoryne()-function.
%
%   The matrices in the cells have the following columns:
%       - dpc
%       - objective value-(min, max, mean, std.)
%       - Euclidean distance of flux vector to nominal flux-(min, max, 
%           mean, std.)
%       - Nominal weighted Euclidean distance (min, max, mean, std.),
%       - Sigma / std. weighted Euclidean distance (min, max, mean, std.),
%       - Time needed to solve an LP problem (min, max, mean)
%       - Iterations of the simplex algorithm (min, max, mean)
    
    if (nargin <= 2)
        % Set defaults
        opts.groupGenerator = @prepareGroupsCoryneTopBM;
        opts.nSamples = 4000;
        opts.sampStr = '4kSamp';
        opts.nBuckets = 30;
        opts.enabledModes = [true, true, true];
        opts.parallel = true;
        opts.dryRun = false;
        opts.noPlots = true;
        opts.figureFormats = [{'.fig'}, {'.tikz'}];
        opts.baseNameStr = 'net_';
    end
    
    if ~isfield(opts, 'groupGenerator') || isempty(opts.groupGenerator)
        opts.groupGenerator = @prepareGroupsCoryneTopBM;
    end

    if ~isfield(opts, 'nBuckets') || (opts.nBuckets <= 0)
        opts.nBuckets = 30;
    end

    if ~isfield(opts, 'nSamples') || (opts.nSamples <= 0)
        opts.nSamples = 4000;
    end
    
    if ~isfield(opts, 'baseNameStr') || isempty(opts.baseNameStr)
        opts.baseNameStr = 'net_';
    else
        % Take care of trailing underscore
        if opts.baseNameStr(end) ~= '_'
            opts.baseNameStr = [opts.baseNameStr, '_'];
        end
    end
    
    if ~isfield(opts, 'sampStr') || isempty(opts.sampStr)
        opts.sampStr = [num2str(floor(opts.nSamples / 1000)) 'kSamp_'];
    else
        % Take care of trailing underscore
        if opts.sampStr(end) ~= '_'
            opts.sampStr = [opts.sampStr, '_'];
        end
    end

    if ~isfield(opts, 'enabledModes') || isempty(opts.enabledModes)
        opts.enabledModes = [true, true, true];
    end
    
    if ~isfield(opts, 'noPlots') || isempty(opts.noPlots)
        opts.noPlots = true;
    end
    
    if ~isfield(opts, 'dryRun') || isempty(opts.dryRun)
        opts.dryRun = false;
    end
    
    if ~isfield(opts, 'parallel') || isempty(opts.parallel)
        opts.parallel = true;
    end

    if ~isfield(opts, 'figureFormats') || isempty(opts.figureFormats)
        opts.figureFormats = [{'.fig'}, {'.tikz'}];
    else
        % Make sure each format has a leading dot
        for i = 1:length(opts.figureFormats)
            if opts.figureFormats{i}(1) ~= '.'
                opts.figureFormats{i} = ['.' opts.figureFormats{i}];
            end
        end
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
            if ~isdir(opts.folders{i}) && opts.enabledModes(i) && ~opts.dryRun
                mkdir(opts.folders{i});
            end
        end
    end
    
    %               -------- vals ---------   -------- diffs ---------
    % Columns: dpc, min, max, mean, std.dev., min, max, mean, std.dev.
    %               --- weightedDiff 1 ----   --- weightedDiff 2 ----
    %               min, max, mean, std.dev., min, max, mean, std.dev.
    %               ---- time ----  - iterations -
    %               min, max, mean, min, max, mean
    metaUn = zeros(length(dpcRange), 23);
    metaNo = zeros(length(dpcRange), 23);
    metaTN = zeros(length(dpcRange), 23);

    modes = {'normal', 'truncnormal', 'uniform'};

    if opts.parallel
        % Open pool for parallel computation
        if matlabpool('size') == 0
            if isfield(opts, 'nWorkers') && ~isempty(opts.nWorkers)
                matlabpool('open', opts.nWorkers);
            else
                matlabpool open;
            end
        end
    end

    firstOptInd = find(model.c ~= 0, 1, 'first');

    tHandle = tic;
    
    groups = opts.groupGenerator(1);

    metaRes.uniform = cell(length(groups), 1);
    metaRes.normal = cell(length(groups), 1);
    metaRes.truncNormal = cell(length(groups), 1);
    
    for k = 1:length(groups)

        for j = 1:length(modes)

            if ~opts.enabledModes(j)
                continue;
            end
            
            for i = 1:length(dpcRange)

                dpc = dpcRange(i);
                groups = opts.groupGenerator(dpc);

                singleCoeff = {groups{k}};
                mode = modes{j};
                
                % Uncomment to select a mode (serial, parallel execution).
                % Parallel execution yields a major speedup.

                if opts.parallel
                    % =============
                    % Parallel mode
                    res = mcGroupsParallel(model, singleCoeff, opts.nSamples, mode);
                else
                    % =============
                    % Serial mode
                    res = mcGroups(model, singleCoeff, opts.nSamples, mode);
                end
                
                if ~opts.dryRun
                    prefix = [opts.folders{j} opts.baseNameStr 'mcG' num2str(k) '_' num2str(dpc) 'pc_' opts.sampStr mode];
                    save([prefix '.mat'], 'res');

                    % ===============
                    % Plot generation

                    if ~opts.noPlots
                        plotHistPercent(res.vals, opts.nBuckets);
                        saveGraph(gcf, [prefix '_objValHist_pc'], opts.figureFormats);

                        plotHistPercent(res.diffs, opts.nBuckets);
                        saveGraph(gcf, [prefix '_2normHist_pc'], opts.figureFormats);

                        plotHistPercent((res.vals - res.origSol(firstOptInd)) / res.origSol(firstOptInd)*100, opts.nBuckets);
                        saveGraph(gcf, [prefix '_objValDeviation_pc'], opts.figureFormats);

                        plotHistPercent((res.diffs) / norm(res.origSol)*100, opts.nBuckets);
                        saveGraph(gcf, [prefix '_2normDeviation_pc'], opts.figureFormats);
                    end
                end
                
                % ================
                % Meta information

                tempMat = [dpc, ...
                        min(res.vals), max(res.vals), mean(res.vals), std(res.vals, 1), ...
                        min(res.diffs), max(res.diffs), mean(res.diffs), std(res.diffs, 1), ...
                        min(res.weightedDiffs(:,1)), max(res.weightedDiffs(:,1)), mean(res.weightedDiffs(:,1)), std(res.weightedDiffs(:,1), 1), ...
                        min(res.weightedDiffs(:,2)), max(res.weightedDiffs(:,2)), mean(res.weightedDiffs(:,2)), std(res.weightedDiffs(:,2), 1), ...
                        min(res.meta(:, 1)), max(res.meta(:, 1)), mean(res.meta(:, 1)), ...
                        min(res.meta(:, 2)), max(res.meta(:, 2)), mean(res.meta(:, 2))];
                if strcmp(mode, 'uniform')
                    metaUn(i, :) = tempMat;
                elseif strcmp(mode, 'normal')
                    metaNo(i, :) = tempMat;
                else
                    metaTN(i, :) = tempMat;
                end
                
            end
            
            if ~opts.dryRun
                if j == 1
                    writeMeta([opts.folders{j} opts.baseNameStr 'mc' num2str(k) '_' opts.sampStr 'Normal.csv'], metaNo);
                elseif j == 2
                    writeMeta([opts.folders{j} opts.baseNameStr 'mc' num2str(k) '_' opts.sampStr 'TruncNormal.csv'], metaTN);
                else
                    writeMeta([opts.folders{j} opts.baseNameStr 'mc' num2str(k) '_' opts.sampStr 'Uniform.csv'], metaUn);
                end
            end
        end
        
        if ~opts.dryRun
            save([opts.folders{j} opts.baseNameStr 'mc' num2str(k) '_' opts.sampStr 'Meta.mat'], 'metaUn', 'metaNo', 'metaTN');
        end
        
        metaRes.uniform{k} = metaUn;
        metaRes.normal{k} = metaNo;
        metaRes.truncNormal{k} = metaTN;
    end
    
    fprintf('Elapsed time: %g sec \n', toc(tHandle));
    
    if opts.parallel
        matlabpool close;
    end
end

function writeMeta(fileName, meta)
    fid = fopen(fileName, 'w');
    fprintf(fid, 'dpc;objMin;objMax;objMean;objStd;diffMin;diffMax;diffMean;diffStd;nwdiffMin;nwdiffMax;nwdiffMean;nwdiffStd;swdiffMin;swdiffMax;swdiffMean;swdiffStd;timeMin;timeMax;timeMean;iterMin;iterMax;iterMean\n');
    fclose(fid);
    dlmwrite(fileName, meta, '-append', 'delimiter', ';');
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
