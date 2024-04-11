function [metaNo, metaTN, metaUn, plugins] = generateDataAndPlots(model, dpcRange, opts)
%GENERATEDATAANDPLOTS Does MC simulations in which the uncertain
%coefficients are randomly disturbed by a given factor (dpcRange) and 
%conducts a FBA for each realisation of the uncertain coefficients.
%
% Results of each run are saved to mat- and csv-files. Histograms are 
% generated and saved as figure and tikz-file if noPlots option is 
% disabled. So by default no plots are generated.
%
% Call: 
%   [metaNo, metaTN, metaUn] = generateDataAndPlots(model, dpcRange, opts)
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
%           Defaults to prepareGroups() with a guessed biomass flux.
%       o nSamples: Number of samples for each run.
%       o baseNameStr: Prefix of saved files. Could be network name.
%           Defaults to 'net'.
%       o typeStr: Type string. Used for file name composition. Defaults to
%           'mcBM'.
%       o sampStr: Sample string. Used for file name composition. Defaults
%           to '[nSamples / 1000]kSamp'.
%           File names: baseNameStr_typeStr_DPCpc_sampStr_mode
%       o nBuckets: Number of buckets for histograms.
%       o enabledModes: Vector of booleans. If an entry is set to false,
%           the MC run with the corresponding perturbation mode will not be
%           performed. Order: [normal, truncnormal, uniform]
%       o noPlots: If true no plots will be saved. Data will still be saved
%           to mat-Files. Default is true.
%       o dryRun: If true data will not be saved to files. Default is
%           false (save data to files and directories).
%       o plugins: A cell array of plugins. A plugin is a function handle 
%           with the following signature: 
%               data = plugin(S, solution, objValue, diff)
%           The result of the plugin will be saved in a cell array.
%           There are no requirements for the return values of a plugin as
%           everything is stored in a cell array. This field is optional.
%       o parallel: If true the computation will be parallelized. The 
%           Parallel Computation Toolbox is required for this feature to
%           work. Default is true (use parallelization).
%       o restartParallelSession: If true the workers will be recreated
%           after each parfor-run. As the glpkmex-Wiki suggests this is
%           necessary to ensure a steady memory consumption (glpkmex in
%           conjunction with parfor causes memory leaks). Default is false.
%       o nWorkers: Number of workers used by the Parallel Computation
%           Toolbox. Optional, leave out to use default configuration.
%       o figureFormats: Cell array with figure formats used by saveas().
%           A figure will be saved in each given format. Use 'tikz' to save
%           a tikz-file. Defaults to [{'fig'}, {'tikz'}]
%
% Returns: Two matrices with meta data about the MC runs.
%   - metaNo: Contains information about the MC runs with Gaussian (normal)
%      distribution.
%   - metaTN: Contains information about the MC runs with truncated normal
%      distribution.
%   - metaUn: Contains information about the MC runs with uniform
%      distribution.
% Each matrix hast the columns dpc, objective value-(min, max, mean, std.),
%   Euclidean distance of flux vector to nominal flux-(min, max, mean,
%   std.), nominal weighted Euclidean distance of flux vector to nominal 
%   flux-(min, max, mean, std.), sigma weighted Euclidean distance of flux 
%   vector to nominal flux-(min, max, mean, std.), time needed to solve an 
%   LP problem (min, max, mean), number of iterations (min, max, mean).
%   - plugins: Structure containing cell arrays with plugin data for each
%       dpc-Value, for each sample, for each plugin (nested cell arrays in 
%       this order).

    if (nargin <= 2)
        % Set defaults
        bmInd = findBiomassFlux(model);
        if isempty(bmInd)
            error('Please specify the groupGenerator.');
        end
        opts.groupGenerator = @(dpc) prepareGroups('mc', model, dpc, bmInd, true);

        opts.nSamples = 10000;
        opts.sampStr = '10kSamp_';
        opts.nBuckets = 30;
        opts.enabledModes = [true, true, true];
        opts.typeStr = 'mcBM_';
        opts.parallel = true;
        opts.dryRun = false;
        opts.noPlots = true;
        opts.plugins = [];
        opts.restartParallelSession = false;
        opts.figureFormats = [{'.fig'}, {'.tikz'}];
        opts.baseNameStr = 'net_';
    end
    
    if ~isfield(opts, 'groupGenerator') || isempty(opts.groupGenerator)
        bmInd = findBiomassFlux(model);
        if isempty(bmInd)
            error('Please specify the groupGenerator.');
        end
        opts.groupGenerator = @(dpc) prepareGroups('mc', model, dpc, bmInd, true);
    end
    
    if ~isfield(opts, 'plugins')
        opts.plugins = [];
    end

    if ~isfield(opts, 'nBuckets') || (opts.nBuckets <= 0)
        opts.nBuckets = 30;
    end

    if ~isfield(opts, 'baseNameStr') || isempty(opts.baseNameStr)
        opts.baseNameStr = 'net_';
    else
        % Take care of trailing underscore
        if opts.baseNameStr(end) ~= '_'
            opts.baseNameStr = [opts.baseNameStr, '_'];
        end
    end
    
    if ~isfield(opts, 'typeStr') || isempty(opts.typeStr)
        opts.typeStr = 'mcBM_';
    else
        % Take care of trailing underscore
        if opts.typeStr(end) ~= '_'
            opts.typeStr = [opts.typeStr, '_'];
        end
    end

    if ~isfield(opts, 'nSamples') || (opts.nSamples <= 0)
        opts.nSamples = 10000;
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
    
    if ~isfield(opts, 'restartParallelSession') || isempty(opts.restartParallelSession)
        opts.restartParallelSession = false;
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
            if (~isempty(opts.folders{i})) && (opts.folders{i}(end) ~= filesep)
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

    plugins.normal = cell(length(dpcRange), 1);
    plugins.truncnormal = cell(length(dpcRange), 1);
    plugins.uniform = cell(length(dpcRange), 1);
    
    % Count runs for status
    nTotalRuns = sum(ones(size(opts.enabledModes)) * length(dpcRange) .* (opts.enabledModes == true));
    nCurrentRun = 1;
    
    if opts.parallel && ~opts.restartParallelSession
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

    for i = 1:length(dpcRange)
        dpc = dpcRange(i);
        
        for j = 1:length(modes)
            
            if ~opts.enabledModes(j)
                continue;
            end
            
            mode = modes{j};
            groups = opts.groupGenerator(dpc);
            
            if opts.parallel
                
                if opts.restartParallelSession
                    if isfield(opts, 'nWorkers') && ~isempty(opts.nWorkers)
                        matlabpool('open', opts.nWorkers);
                    else
                        matlabpool open;
                    end
                end
                
                % =============
                % Parallel mode
                res = mcGroupsParallel(model, groups, opts.nSamples, mode, opts.plugins);

                if opts.restartParallelSession
                    matlabpool close;
                end
            else
                % =============
                % Serial mode
                res = mcGroups(model, groups, opts.nSamples, mode, opts.plugins);
            end

            if ~opts.dryRun
                prefix = [opts.folders{j} opts.baseNameStr opts.typeStr num2str(dpc) 'pc_' opts.sampStr mode];
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
    
                if isfield(res, 'pluginData')
                    plugins.uniform{i} = res.pluginData;
                end
            elseif strcmp(mode, 'normal')
                metaNo(i, :) = tempMat;

                if isfield(res, 'pluginData')
                    plugins.normal{i} = res.pluginData;
                end
            else
                metaTN(i, :) = tempMat;
                
                if isfield(res, 'pluginData')
                    plugins.truncnormal{i} = res.pluginData;
                end
            end
            
            pcBeforeRun = (nCurrentRun - 1) / nTotalRuns * 10;
            pcAfterRun = nCurrentRun / nTotalRuns * 10;
            if floor(pcBeforeRun) < floor(pcAfterRun)
                display(['Status (%): ' num2str(floor(pcAfterRun * 10))]);
            end
            
            nCurrentRun = nCurrentRun + 1;
        end
    end

    if ~opts.dryRun
        if opts.enabledModes(1)
            writeMeta([opts.baseNameStr opts.typeStr opts.sampStr 'Normal.csv'], metaNo);
        end
        if opts.enabledModes(2)
            writeMeta([opts.baseNameStr opts.typeStr opts.sampStr 'TruncNormal.csv'], metaTN);
        end
        if opts.enabledModes(3)
            writeMeta([opts.baseNameStr opts.typeStr opts.sampStr 'Uniform.csv'], metaUn);
        end
        save([opts.baseNameStr opts.typeStr opts.sampStr 'Meta.mat'], 'metaUn', 'metaNo', 'metaTN');
        
        if ~isempty(opts.plugins)
            save([opts.baseNameStr opts.typeStr opts.sampStr 'Plugins.mat'], 'plugins');
        end
    end
    
    fprintf('Elapsed time: %g sec \n', toc(tHandle));
    
    if opts.parallel && ~opts.restartParallelSession
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
