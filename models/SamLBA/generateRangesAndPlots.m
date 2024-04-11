function generateRangesAndPlots(model, dpc, groups, opts)
%GENERATERANGESANDPLOTS Scans a range of one (or two) coefficients and saves the plots.
% This function scans every uncertain coefficient group in a given (relative) range, e.g. 0 - 200%.
% The range is divided by the amount of steps and for each step the corresponding FBA problem
% is solved. The calculated objectiveValue and solution deviation are plotted and saved.
%
% The function also does 2D range scans with two coefficient groups. At the moment the step count
% in 2D mode is limited to 300x300.
%
% Parameters:
%   - model: Model structure.
%   - dpc: Range in percent. Only used to generate filenames. The range and count of 
%		steps of a group is assigned by prepareGroups() in range mode.
%	- groups: Cell array with coefficient groups. See rangeGroups.m for hints.
%   - opts: Optional. Structure with options. Fields are:
%       o folders: Optional. Cell array with names of directories for 1D and 2D ranges.
%           The directory will be used to dump the data and is created if it does not exist. 
%               Order: 1D, 2D.
%       o nSteps: Optional. Number of steps for each range. Used for filenames only. 
%       o baseNameStr: Prefix of saved files. Could be network name.
%           Defaults to 'net'.
%       o stepStr: String. Used for file name composition. Defaults
%           to '[nSteps / 1000]kStep'.
%           File names: baseNameStr_DPCpc_stepStr_GroupNumber
%       o noPlots: If true no plots will be saved. Data will still be saved
%           to mat-Files. Default is true.
%       o dryRun: If true data will not be saved to files. Default is
%           false (save data to files and directories).
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

    if (nargin <= 2)
        % Set defaults
        opts.nSteps = 10000;
        opts.stepStr = '10kStep_';
        opts.parallel = true;
        opts.dryRun = false;
        opts.noPlots = true;
        opts.restartParallelSession = false;
        opts.figureFormats = [{'.fig'}, {'.tikz'}];
        opts.baseNameStr = 'net_';
        opts.enabledModes = [true, true];
    end

    if ~isfield(opts, 'baseNameStr') || isempty(opts.baseNameStr)
        opts.baseNameStr = 'net_';
    else
        % Take care of trailing underscore
        if opts.baseNameStr(end) ~= '_'
            opts.baseNameStr = [opts.baseNameStr, '_'];
        end
    end

    if ~isfield(opts, 'nSteps') || (opts.nSteps <= 0)
        opts.nSteps = 10000;
    end
    
    if ~isfield(opts, 'stepStr') || isempty(opts.stepStr)
        opts.stepStr = [num2str(floor(opts.nSteps / 1000)) 'kSteps_'];
    else
        % Take care of trailing underscore
        if opts.stepStr(end) ~= '_'
            opts.stepStr = [opts.stepStr, '_'];
        end
    end
    
    if ~isfield(opts, 'restartParallelSession') || isempty(opts.restartParallelSession)
        opts.restartParallelSession = false;
    end
    
    if ~isfield(opts, 'enabledModes') || isempty(opts.enabledModes)
        opts.enabledModes = [true, true];
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
        opts.folders = {'', ''};
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
    
    firstOptInd = find(model.c ~= 0, 1, 'first');
    groupPairs = nchoosek(1:length(groups), 2);

    % Count runs for status
    nTotalRuns = length(groups) + size(groupPairs, 1);
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

    tHandle = tic;

    for i = 1:nTotalRuns
        
        if i <= length(groups)
            curGroups = groups(i);
            postfix = num2str(i);
            curFolder = opts.folders{1};

            if ~opts.enabledModes(1)
                continue;
            end
        elseif i > length(groups)
            
            if ~opts.enabledModes(2)
                continue;
            end
            
            postfix = '';
            for j = 1:size(groupPairs, 2)
                ind = groupPairs(i - length(groups), j);
                postfix = [postfix num2str(ind) '-' ];
            end
            postfix = postfix(1:end-1);

            curFolder = opts.folders{2};
            
            % -------- DEBUG
            % Constrain to 300 steps
            
            for j = 1:size(groupPairs, 2)
                ind = groupPairs(i - length(groups), j);
                g = groups{ind};
                
                if isfield(g, 'range')
                    g.range(3) = min(g.range(3), 300);
                end
                
                groups{ind} = g;
            end
            
            curGroups = groups(groupPairs(i - length(groups), :));
        end
        
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
            res = rangeGroups(model, curGroups);

            if opts.restartParallelSession
                matlabpool close;
            end
        else
            % =============
            % Serial mode
            res = rangeGroups(model, curGroups, false);
        end

        if ~opts.dryRun
            prefix = [curFolder opts.baseNameStr num2str(dpc) 'pc_' opts.stepStr postfix];
            save([prefix '.mat'], 'res');

            % ===============
            % Plot generation

            if ~opts.noPlots
                if length(curGroups) == 1
                    
                    if ~isfield(curGroups{1}, 'rawRange')
                        xRange = linspace(curGroups{1}.range(1), curGroups{1}.range(2), curGroups{1}.range(3)) .* 100;
                    else
                        xRange = curGroups{1}.rawRange;
                    end
                    
                    plot(xRange, res.vals);
                    saveGraph(gcf, [prefix '_objValHist_pc'], opts.figureFormats);

                    plot(xRange, res.diffs);
                    saveGraph(gcf, [prefix '_2normHist_pc'], opts.figureFormats);

                    plot(xRange, (res.vals - res.origSol(firstOptInd)) / res.origSol(firstOptInd)*100);
                    saveGraph(gcf, [prefix '_objValDeviation_pc'], opts.figureFormats);

                    plot(xRange, (res.diffs) / norm(res.origSol)*100);
                    saveGraph(gcf, [prefix '_2normDeviation_pc'], opts.figureFormats);
                else
                    visualizeNDArray(curGroups, res.vals);
                    saveGraph(gcf, [prefix '_objValHist_pc'], opts.figureFormats);

                    visualizeNDArray(curGroups, res.diffs);
                    saveGraph(gcf, [prefix '_2normHist_pc'], opts.figureFormats);

                    visualizeNDArray(curGroups, (res.vals - res.origSol(firstOptInd)) / res.origSol(firstOptInd)*100);
                    saveGraph(gcf, [prefix '_objValDeviation_pc'], opts.figureFormats);

                    visualizeNDArray(curGroups, (res.diffs) / norm(res.origSol)*100);
                    saveGraph(gcf, [prefix '_2normDeviation_pc'], opts.figureFormats);
                end
            end
        end

        pcBeforeRun = (nCurrentRun - 1) / nTotalRuns * 10;
        pcAfterRun = nCurrentRun / nTotalRuns * 10;
        if floor(pcBeforeRun) < floor(pcAfterRun)
            display(['Status (%): ' num2str(floor(pcAfterRun * 10))]);
        end

        nCurrentRun = nCurrentRun + 1;
    end
    
    fprintf('Elapsed time: %g sec \n', toc(tHandle));
    
    if opts.parallel && ~opts.restartParallelSession
        matlabpool close;
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
