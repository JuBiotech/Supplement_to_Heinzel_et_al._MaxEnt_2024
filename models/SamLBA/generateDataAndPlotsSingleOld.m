function metaRes = generateDataAndPlotsSingle(model, dpcRange, opts)
%GENERATEDATAANDPLOTSSINGLE Does MC simulations in which the uncertain
%coefficients are randomly disturbed by a given factor (dpcRange) and 
%conducts a FBA for each realisation of the uncertain coefficients.
%
% Instead of perturbing every uncertain coefficient simultaneously, only
% only one coefficient at a time is perturbed and the others are kept
% nominal.
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
%       o sampStr: Sample string. Used for file name composition.
%           File names: coryne_optBM_mcGx_DPCpc_sampStr_mode
%       o nBuckets: Number of buckets for histograms.
%       o enabledModes: Vector of booleans. If an entry is set to false,
%           the MC run with the corresponding perturbation mode will not be
%           performed. Order: [normal, truncnormal, uniform]
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
    
    if (nargin <= 2)
        % Set defaults
        opts.groupGenerator = @prepareGroupsCoryneTopBM;
        opts.nSamples = 4000;
        opts.sampStr = '4kSamp';
        opts.groupGenerator = @prepareGroupsCoryneTopBM;
        opts.nBuckets = 30;
        opts.enabledModes = [true, true, true];
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
    
    % Open pool for parallel computation
    % Comment if parallel execution is disabled!
    if matlabpool('size') == 0
        matlabpool;
    end
    
    groups = opts.groupGenerator(1);
%    groups = prepareGroupsCoryneTopBM(1);

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
%                groups = prepareGroupsCoryneTopBM(dpc);

                singleCoeff = {groups{k}};
                mode = modes{j};
                
                % Uncomment to select a mode (serial, parallel execution).
                % Parallel execution yields a major speedup.

                % =============
                % Parallel mode
                res = mcGroupsParallel(model, singleCoeff, opts.nSamples, mode);

                % =============
                % Serial mode
                %res = mcGroups(model, singleCoeff, 10000, mode);

                prefix = [opts.folders{j} 'coryne_optBM_mcG' num2str(k) '_' num2str(dpc) 'pc_' opts.sampStr '_' mode];
                save([prefix '.mat'], 'res');

                % ===============
                % Plot generation
                %   Comment lines starting with "matlab2tikz" to disable
                %   tikz-file generation.

                plotHistPercent(res.vals, opts.nBuckets);
                saveas(gcf, [prefix '_objValHist_pc.fig']);
                matlab2tikz([prefix '_objValHist_pc.tikz'], 'silent', true, 'height', '\mcHistHeight', 'width', '\mcHistWidth');

                plotHistPercent(res.diffs, opts.nBuckets);
                saveas(gcf, [prefix '_2normHist_pc.fig']);
                matlab2tikz([prefix '_2normHist_pc.tikz'], 'silent', true, 'height', '\mcHistHeight', 'width', '\mcHistWidth');

                plotHistPercent((res.vals - res.origSol(300)) / res.origSol(300)*100, opts.nBuckets);
                saveas(gcf, [prefix '_objValDeviation_pc.fig']);
                matlab2tikz([prefix '_objValDeviation_pc.tikz'], 'silent', true, 'height', '\mcHistHeight', 'width', '\mcHistWidth');

                plotHistPercent((res.diffs) / norm(res.origSol)*100, opts.nBuckets);
                saveas(gcf, [prefix '_2normDeviation_pc.fig']);
                matlab2tikz([prefix '_2normDeviation_pc.tikz'], 'silent', true, 'height', '\mcHistHeight', 'width', '\mcHistWidth');

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
            
            if j == 1
                writeMeta([opts.folders{j} 'coryne_optBM_mc' num2str(k) '_' opts.sampStr '_Normal.csv'], metaNo);
            elseif j == 2
                writeMeta([opts.folders{j} 'coryne_optBM_mc' num2str(k) '_' opts.sampStr '_TruncNormal.csv'], metaTN);
            else
                writeMeta([opts.folders{j} 'coryne_optBM_mc' num2str(k) '_' opts.sampStr '_Uniform.csv'], metaUn);
            end
        end
        
        save([opts.folders{j} 'coryne_optBM_mc' num2str(k) '_' opts.sampStr '_Meta.mat'], 'metaUn', 'metaNo', 'metaTN');

        metaRes.uniform{k} = metaUn;
        metaRes.normal{k} = metaNo;
        metaRes.truncNormal{k} = metaTN;
    end
    
    % Comment if parallel execution is disabled!
    matlabpool close;
        
    fprintf('Needed time: %g sec \n', toc());
end

function writeMeta(fileName, meta)
    fid = fopen(fileName, 'w');
    fprintf(fid, 'dpc;objMin;objMax;objMean;objStd;diffMin;diffMax;diffMean;diffStd\n');
    fclose(fid);
    dlmwrite(fileName, meta, '-append', 'delimiter', ';');
end
