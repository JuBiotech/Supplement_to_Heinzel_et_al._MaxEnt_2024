function generateTradeOffCurves(dpcRange, epsVals, eigenLbs, rhos, model, opts)
%GENERATETRADEOFFCURVES Generates some trade-off curves using the given
%uncertainty levels, epsilon (trade-off parameter) values, lower bounds on
%the eigenvalues of the SDP-matrix and rho (immunization parameter) values.
%
% For each parameter combination of dpcRange (uncertainty level), eigenLbs
% and rhos a trade-off curve is plotted and saved.
%
% Parameters:
%   - dpcRange: Vector with uncertainty levels in percent.
%   - epsVals: Vector with epsilon (trade-off parameter) values.
%   - eigenLbs: Vector with lower bounds on the eigenvalues of the
%       SDP-matrix.
%   - rhos: Vector with rho (immunization parameter) values .
%   - model: Model structure, see mcGroups() for hints.
%   - opts: Structure with further options. Optional.
%       o groupGenerator: Function handle with signature
%           groups = groupGenerator(dpc);
%           Returns a cell array with groups having the uncertainty level
%           dpc.
%       o noPlots: If set to true, no plots will be saved. Defaults to
%           false.
%       o figureFormats: Cell array with figure formats (file extensions)
%           in which the plots are saved. Defaults to '.fig' and '.tikz'.
%       o baseNameStr: Base name of the saved plots and data files.
%           baseNameStr_XXXdpc_XXXomega_XXX_elb.Extension for plots
%           baseNameStr_XXXdpc.mat for data files
%           Defaults to 'ro_'.
%       o folder: Folder to save the files in. Defaults to current folder.
%       o includeConcentrations: Set to true to consider an expanded
%           stoichiometry including the concentrations. This only allows
%           the knots in the network to collect substrate instead of
%           creating it from nothing. Optional, defaults to false.
%
% Saves the following fields in the data file:
%   - groups: Groups cell array.
%   - dpc: The uncertainty level used for the data in the file.
%   - fluxes: Matrix with a header in which each column is a computed
%       robust flux vector. See below for the header.
%   - solution: Matrix with a header in which each column contains the
%       computed objective function value und the worst case stoichiometric
%       error (lambda) in this order. See below for the header.
%
% The header of the fields "solution" and "fluxes" contains the following
% data:
%   Line 1: rho value (immunization parameter)
%   Line 2: eigenLb value (lower bound on the eigenvalues of the
%       SDP-matrix)
%   Line 3: epsilon value (trade-off parameter)

    if (nargin < 5)
        % Set defaults
        bmInd = findBiomassFlux( model );
        opts.groupGenerator = @(dpc) prepareGroups('mc', model, dpc, bmInd, false);
        opts.noPlots = false;
        opts.figureFormats = [{'.fig'}, {'.tikz'}];
        opts.baseNameStr = 'ro_';
        opts.folder = '';
        opts.includeConcentrations = false;
    end

    if ~isfield(opts, 'groupGenerator') || isempty(opts.groupGenerator)
        bmInd = findBiomassFlux( model );
        opts.groupGenerator = @(dpc) prepareGroups('mc', model, dpc, bmInd, false);
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

    if ~isfield(opts, 'noPlots') || isempty(opts.noPlots)
        opts.noPlots = false;
    end

    if ~isfield(opts, 'folder') || isempty(opts.folder)
        opts.folder = '';
    else
        % Make sure the folder has a trailing / (filesep)
        if (~isempty(opts.folder)) && (opts.folder(end) ~= filesep)
            opts.folder = [opts.folder filesep];
        end
                
        % Create directory if it does not exist
        if ~isdir(opts.folder)
            mkdir(opts.folder);
        end
    end
    
    if ~isfield(opts, 'includeConcentrations') || isempty(opts.includeConcentrations)
        opts.includeConcentrations = false;
    end
    
    if ~isfield(opts, 'baseNameStr') || isempty(opts.baseNameStr)
        opts.baseNameStr = 'ro_';
    else
        % Take care of trailing underscore
        if opts.baseNameStr(end) ~= '_'
            opts.baseNameStr = [opts.baseNameStr, '_'];
        end
    end
    
    prefix = [opts.folder opts.baseNameStr];
    
    nTotalCount = length(rhos) * length(dpcRange) * length(eigenLbs);
    counter = 0;
    
    if opts.includeConcentrations
        solver = @roSolveModelWithConcentrations;
    else
        solver = @roSolveModel;
    end
    
    tHandle = tic;
    
    for k = 1:length(dpcRange)
        dpc = dpcRange(k);
        groups = opts.groupGenerator(dpc);

        objVals = zeros(length(epsVals), 1);
        lambdaVals = zeros(size(objVals));
        fluxes = zeros(length(model.c)+3, length(rhos) * length(eigenLbs) * length(epsVals));
        solution = zeros(2+3, length(rhos) * length(eigenLbs) * length(epsVals));
        
        for l = 1:length(rhos)

            for j = 1:length(eigenLbs)

                for i = 1:length(epsVals)
                    [ objVals(i), v, lambdaVals(i) ] = solver( epsVals(i), eigenLbs(j), model, groups, rhos(l) );

                    ind = (l-1) * length(eigenLbs) * length(epsVals) + (j-1) * length(epsVals) + i;
                    fluxes(:, ind) = [rhos(l) ; eigenLbs(j) ; epsVals(i) ; v];
                    solution(:, ind) = [rhos(l) ; eigenLbs(j) ; epsVals(i) ; objVals(i) ; lambdaVals(i)];
                end

                if ~opts.noPlots
                    semilogy(epsVals, objVals);
                    hold on;
                    grid on;
                    semilogy(epsVals, lambdaVals, 'r');
                    hold off;

                    title([num2str(dpc) '% Uncertainty - ' num2str(rhos(l)) ' Omega - ' num2str(eigenLbs(j)) ' Eigen LB']);
                    legend('ObjFun', 'Lambda');
                    saveGraph(gcf, [prefix, num2str(dpc) 'dpc_' num2str(rhos(l)) 'omega_' num2str(eigenLbs(j)) '_elb'], opts.figureFormats);
                end

                counter = counter + 1;
                display(['Status (%): ' num2str(counter / nTotalCount * 100)]);
            end
        end
        save([prefix num2str(dpc) 'dpc.mat'], 'fluxes', 'solution', 'groups', 'dpc');
    end
    
    fprintf('Elapsed time: %g sec \n', toc(tHandle));
    
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
