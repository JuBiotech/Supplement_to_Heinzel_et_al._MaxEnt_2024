function res = generateRobustRanges( model, groups, normPCs, singlePCs, maxIter, angles, nSamples, nSlice, nWorkers )
%GENERATEROBUSTRANGES Calculates robust ranges AKA the robust coefficient
%space for each value pair of (normPC, singlePC) in 1D, 2D and nD.
%
% A robust range is the range in which a coefficent group can vary without
% the solution of the FBA with the uncertain data violating the conditions
%	norm(v_nom - v_sol) <= norm(v_nom) * normPC / 100
%	abs(v_nom - v_sol) <= v_nom * singlePC / 100
% So there are boundaries on each flux vector component and on the flux vector
% as a whole.
%
% Parameters:
%	- model: Model structure, see mcGroups() for hints.
%	- groups: Cell array with group structures. See MCGroups.m for hints.
%	- normPCs: Vector with allowed maximum relative deviation of the flux 
%       vector norm
%	- singlePCs: Vector with allowed maximum relative deviation of each 
%       flux vector component.
%	- maxIter: Maximum iteration count used for bisection search.
%   - angles: Vector with angels (in degree) for the 2D computation.
%   - nSamples: Number of samples for the nD computation.
%   - nSlice: Number of slices for the nD computation, see
%       inverseMCGroupsCoeffEstimationSphere() for details.
%   - nWorkers: Number of worker threads for the nD computation, see
%       inverseMCGroupsCoeffEstimationSphere() for details.
%
% Returns a structure with the following fields:
%   - nWorkers: Number of workers for nD computation
%   - angles: Angle vector for 2D computation
%   - nSamples: Number of samples for nD computation
%   - nSlice: Number of slices for nD computation
%   - maxIter: Maximum iteration count for bisection search
%   - groups: Cell array with groups
%   - runtime: Total runtime of the script in seconds
%   - runs: Cell array of structures with data of each run (i.e. 
%       (normPC, singlePC) pair) with the following fields:
%       o normPC: normPC value
%       o singlePC: singlePC value
%       o runtimes: Vector with runtimes for 1D, 2D and nD computation (in
%           this order)
%       o ranges1D, table1D: Output of calcRobustRanges().
%       o areaTable2D, polys2D: Output of calcRobustRangesDouble().
%       o resSphere: Output of inverseMCGroupsCoeffEstimationSphere().
%       o completeRuntime: Runtime for this run.

    res.runs = cell(length(normPCs), 1);
    res.nWorkers = nWorkers;
    res.angles = angles;
    res.nSamples = nSamples;
    res.nSlice = nSlice;
    res.maxIter = maxIter;
    res.groups = groups;

    tAll = tic;
    for i = 1:length(normPCs)
        
        fprintf('Starting run with %g norm- and %g individual allowed deviation\n', normPCs(i), singlePCs(i));

        set.normPC = normPCs(i);
        set.singlePC = singlePCs(i);
        set.runtimes = zeros(3, 1);
        tRun = tic;

        tMethod = tic;
        [set.ranges1D, set.table1D] = calcRobustRanges( model, groups, normPCs(i), singlePCs(i), maxIter );
        set.runtimes(1) = toc(tMethod);
        fprintf('Finished 1D: %g sec = %g h %g min\n', set.runtimes(3), ...
            floor(set.runtimes(3) / 3600), floor((set.runtimes(3) - floor(set.runtimes(3) / 3600)*3600) / 60));
        
        tMethod = tic;
        [set.areaTable2D, set.polys2D] = calcRobustRangesDouble( model, groups, angles, normPCs(i), singlePCs(i), maxIter, nWorkers );
        set.runtimes(2) = toc(tMethod);
        fprintf('Finished 2D: %g sec = %g h %g min\n', set.runtimes(3), ...
            floor(set.runtimes(3) / 3600), floor((set.runtimes(3) - floor(set.runtimes(3) / 3600)*3600) / 60));
        
        tMethod = tic;
%         set.resSphere = inverseMCGroupsCoeffEstimationSphere( model, groups, nSamples, normPCs(i), singlePCs(i), maxIter, nSlice, true, nWorkers );
        set.runtimes(3) = toc(tMethod);
        fprintf('Finished nD: %g sec = %g h %g min\n', set.runtimes(3), ...
            floor(set.runtimes(3) / 3600), floor((set.runtimes(3) - floor(set.runtimes(3) / 3600)*3600) / 60));
        
        set.completeRuntime = toc(tRun);
        fprintf('Overall time for this run: %g sec = %g h %g min\n', set.completeRuntime, ...
            floor(set.completeRuntime / 3600), floor((set.completeRuntime - floor(set.completeRuntime / 3600)*3600) / 60));
        
        res.runs{i} = set;
    end
    res.runtime = toc(tAll);
    fprintf('Overall time: %g sec = %g h %g min\n', res.runtime, ...
        floor(res.runtime / 3600), floor((res.runtime - floor(res.runtime / 3600)*3600) / 60));
end