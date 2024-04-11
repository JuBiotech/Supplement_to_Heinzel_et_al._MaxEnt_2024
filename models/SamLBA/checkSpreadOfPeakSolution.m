function checkSpreadOfPeakSolution( res, minDiff, maxDiff )
%CHECKSPREADOFPEAKSOLUTION Selects a solution with Euclidean distance in
% the center of [minDiff, maxDiff] and takes the Euclidean distance of all
% solutions in the dataset with diff in [minDiff, maxDiff] to the selected 
% center. Then plots the newly computed diffs and the standard deviations 
% of the fluxes of the selected solutions.
%
% If the selected solutions spread around the center of the [minDiff, maxDiff]
% bin, the histogram of the newly computed diffs should show only mass
% around 0. The standard deviations should be rather small.
%
% Parameters:
%   - res: Dataset of mcGroups() analysis
%   - minDiff: Minimum diff value for flux vector selection
%   - maxDiff: Maximum diff value for flux vector selection

    meanDiff = (maxDiff + minDiff) / 2;
    [sol,val,~] = extractSolutionWithDiff(res, meanDiff);
    
    ind=find(res.diffs > minDiff & res.diffs < maxDiff);
    peakSols = res.sols(:, ind);
    peakVals = res.vals(ind);
    peakDiff = res.diffs(ind);
    
    diff = zeros(size(peakSols, 2), 1);
    stdFlux = zeros(size(peakSols, 1), 1);
    for i = 1:size(peakSols, 1)
        stdFlux(i) = std(peakSols(i, :));
    end
    
    for i = 1:length(diff)
        diff(i) = norm(sol - peakSols(:, i));
    end
    
%     figure;
    hist(diff, 50);
%     figure;
%     hist(peakVals(peakVals>0),50);
%     figure;plot(peakVals,diff,'.r');hold on; plot(val,0,'sk')
%     figure;
%     plot(peakVals,peakDiff,'.r');hold on; plot(val,meanDiff,'sk')    
%     figure;
%     title('Std for each flux');
%     bar(stdFlux);
%     set(gca, 'XTick', 1:length(stdFlux));
%     xlim([0, length(stdFlux)+1]);
%     ylim([0, max(stdFlux)+1]);
    
    fprintf('Min diff: %g,  Max diff: %g,  maxStd: %g\n', min(diff(diff > 0)), max(diff), max(stdFlux));
end

