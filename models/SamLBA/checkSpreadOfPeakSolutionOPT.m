function checkSpreadOfPeakSolutionOPT( res, minVal, maxVal )
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

    meanVal = (maxVal + minVal) / 2;
    sol = extractSolutionWithOpt(res, meanVal);
    
    peakSols = res.sols(:, res.vals > minVal & res.vals < maxVal);

    diff = zeros(size(peakSols, 2), 1);
    for i = 1:length(diff)
        diff(i) = norm(sol - peakSols(:, i));
    end
    
    hist(diff, 50);
    
    fprintf('Mean Value: %g\n', meanVal);
end

function [sol, val, sample] = extractSolutionWithOpt(res, opt)
%EXTRACTSOLUTIONWITHDIFF Extracts a solution with a given diff-value from
%the dataset.
%
% Selectes the flux vector which has the nearest diff-value to the given
% one.
%
% Parameters:
%   - res: Dataset of mcGroups() analysis
%   - opt: objective value
%
% Returns:
%   - sol: Flux vector with the closest diff-value
%   - val: Objective function value of the closest flux vector
%   - sample: Sample used to obtain the flux vector

    [val, index] = min(abs(res.vals - opt));

    sol = res.sols(:, index);
    sample = res.samples(:, index);
    val = res.vals(index);
    diff = res.diffs(index);
    
    fprintf('I selected a sample with diff = %g and val = %g\n', diff, val);
end