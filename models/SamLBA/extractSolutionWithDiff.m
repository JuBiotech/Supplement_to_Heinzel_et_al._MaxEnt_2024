function [sol, val, sample] = extractSolutionWithDiff(res, diff)
%EXTRACTSOLUTIONWITHDIFF Extracts a solution with a given diff-value from
%the dataset.
%
% Selectes the flux vector which has the nearest diff-value to the given
% one.
%
% Parameters:
%   - res: Dataset of mcGroups() analysis
%   - diff: Diff-value
%
% Returns:
%   - sol: Flux vector with the closest diff-value
%   - val: Objective function value of the closest flux vector
%   - sample: Sample used to obtain the flux vector

    [val, index] = min(abs(res.diffs - diff));

    sol = res.sols(:, index);
    sample = res.samples(:, index);
    val = res.vals(index);
    
    fprintf('I selected a sample with diff = %g and val = %g\n', res.diffs(index), val);
end