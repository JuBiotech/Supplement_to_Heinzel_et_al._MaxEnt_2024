function [data, dists] = matchSamplesToAltSolutions( res, altOpts )
% function [data, stat] = matchSamplesToAltSolutions( res, altOpts )
%MATCHSAMPLESTOALTSOLUTIONS Matches the samples to the nearest alternative
%solution.
%   
% Looks at each sample of the given dataset and computes the index of the
% nearest alternative solution given. If the sampled flux vector is the
% null-vector, the index 0 is returned.
%
% Parameters:
%   - res: Dataset of mcGroups() analysis
%   - altOpts: Matrix in which each column denotes an alternative solution
%
% Returns:
%   - data: Matrix in which each row corresponds to a sampled flux vector. 
%       The first column gives the index of the nearest alternative 
%       solution or 0 for the null-vector. The second column contains the
%       Euclidean distance to the nearest alternative solution.
%   - stat: Matrix in which each row corresponds to an alternative
%       solution. The columns contain (in this order): The mean distance to
%       this alternative solution, the standard deviation of the distances,
%       the minimum and maximum distance, the minimum, maximum and mean
%       diff-Values of the sampled flux vectors nearest to this alternative
%       solution

    data = zeros(length(res.diffs), 2);
    dists = zeros(length(res.diffs),size(altOpts, 2));

    for i = 1:length(res.diffs)
        
        for j = 1:size(altOpts, 2)
            dists(i,j) = norm(res.sols(:, i) - altOpts(:, j));
        end
        
        [val, ind] = min(dists(i,:));
        
        if norm(res.sols(:, i)) == 0
            ind = 0;
            val = 0;
        end
        
        data(i, :) = [ind, val];
    end
    
%     stat = zeros(size(altOpts, 2), 7);
%     for i = 1:size(altOpts, 2)
%         
%         inds = data(:, 1) == i;
%         stat(i, :) = [mean(data(inds, 2)), std(data(inds, 2)), min(data(inds, 2)), max(data(inds, 2)), ...
%             min(res.diffs(inds)), max(res.diffs(inds)), mean(res.diffs(inds))];
%         
%     end
    
%     disp(stat);
end

