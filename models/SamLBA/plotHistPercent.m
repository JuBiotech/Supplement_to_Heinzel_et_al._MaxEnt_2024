function plotHistPercent( data, buckets )
%PLOTHISTPERCENT Plots a histogram with a Y-Axis in percent.
%
% Parameters:
%	- data: Data vector for histogram
%	- buckets: Optional count of buckets to use for the histogram
%

    if nargin > 1
        [bin_counts,bin_locations] = hist(data, buckets);
    else
        [bin_counts,bin_locations] = hist(data);
    end

    if nargin > 1
    %	hist(data, buckets);
        bar(bin_locations, 100 * bin_counts ./ sum(bin_counts), 'hist');
    else
    %	hist(data);
        bar(bin_locations, 100 * bin_counts ./ sum(bin_counts), 'hist');
    end

end

