function plotHistFitNormal( data, buckets )
%PLOTHISTFITNORMAL Plots a histogram with a  fitted normal distribution.
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

    bin_width = bin_locations(2) - bin_locations(1);
    hist_area = (bin_width)*100;% * 100;%*(sum(bin_counts));

    if nargin > 1
    %	hist(data, buckets);
        bar(bin_locations, 100 * bin_counts ./ sum(bin_counts), 'hist');
    else
    %	hist(data);
        bar(bin_locations, 100 * bin_counts ./ sum(bin_counts), 'hist');
    end

    mu = mean(data);
    sigma = std(data);

    hold on;

    x = linspace(min(data), max(data), 300);
    plot(x, hist_area .* exp(-((x-mu)/sigma).^2 / 2) / (sqrt(2*pi)*sigma), 'r', 'LineWidth', 1.5);

    hold off;

    fprintf('Fit: Mean = %f and Sigma = %f\n', mu, sigma);

end

