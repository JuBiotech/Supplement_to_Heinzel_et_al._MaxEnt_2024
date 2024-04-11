function saveSingleMetaData( fileName, meta )
%PLOTSINGLEMETADATA Summary of this function goes here
%   Detailed explanation goes here

    n = length(meta);
    sz = size(meta{1});
    data = zeros(sz(1), n * (sz(2)-1) + 1);
    data(:, 1:sz(2)) = meta{1};
    for i = 2:n
        l = (2 + (i-1)*(sz(2)-1));
        u = ((1 + i*(sz(2)-1)));
        mat = meta{i}(:, 2:end);
       data(:, (2 + (i-1)*(sz(2)-1)):((1 + i*(sz(2)-1)))) = mat;
    end
    
    writeMeta(fileName, n, data);
end

function writeMeta(fileName, n, data)
    fid = fopen(fileName, 'w');
    fprintf(fid, 'dpc;G1objMin;G1objMax;G1objMean;G1objStd;G1diffMin;G1diffMax;G1diffMean;G1diffStd');

    for i = 2:n
        fprintf(fid, ';G%1$gobjMin;G%1$gobjMax;G%1$gobjMean;G%1$gobjStd;G%1$gdiffMin;G%1$gdiffMax;G%1$gdiffMean;G%1$gdiffStd', i);
    end
    fprintf(fid, '\n');
    
    fclose(fid);
    dlmwrite(fileName, data, '-append', 'delimiter', ';');
end