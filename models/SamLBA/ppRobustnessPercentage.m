function [ result, resave ] = ppRobustnessPercentage( fileName, res, origSol, normPC, singlePC )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    resave = false;

    diffSol = abs(origSol);
    diffSol(diffSol == 0) = 1;
    diffSol = diffSol * singlePC / 100;
    
    cmpNorm = norm(origSol) * normPC / 100;
    
    result = false(size(res.sols, 2), 1);
    
    for i = 1:size(res.sols, 2)
        
        result(i) = all(abs(res.sols(:,i) - origSol) <= diffSol);
        
        if isfinite(normPC)
            result(i) = result(i) && (norm(res.sols(:, i) - origSol) <= cmpNorm);
        end
    end
    
    maxVals = nan(size(res.samples, 1), 1);
    minVals = nan(size(res.samples, 1), 1);
    
    % Obtain min and max of coefficients
    if any(result)
        for i = 1:size(res.samples, 1)
           maxVals(i) = max(res.samples(i, result));
           minVals(i) = min(res.samples(i, result));
        end
    end
    
    % Parse percentage from fileName
    pcInd = strfind(fileName, 'pc_');
    underScore = find(fliplr(fileName(1:pcInd-1)) == '_', 1, 'first');
    underScore = pcInd - underScore;
    dpc = str2double(fileName(underScore+1:pcInd-1));
    
    % Save results
    result = [{[dpc, sum(result), sum(result) / size(res.sols, 2) * 100]}, ...
        {full([minVals, maxVals, res.origData])}, {fileName}, {result}];
end

