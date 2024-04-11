function [res, resave] = ppWeightedDiffs(fileName, res)

    resave = true;

    wghtOrig = abs(res.origSol);
    wghtOrig(wghtOrig <= 1e-8) = 1;

    wghtSigma = ones(length(res.origSol), 1);
    for i = 1:length(res.origSol)
        compSigma = std(res.sols(i, :));
        if compSigma > 1e-8
            wghtSigma(i) = compSigma;
        end
    end
    
    weightedDiffsA = zeros(length(res.diffs), 1);
    weightedDiffsB = zeros(length(res.diffs), 1);
    
    sols = res.sols;
    
    parfor i = 1:size(res.sols, 2)
%        cur = res.sols(:, i);
        cur = sols(:, i);
        
        weightedDiffsA(i) = sqrt(sum( ((res.origSol - cur) ./ wghtOrig).^2 ));
        weightedDiffsB(i) = sqrt(sum( ((res.origSol - cur) ./ wghtSigma).^2 ));
    end
    res.weightedDiffs = [weightedDiffsA, weightedDiffsB];
    
end