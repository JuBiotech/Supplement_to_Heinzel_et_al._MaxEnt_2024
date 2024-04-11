function [res, resave] = ppCalculateDiffNorms(fileName, res)

    resave = true;
   
    sols = res.sols;
    
    wght = abs(res.origSol) .* (abs(res.origSol) > 1e-8) + (abs(res.origSol) <= 1e-8);
    
    diffs1Norm = zeros(length(res.diffs), 1);
    diffsMaxNorm = zeros(length(res.diffs), 1);
    diffs1NormNomWeight = zeros(length(res.diffs), 1);
    diffs2NormNomWeight = zeros(length(res.diffs), 1);
    diffsMaxNormNomWeight = zeros(length(res.diffs), 1);
%    parfor i = 1:size(res.sols, 2)
    for i = 1:size(res.sols, 2)
        cur = sols(:, i);
        
        diffs1Norm(i) = sum( abs((res.origSol - cur)));
        diffsMaxNorm(i) = max(abs(res.origSol - cur));
        diffs1NormNomWeight(i) = sum( abs((res.origSol - cur) ./ wght));
        diffsMaxNormNomWeight(i) = max( abs((res.origSol - cur) ./ wght));
        diffs2NormNomWeight(i) = norm( (res.origSol - cur) ./ wght);
    end
    res.nomWeightedDiffs = [diffs1NormNomWeight, diffs2NormNomWeight, diffsMaxNormNomWeight];
    res.diffs1Norm = diffs1Norm;
    res.diffsMaxNorm = diffsMaxNorm;
end