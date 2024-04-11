% function [ sols] = compareLooplessFBA( model)
function [ sols] = compareLooplessFBA( model, objValFrac )
%COMPARELOOPLESSFBA Uses the different loop-avoidance methods and puts the
%resulting flux distributions in one matrix.
%
% Parameters:
%   - model: Model structure
%   - objValFrac: Minimal objective function value as described by de Leon
%       et al. 2008.
%
% Returns a matrix in which each but the first row corresponds to a flux.
% The first row contains the objective function value. 
% Columns are (in this order): Original, Norm-1 minimzation, de Leon
% solution (first mode), de Leon solution (second mode), Schellenberger
% solutions.
    lsg=optimizeCbModel(model);
    solStart=lsg.x;
    valStart=lsg.f;
    
    [solNom, valNom, ~] = quickSolveFBA(model, false);
    [solNorm, valNorm, ~] = quickSolveFBA(model, true);

%     [ solLeon1, valLeon1, nActiveRxn1, stat ] = deLeonCycleFreeFBA(model, objValFrac / 100 * valNom, false);
%     [ solLeon2, valLeon2, nActiveRxn2, stat ] = deLeonCycleFreeFBA(model, [], false);
    [ solLeon1, valLeon1, ~, ~ ] = deLeonCycleFreeFBA(model, objValFrac / 100 * valNom, false);
    [ solLeon2, valLeon2, ~, ~ ] = deLeonCycleFreeFBA(model, [], false);

    cobraMod = convertToCOBRAmodel(model);
    
    lsgSchell = optimizeCbModel(cobraMod, 'max', 0, false);
    lsgSchell2 = optimizeCbModel(cobraMod, 'max', 'one', false);
    
%     sols = [valNom, valNorm, lsgSchell.f, lsgSchell2.f; ...
%         solNom, solNorm, lsgSchell.x, lsgSchell2.x];
    sols = [valStart, valNom, valNorm, valLeon1, valLeon2, lsgSchell.f, lsgSchell2.f; ...
        solStart, solNom, solNorm, solLeon1, solLeon2, lsgSchell.x, lsgSchell2.x];
    
end

