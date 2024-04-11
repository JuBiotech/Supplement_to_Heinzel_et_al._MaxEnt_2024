function [ objVal, v, lambda ] = roSolveModel( epsilon, eigenLb, model, groups, rho, quiet )
%ROSOLVEMODEL Solves the supplied ZavRO problem using the CVX package.
%
% Reference: 
% Robust Flux Balance Analysis of Metabolic Networks
% M. Zavlanos et al., 2011 American Control Conference
%
% Parameters:
%   - epsilon: Trade-off parameter epsilon.
%   - eigenLb: Lower bound on the eigenvalues of the SDP-matrix.
%   - model: Model structure. See mcGroups() for hints.
%   - groups: Cell array of groups. See mcGroups() for hints.
%   - rho: Immunization parameter. Optional, defaults to 1.
%   - quiet: If set to true, output of the CVX package is suppressed.
%       Optional, defaults to true.
%
% Returns:
%   - objVal: Objective function value.
%   - v: Flux vector.
%   - lambda: Worst case stoichiometric error.

if (nargin <= 4) || isempty(rho)
    rho = 1;
end

if (nargin <= 5) || isempty(quiet)
    quiet = true;
end

[gS, gB, gCM, gCB, gLb, gUb, gC] = separateGroups(groups);

nFlux = size(model.S, 2);
nMets = size(model.S, 1);
nUncertain = length(gS);
matSize = 1 + nUncertain + nMets;

cvx_begin
    cvx_quiet(quiet);
    
    variable lambda(1);
    variable v(nFlux);
    variable tau(1);
    
    minimize(epsilon * lambda - (1-epsilon) * model.c' * v);

    M = [];
    for i = 1:nUncertain
        uncStruct = gS{i};
        uM = zeros(size(model.S));
        uM(uncStruct.row, uncStruct.column) = uncStruct.sigma .* model.S(uncStruct.row, uncStruct.column);
        
        M = [M, uM * v];
    end

    constraintUncertainInd = zeros(length(model.cB), 1);
    McM = [];
    for i = 1:length(gCM)
        uncStruct = gCM{i};
        uM = zeros(size(model.cM));
        uM(uncStruct.row, uncStruct.column) = uncStruct.sigma .* model.cM(uncStruct.row, uncStruct.column);
        
        McM = [McM, uM * v];
        constraintUncertainInd(uncStruct.row) = i;
    end
    
    subject to
    
        [lambda - rho.^2 * tau, zeros(1, nUncertain), (model.S * v - model.b)'; ...
            zeros(nUncertain, 1), tau * eye(nUncertain), M'; ...
            model.S*v - model.b, M, eye(nMets)] - eigenLb * eye(matSize) == semidefinite(matSize);

        tau >= 0;
        rho.^2 * tau <= lambda;
        
        v <= model.ub;
        v >= model.lb;
        
        if isfield(model, 'cM')
            
            if isempty(gCM)
                model.cM * v <= model.cB;
            else
                for k = 1:size(model.cM, 1)
                    if constraintUncertainInd(k) == 0
                        model.cM(k, :) * v <= model.cB(k);
                    else
                        ind = constraintUncertainInd(k);
                        model.cM(k, :) * v + rho * norm(McM(ind, :)) <= model.cB(k);
                    end
                end
            end
            
        end
cvx_end

objVal = model.c' * v;

end

function [gS, gB, gCM, gCB, gLb, gUb, gC] = separateGroups(groups)
    gS = {};
    gB = {};
    gCM = {};
    gCB = {};
    gLb = {};
    gUb = {};
    gC = {};
    
    for i = 1:length(groups)
        g = groups{i};
        switch lower(g.type)
            case 's'
                gS(end+1) = {g};
            case 'b'
                gB(end+1) = {g};
            case 'cm'
                gCM(end+1) = {g};
            case 'cb'
                gCB(end+1) = {g};
            case 'c'
                gC(end+1) = {g};
            case 'lb'
                gLb(end+1) = {g};
            case 'ub'
                gUb(end+1) = {g};
        end
    end    
end