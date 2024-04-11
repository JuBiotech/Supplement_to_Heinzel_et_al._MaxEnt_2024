function [ sol, objVal, modNew ] = standardROFBA( model, mode, dpc, tol, epsilon, bmInd )
%STANDARDROFBA Constructs and solves the robust counterpart of the FBA problem
%with certain objective function and relaxed steady-state.
%  The problem: max c^T * v  s.t.  S * v = b,  cM * v <= cB,  lb <= v <= ub
%  is relaxed to the problem: max c^T * v  s.t.  b <= S * v <= tol,  cM * v <= cB,  lb <= v <= ub
%  which is then robustified using the box or budget uncertainty set.
%The problem is only robustified against uncertainty in one reaction (e.g. the biomass reaction).
%This function is not usable if you want to have more than one uncertain reaction.
%
% This function invokes solveLPProblem for the solution of the LP problems.
% To control the LP solver, set the global variable LPsolver.
%
% Parameters:
%   - model: Structure with elements c, S, b, lb, ub and optionally cM, cB
%   - mode: 'box' or 'budget'. Specifies the uncertainty set.
%   - dpc: Disturbance in percent
%   - tol: Tolerance for steady-state relaxation. Use 0 to enforce S*v = 0, i.e. no relaxation.
%           Note: The uncertain flux is limited by |v_uncertain| <= tol / (2 * dpc / 100) = tol / (50 * dpc)
%   - epsilon: Used for budget mode, determines the maximum allowed constraint violation probability (in percent).
%           Note: A value below 60.6531% is useless, because then the uncertainty set is the same as the box uncertainty set.
%   - bmInd: Index or name of biomass flux. (Optional)
%
% Returns:
%   - sol: Robust solution
%   - objVal: Robust objective function value
%   - modNew: Robust counterpart LP model with fields c, S, b, cM, cB, lb, ub

    nOriginalFluxLen = length(model.c);
    
    dpc = dpc / 100;
    
    % Check vectors
    if size(model.c, 2) > 1
        model.c = model.c';
    end
    if size(model.lb, 2) > 1
        model.lb = model.lb';
    end
    if size(model.ub, 2) > 1
        model.ub = model.ub';
    end

    % Check constraints
    if ~isfield(model, 'cM')
        model.cM = [];
    end
    if ~isfield(model, 'cB')
        model.cB = [];
    end

    % Combine systems
    A = [-model.S; model.cM];
    b = [-model.b; model.cB];

    if (nargin >= 4) && (isfinite(tol))
        A = [A; model.S];
        b = [b; ones(length(model.b), 1) * tol];
    end
    
    if (nargin < 6) || isempty(bmInd)
        bmInd = findBiomassFlux(model, bmInd);
    end

    inds = A(:, bmInd) ~= 0;
    
    % Extract certain parts
    Anew = A(~inds, :);
    bnew = b(~inds);
    
    if strcmp(mode, 'box')

        inds = find(inds == 1);
        for i = 1:length(inds)
            Anew(end+1, :) = A(inds(i), :) + [zeros(1, bmInd - 1) abs(A(inds(i), bmInd)) * dpc, zeros(1, size(A,2) - bmInd)];
            bnew(end+1) = b(inds(i));
        end

        eqM = [];
        eqB = [];
    elseif strcmp(mode, 'budget')
        
        gamma = sqrt(2 * log(100/epsilon));
        if gamma >= 1
            fprintf('Gamma = %f >= %d, choose epsilon >= %f\n', gamma, 1, exp(-1/2)*100);
        end
        
        % Layout z, w, q, u
        nUnc = sum(inds);
        nCert = size(Anew, 1);
        nFlux = size(A, 2);
        
        Anew = [Anew, zeros(nCert, 4*nUnc)];
        model.lb = [model.lb; -inf(4*nUnc, 1)];
        model.ub = [model.ub; inf(4*nUnc, 1)];
        model.c = [model.c; zeros(4*nUnc, 1)];
        
        % -q_j <= z_j <= q_j
        % -u_j <= w_j <= u_j
        Anew = [Anew; zeros(4*nUnc, nFlux), [-eye(nUnc), zeros(nUnc), -eye(nUnc), zeros(nUnc); ...
            eye(nUnc), zeros(nUnc), -eye(nUnc), zeros(nUnc); ...
            zeros(nUnc), -eye(nUnc), zeros(nUnc), -eye(nUnc); ...
            zeros(nUnc), eye(nUnc), zeros(nUnc), -eye(nUnc)]];
        bnew = [bnew; zeros(4*nUnc, 1)];
        
        % sum q_j + gamma * u_k + (a^0)^T * v <= b^0
        Anew = [Anew; A(inds, :), zeros(nUnc, 2*nUnc), eye(nUnc), gamma*eye(nUnc)];
        bnew = [bnew; b(inds)];
        
        eqM = [zeros(nUnc, nFlux-1), A(inds, bmInd) * dpc, eye(nUnc), eye(nUnc), zeros(nUnc, 2*nUnc)];
        eqB = zeros(nUnc,1);
    end
    
    % Solve
    params.minNorm = false;
    [sol, objVal, status] = solveLPProblem(-1, model.c, Anew, bnew, eqM, eqB, model.lb, model.ub, params);

    % Create model
    modNew.c = model.c;
    modNew.cM = Anew;
    modNew.cB = bnew;
    modNew.S = eqM;
    modNew.b = eqB;
    modNew.lb = model.lb;
    modNew.ub = model.ub;

    modNew.rawSol = sol;
    modNew.rawObjVal = objVal;

    sol = sol(1:nOriginalFluxLen);

    %
%    modNew.ub(isinf(modNew.ub)) = 10*max(modNew.ub(~isinf(modNew.ub)));
%    modNew.lb(isinf(modNew.lb)) = 10*min(modNew.lb(~isinf(modNew.lb)));
end
