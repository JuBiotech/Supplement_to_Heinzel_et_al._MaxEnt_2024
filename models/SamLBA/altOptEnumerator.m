function [ altOpt, rawOpt ] = altOptEnumerator( model, useMosek, tol, tolNonZero, varUbConst )
%ALTOPTENUMERATOR Enumerates all alternative solutions of an FBA problem.
%
% Uses MILP to enumerate all vertices of the optimal face of the flux
% polyhedron.
%
% Reference:
% Recursive MILP model for finding all the alternate optima in LP models for metabolic networks
% S. Lee et al., Computers and Chemical Engineering 24 (2000)
%
% Notice: MOSEK may produce wrong or misleading resutls. Better check with
% POLCO and / or GLPK! GLPK seems to be reliable.
%
% Parameters:
%   - model: Model structure.
%   - useMosek: Optional. Set to true to use the MOSEK optimizer. If set to
%       false, GLPK will be used. Default is false (GLPK).
%   - tol: Optional. Tolerance parameter for the abort condition (objective
%       function value is smaller than the nominal value - tolerance).
%       Defaults to 1e-10. GLPK is sensitive to this parameter.
%   - tolNonZero: Optional. Tolerance for detecting non-zero fluxes.
%       Defaults to 1e-10.
%   - varUbConst: Optional. Global upper bound on the fluxes in the
%       transformed model. Defaults to abs(min(lower)) + abs(max(upperb)) + 1
%       MOSEK is sensitive to this parameter.
%
% Returns:
%   - altOpt: Matrix in which each column is a vertex of the optimal face
%       (i.e. an alternative solution)
%   - rawOpt: Matrix in which each column is a vertex of the optimal face
%       in the transformed model


    % Check constraints
    if ~isfield(model, 'cM')
        model.cM = [];
    end
    if ~isfield(model, 'cB')
        model.cB = [];
    end

    if size(model.c, 2) > 1
        model.c = model.c';
    end
    if size(model.lb, 2) > 1
        model.lb = model.lb';
    end
    if size(model.ub, 2) > 1
        model.ub = model.ub';
    end
   
    if (nargin < 2) || isempty(useMosek)
        useMosek = false;
    end
    
    if (nargin < 3) || isempty(tol)
        tol = 1e-10; % <-- GLPK is sensitive, so you will need to play with this.
    end
    
    if (nargin < 4) || isempty(tolNonZero)
        tolNonZero = 1e-10;
    end
    
    if (nargin < 5) || isempty(varUbConst)
        varUbConst = abs(min(model.lb)) + abs(max(model.ub)) + 1; % <-- You will also need to play with this (especially when using MOSEK).

        % Fun fact: varUbConst influences the number of found alternate optimal
        % solutions. Especially MOSEK is influenced by this parameter.
    end
    
    if useMosek
        lpSolver = @solveMosekLP;
        milpSolver = @solveMosekMILP;
    else
        lpSolver = @solveGlpkLP;
        milpSolver = @solveGlpkMILP;
    end
    
    % Step 1: Transform LP to standard form
    model = transformToStandardFormPaper(model);
    
    % Step 2: Iterate over alternate solutions
    [sol, origValTrans] = lpSolver(model);

    NZ = [{find(sol > tolNonZero)}];
    K = 1;
    
    rawOpt = sol;
    altOpt = [transformSolutionFromStandardFormPaper(model, sol)];
    
    % Add variables y, w
    nVars = length(model.c);
    model.S = [model.S, zeros(size(model.S, 1), 2*nVars)];
    model.c = [model.c; zeros(2 * nVars, 1)];
    model.ub = [2000.*ones(nVars, 1); ones(2*nVars, 1)];
    model.lb = [model.lb; zeros(2*nVars, 1)];

    tHandle = tic;
    
    while (true)
        
        K = K + 1;
        
        lastNZ = NZ{K-1};
        
        ySum = zeros(1, nVars);
        ySum(lastNZ) = 1;
        
        yMatrix = zeros(length(lastNZ), nVars);
        for i = 1:length(lastNZ)
            yMatrix(i, lastNZ(i)) = 1;
        end
        
        % Add constraints
        Ain = [zeros(1, nVars), -ySum, zeros(1, nVars); ...  % Sum y_i >= 1
            eye(nVars), zeros(nVars), -diag(varUbConst .* ones(nVars, 1)); ... % z_i <= ub * w_i
            zeros(nVars), eye(nVars), eye(nVars)]; % y_i + w_i <= 1

        bin = [-1; zeros(nVars, 1); ones(nVars, 1)];

        for i = 1:K-1
            wSum = zeros(1, nVars);
            wSum(NZ{i}) = 1;
            Ain = [Ain; zeros(1, 2*nVars), wSum];
            bin = [bin ; length(NZ{i})-1];
        end
        
        
        % Now solve
        [sol,val,solverStat] = milpSolver(model, Ain, bin, nVars+1:3*nVars);
        if (~solverStat) || (val < origValTrans-tol)
            break;
        end
        
%        fprintf('Original %8.6f  vs. %8.6f  = Diff %8.6f\n', origValTrans, val, abs(origValTrans - val));
        
        NZ{K} = find(sol(1:nVars) > tolNonZero);
        altOpt = [altOpt, transformSolutionFromStandardFormPaper(model, sol(1:nVars))];
        rawOpt = [rawOpt, sol(1:nVars)];
    end
    
    fprintf('Elapsed time: %g sec \n', toc(tHandle));
    fprintf('Found %g alternate solutions with objective function value %g \n', size(altOpt, 2), origValTrans);
end

function model = transformToStandardFormPaper(model)

    % Eliminate fixed variables
    inds = find(model.lb == model.ub);
    model.fixedInds = inds;
    model.nonFixedInds = 1:length(model.c);
    model.nonFixedInds(inds) = [];
    
    model.fixedVarVals = model.lb(inds);
    if isfield(model, 'S') && ~isempty(model.S)
        model.b = model.b - model.S(:, inds) * model.fixedVarVals;
        model.S(:, inds) = [];
    end
    
    if isfield(model, 'cM') && ~isempty(model.cM)
        model.cB = model.cB - model.cM(:, inds) * model.fixedVarVals;
        model.cM(:, inds) = [];
    end
    
    model.c(inds) = [];
    model.lb(inds) = [];
    model.ub(inds) = [];


% Final variable layout: s_L, s_U, s
    nVars = length(model.c);

    % Introduce slack variables for inequalities
    if ~isempty(model.cM)
        nSlack = size(model.cM, 1);

        if ~isempty(model.S)
            model.b = [model.b - model.S * model.lb; model.cB - model.cM * model.lb; model.ub - model.lb];
            model.S = [model.S, zeros(size(model.S, 1), nVars + nSlack); ...
                model.cM, zeros(size(model.cM, 1), nVars) eye(nSlack); ...
                eye(nVars) eye(nVars) zeros(nVars, nSlack)];
        else
            model.b = [model.cB - model.cM * model.lb; model.ub - model.lb];
            model.S = [model.cM, zeros(size(model.cM, 1), nVars) eye(nSlack); ...
                eye(nVars) eye(nVars) zeros(nVars, nSlack)];
        end
        
        model.c = [model.c ; zeros(nVars + nSlack, 1)];

        model.origLb = model.lb;
        model.nSlack = nSlack;
        model.nVars = nVars;
        model = rmfield(model, 'cM');
        model = rmfield(model, 'cB');        

        model.lb = zeros(nVars * 2 + nSlack, 1);
    else
        model.b = [model.b - model.S * model.lb; model.ub - model.lb];

        model.S = [model.S, zeros(size(model.S, 1), nVars); ...
            eye(nVars) eye(nVars)];

        model.c = [model.c ; zeros(nVars, 1)];

        model.origLb = model.lb;
        model.nVars = nVars;

        model.lb = zeros(nVars * 2, 1);
    end
    
    model.varUb = model.ub;
    model.ub = [];
end

function transSol = transformSolutionFromStandardFormPaper(model, sol)
    transSol = zeros(length(model.fixedVarVals) + length(model.nonFixedInds), 1);
    transSol(model.nonFixedInds) = sol(1:model.nVars) + model.origLb;
    transSol(model.fixedInds) = model.fixedVarVals;
end


function model = transformToStandardForm(model)

% Final variable layout: normal vars, slack vars, negative split vars

% varType:
%   - v: Normal variable
%   - s: Slack variable
%   - n: Negative split variable
%   - r: Sign reversed

% eqType:
%   - e: Equation (S matrix)
%   - l: Inequality (<=, constraint Matrix cM)

    model.varType = repmat('v', length(model.c), 1);
    model.eqType = repmat('e', length(model.b), 1);
    model.refEq = zeros(length(model.b), 1);
    model.refVar = zeros(length(model.c), 1);

    % Introduce slack variables for inequalities
    if ~isempty(model.cM)
        nSlacks = size(model.cM, 1);
        model.S = [model.S zeros(size(model.S, 1), nSlacks); model.cM eye(nSlacks) ];
        model.b = [model.b ; model.cB];
        
        model.c = [model.c ; zeros(nSlacks, 1)];
        model.lb = [model.lb; zeros(nSlacks, 1)];
%        model.ub = [model.ub; inf(nSlacks, 1)];
        model.ub = [model.ub; 1e3 .* ones(nSlacks, 1)];
        
        model = rmfield(model, 'cM');
        model = rmfield(model, 'cB');
        
        model.varType = [model.varType; repmat('s', nSlacks, 1)];
        model.eqType = [model.eqType; repmat('l', nSlacks, 1)];
        model.refEq = [model.refEq; zeros(nSlacks, 1)];
        model.refVar = [model.refVar; zeros(nSlacks, 1)];
    end

    % Transform all variables to non-negative (split)

    % a) Variable has pure negative range (reverse sign)
    negInds = find((model.lb < 0) & ( model.ub <= 0));
    for i = 1:length(negInds)
        ni = negInds(i);
            
        model.S(:, ni) = -model.S(:, ni);

        u = model.ub(ni);
        model.ub(ni) = -model.lb(ni);
        model.lb(ni) = -u;

        model.c(ni) = -model.c(ni);

        model.varType(ni) = 'r';
    end
    
    % b) Add negative variable (Split)
    negInds = find((model.lb < 0) & ( model.ub > 0));
    for i = 1:length(negInds)
        ni = negInds(i);

        model.S = [model.S, -model.S(:, ni) ];

        model.ub = [model.ub; -model.lb(ni)];
        model.lb = [model.lb; 0];
        model.lb(ni) = 0;

        model.c = [model.c; -model.c(ni)];

        model.varType = [model.varType; 'n'];
        model.refVar = [model.refVar; ni];
    end
    
    % Fix offsets
    offInds = find(model.lb > 0);
    model.varOffset = zeros(length(model.c), 1);
    for i = 1:length(offInds)
        oi = offInds(i);
        offset = model.lb(oi);
        model.lb(oi) = 0;
        model.ub(oi) = model.ub(oi) - offset;
        model.varOffset(oi) = offset;
        model.b = model.b - model.S(:, oi) .* offset;
    end
    
    model.varUb = model.ub;
end

function transSol = transformSolutionFromStandardForm(model, sol)
    lastInd = find(model.varType == 'v', 1, 'last' );
    revInds = find(model.varType == 'r');
    firstInd = find(model.varType == 'n', 1, 'first' );
    
    sol = sol + model.varOffset;
    
    transSol = sol(1:lastInd);
    transSol(revInds) = -transSol(revInds);
    refs = model.refVar(firstInd:end);
    transSol(refs) = transSol(refs) - sol(firstInd:end);
end

function [sol, val] = solveGlpkLP(model)
    params.minNorm = true;
    params.dual = 2;
    params.lpsolver = 3;
    [sol, val] = solveLPProblem(-1, model.c, [], [], model.S, model.b, model.lb, model.ub, params);
end

function [sol, val] = solveMosekLP(model)
%MSK_OPTIMIZER_FREE_SIMPLEX    

    prob.c        = model.c;
    prob.a        = sparse(model.S);
    prob.blc      = model.b;
    prob.buc      = model.b;
    prob.blx      = model.lb;
    prob.bux      = model.ub;

    % Specify parameters
    params.MSK_IPAR_OPTIMIZER = 7; %'MSK_OPTIMIZER_FREE_SIMPLEX';
    params.MSK_IPAR_LOG = 0;
    
    % Optimize the problem.
    [r,res] = mosekopt('maximize', prob, params);

    try 
      % Display the optimal solution.
      val = res.sol.bas.pobjval;
      sol = res.sol.bas.xx;
      stat = (r == 0);
    catch
      fprintf('MSKERROR: Could not get solution')
      stat = false;
    end
    
%    [sol, val] = linprog(-model.c, [], [], model.S, model.b, model.lb, model.ub);
%    val = -val;
end

function [sol, val, stat] = solveGlpkMILP(model, Ain, bin, intVars)
    % Now solve
    bigA = [model.S; Ain];
    bigB = [model.b; bin];
    bigB = full(bigB);
    csense = [repmat('S', size(model.S, 1), 1); repmat('U', size(Ain, 1), 1)];
    varType = repmat('C', size(bigA, 2), 1);
    varType(intVars) = 'B';

    glpkParams.dual = 2;
    glpkParams.lpsolver = 1; % 3 = Exact arithmetic seems not to work with MILP
%    glpkParams.msglev = 3;
    
    [sol,val,solverStat] = glpk(model.c, bigA, bigB, model.lb, model.ub, csense, varType, -1, glpkParams);
    stat = (solverStat == 5);
end

function [sol, val, stat] = solveMosekMILP(model, Ain, bin, intVars)
    prob.c        = model.c;
    prob.a        = sparse([model.S; Ain]);
    prob.blc      = [model.b; -inf(length(bin), 1)];
    prob.buc      = [model.b; bin];
    prob.blx      = model.lb;
    prob.bux      = model.ub;

    % Specify indexes of variables that are integer
    % constrained.

    prob.ints.sub = intVars;

    % Specify parameters
    params.MSK_IPAR_OPTIMIZER = 7; %'MSK_OPTIMIZER_FREE_SIMPLEX';
    params.MSK_IPAR_MIO_ROOT_OPTIMIZER = 7; %'MSK_OPTIMIZER_FREE_SIMPLEX';
%    params.MSK_IPAR_INFEAS_REPORT_AUTO = 1;
%    params.MSK_IPAR_INFEAS_REPORT_LEVEL = 100;
    params.MSK_IPAR_LOG = 0;
    
    % Optimize the problem.
    [r,res] = mosekopt('maximize', prob, params);

    try 
      % Display the optimal solution.
      val = res.sol.int.pobjval;
      sol = res.sol.int.xx;
      stat = (r == 0);
    catch
      fprintf('MSKERROR: Could not get solution')
      stat = false;
    end
end
