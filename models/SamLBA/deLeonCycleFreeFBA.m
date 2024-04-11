function [ sol, val, nActiveRxn, stat ] = deLeonCycleFreeFBA( model, minObjVal, useMosek )
%DELEONCYCLEFREEFBA Employs a method proposed by de Leon et al. 2008 which
%first determines the minimum number of active reactions to achieve a given
%minimum objective value. Then restricts the network to this number of
%active reactions and optimizes for the original objective.
%
% This function also supports another operation mode. It solves the
% original FBA problem and then minimizes the amount of active reactions
% while maintaining the optimal original objective value.
%
% Parameters:
%   - model: Model structure
%   - minObjVal: Optional. If set, uses the first (de Leon) operation mode.
%       It is then specified as pecentage (0 - 100) of the original
%       objective function value. If left out or empty, the function
%       switches to the second operation mode.
%   - useMosek: Optional. Set to true to use the MOSEK optimizers. Defaults
%       to false
%
% Returns:
%   - sol: Flux distribution
%   - val: Objective function value
%   - nActiveRxn: Number of active reactions
%   - stat: Solution status

    % Check constraints
    if ~isfield(model, 'cM')
        model.cM = [];
    end
    if ~isfield(model, 'cB')
        model.cB = [];
    end
    
    if (nargin <= 2) || isempty(useMosek)
        useMosek = false;
    end
    
    if useMosek
        lpSolver = @solveMosekLP;
        milpSolver = @solveMosekMILP;
    else
        lpSolver = @solveGlpkLP;
        milpSolver = @solveGlpkMILP;
    end
    
    % Identify exchange reactions
    indExch = [];
    for i = 1:size(model.S, 2)
        if sum(model.S(:, i) ~= 0) == 1
            indExch = [indExch; i ];
        end
    end
    
    indInt = 1:size(model.S, 2)';
    indInt(indExch) = [];
    
    nExch = length(indExch);
    nInt = length(indInt);
    
    % Enforce a minimum objective value and obtain the minimum amount of
    % active reactions to maintain this condition
    
    intMatrix = zeros(nInt, nInt + nExch);
    for i = 1:nInt
        intMatrix(i, indInt(i)) = 1;
    end
    
    if ~isempty(minObjVal)
        [sol, nActiveRxn, stat] = solveForMinimalActiveReations(model, minObjVal, intMatrix, nInt, nExch, indInt, milpSolver);

        % Now solve nominal problem with constraint on the number of active
        % reactions

        [sol, val, stat] = solveNormalProblemWithGivenActiveRxn(model, abs(nActiveRxn), intMatrix, nInt, nExch, indInt, milpSolver);
    else
        [ sol, val, status ] = quickSolveFBA(model, false);
        [sol, nActiveRxn, stat] = solveForMinimalActiveReations(model, val, intMatrix, nInt, nExch, indInt, milpSolver);
    end
    
    nActiveRxn = abs(nActiveRxn);
end

function [sol, nActiveRxn, stat] = solveForMinimalActiveReations(model, minObjVal, intMatrix, nInt, nExch, indInt, milpSolver)
    % Layout: v, y
    minReactModel.S = [model.S, zeros(size(model.S, 1), nInt)];
    minReactModel.b = model.b;
    minReactModel.lb = [model.lb; zeros(nInt, 1)];
    minReactModel.ub = [model.ub; ones(nInt, 1)];
    minReactModel.cM = [model.cM, zeros(size(model.cM, 1), nInt); ...
                        intMatrix, -diag(model.ub(indInt)); ...
                        -intMatrix, diag(model.lb(indInt)); ...
                        -model.c', zeros(1, nInt)];
    minReactModel.cB = [model.cB; zeros(2*nInt, 1); -minObjVal];
    minReactModel.c = [zeros(nExch + nInt, 1); -ones(nInt, 1)];
    
    [sol, nActiveRxn, stat] = milpSolver(minReactModel, minReactModel.cM, minReactModel.cB, (nInt+nExch+1):(nInt+nExch+nInt));
    sol = sol(1:(nInt+nExch));
end

function [sol, val, stat] = solveNormalProblemWithGivenActiveRxn(model, nActiveRxn, intMatrix, nInt, nExch, indInt, milpSolver)
    % Layout: v, y
    nomModel.S = [model.S, zeros(size(model.S, 1), nInt)];
    nomModel.b = model.b;
    nomModel.lb = [model.lb; zeros(nInt, 1)];
    nomModel.ub = [model.ub; ones(nInt, 1)];
    nomModel.cM = [model.cM, zeros(size(model.cM, 1), nInt); ...
                        intMatrix, -diag(model.ub(indInt)); ...
                        -intMatrix, diag(model.lb(indInt)); ...
                        zeros(1, nInt + nExch), ones(1, nInt)];
    nomModel.cB = [model.cB; zeros(2*nInt, 1); nActiveRxn];
    nomModel.c = [model.c; zeros(nInt, 1)];
    
    [sol, val, stat] = milpSolver(nomModel, nomModel.cM, nomModel.cB, (nInt+nExch+1):(nInt+nExch+nInt));
    sol = sol(1:(nInt+nExch));
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
