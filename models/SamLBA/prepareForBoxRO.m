function [f, iqA, iqB, eqA, eqB, lower, upper, solIndices] = ...
    prepareForBoxRO(U, f, S, reactions, pools, cons, consB, lower, upper )
%PREPAREFORBBOXRO prepares the given LP for robust optimization using
%box uncertainty.
%
% Parameters:
%   - U: Uncertainty structure (see parseUncertainty.m for help)
%   - f: Objective function as string or vector
%   - S: Stoichiometric matrix
%   - reactions: Cell array of reactions
%   - pools: Cell array of pools
%   - cons: Constraint matrix
%   - consB: Constraint right hand side
%   - lower: Lower bound of fluxes
%   - upper: Upper bound of fluxes
%
% Returns:
%   - f: Objective function
%   - iqA: Inequality system matrix
%   - iqB: Inequality right hand side
%   - eqA: Equality system matrix
%   - eqB: Equality right hand side
%   - lower: Lower flux bounds
%   - upper: Upper flux bounds
%   - solIndices: Indices of the solution vector containing the actual
%       solution (ignores auxiliary variables)

    iqA = [];
    iqB = [];
    mapping = [];

    if isfield(U, 'stoichiometry') && ~isempty(U.stoichiometry)
        % Turn equations S * x = 0 into inequalities S * x <= 0, -S*x <= 0
        rowUncertaintyS = summarizeRowUncertainty(U.stoichiometry, length(pools), @(sto) getIndexOfElement(pools, sto.pool));
        [eqA, eqB, iqA, iqB, mapping] = processUncertainStoichiometry(rowUncertaintyS, S, mapping);
    else
        eqA = S;
        eqB = zeros(size(S, 1), 1);

        rowUncertaintyS = [];
    end

    % Treat constraints
    nIqsFromS = size(iqA, 1);
    if isfield(U, 'constreq') && ~isempty(U.constreq)
        rowUncertaintyC = summarizeRowUncertainty(U.constreq, length(consB), @(sto) sto.index);
        mapping = processUncertainConstraintsMatrix(rowUncertaintyC, mapping, nIqsFromS);
    else
        rowUncertaintyC = cell(length(consB), 1);
    end
    iqA = [iqA; cons];

    if isfield(U, 'constrrhs') && ~isempty(U.constrrhs)
        [mapping, rowUncertaintyC] = processUncertainConstraintsRHS(U.constrrhs, rowUncertaintyC, mapping, nIqsFromS);
    end
    iqB = [iqB; consB];
    solIndices = 1:length(f);

    % Process objective function
    if isfield(U, 'target') && ~isempty(U.target)
        [eqA, eqB, iqA, iqB, lower, upper, mapping] = processUncertainObjective(f, eqA, eqB, iqA, iqB, lower, upper, mapping);
        f = [zeros(size(eqA, 2) - 1, 1); 1];
    else
        U.target = [];
    end
    
    if isempty(mapping)
		lower = [lower; -inf * ones(size(iqA, 2) - length(lower), 1)];
		upper = [upper; inf * ones(size(iqA, 2) - length(upper), 1)];
        return;
    end
    
    [iqA, iqB, eqA, eqB, f] = immunizeInequalities(iqA, iqB, eqA, eqB, f, mapping, rowUncertaintyS, rowUncertaintyC, U.target, reactions);
	lower = [lower; -inf * ones(size(iqA, 2) - length(lower), 1)];
	upper = [upper; inf * ones(size(iqA, 2) - length(upper), 1)];
end

function [iqA, iqB, eqA, eqB, f] = immunizeInequalities(iqA, iqB, eqA, eqB, f, mapping, rowUncertaintyS, rowUncertaintyC, uncertainObjective, reactions)
%IMMUNIZEINEQUALITIES immunizes all inequalities against uncertainty.
%
% Parameters:
%   - iqA: Inequality system matrix
%   - iqB: Inequality right hand side
%   - eqA: Equality system matrix
%   - eqB: Equality right hand side
%   - f: Objective function
%   - mapping: Mapping matrix in which each row represents a mapping of the
%       inequality system to the equation or inequality which produced the
%       inequality. First column is the index of the sub-system and second
%       column is the type (1 = stoichiometry, 2 = constraints, 3 =
%       objective)
%   - rowUncertaintyS: Cell array with the uncertainty structures of each
%       row in the equation system. Empty cell means, the corresponding row
%       of the stoichiometric matrix is certain
%   - rowUncertaintyC: Same as rowUncertaintyS for constraint system
%   - uncertainObjective: Cell array with objective function uncertainty
%   - reactions: Cell array with reaction names
%
% Returns:
%   - iqA: Generated inequality system
%   - iqB: Generated inequality right hand side
%   - eqA: Generated equation system
%   - eqB: Generated equation right hand side
%   - f: Augmented objective function

    origSize = length(f);
    for i = 1:size(mapping, 1)
        % Read type
        if (mapping(i, 2) == 1) || (mapping(i, 2) == 2)
            if mapping(i, 2) == 1
                % Stoichiometry
                uList = rowUncertaintyS{mapping(i, 1)};
            else
                % Constraints
                uList = rowUncertaintyC{mapping(i, 1)};
            end
            
            [iqA, iqB, eqA, eqB, f] = immunizeInequality(iqA, iqB, eqA, eqB, i, f, uList, origSize, @(sto) getIndexOfElement(reactions, sto.reaction));
        elseif mapping(i, 2) == 3
            [iqA, iqB, eqA, eqB, f] = immunizeInequality(iqA, iqB, eqA, eqB, i, f, uncertainObjective, origSize, @(obj) obj.index);
        end
    end
    
    % Delete immunized rows
    inds = mapping(:, 1) ~= 0;
    iqA(inds, :) = [];
    iqB(inds) = [];
end

function [iqA, iqB, eqA, eqB, f] = immunizeInequality(iqA, iqB, eqA, eqB, iqRow, f, uncertainty, origSize, indexSelector)
%IMMUNIZEINEQUALITY immunizes a single inequality against uncertainty.
%
% Parameters:
%   - iqA: Inequality system matrix
%   - iqB: Inequality right hand side
%   - eqA: Equality system matrix
%   - eqB: Equality right hand side
%   - iqRow: Row index in the inequality system matrix
%   - f: Objective function
%   - uncertainty: Cell array with uncertainty structures of the given row
%   - origSize: Original count of the solution variables without
%       augmentations
%   - indexSelector: Function converting an uncertainty structure from the
%       given cell array to a column index in the uncertainty structure
%
% Returns:
%   - iqA: Generated inequality system
%   - iqB: Generated inequality right hand side
%   - eqA: Generated equation system
%   - eqB: Generated equation right hand side
%   - f: Augmented objective function
%
% Remarks:
%   The variables of the systems will be augmented by [z, w, u, v] where u
%   and v are auxiliary variables for abs(z) and abs(w) respectively.
%   abs(x) is transformed into -x <= u and x <= u, thus abs(x) = u.
%   max abs(w_i) <= b is transformed into v_i <= b for all i 
%   (Note: abs(w_i) = v_i). This is expanded into n inequalities like:
%       foo + v_1 <= b, foo + v_2 <= b, ..., foo + v_n <= b
%   
%   Augmentation: +n variables, +2n+1 inequalities
%       where n is the count of uncertain coefficients in the inequality

    n = length(uncertainty);
    origCoeff = iqA(iqRow, 1:origSize);
    origRHS = iqB(iqRow);
    
    rightHandSide = zeros(n, origSize);
    for j = 1:n
        if isstruct(uncertainty)
            element = uncertainty;
        else
            element = uncertainty{j};
        end
                
        ind = indexSelector(element);
        rightHandSide(j, ind) = abs(diff(element.bounds)) / 2;
    end

    % Enhance systems
    % Insert v - u <= 0, -v - u <= 0 and s^T * v + sum u <= 0
    if (size(eqA, 2) == origSize) || (size(iqA, 2) == origSize)
        iqA = [iqA, zeros(size(iqA, 1), n); ...
            rightHandSide, -eye(n); ...
            -rightHandSide, -eye(n); ...
            origCoeff, ones(1, n)];
    else
        padding = zeros(n, size(iqA, 2) - origSize);
        iqA = [iqA, zeros(size(iqA, 1), n); ...
            rightHandSide, padding, -eye(n); ...
            -rightHandSide, padding, -eye(n); ...
            origCoeff, padding, ones(1, n)];
    end
        
    iqB = [iqB; zeros(2*n, 1); origRHS];
        
    f = [f; zeros(n, 1)];
    eqA = [eqA, zeros(size(eqA, 1), n)];
end

function [eqA, eqB, iqA, iqB, lower, upper, mapping] = processUncertainObjective(f, eqA, eqB, iqA, iqB, lower, upper, mapping)
%PROCESSUNCERTAINOBJECTIVE processes uncertainty in the objective function
%and converts the objective function to an inequality augmenting the whole
%problem by a auxiliary variable (if necessary).
%
% Parameters:
%   - f: Objective function vector
%   - eqA: Equation system matrix
%   - eqB: Equation system right hand side
%   - iqA: Inequality sytem matrix
%   - iqB: Inequality system right hand side
%   - lower: Lower bound of fluxes
%   - upper: Upper bound of fluxes
%   - mapping: Mapping matrix in which each row corresponds to the same row
%       in the global inequality system. First column gives index of the
%       uncertain subsystem matrix row responsible for creation of the
%       inequality. Second column is the subsystem type (1 = Stoichiometry,
%       2 = Constraints, 3 = Objective)
%
% Returns:
%   - eqA: Edited equation system matrix
%   - eqB: Edited equation system right hand side
%   - iqA: Edited inequality system matrix
%   - iqB: Edited inequality system right hand side
%   - lower: Edited lower bound of fluxes
%   - upper: Edited upper bound of fluxes
%   - mapping: Edited mapping matrix

    % Augment all equations by parameter t
    eqA = [eqA, zeros(size(eqA, 1), 1)];
    iqA = [iqA, zeros(size(iqA, 1), 1)];
    lower(end + 1) = 0;
    upper(end + 1) = inf;

    % Add objective function as last row of inequalities
    iqA = [iqA; -f', 1];
    iqB = [iqB; 0];
    
    mapping = [mapping; length(iqB), 3];
end

function [mapping] = processUncertainConstraintsMatrix(constreq, mapping, offset)
%PROCESSUNCERTAINTYCONSTRAINTSMATRIX processes uncertainty in the
%constraint inequality system matrix.
%
% Parameters: 
%   - constreq: Cell array of uncertainty structures in which each row
%       corresponds to a row in the underlying constraint matrix
%   - mapping: Mapping matrix mapping each row in the global inequality
%       system to its generation cause.
%   - offset: Index offset, number of inequalities generated by previous
%       functions.
%
% Returns:
%   - mapping: Edited mapping

    for i = 1:length(constreq)
        if isempty(constreq{i})
            % Skip empty lines
            continue;
        end
        
        index = offset + i;
        
        if size(mapping, 1) < index
            mapping = [mapping; zeros(index - size(mapping, 1), 2)];
        end
        mapping(index, :) = [i, 2];
    end
end

function [mapping, rowUncertaintyC] = processUncertainConstraintsRHS(constrrhs, rowUncertaintyC, mapping, offset)
%PROCESSUNCERTAINCONSTRAINTS processes uncertainty in the constraint right
%hand side, taking the minimum of the uncertatinty bounds.
%
% Parameters:
%   - constrrhs: Uncertainties
%   - consB: Constraint right hand side
%   - mapping: Mapping matrix mapping each row in the global inequality
%       system to its generation cause.
%   - offset: Index offset, number of inequalities generated by previous
%       functions (Attention: Should be the same as in
%       processUncertainConstraintsMatrix).
%
% Returns: Updated mapping and edited constraint right hand side

    for i = 1:length(constrrhs)
        cr = constrrhs{i};
        constrrhs{i}.coefftype = 'offset';
        rowUncertaintyC{cr.index} = [rowUncertaintyC{cr.index}, constrrhs{i}];
        if size(mapping, 1) < cr.index + offset
            mapping = [mapping; zeros(cr.index + offset - size(mapping, 1), 2)];
        end
        mapping(cr.index + offset, :) = [cr.index, 2]; 
    end
end

function [eqA, eqB, ineqA, ineqB, mapping] = processUncertainStoichiometry(uncertain, S, mapping)
%PROCESSUNCERTAINSTOICHIOMETRY processes uncertainty in the stoichiometry.
%This function transforms every uncertain row in the equation system into
%two inequalities.
%
% Parameters:
%   - uncertain: Summarized uncertainty cell array from
%       summarizeRowUncertainty()
%   - S: Stoichiometric matrix
%   - mapping: Mapping matrix mapping each row of the inequality structure
%       to its generation cause (uncertain stoichiometric row, uncertain
%       constraint row)
%
% Returns:
%   - eqA: Equation system matrix
%   - eqB: Equation right hand side
%   - ineqA: Inequality system created by uncertain rows
%   - ineqB: Inequality right hand side
%   - mapping: Edited mapping matrix

    ineqA = [];
    ineqB = [];
    indChanged = [];
    for i = 1:length(uncertain)
        if isempty(uncertain{i})
            % Skip empty lines
            continue;
        end
        
        sto = uncertain{i};
		if iscell(sto)
			sto = sto{1};
		end
        
        ineqA = [ineqA; S(i, :); -S(i, :)];
        ineqB = [ineqB; sto.rhsTolerance(2); -sto.rhsTolerance(1)];
        
        indChanged = [indChanged; i];
        mapping = [mapping; i, 1; i, 1];
    end

    % Add equations
    indsEq = 1:size(S, 1);
    indsEq(indChanged) = [];
    eqA = S(indsEq, :);
    eqB = zeros(length(indsEq), 1);
end
