function [C,b] = parseConstraints(const, names)
%PARSECONSTRAINTS parses the given constraints using the given variables.
%The contstraints are given as string with ';' terminating a constraint.
%The output is a system C*x <= b equivalent to the given constraints.
%
% String format: '12.7v - 2e-1 u_3 >= 3; 2x >= 3' using the variable names 
%   v, u_3, x.
%
% Parameters:
%   - const: Constraint string. Or structure with the fields C, b.
%   - names: Cell array with variable names
%
% Returns:
%   System C*x <= b equivalent to the given constraints

    C = zeros(1, length(names));
    b = zeros(1,1);
    
    if isempty(const)
        C = [];
        b = [];
        return;
    end
    
    if isstruct(const) && isfield(const, 'C') && isfield(const, 'b')
        C = const.C;
        b = const.b;
        return;
    end
    
    % Split at ;
    cos = regexp(const, ';', 'split');
    for i = 1:length(cos)
        [temp, bAdd] = parseSingleConstraint(cos{i}, names);
        b(i:i+size(temp,1)-1,1) = bAdd;
        C(i:i+size(temp,1)-1,:) = temp;
    end
    
end

function [C, b] = parseSingleConstraint(target, names)
% Parses a single constraint using a finite state machine
    
    cos = regexp(target, '=', 'split');
    [leftC, leftB] = parseConstraintSide(cos{1}, names);
    [rightC, rightB] = parseConstraintSide(cos{2}, names);
    
    if (cos{1}(end) == '<') || (cos{1}(end) == '>')
        if (cos{1}(end) == '>')
            b = -rightB + leftB;
            C = rightC - leftC;
        else
            b = -leftB + rightB;
            C = leftC - rightC;
        end
    else
        % Equality
        b = [-leftB + rightB; leftB - rightB];
        C = [leftC - rightC; -leftC + rightC];
    end
end

function [v, b] = parseConstraintSide(target, names)
% Parses a single constraint using a finite state machine
    
    atoms = [];
    sgns = [];
    cos = regexp(target, '+', 'split');
    for i = 1:length(cos)
        temp = regexp(cos{i}, '-', 'split');
        atoms = [atoms, temp];
        sgns = [sgns, 1, -ones(1, length(temp)-1)];
    end
    
    b = 0;
    v = zeros(1, length(names));
    
    for i = 1:length(atoms)
        [val, ind] = parseValueAndIndex(atoms{i}, names);
        if ind == -1
            b = b + sgns(i) * val;
        else
            v(ind) = v(ind) + sgns(i) * val;
        end
    end
end

function [val, ind] = parseValueAndIndex(str, names)
    i = length(str);
    
    while i >= 1
        if isspace(str(i)) || (str(i) == '<') || (str(i) == '>')
            i = i - 1;
            continue
        end
        break;
    end
    str = str(1:i);

    start = find(isalpha(str), 1, 'first' );
    if isempty(start)
        ind = -1;
        val = str2double(str);
        return;
    end
    
    sec = find(isalpha(str(start+1:end)), 1, 'first' );
    
    if 1 ~= sec
        start = sec;
    end
    
    ind = find(strcmpi(str(start:end),names));
    if ~all(isspace(str(1:start-1)))
        val = str2double(str(1:start-1));
    else
        val = 1;
    end
        
end

function y = isalpha(x)
    y = ( (x>='a' & x<='z') | (x>='A' & x<='Z') | (x=='_') );
end

function y = isnumber(x)
    y = ( (x>='0' & x<='9') );
end
