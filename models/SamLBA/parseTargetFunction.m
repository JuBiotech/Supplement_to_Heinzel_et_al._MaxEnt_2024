function t = parseTargetFunction(target, names)
%PARSETARGETFUNCTION parses the given target function using the given
%variables.
%The output is a vector with coefficients of the given variables (order is 
%relevant).
%
% String format: '12.7v - 2e-1 u_3' using the variable names v, u_3.
%
% Parameters:
%   - target: Constraint string.
%   - names: Cell array with variable names
%
% Returns:
%   Vector with the given constants in the order given by names.

    t = zeros(length(names), 1);
    indices = 1:length(names);
    target = [target '!'];

    state = 1;
    inStart = -1;
    count = 0;
    num = NaN;
    i = 1;
    sign = 1;
    
    while i <= length(target)
        if state == 2
            if (target(i-1) == 'e' || target(i-1) == 'E') && (target(i) == '-')
                % Skip over the following sign in the exponent
                i = i+1;
                continue;
            end
            if target(i) == '.'
                % Skip over dot
                i = i+1;
                continue;
            end
            if (~isnumber(target(i)) && (target(i) ~= 'e') && ...
                    (target(i) ~= 'E')) || (((target(i) == 'e') || ...
                    (target(i) == 'E')) && count == 1)
                
                if (target(i-1) == 'e') || (target(i-1) == 'E')
                    % Rewind to the 'e'
                    i = i-1;
                end
                
                % Number end - read it and return to standard state (1)
                num = sign * str2double(target(inStart:(i-1)));
                state = 1;
                                                
                % Do not increment the counter as we do not want to
                % neglect the just read character
                count = 0;
                continue;
            end
            if (target(i) == 'e') || (target(i) == 'E')
                % Increment e-counter
                count = 1;
            end
            % Consume number
        end
        if state == 3
            if ~isalpha(target(i)) && ~isnumber(target(i))
                % Set coefficient in target function
                var = target(inStart:(i-1));
                index = indices(strcmp(var, names));
                
                if isnan([num])
                    num = sign;
                end
                t(index) = num;
                num = NaN;
                
                state = 1;
                continue;
            end
        end
        
        if state == 1
            if (target(i) == '*') || isspace(target(i))
                % Ignore * character and whitespace
                i = i+1;
                continue;
            end
            
            inStart = i;
            
            if (target(i) == '+') || (target(i) == '-') || isnumber(target(i))
                % Detected number - switch to number reading state (2)
                state = 2;
                sign = 1;
                if (target(i) == '-')
                    sign = -1;
                    inStart = i+1;
                end
            end
            if isalpha(target(i))
                % Detected string character - switch to variable reading
                % state (3)
                state = 3;
            end
        end
        
        % Move along the string
        i = i+1;
    end
end

function y = isalpha(x)
    y = all( (x>='a' & x<='z') | (x>='A' & x<='Z') | (x=='_') );
end

function y = isnumber(x)
    y = all( (x>='0' & x<='9') );
end
