function [ cobraModel ] = convertToCOBRAmodel( model )
%CONVERTTOCOBRAMODEL Converts the given model to a COBRA toolbox model.
%
% Parameters:
%   - model: Model to convert. Structure with fields:
%		o S: Stoichiometric matrix
%		o c: Target function vector
%		o b: Right hand side of the equation S * v = b. 
%			(Equals zero for FBA => Steady State)
%       o cM: Constraint matrix (optional)
%       o cB: Right hand side of the inequality cM * v <= cB (optional)
%		o lb: Lower bound on the fluxes
%		o ub: Upper bound on the fluxes
%
% Returns: A valid COBRA toolbox model

    cobraModel = model;
    
    if isfield(cobraModel, 'csense')
        cobraModel = rmfield(cobraModel, 'csense');
    end
    
    % Check constraints
    if isfield(model, 'cM') && isfield(model, 'cB')
        cobraModel.S = [cobraModel.S ; cobraModel.cM];
        cobraModel.b = [cobraModel.b ; cobraModel.cB];
        
        if length(cobraModel.mets) < size(cobraModel.S, 1)
            
            i = 1;
            while(length(cobraModel.mets) < size(cobraModel.S, 1))
                cobraModel.mets = [cobraModel.mets, {['ineq' num2str(i)]}];
                cobraModel.metNames = [cobraModel.metNames; {['ineq' num2str(i)]}];
                i = i + 1;
            end
            
        end
        
        cobraModel.csense = [repmat('E', size(model.S, 1), 1); repmat('L', size(model.cM, 1), 1)];
        
        cobraModel = rmfield(cobraModel, 'cM');
        cobraModel = rmfield(cobraModel, 'cB');
    end

end

