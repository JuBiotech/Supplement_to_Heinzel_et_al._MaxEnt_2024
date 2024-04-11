function plotTradeOffCurve(epsVals, eigenLb, model, groups, includeConcentration)
%PLOTTRADEOFFCURVE Plots a single trade-off curve varying the trade-off
%parameter epsilon and plotting objective function value and worst-case
%stoichiometric error lambda against epsilon.
%
% Parameters:
%   - epsVals: Trade-off parameter (epsilon) values.
%   - eigenLb: Lower bound on the eigenvalues of the SDP-matrix.
%   - model: Model structure. See mcGroups() for hints.
%   - groups: Cell array of groups. See mcGroups() for hints.
%   - includeConcentration: Optional. If set to true uses
%       roSolveModelWithConcentrations() to solve the RO problems. Defaults
%       to false.

    objVals = zeros(length(epsVals), 1);
    lambdaVals = zeros(length(epsVals), 1);

    solver = @roSolveModel;
    if (nargin >= 5) && ~isempty(includeConcentration) && includeConcentration
        solver = @roSolveModelWithConcentrations;
    end
    
    for i = 1:length(epsVals)

        [ objVals(i), blubb, lambdaVals(i) ] = solver( epsVals(i), eigenLb, model, groups );

    end

    semilogy(epsVals, objVals);
    hold on;
    grid on;
    semilogy(epsVals, lambdaVals, 'r');
    hold off;

    legend('Biomass', 'Lambda');

end