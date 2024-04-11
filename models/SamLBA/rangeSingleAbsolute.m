function [ res ] = rangeSingleAbsolute( model, row, column, range, steps )
%RANGESINGLEABSOLUTE Samples a given parameter range and records the
%results.
% The range is linearly sampled with the given amount of intermediate steps.
% The specified coefficient is then changed to the sampled value and a FBA 
% is conducted.
%
% Example: 
%	range = [-3, 1], steps = 4
%	Values: -3, -2, -1, 0, 1
%	So actually there are 5 samples.
%
% Parameters:
%	- model: Model structure
%	- row: Row of coefficient in stoichiometric matrix
%	- column: Column of coefficient in stoichiometric matrix
%	- range: [min, max] of coefficient
%	- steps: Number of steps
%
% Returns:
%	Structure with fields:
%	- diffs: Distance of original solution to sampled solution in Euclidian 
%		metric.
%	- vals: Objective function values
%	- solutions: Matrix in which each column is a solution to the FBA
%	- coeffVals: Sampled coefficient values

	S = model.S;
	
	diffs = zeros(steps+1, 1);
	sol = zeros(length(model.lb), steps+1);
	val = zeros(steps+1, 1);
	coeffVal = zeros(steps+1, 1);
	
	origSol = solveLPProblem(-1, model.c, [], [], S, model.b, model.lb, model.ub);
	
	delta = (range(2)-range(1));
	c = delta / steps;

	for i = 1:steps+1
		S_neu = S;
		S_neu(row, column) = range(1) + (i-1)*c;
				
		[sol(:, i), val(i)] = solveLPProblem(-1, model.c, [], [], S_neu, model.b, model.lb, model.ub);
		diffs(i) = norm(origSol - sol(:, i));
		coeffVal(i) = range(1) + (i-1)*c;
	end

	res.diffs = diffs;
	res.vals = val;
	res.solutions = sol;
	res.coeffVals = coeffVal;
end

