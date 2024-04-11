function [ res ] = rangeGroupRelative( model, group, delta, steps )
%RANGEGROUPRELATIVE Samples a given parameter range and records the
%results.
% The range is linearly sampled with the given amount of intermediate steps.
% At each step the current value of the range is taken and (1 + value) is
% multiplied onto the stoichiometric coefficients specified by the group.
% Then a FBA is conducted and the results are recorded.
%
% Example:
%	delta = 0.1, steps = 10
%		Factors: 0.9, 0.92, ... 0.98, 1.0, 1.02, ..., 1.08, 1.1
%		So there are actually 11 samples created. The 6th sample is equal
%		to the unperturbed problem.
%
% Parameters:
%	- model: Model structure
%	- group: row and column of coefficients in stoichiometric matrix.
%		Structure with fields row and column. Fields may be vectors.
%	- delta: Relative deviation in one direction. 0.1 means that the
%		coefficient will be varied from -10% to 10% of its value.
%	- steps: Number of steps
%
% Returns:
%	Structure with fields:
%	- diffs: Distance of original solution to sampled solution in Euclidian 
%		metric.
%	- vals: Objective function values
%	- solutions: Matrix in which each column is a solution to the FBA

	S = model.S;
	
	diffs = zeros(steps+1, 1);
	sol = zeros(length(model.lb), steps+1);
	val = zeros(steps+1, 1);
	
	origSol = solveLPProblem(-1, model.c, [], [], S, model.b, model.lb, model.ub);
	c = 2*delta / steps;

	for i = 1:steps+1
		S_neu = S;
		S_neu(group.row, group.column) = (1-delta + (i-1)*c) * S_neu(group.row, group.column);
		
		[sol(:, i), val(i)] = solveLPProblem(-1, model.c, [], [], S_neu, model.b, model.lb, model.ub);
		diffs(i) = norm(origSol - sol(:, i));
	end

	res.diffs = diffs;
	res.vals = val;
	res.solutions = sol;
end

