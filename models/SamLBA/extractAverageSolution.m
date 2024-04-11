function sol = extractAverageSolution(res, funcLower, funcUpper)
%EXTRACTAVERAGESOLUTION extracts a solution from the samples by taking the average of each solution vector component
%whose objective function value lies in the specified bounds (funcLower, funcUpper).
%
% Parameters:
%	- res: Results structure from Monte-Carlo simulations. Should contain the following fields:
%		o sols: Matrix with solution vectors (each column is a solution)
%		o vals: Vector with objective function values
%
% Returns: A solution vector

	inds = (res.diffs >= funcLower) & (res.diffs <= funcUpper);
	sol = zeros(size(res.sols, 1), 1);
	for i = 1:length(sol)
		sol(i) = mean(res.sols(i, inds));
	end

end
