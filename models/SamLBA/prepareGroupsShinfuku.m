function groups = prepareGroupsCoryneTopBM( sigmaPC )
%PREPAREGROUPSCORYNE Prepares the groups for the Coryne Glutamicum network.
% The groups created by this function preserve the correlation of the
% coefficients of the biomass reaction.
%
% Note: Only the top level biomass flux / equation will be perturbed.
%
% Correlated coefficients will be put together in one group.
% For description of the group structure see mcGroups.m .
%
% Parameters:
%	- sigmaPC: Standard deviation of the biomass coefficients in percent(!)
%
% Returns a cell array of group structures.

	%%%%
	% Biomass groups (42 cells)
	% In Coryne Glutamicum reaction 110 is biomass.
	col = 110;
	corrs = {[11, 21, 22, 29], [23, 33], [68, 81], [138, 139, 140], [147, 160], [148, 149], [158, 159]};
	uncorr = [1, 4, 5, 12, 26, 27, 31, 36, 37, 49, 61, 62, 72, 91, 93, 110, 122, 125, 137, 141, 142, 143, 144, 145, 146, 150, 151, 152, 153, 154, 155, 156, 157, 161, 162];

	sigmaPC = sigmaPC / 100;
	groups = cell(42, 1);	% Allocate memory for the groups. There are 42 groups
	
	% Save correlated
	g.column = col;
	g.relative = true;
	g.mean = 1;
	g.sigma = sigmaPC;
	
	for i = 1:length(corrs)
		g.row = corrs{i};
		groups{i} = g;
	end

	% Save uncorrelated
	for i = 1:length(uncorr)
		g.row = uncorr(i);
		groups{length(corrs) + i} = g;
	end
	
end