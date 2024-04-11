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
	% Biomass groups (7 cells)
	% In Coryne Glutamicum reaction 300 is biomass.
	col = 300;
	corrs = {[3, 5, 10], [304, 305]};
	uncorr = [300, 301, 302, 303, 306];

	sigmaPC = sigmaPC / 100;
	groups = cell(7, 1);	% Allocate memory for the groups. There are 7 
							% groups
	
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

