function groups = prepareGroupsEcoli( sigmaPC )
%PREPAREGROUPSECOLI Prepares the groups for the E. coli network.
% The groups created by this function preserve the correlation of the
% coefficients of the biomass reaction.
%
% Correlated coefficients will be put together in one group.
% For description of the group structure see mcGroups.m .
%
% Parameters:
%	- sigmaPC: Standard deviation of the biomass coefficients in percent(!)
%
% Returns a cell array of group structures.

	%%%%
	% Biomass groups (16 cells)
	% In E. coli reaction 13 is biomass.
	col = 13;
	corrs = {[13, 17, 41, 43, 60], [10, 21], [50, 51], [52, 53]};
	uncorr = [3, 14, 23, 26, 33, 34, 36, 38, 58, 59, 62, 66];

	sigmaPC = sigmaPC / 100;
	groups = cell(18, 1);	% Allocate memory for the groups. There are 18 
							% groups (16 biomass, 2 separate: O_2, H)
	
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
	
	%%%%
	% EX_h(e) group (1 cell)
	g.column = 31;
	g.row = 44;
	g.relative = false;
	g.bounds = [-3, 1];
	g.isInteger = true;
	g.nonZero = true;
	groups{17} = g;

	%%%%
	% EX_O2(e) group (1 cell)
	g.column = 36;
	g.row = 57;
	g.relative = false;
	g.bounds = [-3, 1];
	g.isInteger = true;
	g.nonZero = true;
	groups{18} = g;
	
end

