function groups = prepareGroupsCglut( sigmaPC )
%PREPAREGROUPSCGLUT Summary of this function goes here
%   Detailed explanation goes here

	%%%%
	% Biomass groups (6 cells)
	col = 300;
	corrs = {[3, 5, 10], [304, 305]};
	uncorr = [300, 301, 302, 303, 306];
	sigmaPC = sigmaPC / 100;
	groups = cell(7, 1);
	
	% Correlated
	g.column = col;
	g.relative = true;
	g.mean = 1;
	g.sigma = sigmaPC;
	for i = 1:length(corrs)
		g.row = corrs{i};
		groups{i} = g;
    end
	% Uncorrelated
	for i = 1:length(uncorr)
		g.row = uncorr(i);
		groups{length(corrs) + i} = g;
    end
	