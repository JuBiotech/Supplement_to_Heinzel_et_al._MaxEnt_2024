function groups = prepareGroupsCglutSoxPhos3( sigmaPC )
%PREPAREGROUPSCGLUT Summary of this function goes here
%   Detailed explanation goes here

	%%%%
	% Biomass groups (4 cells)
	col = 523;
	corrs = {[1,40]};%{[1,2]};
    
	sigmaPC = sigmaPC / 100;
	groups = cell(1, 1); % 19 bzw 21
	
	% Correlated
	g.column = col;
	g.relative = true;
	g.mean = 1;
	g.sigma = sigmaPC;
	for i = 1:length(corrs)
		g.row = corrs{i};
		groups{i} = g;
    end