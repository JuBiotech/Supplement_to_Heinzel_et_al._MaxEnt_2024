function groups = prepareGroupsCglutSoxPhos( sigmaPC )
%PREPAREGROUPSCGLUT Summary of this function goes here
%   Detailed explanation goes here
    
    %%%
	% Biomass groups (6 cells)
	col = [44, 523];%44
	corrs = {[1, 2]};

	sigmaPC = sigmaPC / 100;
	groups = cell(1, 1);
	
	% Correlated
	
	g.relative = true;
	g.mean = 1;
	g.sigma = sigmaPC;
    for j = 1:length(col)
        g.column = col(j);
        g.row = corrs{1};
		groups{j} = g;
    end