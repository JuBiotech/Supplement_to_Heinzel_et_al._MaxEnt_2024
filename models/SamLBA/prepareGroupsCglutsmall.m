function groups = prepareGroupsCglutsmall( sigmaPC )
%PREPAREGROUPSCGLUT Summary of this function goes here
%   Detailed explanation goes here

	%%%%
	% Biomass groups (6 cells)
	col = 44;
	corrs = {[11, 20, 22], [30, 31], [6, 10, 14], [38, 40], [52, 53], [55, 57, 58], [4, 62, 63]};
    uncorr = [7, 16, 36, 44, 46];
% 	corrs = {[11, 20, 22], [4, 53], [10, 14]};
% 	uncorr = [6 , 7, 16, 30, 31, 36, 38, 40, 44, 46, 52, 55, 57, 58, 62, 63];
%     corrs = {[11, 20, 22]};
% 	uncorr = [4, 6, 7, 10, 14, 16, 30, 31, 36, 38, 40, 44, 46, 52, 53, 55, 57, 58, 62, 63];    

	sigmaPC = sigmaPC / 100;
	groups = cell(12, 1); % 19 bzw 21
	
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