function ecoliGeneratePlotsForDisturbance( model, dpc, mode )
%ECOLIGENERATEPLOTSFORDISTURBANCE Summary of this function goes here
%   Detailed explanation goes here

	if nargin == 2
		mode = 'normal';
	end

	groups = prepareGroupsEcoli(dpc);
	res = mcGroups(model, groups, 4000, mode);
	
	prefix = ['ecoli_optBM_mcBM_' num2str(dpc) 'pc_4000samples'];
	save([prefix '.mat'], 'res');

	plotHistFitNormal((res.vals - res.origSol(13)) / res.origSol(13)*100, 20);
	saveas(gcf, [prefix '_objValDeviation_pc.fig']);
	matlab2tikz([prefix '_objValDeviation_pc.tikz']);

	plotHistPercent((res.diffs) / norm(res.origSol)*100, 20);
	saveas(gcf, [prefix '_2normDeviation_pc.fig']);
	matlab2tikz([prefix '_2normDeviation_pc.tikz']);	
end

