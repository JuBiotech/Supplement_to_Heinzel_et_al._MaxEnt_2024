function mets = getMetaboliteHeadingsOfGroups( model, gp )
%GETMETABOLITEHEADINGSOFGROUPS Returns the metabolite names of the given
%groups. It defaults to the first metabolite (row in stoichiometry).
%If no metabolite names are available (group does not affect stoichiometry)
%it returns 'Group [RowNumber]'.
%
% Parameters:
%   - model: Model structure.
%   - gp: Cell array of group structures.
%
% Returns a cell array with metabolite headings for each given group.

    mets = [];
	for i = 1:length(gp)
		if lower(gp{i}.type) == 's'
		   	mets = [mets, model.mets(gp{i}.row(1))];
		else
			for j = 1:length(gp{i}.row)
    			mets = [mets, {['Group ' gp{i}.type ' ' num2str(gp{i}.row(j))]}];
            end
		end
	end
	
	if size(mets, 1) > 1
		mets = mets';
	end
end
