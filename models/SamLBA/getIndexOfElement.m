function index = getIndexOfElement(vals, e)
%GETINDEXOFELEMENT determines the index of the value e in the string cell
%array vals.
%
% Parameters:
%   - vals: Collection of values
%   - e: Element to search for
%
% Returns:
%   - index: Index of the element in the collection

	if iscell(vals) && isnumeric(e)
		index = e;
		return;
	end

    inds = 1:length(vals);
    index = inds(strcmp(vals, e));
end
