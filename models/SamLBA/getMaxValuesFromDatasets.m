function maxVals = getMaxValuesFromDatasets( folder, field, normalizer )
%GETMAXVALUESFROMDATASETS Postprocesses a given set of data files (.mat) and extracts the maximum absolute values.
% The given folder is scanned for .mat files. Each mat-file is loaded and the given field of the res structure
% is processed. For this the given normalizer function is used, which has the following signature:
%   normalizer(res.<field>, res.origSol)
% For each data set (mat-file) the max( abs( normalizer(...) ) ) is recorded.
%
% Parameters:
%	- folder: Folder to look for mat-files.
%	- field: String. Name of the field of the res-structure to process.
%	- normalizer: Function handle.
%
% Returns a vector with the maximum absolute values. 
%	The order corresponds to the order of found mat-files by the what()-function.

	files = what(folder);
	maxVals = zeros(length(files.mat), 1);
	
	for i = 1:length(files.mat)
		f = [folder '/' files.mat{i}];
		res = load(f);

		if ~isfield(res, 'res')
			continue;
		end
		disp(f);
		
		vals = res.res.(field);
		
		maxVals(i) = max(abs(normalizer(vals, res.res.origSol)));
	end

end

