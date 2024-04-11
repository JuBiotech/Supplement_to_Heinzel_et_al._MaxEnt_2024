function results = postProcessExistingDatasets(folder, postProcessor, returnGlobalResults)
%POSTPROCESSEXISTINGDATASETS Postprocesses existing MC datasets.
%
% Applies the postProcessor to each MC dataset (mat-File) in a directory
% and saves the results back or collects the results for each dataset and
% returns it.
%
% Parameters:
%   - folder: Folder to look for mat-Files. Each mat-File containing a
%       "res" structure will be processed. The other files are ignored.
%   - postProcessor: Function handle which is called for every dataset.
%       Signature is [res, resave] = postProcessor(fileName, res);
%       res is the loaded dataset which will be rewritten to the file if
%       resave is true. If resave is false, the file will not be touched.
%   - returnGlobalResults: Optional. Set to true if the results of the
%       postProcessor shall be collected and returned. This is independent
%       from the resave option above.
%
% Returns (if returnGlobalResults is enabled) a cell array with the results
% for each processed file (i.e. the returned res variable from
% postProcessor).

    if (nargin < 3) || isempty(returnGlobalResults)
        returnGlobalResults = false;
    end

	files = what(folder);
    
    if returnGlobalResults
        results = cell(length(files.mat),1);
    else
        results = [];
    end
    
	for i = 1:length(files.mat)
		f = [folder '/' files.mat{i}];
		res = load(f);

		if ~isfield(res, 'res')
			continue;
		end
		disp(f);
		res = res.res;
        
		[res, resave] = postProcessor(files.mat{i}, res);
                		
		if resave
            save(f, 'res');
        end
    
        if returnGlobalResults
            results{i} = res;
        end
	end

end