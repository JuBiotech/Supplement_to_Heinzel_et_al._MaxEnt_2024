function rowUncertainty = summarizeRowUncertainty(uncertain, len, determineRow)
%SUMMARIZEROWUNCERTAINTY creates a row uncertainty mapping which is a cell
%array in which each row corresponds to a row in the underlying system
%(stoichiometry or constraints) and summarizes the uncertainty belonging to
%a given row in a system.
%Result is a cell array in which each cell contains an array of uncertainty
%structures relating to the row in the underlying system.
%
% Parameters:
%   - uncertain: Cell array of uncertainty structures
%   - len: Row count of the underlying system matrix
%   - determineRow: Function which determines the row index in the
%       underlying system matrix of the given uncertainty structure
%       Usage: index = determineRow(uncertaintyStructure);
%
% Returns:
%   - rowUncertainty: Cell array in which each cell corresponds to the row
%       with the same index in the underlying system matrix. The cells
%       contain a cell array of uncertainty structures affecting the given
%       row.

    rowUncertainty = cell(len, 1);
    for i = 1:length(uncertain)
        sto = uncertain{i};
        p = determineRow(sto);
        
        if isempty(rowUncertainty{p})
            rowUncertainty{p} = { sto };
        else
            if strcmp(sto.reaction, '#row')
                % Entire row is unertain
                rowUncertainty{p} = sto;
            else
                if isempty(rowUncertainty{p})
                    rowUncertainty{p} = [rowUncertainty{p}, { sto }];
                else
                    temp = rowUncertainty{p,1};
                    if iscell(temp)
                        temp = temp{1};
                    end
                    if ~strcmp(temp.reaction, '#row')
                        rowUncertainty{p} = [rowUncertainty{p}, { sto }];
                    end
                end
            end
        end
    end
end