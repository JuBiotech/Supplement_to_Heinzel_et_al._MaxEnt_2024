function [poolNames, reationNames] = getPoolsAndReactions(file)
%GETPOOLSANDREACTIONS reads pool and reaction names from an HDF5-file.

    pools = hdf5read(file,'/stoichiometry/pools');
    reactions = hdf5read(file,'/stoichiometry/reactions');

    poolNames = cell(length(pools),1);
    reationNames = cell(length(reactions),1);

    % Get pool and reaction names
    for i=1:length(pools)
        poolNames{i} = pools(i).Data;
    end
    for i=1:length(reactions)
        reationNames{i} = reactions(i).Data;
    end
end