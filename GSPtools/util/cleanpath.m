% Add folder and subdirectories to Matlab's search path, but ignore
% Subversion litter. Ay, Subversion.
function cleanpath(folder)
    addpath(genpath(folder));
    s = matlabpath;
    indices = [0 strfind(s, pathsep)];
    for i = 1:length(indices)-1
        t = s(indices(i)+1:indices(i+1)-1);
        if findstr(t, '.svn'); rmpath(t); end
    end
end