% Add folder and subdirectories to Matlab's search path, but ignore
% Subversion litter. Ay, Subversion.
function cleanpath(folder)
    addpath(genpath(folder));
    s = matlabpath;
    p = 1;

    while true
        t = strtok(s(p:end), pathsep);
        if findstr(t, '.svn')
            rmpath(t);
        end
        p = p + length(t) + 1;
        if isempty(strfind(s(p:end),';')) break, end;
    end
end