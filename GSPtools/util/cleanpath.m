% Add folder and subdirectories to Matlab's search path, but ignore
% Subversion litter. Ay, Subversion. AY.
function cleanpath(folder)
    warning off MATLAB:rmpath:DirNotFound;

    % What's the old path?
    old = which('classify.m');
    if length(old) > 0
        % -1 to remove the pathsep.
        old(strfind(old, 'unities')-1:end) = [];

        % Remove the old path.
        s = path;
        indices = [0 strfind(s, pathsep)];
        for i = 1:length(indices)-1
            idx = indices(i)+1:indices(i+1);
            t = s(idx);
            if strfind(t, old); rmpath(t); end;
        end
    end

    % Add the new path.
    s = genpath(folder);
    indices = [0 strfind(s, pathsep)];
    for i = 1:length(indices)-1
        idx = indices(i)+1:indices(i+1);
        t = s(idx);

        if strfind(t, 'pearls\Tests'); s(idx) = '`'; end
        if strfind(t, '.svn'); s(idx) = '`'; end
        if strfind(t, '.hg'); s(idx) = '`'; end
    end
    addpath(s(find(s ~= '`')));
end
