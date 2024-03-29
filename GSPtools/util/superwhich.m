% A which() that works automatically if the file already exists.
function [result] = superwhich(file)
    if strfind(file, 'http://')
        result = file;
        return;
    end

    result = which(file);
    if result, return; end;

    existence = exist(file);
    if existence > 0
        result = file;
        return;
    end

    error(sprintf('Could not find file %s', file));
end
