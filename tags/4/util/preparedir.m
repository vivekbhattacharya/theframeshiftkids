function preparedir(folder)
    if isdir(folder), rmdir(folder, 's'); end;
    mkdir(folder);
end