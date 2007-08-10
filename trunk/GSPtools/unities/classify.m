function classify(folder, subfolder, crusade)
    subdir = fullfile(folder, subfolder);
    try
        mkdir(subdir);
    catch
        disp(sprintf('Folder %s does not exist', folder));
        return;
    end
    
    function helper(path)
        [displacement, n] = walrus_surprise(path);
        [folder, file] = fileparts(path);
        
        image = fullfile(subdir, [file '.png']);
        crusade(displacement, n, file, image);
    end
    
    if isdir(folder)
        d = [dir([folder '/*.txt']); dir([folder '/*.fasta'])];
        for i = 1:length(d)
            file = fullfile(folder, d(i).name);
            helper(file);
        end
    else
        % `folder` is actually a file.
        helper(which(folder));
    end
end