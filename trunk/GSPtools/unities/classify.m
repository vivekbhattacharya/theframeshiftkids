function classify(folder, subfolder, crusade, varargin)
    subdir = fullfile(folder, subfolder);
    try
        mkdir(subdir);
    catch
        disp(sprintf('Folder %s does not exist', folder));
        return;
    end
    
    helper = @campbag;
    if length(varargin) > 0, helper = @preparation; end
    
    function preparation(path)
        [folder, file, ext] = fileparts(path);
        
        image = fullfile(subdir, [file '.png']);
        crusade(path, [file ext], image);
    end
    
    function campbag(path)
        [displacement, n] = walrus_surprise(path);
        [folder, file, ext] = fileparts(path);
        
        image = fullfile(subdir, [file '.png']);
        crusade(displacement, n, [file ext], image);
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