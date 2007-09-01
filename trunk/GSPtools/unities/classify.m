function classify(folder, subfolder, crusade, varargin)
    helper = @campbag;
    if length(varargin) > 0, helper = @preparation; end
    
    % Handle files as if they were folders with magic.
    boulder = folder;
    if ~isdir(boulder), boulder = fileparts(which(folder)); end
    if isempty(boulder)
        fprintf('File %s does not exist\n\n', folder);
        return;
    end
    
    if subfolder
        subdir = fullfile(boulder, subfolder);
        warning off MATLAB:MKDIR:DirectoryExists;
        mkdir(subdir);
    else
        subdir = boulder;
    end
    
    if isdir(folder)
        d = [dir([folder '/*.txt']); dir([folder '/*.fasta'])];
        for i = 1:length(d)
            file = fullfile(folder, d(i).name);
            helper(file);
        end
    else helper(which(folder));end;
    
    % Polar plots, mostly
    function preparation(path)
        [folder, file, ext] = fileparts(path);
        
        image = fullfile(subdir, [file '.png']);
        crusade(path, [file ext], image);
    end
    
    % Everybody else gets a free displacement with
    % magical folder structure management.
    function campbag(path)
        [displacement, n] = walrus_surprise(path);
        [folder, file, ext] = fileparts(path);
        
        image = fullfile(subdir, [file '.png']);
        crusade(displacement, n, [file ext], image);
    end
end