% Meta-meta-meta-metafunction for unities.
function classify(place, subdir, fs, crusade)
    config;
    global Config;

    % Special-case single-file `place`s.
    if ~isdir(place)
        file   = superwhich(place);
        place  = fileparts(file);

        if subdir
            subdir = prepare_subdir(place, subdir);
        end

        helper(superwhich(file));
        return;
    end

    % Qualify subdir.
    if subdir
        subdir = prepare_subdir(place, subdir);
    end

    files = [dir([place '/*.txt']); dir([place '/*.fasta'])];
    for i = 1:length(files)
        file = fullfile(place, files(i).name);
        helper(file);
    end

    % Everybody else gets a free displacement with
    % magical folder structure management.
    function helper(path)
        m = Config.model(path, fs);
        [folder, file, ext] = fileparts(path);

        if subdir
            png = fullfile(subdir, [file '.png']);
            crusade(m, [file ext], png);
        else
            crusade(m, [file ext]);
        end
    end
end

% Create a subdirectory under a parent directory if necessary. Returns
% the full path of the subdirectory.
function subdir = prepare_subdir(dir, subdir)
    subdir = fullfile(dir, subdir);
    warning off MATLAB:MKDIR:DirectoryExists;
    mkdir(subdir);
end