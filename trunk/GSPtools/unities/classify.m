% Meta-meta-meta-metafunction for unities.
function classify(place, subdir, fs, crusade)
    config;
    global Config;

    % Handle files as if they were folders with magic.
    if ~isdir(place)
        if subdir
            error('classify: no subdir allowed when passed a file');
        end

        file = superwhich(place);
        helper(superwhich(file));
        return;
    end

    adir = place;

    % Qualify subdir.
    if subdir
        subdir = fullfile(adir, subdir);
        warning off MATLAB:MKDIR:DirectoryExists;
        mkdir(subdir);
    end

    files = [dir([adir '/*.txt']); dir([adir '/*.fasta'])];
    for i = 1:length(files)
        file = fullfile(adir, files(i).name);
        helper(file);
    end

    % Everybody else gets a free displacement with
    % magical folder structure management.
    function helper(path)
        m = Config.model(path, fs);
        [folder, file, ext] = fileparts(path);

        if subdir
            image = fullfile(subdir, [file '.png']);
            crusade(m, [file ext], image);
        else
            crusade(m, [file ext]);
        end
    end
end
