function thrushbaby(file, work_folder, times)
    global beached_whale; beached_whale = 1;
    start = 0;

    [folder, name, ext] = fileparts(which(file));
    contents = struct('name', [name ext]);
    codon = '_';
    while start > -1
        [file, s] = runner(folder, times, codon, [contents]);
        
        % J:\work\0, J:\work\25, J:\work\200, etc.
        folder = fullfile(work_folder, num2str(start));
        if isdir(folder), rmdir(folder, 's'); end;
        disp(['New folder: ' folder]);
        mkdir(folder);
        
        disp(['Now running Perl on: ' file]);
        a = pearl('starling.pl', [num2str(start) ' "' file '" "' folder '" ' stringify(s)]);
        disp(['Perl says: ' a]);
        a = eval(a);
        start = a{1}; codon = a{2};
        disp(['Next codon target: ' codon sprintf('\n')]);
        contents = [dir(fullfile(folder, '/*.txt')); dir(fullfile(folder, '*.fasta'))];
    end
end

function [best_name, best_struct] = runner(folder, times, codon, d)
    best_time = times*2;
    best_struct = struct();
    global ants termites;
    for i = 1:length(d)
        name = d(i).name;
        
        % Print the data.
        [Signal, S] = get_signal([folder '/' name]);
        [Mag, Phase, n] = cumm_mag_phase(Signal);
        Dvec = diff_vectors(Mag, Phase, n);
        time = 0; s = struct();
        for i=1:times
            displacement(S(13:end), Phase, n, Dvec, {}, {});
            insects = horzcat(ants, termites);
            s = helper(s, insects);
            if strmatch(codon, insects)
                time = time + 1;
            end
        end
        disp([name ': ' num2str(time)]);
        if time < best_time
            best_time = time;
            best_name = name;
            best_struct = s;
            if time == 0, break; end;
        end
    end
    best_name = fullfile(folder, best_name);
    disp(['I choose ' best_name]);
end

function [s] = helper(s, insect)
    for i=1:length(insect)
        key = strrep(insect{i}, ',', '_');
        if isfield(s, key), s.(key) = s.(key) + 1;
        else s = setfield(s, key, 0);
        end
    end
end

function [sorrow] = stringify(s)
    f = fieldnames(s);
    sorrow = [];
    for i=1:length(f)
        key = f{i}; value = s.(key);
        if value ~= 0, sorrow = [sorrow key '_' num2str(value) ' ']; end
    end
end