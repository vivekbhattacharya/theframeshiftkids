% Optimizes a gene sequence to produce one
% with no frameshifts, theoretically anyway.
% Arguments: file, work folder, number of iterations
function thrushbaby(file, work_folder, times)
    start = 0; codon = '_';
    
    disp(['I''m about to obliterate ' work_folder '. Proceed?']); pause;
    preparedir(work_folder);
    
    folder = fullfile(work_folder, num2str(start)); preparedir(folder);
    copyfile(which(file), folder);
    while start > -1
        folder = fullfile(work_folder, num2str(start));
        disp(['Next codon target: ' codon sprintf('\n')]);
        
        [file, s] = runner(folder, times, codon);
        
        disp(['Now running Perl on: ' file]);
        a = pearl('starling.pl', sprintf('%g "%s" "%s" %s', start, file, work_folder, stringify(s)));
        
        disp(['Perl says: ' a]);
        a = eval(a); start = a{1}; codon = a{2};
    end
    disp(sprintf('\nThe final file is %s\n', file));
end

function [best_name, best_struct] = runner(folder, times, codon)
    d = [dir(fullfile(folder, '/*.txt')); dir(fullfile(folder, '*.fasta'))];
    best_time = times*2;
    best_struct = struct();
    global ants termites;
    for i = 1:length(d)
        name = d(i).name;
        
        % Print the data.
        [S, n, Dvec] = walrus_surprise(fullfile(folder, name));
        time = 0; s = struct();
        for i=1:times
            displacement(S(13:end), n, Dvec, {}, {});
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
end

% Update the tabulation struct `s` with
% counts of codon frameshifts
function [s] = helper(s, insect)
    for i=1:length(insect)
        key = strrep(insect{i}, ',', '_');
        if ~isfield(s, key), s = setfield(s, key, 0); end;
    end
end

% Parse the tabulation struct `s` into
% something usable for starling.pl
function [sorrow] = stringify(s)
    f = fieldnames(s);
    sorrow = '';
    for i=1:length(f)
        sorrow = [sorrow f{i} ' '];
    end
end