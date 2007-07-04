function thrushbaby(file)
    [Signal, S] = get_signal(file);
    [Mag, Phase, numcodons] = cumm_mag_phase(Signal);
    Dvec = diff_vectors(Mag, Phase, numcodons);
    
    global beached_whale;
    global ants termites;
    beached_whale = 1;
    
    start = 0;
    %while start > -1
        s = struct();
        for i=1:15
            [x,diffx] = displacement(S(13:end),Phase,numcodons,Dvec,{},{});
            s = helper(s, ants);
            s = helper(s, termites);
        end
        folder = 'J:\Frameshift\tb1';
        a = pearl('starling.pl', [num2str(start) ' "' which(file) '" "' folder '" ' stringify(s)]);
        a = eval(a)
        start = a{1};
        codon = a{2};
        runner(folder, 10, codon);
        pause;
    %end
end

function [best_name] = runner(folder, limit, codon)
    d = [dir([folder '/*.txt']); dir([folder '/*.fasta'])];
    global beached_whale; beached_whale = 1;
    
    best_name = ''; best_time = limit*2;
    global ants termites;
    for i = 1:length(d)
        name = d(i).name;
        
        % Print the data.
        [Signal, S] = get_signal([folder '/' name]);
        [Mag, Phase, n] = cumm_mag_phase(Signal);
        Dvec = diff_vectors(Mag, Phase, n);
        time = 0;
        for i=1:limit
            displacement(S(13:end), Phase, n, Dvec, {}, {});
            insects = horzcat(ants, termites);
            if strmatch(codon, insects)
                time = time + 1;
            end
        end
        disp([name ': ' num2str(time)]);
        if time < best_time
            best_time = time;
            best_name = name;
            if time == 0, break; end;
        end
    end
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