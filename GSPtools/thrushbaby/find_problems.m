function find_problems(file)
    [Signal, S] = get_signal(file);
    [Mag, Phase, numcodons] = cumm_mag_phase(Signal);
    Dvec = diff_vectors(Mag, Phase, numcodons);
    
    global beached_whale;
    global ants termites;
    beached_whale = 1;
    
    start = 0;
    while start > -1
        s = struct();
        for i=1:50
            [x,diffx] = displacement(S(13:end),Phase,numcodons,Dvec,{},{});
            s = helper(s, ants);
            s = helper(s, termites);
        end
        %stringify(s)
        pearl('jovial.pl', [num2str(start) ' "' which(file) '" J:\Frameshift\tb1 ' stringify(s)])
        impunity('J:\Frameshift\tb1', {}, {}, 100);
        pause;
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