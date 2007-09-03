function pokegod(folder)
    load TAV.mat;
    load Codons.mat;
    
    Pool_n = 64;
    gen = generation_zero(TAV, Pool_n);
    times = 10;
    optimal = [];
    for N = 1:times
        yields = zeros(1, Pool_n);
        for row_n = 1:Pool_n
            taz = gen(row_n, :);
            yields(row_n) = get_yield(folder, get_travel(Names, taz));
        end
        [gen, yields] = sort_and_kill(gen, yields);
        optimal = gen(1, :);
        
        tmp = zeros(Pool_n, 64);
        for row_n = 1:Pool_n
            tmp(row_n) = chooser(gen, yields);
        end
    end
end

function [value] = get_yield(folder, travails)
    global Travel;
    Travel = travails;

    good = 0;
    bad = 0;
    times = 10;
    classify(folder, '', @helper);
    function helper(d, n, filename, image)
        global shoals sands;
        for i = 1:times, d([]); end;
        yield = shoals/sands;
        
        disp(filename);
        if strfind(filename, '101')
            good = good + yield;
        elseif strfind(filename, '105')
            good = good + yield;
        else
            bad = bad + yield;
        end
    end
    
    % (good/2) / (bad / 6)
    value = good*3/bad;
end