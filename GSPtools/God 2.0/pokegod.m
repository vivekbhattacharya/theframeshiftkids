% Run the genetic algorithm for fudging tRNA availabilites in order to
% obtain the least deviation possible on a displacement plot, thereby
% improving protein efficiency magically. It takes a folder of
% "calibration" genes, which for now had better be bGH.
function pokegod(folder)
    load TAV.mat;
    load Codons.mat;

    % How big should the gene pool be?
    % How many generations should we have?
    % How big of a sample size for deviation?
    Pool_n = 4;
    times = 1;
    sample_size = 1;
    
    % Save optimal row because we won't have it after the for loop
    % completes looping.
    gen = generation_zero(TAV, Pool_n);
    optimal = [];
    for N = 1:times
        fprintf('Creating generation %g\n', N);
        yields = zeros(1, Pool_n);
        for row_n = 1:Pool_n
            taz = get_travel(Names, gen(row_n, :));
            yields(row_n) = get_yield(folder, taz, sample_size);
            fprintf(': Created person %g with yield %g\n', row_n, yields(row_n));
        end
        
        % Use yields as weights to determine who lives in the next
        % generation and who dies.
        [gen, yields] = sort_and_kill(gen, yields);
        optimal = gen(1, :);
        
        % Spawn the next generation.
        tmp = zeros(Pool_n, 64);
        for row_n = 1:Pool_n
            tmp(row_n, :) = chooser(gen, yields);
        end
        
        % Throw away the parents.
        gen = tmp;
    end
end

function [value] = get_yield(folder, travails, times)
    global Travel;
    Travel = travails;

    good = 0; bad = 0;
    classify(folder, '', @helper);
    function helper(d, n, filename, image)
        global shoals sands;
        for i = 1:times, d([]); end;
        yield = shoals/sands;
        
        % OR hack
        if [strfind(filename, '101') strfind(filename, '105')]
            good = good + yield;
        else
            bad = bad + yield;
        end
    end

    % There are six sequences besides 101 and 105.
    % (good/2) / (bad / 6)
    value = good*3/bad;
end