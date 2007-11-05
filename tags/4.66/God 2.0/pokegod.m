% Run the genetic algorithm for fudging tRNA availabilites in order to
% obtain the least deviation possible on a displacement plot, thereby
% improving protein efficiency magically. It takes a folder of
% "calibration" genes, which for now had better be bGH.
function [gen] = pokegod(folder)
    load TAV.mat;
    load Codons.mat;

    % How big should the gene pool be?
    % How many generations should we have?
    % How big of a sample size for deviation?
    % How many nucleotides of a person to mutate per spawn?
    pool_n = 65;
    times = 100;
    sample_n = 50;
    radiation_n = 6;
    
    % Save optimal row because we won't have it after the for loop
    % completes looping.
    gen = generation_zero(TAV, pool_n);
    optimal = [];
    for N = 1:times
        fprintf('Creating generation %g\n', N);
        yields = zeros(1, pool_n);
        for row_n = 1:pool_n;
            taz = get_travel(Names, gen(row_n, :));
            yields(row_n) = get_yield(folder, taz, sample_n);
            fprintf(': Created person %g with yield %g\n', row_n, yields(row_n));
        end
        
        % Use yields as weights to determine who lives in the next
        % generation and who dies.
        [gen, yields] = sort_and_kill(gen, yields);
        optimal = gen(1, :);
        
        % Spawn the next generation.
        tmp = zeros(pool_n, 64);
        for row_n = 1:pool_n
            tmp(row_n, :) = spawn(gen, yields, radiation_n);
        end
        
        % Throw away the parents.
        gen = tmp;
        % Let Ctrl-C work.
        myfile = ['banana' num2str(N) '.mat'];
        save(fullfile(folder, myfile), 'optimal');
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
        if [strfind(filename, '101') strfind(filename, '105') strfind(filename, '112') strfind(filename, '115')]
            good = good + yield;
        else
            bad = bad + yield;
        end
    end

    % There are six sequences besides 101 and 105.
    % (good/2) / (bad / 6)
    value = good/bad;
end