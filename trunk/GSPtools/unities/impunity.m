% ----------------------------------------------------------------
% Impunity.m calculates the yield for all the genes in a given
% work folder, storing those results to a results.txt file in the
% Matlab workspace and also printing them to the standard output.
% It needs the predicted frameshifts in order to calculate yield
% and also a sample size.
%
% Usage: impunity('c:\work folder', [25], 100)
% ----------------------------------------------------------------
function impunity(folder, fshifts, limit)
    d = 0; singleton = 0;
    if isdir(folder)
        d = [dir([folder '/*.txt']); dir([folder '/*.fasta'])];
    else
        d = which(folder);
        singleton = 1;
    end
    
    function [yields] = find_yield(file)
        % ----------------------------------------------------------------
        % The meat of `impunity`, this is `megaunity` without the graphs
        % or the infinite loop, instead capped at `limit`. In addition,
        % it returns the yield instead of displaying it.
        % ----------------------------------------------------------------
        displacement = walrus_surprise(file);
        global shoals sands;
        yields = zeros(1, limit);
        for i=1:limit
            displacement(fshifts);
            yields(i) = shoals/sands;
            shoals = 0; sands = 0;
        end
    end
    
    if singleton
        yield = find_yield(d);
        disp([d ': ' num2str(yield)]);
        fprintf('\n');
        return;
    end
    
    for i = 1:length(d)
        name = d(i).name;
        yields = find_yield(fullfile(folder, name));
        fprintf('%s: %g +/- %g\n', name, mean(yields), std(yields));
    end
end