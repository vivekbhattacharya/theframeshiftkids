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
    
    function [yield] = find_yield(file)
        % ----------------------------------------------------------------
        % The meat of `impunity`, this is `megaunity` without the graphs
        % or the infinite loop, instead capped at `limit`. In addition,
        % it returns the yield instead of displaying it.
        % ----------------------------------------------------------------
        displacement = walrus_surprise(file);        
        global shoals sands;
        for i=1:limit, displacement(fshifts); end;
        yield = shoals/sands;
    end
    
    if singleton
        yield = find_yield(d);
        disp([d ': ' num2str(yield)]);
        fprintf('\n');
        return;
    end
    
    % Print the header.
    fid = fopen('results.txt', 'a+');
    fwrite(fid, ['#' num2str(limit) ' iterations']);
    fprintf(fid, '\n');
    
    for i = 1:length(d)
        name = d(i).name;
        yield = find_yield(fullfile(folder, name));
        
        % Print the data.
        disp([name ': ' num2str(yield)]);
        fwrite(fid, [name ': ' num2str(yield)]);
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n\n');
    fclose(fid);
end