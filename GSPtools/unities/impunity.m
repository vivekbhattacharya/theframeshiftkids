% ----------------------------------------------------------------
% Stores yields to a results.txt file for perusal. `limit`
% specifies the number of iterations while `folder` is the folder
% path where `cornerstone.pl` dumped all its proteins. It can be
% any folder of txt files with protein sequences.
%
% Results file will ultimately end up in Matlab's work folder,
% appended.
%
% Usage: impunity(work folder, list of genes that should
%   frameshift +1 a la megaunity, frameshifters -1, number of 
%   iterations to run the model before recording yield)
% ----------------------------------------------------------------
function impunity(folder, fshifts, bshifts, limit)
    d = [dir([folder '/*.txt']); dir([folder '/*.fasta'])];
    
    function [yield] = find_yield(file)
        % ----------------------------------------------------------------
        % The meat of `impunity`, this is `megaunity` without the graphs
        % or the infinite loop, instead capped at `limit`. In addition,
        % it returns the yield instead of displaying it.
        % ----------------------------------------------------------------
        [S, n, Dvec] = walrus_surprise(file);        
        global shoals sands;
        for i=1:limit
            x = displacement(S(13:end), n, Dvec, fshifts, bshifts);
        end
        yield = shoals/sands;
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