% Finds the minimum codon changes necessary given
% a lower yield sequence and a newer yield sequence
% to modify the worse sequence into a sequence of
% the same yield of the newer sequence. Very compatible
% with ThrushBaby. Also, I need the number of iterations.
%
% backtrack('rpoS.txt', 'rpoSshiny.txt', 'C:\Work folder', 327)
%
% DID YOU KNOW? "Hansel and Gretel" is AT 327A.
function backtrack(old, new, folder, times)
    disp(['I''m about to obliterate ' folder '. Proceed?']); pause;
    preparedir(folder);
    
    global beached_whale; beached_whale = 1;
    
    start = 0;
    new = which(new);
    while (start > -1)
        disp('');
        mold = pearl('hansel.pl', sprintf('"%s" "%s" %s "%s"', which(old), new, num2str(start), folder));
        mold = eval(mold);
        start = mold{1}; start_folder = mold{2};
        
        disp(['Next codon target: ' num2str(start)]);
        new = runner(start_folder, times);
        disp(['         I choose ' new]);
    end
end

function [best_name] = runner(folder, times)
    d = [dir(fullfile(folder, '/*.txt')); dir(fullfile(folder, '*.fasta'))];
    best_yield = -1;
    for i = 1:length(d)
        name = d(i).name;
        
        % Print the data.
        [S, n, Dvec] = walrus_surprise(fullfile(folder, name));
        % Do this here after walrus_surprise clears the variables.
        global sands shoals;
        for i=1:times
            displacement(S(13:end), n, Dvec, {}, {});
        end
        yield = shoals/sands;
        disp(['         ' name ': ' num2str(yield)]);
        if yield > best_yield
            best_yield = yield;
            best_name = name;
        end
    end
    best_name = fullfile(folder, best_name);
end