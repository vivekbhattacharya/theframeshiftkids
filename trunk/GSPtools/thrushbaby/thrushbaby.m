% Optimizes a gene sequence to produce one
% with no frameshifts, theoretically anyway.
% Arguments: file, work folder, number of iterations
function thrushbaby(file, work_folder, times)
    codon = '_';
    
    disp(['I''m about to obliterate ' work_folder '. Proceed?']); pause;
    preparedir(work_folder);
    
    folder = fullfile(work_folder, '0'); preparedir(folder);
    copyfile(which(file), folder);
    current = 0;
    while 1
        folder = fullfile(work_folder, num2str(current));
        [yield, file, current] = runner(folder, times, codon, current);
        if current == -1, disp('Finished'); break; end;
        disp(['Best yield: ' num2str(yield)]);
        disp(['Next codon target: ' num2str(current) sprintf('\n')]);
        
        disp(['Now running Perl on: ' file]);
        pearl('starling.pl', sprintf('"%s" "%s" %g', file, work_folder, current));
    end
    disp(sprintf('\nThe final file is %s\n', file));
end

function [best_yield, best_name, next] = runner(folder, times, codon, last)
    best_yield = 300;
    classify(folder, '', @helper);
    
    function helper(displacement, n, title, image)
        for i = 1:times, displacement([]); end;
        global shoals sands;
        yield = shoals/sands;

        if yield < best_yield
            disp([title ': ' num2str(yield)]);
            best_yield = yield;
            best_name = title;
        end
    end
    best_name = fullfile(folder, best_name);
    d = walrus_surprise(best_name);
    railyard = abs(d([])); next = -1;
    for i = last+1:length(railyard)
        if railyard(i) > 0.138, next = i; break; end;
    end
end 