% Iterate over the cell array and print the results
% sanely. If only Matlab had a join function like
% every other modern language in the world.
function disp_shifts()
    global ants termites anthill;
    if length(termites)
        if termites{end} == 0, return; end;
    end
    
    acc = '> ';
    for i = 1:length(anthill)
        acc = [acc sprintf('%s %g; ', ants{i}, anthill(i))];
    end

    acc = [acc sprintf('\n< ')];
    for i = 1:length(termites)
        acc = [acc sprintf('%s; ', termites{i})];
    end
    acc = [acc sprintf('\n')];
    disp(acc);
end