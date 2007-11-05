% Iterate over the cell array and print the results
% sanely. If only Matlab had a join function like
% every other modern language in the world.
function disp_shifts()
    global ants termites anthill;
    if length(termites)
        if termites{end} == 0, return; end;
    end
    
    fprintf('> ');
    for i = 1:length(anthill)
        fprintf('%s %g; ', ants{i}, anthill(i));
    end

    fprintf('\n< ');
    for i = 1:length(termites)
        fprintf('%s; ', termites{i});
    end
    fprintf('\n');
end