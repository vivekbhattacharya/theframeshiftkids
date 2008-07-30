% Iterate over the cell array and print the results
% sanely. If only Matlab had a join function like
% every other modern language in the world.
function disp_shifts()
    global store;
    if length(store.termites)
        if store.termites{end} == 0, return; end;
    end

    acc = '> ';
    for i = 1:length(store.anthill)
        acc = [acc sprintf('%s %g; ', store.ants{i}, store.anthill(i))];
    end

    acc = [acc sprintf('\n< ')];
    for i = 1:length(store.termites)
        acc = [acc sprintf('%s; ', store.termites{i})];
    end
    acc = [acc sprintf('\n')];
    disp(acc);
    clear global store;
end
