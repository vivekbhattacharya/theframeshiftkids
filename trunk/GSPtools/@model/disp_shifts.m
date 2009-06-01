% Displays frameshifts from loop() to standard output.
function disp_shifts(model)
    acc = '> ';
    for i = 1:length(model.wal)
        acc = [acc sprintf('%s; ', model.wal{i})];
    end

    acc = [acc sprintf('\n< ')];
    for i = 1:length(model.rus)
        acc = [acc sprintf('%s; ', model.rus{i})];
    end
    acc = [acc sprintf('\n')];
    disp(acc);
end
