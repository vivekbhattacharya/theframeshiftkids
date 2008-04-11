% Given the folder for bGH sequences and a matrix of TAVs, I sort the
% TAVs by optimal separation between the 101 and the 105 yields.
function paradiselost(folder, tavs)
    [rows cols] = size(tavs);
    yields = zeros(rows, 2);

    load Codons.mat;
    for i = 1:rows
        yields(i, 1) = i;
        travel = get_travel(Names, tavs(i, :));
        yields(i, 2) = get_yield(folder, travel, 10);
        fprintf('%g: %g\n', i, yields(i, 2));
    end
    sortrows(yields, -2)
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
        if [strfind(filename, '101') strfind(filename, '105')]
            good = good + yield;
        else
            bad = bad + yield;
        end
    end

    % There are six sequences besides 101 and 105.
    % (good/2) / (bad / 6)
    value = good*3/bad;
end
