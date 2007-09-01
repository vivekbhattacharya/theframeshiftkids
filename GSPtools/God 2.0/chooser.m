% Randomly select two TAV rows from the matrix, biasing toward
% those with lower deviation.
function [tav1, tav2, wt1, wt2] = chooser(tavmatrix, weights)
    [rows columns] = size(tavmatrix);

    % Simple distribution where the best TAV gets `rows` chips,
    % the next gets `rows - 1` chips, and so forth.
    toplimit = rows*(rows+1)/2;
    r1 = selectnum(ceil(toplimit * rand));
    r2 = selectnum(ceil(toplimit * rand));

    tav1 = tavmatrix(r1,:);
    tav2 = tavmatrix(r2,:);
    wt1 = weights(r1);
    wt2 = weights(r2);
end

function num = selectnum(init)
    % Solve init = n*(n-1)/2
    r = roots([.5 .5 -init])';
    num = max(ceil(r));
end