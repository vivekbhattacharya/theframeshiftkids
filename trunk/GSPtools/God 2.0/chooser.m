% Randomly select two TAV rows from the matrix, biasing toward
% those with lower deviation.
function [child] = chooser(tavmatrix, weights)
    [rows columns] = size(tavmatrix);

    % Simple distribution where the best TAV gets `rows` chips,
    % the next gets `rows - 1` chips, and so forth.
    toplimit = rows*(rows+1)/2;
    r1 = selectnum(ceil(toplimit * rand));
    r2 = selectnum(ceil(toplimit * rand));

    tav1 = tavmatrix(r1,:);
    tav2 = tavmatrix(r2,:);
    
    child = mean([tav1; tav2]);
end

function num = selectnum(init)
    % Solve init = n*(n+1)/2 => n^2 + n - 2*init = 0
    b = 1; c = -2*init;
    q = sqrt(b^2 - 4*c);
    r = ([q -q] - b)/2;
    num = max(ceil(r));
end