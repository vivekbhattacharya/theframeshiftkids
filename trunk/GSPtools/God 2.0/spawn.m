% Randomly select two TAV rows from the matrix, biasing toward those
% with lower deviation. Then mutate mutt_n number of times
% conseratively, moreso than generation_zero.
%
% Precondition: The matrix of TAVs must be sorted in descending
% order of corresponding yield. The yields must be sorted in
% descending order, period. See sort_and_kill.
function [child] = spawn(tavmatrix, weights, mutt_n)
    [rows columns] = size(tavmatrix);

    % Simple distribution where the best TAV gets `rows` chips,
    % the next gets `rows - 1` chips, and so forth.
    toplimit = rows*(rows+1)/2;
    r1 = selectnum(ceil(toplimit * rand));
    r2 = selectnum(ceil(toplimit * rand));

    % Start by assuming tav1 has the lower yield.  If not, fix
    % accordingly. (r1 has a lower yield than r2 iff r1 > r2 because
    % we sorted in descending order.)
    C = 0.7;
    if r1 < r2, C = 1 - C; end;

    tav1 = tavmatrix(r1,:);
    tav2 = tavmatrix(r2,:);
    % Weighted mean
    child = C * tav1 + (1 - C) * tav2;

    % See generation_zero. Mutate within interval [0.75, 1.25] =
    % 0.75 + 0.5*[0, 1].
    changes = rand_int(columns, mutt_n);
    r = 0.75 + 0.5 * rand(1, mutt_n);
    child(changes) = child(changes) .* r;
end

function num = selectnum(init)
    % Solve init = n*(n+1)/2 => n^2 + n - 2*init = 0
    b = 1; c = -2*init;
    q = sqrt(b^2 - 4*c);
    r = ([q -q] - b)/2;
    num = max(ceil(r));
end
