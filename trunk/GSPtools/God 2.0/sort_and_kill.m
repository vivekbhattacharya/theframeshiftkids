% Kill half the population and sort the rest by weights.                                                                     
function [tav, weights] = sort_and_kill(tav, weights);
    [rows cols] = size(tav);
    weight_col = cols + 1;
    cutoff = ceil(rows/2);

    % Schwartzian transform
    tav(:, weight_col) = weights;
    % Sort rows by weights
    weights = sort(weights, 'descend');
    % Negative to sort descending
    tav = sortrows(tav, -weight_col)
    % Kill half the population.
    tav(1:cutoff, :) = [];
    % And delete their weights.
    weights(1:cutoff) = [];
    % Undo transform.
    tav(:, weight_col) = [];
end