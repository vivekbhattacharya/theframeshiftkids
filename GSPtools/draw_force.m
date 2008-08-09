function draw_force(aforce)
    max = length(aforce);
    dx = get_force(2:0.1:max, aforce);
    plot(2:0.1:max, dx);
end

function [dx] = get_force(n_base, aforce)
    n_base_index = round(n_base);

    previous = aforce(n_base_index - 1)
    next = aforce(n_base_index);

    % Slope is rise over run, and run is one here. The code condenses the
    % following three lines into one line.

    % slope = next - previous;
    % intercept = next - slope * (n_base_index + 0.5);
    % dx = slope * n_base + intercept;

    dx = (next - previous) .* (n_base - n_base_index - 0.5) + next;
end
