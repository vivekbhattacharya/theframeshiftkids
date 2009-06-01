% Return the yield, whatever that means at the current time, of the
% current model after calling loop().
function [yield] = yield(model, xs)
    global Config;
    fs = model.fs;

    % 0 = deviation, 1 = error-free rate.
    if Config.yield == 0
        % Generate the ideal displacement.
        theory = zeros(1, length(xs));
        theory(fs + 1) = 2;
        theory = cumsum(theory);

        yield = mean((xs - theory) .^ 2);
        yield = sqrt(yield);

    elseif Config.yield == 1
        yield = 0;

        % Backframeshifts are always an error.
        if length(model.rus) > 0
            return;
        elseif size(fs) == size(model.nut)
            if isempty(fs) || all(fs == model.nut)
                yield = 1;
            end
        end
    end
end
