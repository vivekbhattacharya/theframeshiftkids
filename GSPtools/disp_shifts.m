% Iterate over the cell array and print the results
% sanely. If only Matlab had a join function like
% every other modern language in the world.
function disp_shifts()
    function beached_whale(insects)
        if length(insects) ~= length({})
            for i=1:length(insects), fprintf([insects{i} '; ']); end;
        end
        fprintf('\n');
    end
    global ants termites;
    fprintf('> '); beached_whale(ants);
    fprintf('< '); beached_whale(termites);
    fprintf('\n');
end