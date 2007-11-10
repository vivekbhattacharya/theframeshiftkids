% ----------------------------------------------------------------
% Impunity.m calculates the yield for all the genes in a given
% work folder, storing those results to a results.txt file in the
% Matlab workspace and also printing them to the standard output.
% It needs the predicted frameshifts in order to calculate yield
% and also a sample size.
%
% Usage: impunity('c:\work folder', [25], 100)
% ----------------------------------------------------------------
function impunity(folder, fshifts, limit)
    classify(folder, 'impunity', @helper);
    function helper(displacement, n, file, image)
        global shoals sands;
        yields = zeros(1, limit);
        for i=1:limit
            displacement(fshifts);
            yields(i) = shoals/sands;
            shoals = 0; sands = 0;
        end
        fprintf('%s: %g +/- %g\n', file, mean(yields), std(yields));
    end
end