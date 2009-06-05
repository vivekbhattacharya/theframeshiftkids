% yields calculates the yield for all the genes in a given work
% folder, storing those results to a results.txt file in the Matlab
% workspace and also printing them to the standard output. It needs
% the predicted frameshifts in order to calculate yield and also a
% sample size.
%
% Usage: yields('c:\work folder', 25, 100)
function yields(folder, fshift, limit)
    classify(folder, [], fshift, @helper);

    function helper(model, file)
        yields = zeros(1, limit);

        for i = 1:limit
            [model, x]  = displacement(model);
            yields(i)   = yield(model, x);
        end

        fprintf('%s: %g (std:%g)\n', file, mean(yields), std(yields));
    end
end
