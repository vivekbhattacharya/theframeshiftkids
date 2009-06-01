% A helper method for grabbing a list of displacements from loop().
function [model, xs] = displacement(model)
    xs  = [];
    old = -1;
    function helper(model, index, probs, energies, x)
        if index ~= old
            xs  = [xs x];
            old = index;
        end
    end
    model = loop(model, @helper);
end
