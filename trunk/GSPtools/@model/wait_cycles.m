% Calculate normalized wait cycles vector (one for each {-1, 0, +1}
% frame) for a given nucleotide base.
function [loops] = wait_cycles(model, base)
    x = (base + 0):(base + 2);
    a = model.seq(x);
    b = model.seq(x + 1);
    c = model.seq(x + 2);

    loops = [model.travel.(a) model.travel.(b) model.travel.(c)];
    loops = 2 .^ (1 ./ ceil(loops));
    loops = loops ./ (loops - 1);
end
