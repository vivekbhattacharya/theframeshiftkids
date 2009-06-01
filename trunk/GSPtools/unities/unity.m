% Plots the ribosomal displacement over codons given the file name of
% an mRNA sequence.
%
% > unity('prfB.txt')
function unity(file, varargin)
    shifts = [];
    if length(varargin) > 0
        shifts = varargin{1};
    end

    m = model(file, shifts);
    [m, x] = displacement(m);
    disp_shifts(m);
    disp(sprintf('Yield: %g', yield(m, x)));

    figure;
    plot(x, 'LineWidth', 2);
    axis([1 length(x) min(0, min(x)) max(3, max(x))]);
    grid;
    xlabel('Codon Number');
    ylabel('Displacement');
end
