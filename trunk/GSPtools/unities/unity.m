% Plots the ribosomal displacement over codons given the file name of
% an mRNA sequence.
%
% > unity('prfB.txt')
function unity(file, varargin)
    config;
    global Config;

    shifts = [];
    if length(varargin) > 0
        shifts = varargin{1};
    end

    m = Config.model(file, shifts);
    [m, x] = displacement(m);
    disp_shifts(m);
    fprintf('Yield: %g\n', yield(m, x));

    plot_displacement(x);
end
