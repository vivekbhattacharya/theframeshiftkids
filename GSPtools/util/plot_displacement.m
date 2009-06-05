% Plots displacement in a nice graph.
function plot_displacement(x, varargin)
    narg = nargin - 1;

    switch(narg)
      case 0
        figure;
      case 1
        figure(varargin{1});
      otherwise
        error('plot_displacement: too many arguments');
    end

    plot(x, 'LineWidth', 2);
    axis([1 length(x) min(0, min(x)) max(3, max(x))]);
    grid;
    xlabel('Codon Number');
    ylabel('Displacement');
end
