% This function repeatedly measures the yield of a given sequence as
% defined by Config.yield by acting as if unity is being called
% multiple times.
%
% USAGE:
%   megaunity('rpoS.txt');
%   megaunity('prfB.txt', [25])
%   % 300 iterations
%   megaunity('prfB.txt', [25], 300);
%   % Display the graph each time.
%   megaunity('prfB.txt', [25], 300, 'graph');
function megaunity(file, varargin)
    config;
    global Config;

    x = nargin;
    shift = [];
    limit = inf;
    quiet = 1;

    if x >= 1, shift = varargin{1}; end;
    if x >= 2, limit = varargin{2}; end;
    if x >= 3, quiet = 0; end;

    m = Config.model(file, shift);
    yields = [];

    for i = 1:limit
        [bye, x] = displacement(m);
        yields = [yields yield(bye, x)];
        disp_shifts(bye);

        fprintf('Yield: %g (mean:%g) (std:%g) (n:%g)\n\n', ...
                yields(i), mean(yields), std(yields), i);

        if quiet
            continue;
        end

        plot_displacement(x, 2000);
    end
end
