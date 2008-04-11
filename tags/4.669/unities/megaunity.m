function megaunity(file, varargin)
% --------------------------------------------------------------
% This function is used to simulate repeated calls to unity
% in order to assess the "yield" of a given sequence (file)
% based on the deviation of the displacement from the predicted
% frame, which is usually 0.
%
%
% USAGE:
%   megaunity('rpoS.txt');
%   megaunity('prfB.txt', {'uga,25'});
%   % 300 iterations
%   megaunity('prfB.txt', {'uga,25'}, 300);
%
%   % Display the graph each time.
%   megaunity('prfB.txt', {'uga,25'}, 300, 'graph');
% --------------------------------------------------------------

displacement = walrus_surprise(file);
global shoals sands;

x = length(varargin);
fshifts = []; limit = inf; quiet = 1;

if x >= 1, fshifts = varargin{1}; end;
if x >= 2, limit = varargin{2}; end;
if x >= 3, quiet = 0; end;

for i = 1:limit
    x = displacement(fshifts);
    disp_shifts;
    fprintf('Yield: %g (%g)\n\n', shoals/sands, sands);

    if quiet, continue; end;
    h = figure(1);
        plot(1:length(x), x);
        axis([xlim min(0, min(x)) max(3, max(x))]);
        grid; xlabel('Codon'); ylabel('Displacement');
end
