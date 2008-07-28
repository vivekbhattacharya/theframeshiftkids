% Estimates species angle for a single gene. Connect to a for-loop or
% classify to run on a directory.
%
% Example: angle = species_angle('rpoS.txt')
%
% Example:
%  dirs = dir('c:/genes');
%  dirs(1:2) = []; % Remove . and ..
%  b = zeros(1, length(dirs));
%  for i = 1:length(dirs)
%    file = fullfile('c:/genes', dir(i).name);
%    b(i) = species_angle(file);
%  end
function [angle] = species_angle(gene)
    config;
    signal = get_signal(gene);

    N = length(signal);
    x = ones(N, 1);
    y = (0:N-1)';
    regressor = [x sin(2*pi*y/3) cos(2*pi*y/3)];
    [b, b_int, r] = regress(signal', regressor);

    angle = 180/pi*atan2(b(3),b(2));
end
