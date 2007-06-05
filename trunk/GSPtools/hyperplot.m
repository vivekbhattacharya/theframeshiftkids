% ------------------------------------------------
% Plots 5 iterations simultaneously on one graph
% and saves it to a file in a work folder
% 
% Usage: hyperplot('C:\folder\prfb.txt',
%   'pictures', 5)
% ------------------------------------------------
function hyperplot(file, subfolder, limit, mode)
%%%%% Same old, same old
[Signal, S] = get_signal(file);
global TAV Names;
load TAV.mat; load Codons.mat;

[Mag, Phase, numcodons] = calc_cumm_mag_phase(Signal);
[Dvec] = diff_vectors(Mag, Phase, numcodons);

% Disable verbosity with beached_whale
global beached_whale; beached_whale = 1;
%%%%% Same old ends right now

if strcmp(mode, 'superimpose')
    plots = cell(1, limit*2);
    for i=1:limit
        [theta,x,diffx] = displacement(S(13:end),Phase,numcodons,Dvec,{});
        plots(i*2 - 1) = {1:length(x)};
        plots(i*2) = {x};
    end

    [folder, file, ext] = fileparts(file);
    h = figure(1); set(h, 'Renderer', 'OpenGL');
        set(h, 'Visible', 'off');
        plot(0,0); plot(plots{:}); title(file);
            grid; xlabel('Codon Number'); ylabel('x(k)');
    saveas(h, fullfile(folder, subfolder, [file '.png']), 'png');
elseif strcmp(mode, 'errorbars')
    a = cell(1, limit);
    for i=1:limit
        [theta,x,diffx] = displacement(S(13:end),Phase,numcodons,Dvec,{});
        a(i) = {x};
    end
	upper = size(a{1}, 2) * 1.5;
    b = zeros(1, upper); c = zeros(1, upper);
    for i=1:limit % Which array?
        y = a{i};
        for j=1:size(y,2) % Which element?
			b(j) = b(j) + y(j);
			c(j) = c(j) + 1;
        end
    end
	avg = b ./ c;
	h = figure(1); set(h, 'Renderer', 'OpenGL');
	plot(0,0); plot(1:length(x), avg(1:length(x))); title(file);
		grid; xlabel('Codon Number'); ylabel('x(k)');
end