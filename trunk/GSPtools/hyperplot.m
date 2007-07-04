% ------------------------------------------------
% Plots 5 iterations simultaneously on one graph
% and saves it to a file in a work folder
% 
% Usage: hyperplot('C:\folder\prfb.txt',
%   'pictures', 5, 'superimpose')
% ------------------------------------------------
function hyperplot(file, subfolder, limit, mode)
[Signal, S] = get_signal(file);
[Mag, Phase, n] = cumm_mag_phase(Signal);
Dvec = diff_vectors(Mag, Phase, n);

% Disable verbosity with beached_whale
global beached_whale; beached_whale = 1;
if strcmp(mode, 'superimpose')
    plots = cell(1, limit*2);
    for i=1:limit
        [x,diffx] = displacement(S(13:end),Phase,n,Dvec,{},{});
        plots(i*2 - 1) = {1:length(x)};
        plots(i*2) = {x};
    end
    
    top = 3; bottom = -3;
    if max(x) > 3, top = ceil(max(x)); end;
    if min(x) < -3, bottom = floor(min(x)); end;
    
    max_dom = 50*(ceil(length(x)/50));

    [folder, file, ext] = fileparts(file);
    h = figure(1); set(h, 'Renderer', 'OpenGL');
        set(h, 'Visible', 'off');
        plot(0,0); plot(plots{:}); title(file);
            grid; xlabel('Codon Number'); ylabel('x(k)');
            axis([0, max_dom, bottom, top]);
    saveas(h, fullfile(folder, subfolder, [file '.png']), 'png');
elseif strcmp(mode, 'errorbars')
	upper = ceil(n);
    a = zeros(limit, upper);
    for i=1:limit
        [x,diffx] = displacement(S(13:end),Phase,n,Dvec,{},{});
		maxiderm = size(x, 2);
        for j=1:upper
			if j <= maxiderm, a(i,j) = x(j);
			else a(i,j) = -3.14;
			end;
		end
    end
	[folder, file, ext] = fileparts(file);
	[avg, stddev] = exmean(a);
	h = figure(1); set(h, 'Renderer', 'OpenGL');
        set(h, 'Visible', 'off');
        plot(0,0); errorbar(1:upper, avg, stddev); title(file);
        grid; xlabel('Codon Number'); ylabel('x(k)');
    saveas(h, fullfile(folder, subfolder, [file '.png']), 'png');
end