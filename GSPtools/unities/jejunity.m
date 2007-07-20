% ------------------------------------------------
% Plots errorbars for displacmenet and saves it
% to a file for all files in a given folder where
% each file contains a gene sequence
% 
% Usage: jejunity('C:\folder', 10)
% ------------------------------------------------
function jejunity(folder, limit)
    d = [dir([folder '/*.txt']); dir([folder '/*.fasta'])];
    subfolder = 'carnivale';
    mkdir(fullfile(folder, subfolder));

    for i = 1:length(d)
       disp(d(i).name);
       helper(fullfile(folder, d(i).name), subfolder, limit);
    end
end

function helper(file, subfolder, limit, mode)
    [S, n, Dvec] = walrus_surprise(file);

	upper = ceil(n);
    a = zeros(limit, upper);
    for i=1:limit
        x = displacement(S(13:end),n,Dvec,{},{});
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