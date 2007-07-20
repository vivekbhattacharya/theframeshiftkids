% ------------------------------------------------
% Superimposes 5 iterations and saves that plot
% to a file for all files in a given folder where
% each file contains a gene sequence
% 
% Usage: opportunity('C:\folder', 10)
% ------------------------------------------------
function opportunity(folder, limit)
    d = [dir([folder '/*.txt']); dir([folder '/*.fasta'])];
    subfolder = 'thepictureshow';
    mkdir(fullfile(folder, subfolder));

    for i = 1:length(d)
       disp(d(i).name);
       helper(fullfile(folder, d(i).name), subfolder, limit);
    end
end

function helper(file, subfolder, limit)
    [S, n, Dvec] = walrus_surprise(file);

    plots = cell(1, limit*2);
    for i=1:limit
        x = displacement(S(13:end),n,Dvec,{},{});
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
end