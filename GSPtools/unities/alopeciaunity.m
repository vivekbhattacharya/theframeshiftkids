% ------------------------------------------------
% Polar points become saved inside your folder.
% 
% Usage: alopeciaunity('C:\folder')
% ------------------------------------------------
function alopeciaunity(folder)
    d = [dir([folder '/*.txt']); dir([folder '/*.fasta'])];
    subfolder = 'merryround';
    mkdir(fullfile(folder, subfolder));
    
    for i = 1:length(d)
        file = d(i).name; disp(file);
        h = figure(2); set(h, 'Renderer', 'OpenGL');
            set(h, 'Visible', 'off');
            walrus_surprise(fullfile(folder, file), 'polar');
            set(h, 'Visible', 'off');
        saveas(h, fullfile(folder, subfolder, [file '.png']), 'png');
    end
    disp(' ');
end