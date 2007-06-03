function opportunity(folder)
d = dir([folder '\*.txt']); % Ignore .fasta files laying around
subfolder = 'thepictureshow';
mkdir(fullfile(folder, subfolder));

for i = 1:length(d)
   disp(['---------- [' d(i).name '] ----------']);
   hyperplot(fullfile(folder, d(i).name), subfolder, 5);
end

% ------------------------------------------------
% Plots 5 iterations simultaneously on one graph
% and saves it to a file in a work folder
% ------------------------------------------------
function hyperplot(file, subfolder, limit)

%%% Same old, same old
[Signal, S] = get_signal(file);
global TAV Names;
load TAV.mat; load Codons.mat;

[Mag, Phase, numcodons] = calc_cumm_mag_phase(Signal);
[Dvec] = diff_vectors(Mag, Phase, numcodons);
global beached_whale;
    % Disable verbosity with beached_whale
    beached_whale = 1;
%%% Same old ends right now

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