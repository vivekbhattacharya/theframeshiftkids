% ----------------------------------------------------------------
% Stores yields to a results.txt file for perusal. `limit`
% specifies the number of iterations while `folder` is the folder
% path where `cornerstone.pl` dumped all its proteins. It can be
% any folder of txt files with protein sequences.
%
% Results file will ultimately end up in Matlab's work folder,
% appended.
%
% Usage: impunity(work folder, list of genes that should
%   frameshift a la megaunity, number of iterations to run
%   the model before recording yield)
% ----------------------------------------------------------------
function impunity(folder, frameshifters, limit)
d = dir([folder '\*.txt']); % Ignore .fasta files laying around

%% Print the header. %%
fid = fopen('results.txt', 'a+'); % Append
fwrite(fid, ['------ [' num2str(limit) ' iterations] ------']);
fprintf(fid, '\n');

for i = 1:length(d)
   if(strmatch(d(i).name, ['. ';'..'])), continue; end;
   
   % Num2str is required by fwrite.
   disp(['---------- [' d(i).name '] ----------']);
   yield = num2str(find_yield(d(i).name, frameshifters, limit));
   
   %% Print the data. %%
   fwrite(fid, [d(i).name '          ' yield]);
   fprintf(fid, '\n');
end
fprintf(fid, '\n\n\n');
fclose(fid);

% ----------------------------------------------------------------
% The meat of `impunity`, this is `megaunity` without the graphs
% or the infinite loop, instead capped at `limit`. In addition,
% it returns the yield instead of displaying it.
% ----------------------------------------------------------------
function [yield] = find_yield(file, frameshift_genes, limit)
[Signal, S] = get_signal(file);

global TAV Names;
load TAV.mat; load Codons.mat;

[Mag, Phase, numcodons] = calc_cumm_mag_phase(Signal);
[Dvec] = diff_vectors(Mag, Phase, numcodons);
global shoals sands beached_whale;
    % Disable verbosity with beached_whale
    shoals = 0; sands = 0; beached_whale = 1;
for i=1:limit
    [theta,x,diffx] = displacement(S(13:end),Phase,numcodons,Dvec,frameshift_genes);
end

yield = shoals/sands;