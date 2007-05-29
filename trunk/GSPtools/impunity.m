% ----------------------------------------------------------------
% Stores yields to a results.txt file for perusal. `limit`
% specifies the number of iterations while `folder` is the folder
% path where `cornerstone.pl` dumped all its proteins. It can be
% any folder of txt files with protein sequences.
%
% Results file will ultimately end up in Matlab's work folder,
% appended.
% ----------------------------------------------------------------
function impunity(folder, frameshift_genes, limit)
d = dir([folder '\*.txt']); % Ignore .fasta files laying around
fid = fopen('results.txt', 'a+'); % Append
fwrite(fid, ['------ [' limit ' iterations] ------']);

for i = 1:length(d)
   if(strmatch(d(i).name, ['. ';'..'])), continue; end;
   
   % Num2str is required by fwrite.
   disp(['---------- [' d(i).name '] ----------']);
   yield = num2str(find_yield(d(i).name, frameshift_genes, limit));
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
global shoals sands; shoals = 0; sands = 0;
for i=1:limit
    [theta,x,diffx] = displacement(S(13:end),Phase,numcodons,Dvec,frameshift_genes);
    fprintf('\n');
end

yield = shoals/sands;