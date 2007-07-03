% ------------------------------------------------
% Superimposes 5 iterations and saves that plot
% to a file for all files in a given folder where
% each file contains a gene sequence
% 
% Usage: opportunity('C:\folder')
% ------------------------------------------------
function opportunity(folder, limit)
d = [dir([folder '/*.txt']); dir([folder '/*.fasta'])];
subfolder = 'thepictureshow';
mkdir(fullfile(folder, subfolder));

for i = 1:length(d)
   disp(['---------- [' d(i).name '] ----------']);
   hyperplot(fullfile(folder, d(i).name), subfolder, limit, 'superimpose');
end