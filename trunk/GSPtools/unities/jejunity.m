% ------------------------------------------------
% Plots errorbars for displacmenet and saves it
% to a file for all files in a given folder where
% each file contains a gene sequence
% 
% Usage: jejunity('C:\folder')
% ------------------------------------------------
function jejunity(folder, limit)
d = [dir([folder '/*.txt']); dir([folder '/*.fasta'])];
subfolder = 'carnivale';
mkdir(fullfile(folder, subfolder));

for i = 1:length(d)
   disp(d(i).name);
   hyperplot(fullfile(folder, d(i).name), subfolder, limit, 'errorbars');
end