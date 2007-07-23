function megaunity(file, fshifts, bshifts)
% --------------------------------------------------------------
% This function is used to simulate repeated calls to unity
% in order to assess the "yield" of a given sequence (file)
% based on whether or not the stochastic model frameshifts
% completely correctly. --Frameshift Kids
%
%
% USAGE:
% function [] = megaunity(file, str1)
% file = string with the name of the text file, with extension,
%        that includes the leader
% frameshift_genes = string with list of desired frameshifts
%        ex. ['uga,25';'cau,73'] or '[]'
% --------------------------------------------------------------

displacement = walrus_surprise(file, 'polar');
global shoals sands;
while 1
    x = displacement(fshifts, bshifts);
    disp_shifts;

    h = figure(1); set(h, 'Renderer', 'OpenGL');
        plot(0,0); plot(1:length(x), x);
        axis([1 length(x) min(0,min(x)) max(3,max(x))]);
        grid; xlabel('Codon Number'); ylabel('x(k)');    
    disp(sprintf('Yield so far: %g (%g)', shoals/sands, sands));
    fprintf('\n');
end