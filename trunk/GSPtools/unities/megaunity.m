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

% These files need to be in the include path or working directory.
% See code.google.com website for copies.
[S, n, Dvec] = walrus_surprise(file);
global shoals sands;
while 1
    [x,diffx] = displacement(S(13:end),n,Dvec,fshifts,bshifts);

    cp = 0;
    h = figure(1); set(h, 'Renderer', 'OpenGL');
        plot(0,0);plot(1+cp:length(x)+cp, x);
        axis([1 length(x)+cp min(0,min(x)) max(3,max(x))]);
        grid; xlabel('Codon Number'); ylabel('x(k)');    
    disp(['Yield so far: ' num2str(shoals/sands) ' (' num2str(sands) ')']);
    fprintf('\n');
end