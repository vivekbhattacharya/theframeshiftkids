function unity(file)
% --------------------------------------------------------------
% This function is used to unify the elements in this toolbox so
% that an mRNA sequence in a text file with a 12-character
% leader sequence as a string can output a displacement plot 
% and a polar plot using calcmpx.
%
% It was created by The Frameshift Kids to complement Dr. Lalit
% Ponnala's toolbox.
%
% USAGE:
% function [] = unity(file)
% file = string with the name of the text file, with extension,
%        that includes the leader
% --------------------------------------------------------------

% These files need to be in the include path or working directory.
% See code.google.com website for copies.
[S, n, Dvec] = walrus_surprise(file);
x = displacement(S(13:end),n,Dvec,{},{});

disp_shifts;

figure(1); plot(0,0);plot(1:length(x), x);
    axis([1 length(x) min(0,min(x)) max(3,max(x))]);
    grid; xlabel('Codon Number'); ylabel('x(k)');