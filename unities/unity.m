function unity(file, varargin)
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

displacement = walrus_surprise(file, 'polar');
global shoals;

fshifts = [];
if length(varargin) > 0, fshifts = varargin{1}; end
x = displacement(fshifts);

disp_shifts;
disp(sprintf('Yield: %g', shoals));

figure(1); plot(0,0); plot(1:length(x), x);
    axis([1 length(x) min(0, min(x)) max(3, max(x))]);
    grid; xlabel('Codon Number'); ylabel('x(k)');