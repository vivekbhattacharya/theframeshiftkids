% --------------------------------------------------------------
% This function will turn Recode or Genbank gene sequences
% (capital letters, DNA) with line numbers removed into a
% toolbox-compatible format (FASTA).
%
% It was created by The Frameshift Kids to compliment Dr. Lalit
% Ponnala's toolbox.
%
% USAGE:
% function [] = tofasta(file)
% file = file with a Recode/Genbank sequence
% --------------------------------------------------------------

function return_file = tofasta(file)

file = which(file);
if(isempty(file)); error(['Cannot find ``' file '" in your path']); end;

% Load prfb, and generate a fasta file with column width of 60.
S = getseq(file);
S = num2char(char2num(S, 1, 0), 0, 1);
write2fasta(file, S, 'candy', 60);

return_file = file;
disp(['Now please edit `' file '` to remove the >candy.']);