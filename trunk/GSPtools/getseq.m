% Usage: S = getseq(file, fasta)
%
% File may contain a FASTA header or spaces. It may either be
% DNA or RNA. Function returns the sequence in mRNA format,
% lowercase and otherwise sanitized perfection. In addition,

function S = getseq(file)
% Matlab does not sync well with system, so I aid
% its module-finding skills w/r/t Smooth.pm.
S = pearl('getseq.pl', ['"' file '"']);
