function fabio(freq, names, id)
% Vivek doesn't want you to use this program.
%
% Given a list of frequencies generated from
% the third link on the search for codon bias,
% Codons.mat in a text file, and the file name
% for the .mat file fabio will generate, fabio
% will do as you wish.
dire = fileparts(which('Smooth.pm'));
[S, S] = system(sprintf('perl -I"%s" "%s" "%s" "%s"', dire, which('fabio.pl'), freq, names));
TAV = str2num(S);
save(id, 'TAV');