% Vivek doesn't want you to use this program

function fabio(freq, names, id)
dire = fileparts(which('Smooth.pm'));
[S, S] = system(sprintf('perl -I"%s" "%s" "%s" "%s"', dire, which('fabio.pl'), freq, names));
TAV = str2num(S)
%save(id, 'TAV');