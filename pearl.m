function S = pearl(file, args)
% Matlab does not sync well with system, so I aid
% its module-finding skills w/r/t Smooth.pm.
dire = fileparts(which('Smooth.pm'));
command = sprintf('perl -I"%s" "%s" %s', dire, which(file), args);
[status, S] = system(command);
