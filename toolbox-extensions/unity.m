% --------------------------------------------------------------
% This function is used to unify the elements in this toolbox so
% that an mRNA sequence in a text file with a 12-character
% leader sequence as a string can output a displacement plot 
% and a polar plot using calcmpx.
%
% It was created by Hao Lian, Vivek Bhattacharya, and Daniel
% Vitek to compliment Dr. Lalit Ponnala's toolbox.  It uses a
% modified version of Dr. Josh Starmer's free2bind program's
% free_scan.pl.  In the for loop at line 490, just add `&& 0'
% (without quotation marks, of course) to the conditions. This
% ensures that the output file will only be one line long.
%
% USAGE:
% function [] = unity(file)
% file = string with the name of the text file, with extension,
%        that includes the leader
% --------------------------------------------------------------

function [] = unity(file)

file_not_there = isempty(which(file));
if(file_not_there)
  error(['Cannot find ' file ' in your path']);
end

% These files need to be in the include path or working directory (GSPdemos).
load TAV.mat; load Codons.mat;

% Load prfb and convert to lowercase whatever format.
S = getseq(file);

% Generate a fasta file for free2bind.
Fasta = strcat(file, '.fasta');
write2fasta(Fasta, S, 'prfb', 60);

% Run free2bind. Make sure to include its directory into -I.
Include = sprintf('-I"%s"', fileparts(which('FreeAlign.pm')));
[status, Signal] = system(sprintf('perl %s %s -r -e -q -p FREIER auuccuccacuag %s', Include, which('free_scan.pl'), Fasta));

% Simulate load() on a string instead of a file.
Signal = str2num(Signal);

if(isempty(Signal) && isempty(which('perl.exe')))
    error('I cannot pull signals. Ensure `perl.exe` is in your path (`dos(''path'')`).');
end

% Magic numbers.
cp = 0; phi_sp=-13*(pi/180); initialx = 0.001; C1 = 0.005; C2 = initialx; Nstop=1000; spc=1;

% For debugging.
% disp(length(Signal) - length(S));

%%% demo4 code %%%
% Run it through the model
[mag,theta,x] = calcmpx(S(13:end),Signal,phi_sp,Names,TAV,C1,C2,Nstop,spc);

% Polar plot
figure(1); polar(theta,mag);

% Displacement plot
figure(2); plot(1+cp:length(x)+cp, x);
axis([1 length(x)+cp min(-4,min(x)) max(4,max(x))]);
grid; xlabel('Codon Number'); ylabel('x(k)');

% Newline.
fprintf('\n');