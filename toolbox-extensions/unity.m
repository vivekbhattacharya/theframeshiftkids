                                                                     
                                                                     
                                                                     
                                             
% --------------------------------------------------------------
% This function is used to unify the elements in this toolbox so
% that an mRNA sequence in a text file along with a 12-character
% leader sequence as a string.  It outputs a displacement plot 
% and a polar plot using calcmpx.
%
% It was created by Hao Lian, Vivek Bhattacharya, and Daniel
% Vitek to compliment Dr. Lalit Ponnala's toolbox.  It uses a
% modified version of Dr. Josh Starmer's free2bind program's
% free_scan.pl.  In the for loop at line 490, just add '&& 0'
% (without the quotations, of course) to the conditions. This
% ensures that the output file will only be one line long.
%
% USAGE:
% function [] = unity(file,leader)
% file = string with the name of the text file, with extension
% leader = 12-character string of the mRNA leader
% --------------------------------------------------------------

function [] = unity(file, leader)

if(isempty(which(file))
  message = strcat('Cannot find ',file,' in your path');
  error(message);
end

% These files need to be in the include path or working directory (GSPdemos).
load TAV.mat; load Codons.mat;

% Load prfb and convert to lowercase whatever format.
S = getseq(file);
% try
%     S = num2char(char2num(S, 1, 0), 0, 1);
% catch
%     error('I cannot run char2num. Ensure your current working directory is correct (e.g. K:).');
% end

% Generate a fasta file for free2bind.
Fasta = strcat(file, '.fasta');
write2fasta(Fasta, [S], 'prfb', 60);

% Run free2bind. Make sure to include its directory into -I.
Include = sprintf('-I"%s"', fileparts(which('FreeAlign.pm')));
[status, Signal] = system(sprintf('perl %s %s -r -e -q -p FREIER auuccuccacuag %s', Include, which('free_scan.pl'), Fasta));
% fid = fopen('signals.txt', 'w')
% fwrite(fid, Signal);
% fclose(fid);

% Simulate load() on a string instead of a file.
Signal = str2num(Signal);

if(isempty(Signal) && isempty(which('perl.exe')))
    E = 'I cannot pull signals. Ensure `perl.exe` is in your path (`dos(''path'')`).';
    error(E);
end

% Magic numbers.
cp = 0; phi_sp=-13*(pi/180); initialx = 0.001; C1 = 0.005; C2 = initialx; Nstop=1000; spc=1;

% For debugging.
% disp(length(Signal) - length(S));

%%% demo4 code %%%
% Run it through the model
[mag,theta,x] = calcmpx(S(length(leader)+1:end),Signal,phi_sp,Names,TAV,C1,C2,Nstop,spc);

% Polar plot
figure(1); polar(theta,mag);

% Displacement plot
figure(2); plot(1+cp:length(x)+cp, x);
axis([1 length(x)+cp min(-4,min(x)) max(4,max(x))]);
grid; xlabel('Codon Number'); ylabel('x(k)');

% Newline.
fprintf('\n');