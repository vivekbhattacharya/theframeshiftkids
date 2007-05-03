function unity(file)
% --------------------------------------------------------------
% This function is used to unify the elements in this toolbox so
% that an mRNA sequence in a text file with a 12-character
% leader sequence as a string can output a displacement plot 
% and a polar plot using calcmpx.
%
% It was created by The Frameshift Kids to complement Dr. Lalit
% Ponnala's toolbox.  It uses a modified version of Dr. Josh
% Starmer's free2bind program's free_scan.pl.  In the for loop
% at line 490, just add `&& 0' (without quotation marks, of
% course) to the conditions. This ensures that the output file
% will only be one line long.
%
% USAGE:
% function [] = unity(file)
% file = string with the name of the text file, with extension,
%        that includes the leader
% --------------------------------------------------------------

file = which(file);
if(isempty(file)); error(['Cannot find ``' file '" in your path']); end;

% Load prfb, and generate a fasta file with column width of 60.
S = getseq(file); Fasta = [file, '.fasta'];
write2fasta(Fasta, S, 'prfb', 60);

% Run free2bind. Make sure to include its directory into -I.
if(isempty('perl.exe')); error(['I cannot find Perl 5.6 or above.' ...
   'Please download it to the VCL C: drive.']); end;
   Include = fileparts(which('FreeAlign.pm'));
   Template = 'perl.exe -I"%s" "%s" -r -e -q -p FREIER auuccuccacuag "%s"';
   Command = sprintf(Template, Include, which('free_scan.pl'), Fasta);
[status, Signal] = dos(Command);

% Simulate load() on a string instead of a file.
Signal = str2num(Signal);
if(isempty(Signal))
    ensure = 'I cannot pull signals. Ensure `perl.exe` is outputting the rite stuff.';
    ensure = [ensure '\nAlso, ensure Perl is of version 5.6 or above.'];
    error(ensure);
end

% Magic numbers.


%%% demo4 code %%%
% These files need to be in the include path or working directory
% (GSPdemos).
global TAV Names;
load TAV.mat; load Codons.mat;

[Mag, Phase, numcodons] = calc_cumm_mag_phase(Signal);
[Dvec] = diff_vectors(Mag, Phase, numcodons);
[theta,x] = displacement(S(13:end),1000,1,Phase,numcodons,Dvec);

for k=1:length(theta) % No negative values
    if theta(k)<0; theta(k)=theta(k)+(2*pi); end
end
for k=1:length(x)-1
    diffx(k)=x(k+1)-x(k);
end

cp = 0;
figure(1);
	subplot(211);plot(0,0);plot(1+cp:length(x)+cp, x);
    	axis([1 length(x)+cp min(0,min(x)) max(3,max(x))]);
    	grid; xlabel('Codon Number'); ylabel('x(k)');    
    subplot(212); plot(0,0);plot(1:length(diffx),diffx); xlabel('Codon number'); ylabel('Force on ribosome'); title('Plot of "force", i.e. incremental displacement');	
% figure(2);
%     subplot(211); plot(0,0);polar(theta,mag);
%     subplot(212);plot(0,0);plot(1:length(theta),(180/pi)*theta); xlabel('Codon number'); ylabel('Phase angle (degrees)'); title('Plot of cumulative phase');

% The toolbox code in demo4 need a newline.
fprintf('\n');