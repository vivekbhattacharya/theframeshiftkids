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

[Signal, S] = get_signal(file);

% These files need to be in the include path or working directory
% (GSPdemos).
global TAV Names;
load TAV.mat; load Codons.mat;

[Mag, Phase, numcodons] = cumm_mag_phase(Signal);
[Dvec] = diff_vectors(Mag, Phase, numcodons);
[theta,x,diffx] = displacement(S(13:end),Phase,numcodons,Dvec,{},{});

cp = 0;
figure(1);
	subplot(211);plot(0,0);plot(1+cp:length(x)+cp, x);
    	axis([1 length(x)+cp min(0,min(x)) max(3,max(x))]);
    	grid; xlabel('Codon Number'); ylabel('x(k)');    
    subplot(212); plot(0,0);plot(1:length(diffx),diffx);
        xlabel('Codon number'); ylabel('Force on ribosome');
        title('Plot of "force", i.e. incremental displacement');	
% figure(2);
%     subplot(211); plot(0,0);polar(theta,mag);
%     subplot(212);plot(0,0);plot(1:length(theta),(180/pi)*theta);
%     xlabel('Codon number'); ylabel('Phase angle (degrees)');
%     title('Plot of cumulative phase');    

% The toolbox code in demo4 need a newline.
fprintf('\n');