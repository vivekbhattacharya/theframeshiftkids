function megaunity(file, frameshift_genes_forward, frameshift_genes_back)
% --------------------------------------------------------------
% This function is used to simulate repeated calls to unity
% in order to assess the "yield" of a given sequence (file)
% based on whether or not the stochastic model frameshifts
% completely correctly. --Frameshift Kids
%
%
% USAGE:
% function [] = megaunity(file, str1)
% file = string with the name of the text file, with extension,
%        that includes the leader
% frameshift_genes = string with list of desired frameshifts
%        ex. ['uga,25';'cau,73'] or '[]'
% --------------------------------------------------------------

[Signal, S] = get_signal(file);

% These files need to be in the include path or working directory
% (GSPdemos).
global TAV Names;
load TAV.mat; load Codons.mat;

[Mag, Phase, numcodons] = calc_cumm_mag_phase(Signal);
[Dvec] = diff_vectors(Mag, Phase, numcodons);
global shoals sands; shoals = 0; sands = 0;
while 1
    [theta,x,diffx] = displacement(S(13:end),Phase,numcodons,Dvec,frameshift_genes_forward, frameshift_genes_back);

    cp = 0;
    h = figure(1); set(h, 'Renderer', 'OpenGL');
        subplot(211);plot(0,0);plot(1+cp:length(x)+cp, x);
            axis([1 length(x)+cp min(0,min(x)) max(3,max(x))]);
            grid; xlabel('Codon Number'); ylabel('x(k)');    
        subplot(212); plot(0,0);plot(1:length(diffx),diffx);
            xlabel('Codon number'); ylabel('Force on ribosome');
            title('Plot of "force", i.e. incremental displacement');
    
    fprintf('\n');
    disp(['Yield so far: ' num2str(shoals/sands) ' (' num2str(sands) ')']);
    fprintf('\n');
end