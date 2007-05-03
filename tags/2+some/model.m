% Performs model-related tasks
clear; clc;

% NOTES:
% 16S tail = 'auuccuccacuag' 
% Use RECODE/GENBANK prfB sequence


% % ----------PART 1----------
% % % ----Link Genes----
% % load LinkGeneDataFile.mat; LINKSIGFILE = 'LinkGeneSignals.txt'; sigpreset = 0; sigpostset = 0;
% % fprintf('\nChecking if Link signals are written-out correctly...');
% % I = checkwrsigs(LinkGeneData,LINKSIGFILE,sigpreset,sigpostset);
% % if ~isempty(I)
% %     fprintf('\nLink signals having the following indices are not written-out correctly:\n'); disp(I); 
% %     error('\n!!Error in writing-out Link signals!!');
% % end
% % fprintf('\nPicking good genes among Link genes...');
% % minL = 600; alpha = 0.05; showpgram = 0;
% % GG = pickgg(LinkGeneData,LINKSIGFILE,sigpreset,sigpostset,minL,alpha,showpgram); 
% % save LinkGG.mat GG;
% %
% % fprintf('\nPlotting histograms...');
% % LineNum = GG(:,1); Length = GG(:,2); mu = GG(:,3); Amp = GG(:,4); Phase = GG(:,5); SNR = GG(:,6);
% % fprintf('\nPlotting histogram of phase...'); 
% % figure; hist(Phase,100); xlabel('Phase (degrees)'); ylabel('Number of genes');
% % fprintf('\n [mean(Phase), var(Phase), std(Phase)] = [%f\t %f\t %f]',mean(Phase),var(Phase),std(Phase)); 
% % fprintf('\nPlotting histogram of SNR...');
% % figure; hist(SNR,100); xlabel('SNR (dB)'); ylabel('Number of genes');
% % fprintf('\n [mean(SNR), var(SNR), std(SNR)] = [%f\t %f\t %f]',mean(SNR),var(SNR),std(SNR)); 
%
% % ----Long Genes----
% load LongGeneDataFile.mat; LONGSIGFILE = 'LongGeneSignals.txt'; sigpreset = 0; sigpostset = 0;
% fprintf('\nChecking if Long signals are written-out correctly...');
% I = checkwrsigs(LongGeneData,LONGSIGFILE,sigpreset,sigpostset);
% if ~isempty(I)
%     fprintf('\nLong signals having the following indices are not written-out correctly:\n'); disp(I);
%     error('\n!!Error in writing-out Long signals!!'); 
% end
% fprintf('\nPicking good genes among Long genes...');
% minL = 600; alpha = 0.05; showpgram = 0;
% GG = pickgg(LongGeneData,LONGSIGFILE,sigpreset,sigpostset,minL,alpha,showpgram);
% save LongGG.mat GG;
%
% fprintf('\nPlotting histograms...');
% LineNum = GG(:,1); Length = GG(:,2); mu = GG(:,3); Amp = GG(:,4); Phase = GG(:,5); SNR = GG(:,6);
% fprintf('\nPlotting histogram of phase...');
% figure; hist(Phase,100); xlabel('Phase (degrees)'); ylabel('Number of genes');
% fprintf('\n [mean(Phase), var(Phase), std(Phase)] = [%f\t %f\t %f]',mean(Phase),var(Phase),std(Phase));
% fprintf('\nPlotting histogram of SNR...'); 
% figure; hist(SNR,100); xlabel('SNR (dB)'); ylabel('Number of genes');
% fprintf('\n [mean(SNR), var(SNR), std(SNR)] = [%f\t %f\t %f]',mean(SNR),var(SNR),std(SNR));


% % ----------PART 2---------- 
% % Estimate tRNA availability for each codon
% % ****Processes 180,000 codons per minute****
% SEQFILE='CertSeqs.txt'; Case=0; Type=0; seqpreset=0; seqpostset=0;
% NamesFile='Codons.mat'; CARFile=' CAR.mat';
% TAVFile='TAV.mat'; tRNAFile='tRNAs.txt';
% calctav(SEQFILE,Case,Type,seqpreset,seqpostset,NamesFile,CARFile,TAVFile,tRNAFile);


% ----------PART 3----------

% ----prfB---- 
fprintf('\n\nAnalyzing prfB...');

% % Use RECODE's prfB sequence and signal
% Seq=readfasta('prfB_seq_RECODE.fasta');
% Signal=load('prfB_signal_RECODE.txt');
% % Signal=load('prfB_FREIER_signal_RECODE.txt'); 
% cp=4; % codon preset

% Use GENBANK's prfB sequence and signal
Seq=readfasta('prfB_seq_GENBANK.fasta'); seqpreset=12; seqpostset=0;
% Seq=readfasta('temp_prfB.fasta'); seqpreset=12; seqpostset=0; 
Seq=Seq(seqpreset+1:end-seqpostset);
% Signal=load('prfB_signal_GENBANK.txt');
Signal=load('prfB_FREIER_signal_GENBANK.txt');
% Signal=load('temp_prfB_Signal.txt');
cp=0; % codon preset 

% Set C1, initialx, Nstop and phi_sp (keep changing these until the state
% transition looks good)
% spc=0, phi_sp=-13*(pi/180); spc=1, phi_sp=-13*(pi/180); spc=2, phi_sp=-50*(pi/180);
load TAV.mat; load Codons.mat;
phi_sp=-13*(pi/180); initialx = 0.001; C1 = 0.005; C2 = initialx; Nstop=1000; spc=1;
% Run it through the model
[mag,theta,x] = calcmpx(Seq(3*cp+1:end),Signal,phi_sp,Names,TAV,C1,C2,Nstop,spc);
figure(2); polar(theta,mag); 
figure(3); plot(1+cp:length(x)+cp, x); axis([1 length(x)+cp min(-4,min(x)) max(4,max(x))]); grid; xlabel('Codon number k'); ylabel('Displacement x(k)');
figure(4); plot(1:length(theta),(180/pi)*theta);

% % ----aceF----
% fprintf('\n\nAnalyzing aceF...');
%
% % Use GENBANK's aceF sequence and signal
% load LinkGeneDataFile.mat ; SEQFILE='LinkGeneSeqs.txt'; seqpreset=0; seqpostset=0;
% SIGFILE='LinkGeneSignals.txt'; sigpreset=0; sigpostset=0;
% [Seq,Signal] = getseqsig1('aceF',LinkGeneData,SEQFILE,seqpreset,seqpostset,SIGFILE,sigpreset,sigpostset); 
% cp=0; % codon preset
%
% % Set C1, initialx, Nstop and phi_sp (keep changing these until the state
% % transition looks good)
% load TAV.mat; load Codons.mat;
% phi_sp=-13*(pi/180); initialx = 0.001; C1 = 0.005; C2 = initialx; Nstop=1000; spc=0;
% % Run it through the model
% [mag,theta,x] = fcalcmpx(Seq(3*cp+1:end),Signal,phi_sp,Names,TAV,C1,C2,Nstop,spc);
% figure(1); polar(theta,mag);
% figure(2); plot(1+cp:length(x)+cp, x); axis([1 length(x)+cp min(-4,min(x)) max(4,max(x))]); grid; xlabel('Codon number k'); ylabel('Displacement x(k)');

% % ----Link Genes----
% fprintf('\n\nAnalyzing Link Genes...');
% LINKSEQFILE='LinkGeneSeqs.txt'; seqpreset=0; seqpostset=0; LINKSIGFILE='LinkGeneSignals.txt '; sigpreset=0; sigpostset=0;
% load LinkGeneDataFile.mat; load TAV.mat; load Codons.mat; SelInd=[1:size(LinkGeneData,1)]';
% phi_sp=-13*(pi/180); initialx = 0.001; C1 = 0.005; C2 = initialx; Nstop=1000;
% % Run function to evaluate results
% spc=0; display=1;
% [Npass,Ntotal] = fsetfsparams(LINKSEQFILE,seqpreset,seqpostset,LINKSIGFILE,sigpreset,sigpostset,SelInd,LinkGeneData,phi_sp,Names,TAV,C1,C2,Nstop,spc,display); 
% fprintf('\nNumber of genes analyzed = %d, passed = %d',Ntotal,Npass);

% % ----Link Genes that don't work----
% fprintf('\n\nAnalyzing Link Genes that dont work...');
% LINKSEQFILE='LinkGeneSeqsFail.txt'; seqpreset=0; seqpostset=0; LINKSIGFILE='LinkGeneSignalsFail.txt'; sigpreset=0; sigpostset=0;
% load LinkGeneDataFailFile.mat; load TAV.mat; load Codons.mat; SelInd=[1:size(LinkGeneDataFail,1)]';
% phi_sp=-13*(pi/180); initialx = 0.001; C1 = 0.005; C2 = initialx; Nstop=1000;
% % Run function to evaluate results
% spc=0; display=1;
% [Npass,Ntotal] = fsetfsparams(LINKSEQFILE,seqpreset,seqpostset,LINKSIGFILE,sigpreset,sigpostset,SelInd,LinkGeneDataFail,phi_sp,Names,TAV,C1,C2,Nstop,spc,display); 
% fprintf('\nNumber of genes analyzed = %d, passed = %d',Ntotal,Npass);

% % ----Long Genes----
% fprintf('\n\nAnalyzing Long Genes...');
% LONGSEQFILE='LongGeneSeqs.txt'; seqpreset=0; seqpostset=0; LONGSIGFILE=' LongGeneSignals.txt'; sigpreset=0; sigpostset=0;
% load LongGeneDataFile.mat; load TAV.mat; load Codons.mat; SelInd=[1:size(LongGeneData,1)]';
% phi_sp=-13*(pi/180); initialx = 0.001; C1 = 0.005; C2 = initialx; Nstop=1000; 
% % Run function to evaluate results
% spc=0; display=0;
% [Npass,Ntotal] = fsetfsparams(LONGSEQFILE,seqpreset,seqpostset,LONGSIGFILE,sigpreset,sigpostset,SelInd,LongGeneData,phi_sp,Names,TAV,C1,C2,Nstop,spc,display); 
% fprintf('\nNumber of genes analyzed = %d, passed = %d',Ntotal,Npass);

% Draw figure repeatedly with pausing
% Using EMBC2006 data

% ---- To do ----
% make the tail move one base forward after the frameshift.... should it? 
% Only the codon chosen for wait-time calculation differs, not the
% tail-mRNA alignment for calculating free energy

clear; clc;

load prfB_signal.mat; % From EMBC2006/Code
load prfB_seq.mat; seqpreset=12; mark_point=26; 
load EcoliTAV.mat; load EcoliCodons.mat; % From EMBC2006/Code

% Set parameters
phi_sp = -20*(pi/180); initialx = 0.05; C1 = 0.005; C2 = initialx; Nstop=1000;
% Run it through the model
[mag,theta,x,Dvec,Nloop] = calcmpx4(seq(seqpreset+1:end),signal,phi_sp,Names,TAV,C1,C2,Nstop); 

figure(1);
for k=1:length(x)-1
    % ---- Sequence alignment ----
    mRNA = seq(3*k:min(3*k+18,length(seq))); tail = 'auuccuccacuag';
    subplot(221); plot(0,0);
    text(-0.75,0.5,'16Stail-mRNA alignment'); text(- 0.75,0.4,'------------------------------------');
    text(-0.5,-0.35,tail); text(-0.5,-0.5,mRNA);
    set(gca,'xtick',[]); set(gca,'ytick',[]);

    % ---- Table of values ----
    subplot(222); plot(0,0); 
    text(-0.5,0.75,strcat('k=',num2str(k))); text(-0.5,0.5,strcat('\theta_{k}=',num2str((180/pi)*theta(k))));
    text(-0.5,0.25,strcat('Mag(D)=',num2str(Dvec(k,1)))); text(-0.5,0,strcat('Phase(D)=',num2str((180/pi)*Dvec(k,2)))); 
    text(-0.5,-0.25,strcat('Nloop=',num2str(Nloop(k)))); text(-0.5,-0.5,strcat('x=',num2str(x(k))));
    set(gca,'xtick',[]); set(gca,'ytick',[]);

    % ---- Polar plot ----
    subplot(223); polar(theta(1:k),mag(1:k)); title('Polar plot');
    % axis([0 mag(end) 0 mag(end)]);
    hold on; if ge(k,mark_point), polar(theta(mark_point), mag(mark_point), 'r*'), end; hold off; 

    % ---- Displacement plot ----
    subplot(224); plot(1:k, x(1:k)); title('Displacement plot'); grid on; xlabel('Codon Number'); ylabel('x(k)');
    axis([1 length(x) min(-4,min(x)) max(4,max(x))]);

    % ---- Choose the appropriate pause option ----
    if k>20, pause, else, pause(0.001), end
    % keyboard; 
    % pause;
    % pause(0.001);
end