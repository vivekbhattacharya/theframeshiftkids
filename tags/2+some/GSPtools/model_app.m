% Produces additional plots
clear; clc;

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
[mag,theta,inst_theta,x,Dvec] = fcalcmpx(Seq(3*cp+1:end),Signal,phi_sp,Names,TAV,C1,C2,Nstop,spc); 
figure(1); polar(theta,mag); title('Polar plot'); 
figure(2); plot(1+cp:length(x)+cp, x); axis([1 length(x)+cp min(-4,min(x)) max(4,max(x))]); grid; xlabel('Codon number k'); ylabel('Displacement x(k)'); title('Displacement plot'); 
for k=1:length(theta)
    if theta(k)<0
        theta(k)=theta(k)+(2*pi);
    end
end
figure(3); plot(1:length(theta),(180/pi)*theta); xlabel('Codon number'); ylabel('Phase angle (degrees)'); title('Plot of cumulative phase'); 
figure(4); plot(1:length(inst_theta),inst_theta); xlabel('Codon number'); ylabel('Phase angle (degrees)'); title('Plot of instantaneous phase'); 
figure(5); plot(1:size(Dvec,1),Dvec(:,1)); xlabel('Codon number'); title('Magnitude of the differential vector'); 

for k=1:length(x)-1
    diffx(k)=x(k+1)-x(k);
end
figure(6); plot(1:length(diffx),diffx); xlabel('Codon number'); title('Plot of "force", i.e. incremental displacement'); 
