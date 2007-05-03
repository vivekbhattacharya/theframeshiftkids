% slam: "Function to calculate magnitude, phase, displacement"
% Calculates the cummulative magnitude, phase and the
% displacement vector based on the frameshift model. 
% Aggregates all useful and desirable features into one function. 
% 
% 1. Allows for spacing between the tail and the A-site. 
% 2. Starts calculating wait-time from codon #2 onwards (since the start codon
%    is locked into the P-site at translation initiation)
% 3. Returns the differential vectors (Dvec) and the number of wait-cycles
%    for each codon (Nloop). 
% 
% USAGE: [Mag,Phase,InstPhase,x,Dvec,Nloop] = slam(seq,signal,phi_sp,C1,C2)
% seq = string containing the RNA sequence in lower case
% signal = row vector containing the free energy signal
% phi_sp = species-specific phase angle, typically the mean of all the
%          individual signal phase angles
% C1 = step-size of the incremental displacement
% C2 = estimate of initial displacement
% 
% Note: 
% 1. This function does not check if the signal length is a multiple of 3
% 2. Code for calculating displacement is taken from xmodel7.m
%
% Return values:
% Mag, Phase are arrays of data for the polar plot.
% x is the displacement array.

function [Mag,Phase,InstPhase,x,Dvec,Nloop] = slam(seq,signal,phi_sp,C1,C2)
if length(seq) ~= length(signal)
    error('Sequence and signal are NOT of the same length!');
end

% ------------------------------------------------------------
% CALCULATE CUMMULATIVE MAGNITUDE AND PHASE
% ------------------------------------------------------------    
[Mag, Phase, Err] = cumm_mag_phase(signal);
I1 = find(Err(1,:)); I2 = find(Err(2,:));
if ~isempty(I1)
    fprintf(1,'\nMagnitude negative at %d indices:',length(I1));
end
if ~isempty(I2)
    fprintf(1,'\nEquations not satisfied at %d indices:',length(I2));
end
fprintf(1,'\nNumber of codons: %d',length(Mag));
numcodons = length(Mag);

% ------------------------------------------------------------
% CALCULATE DIFFERENTIAL VECTORS
% ------------------------------------------------------------    
L = 3; % L must be odd
P = 1; % Order of the polynomial
for k=1:numcodons
    x = min(max(1,k-1),numcodons-2);
    index = (x:x+2);
    polyMag = polyfit(1:L,Mag(index),P);
    polyPhase = polyfit(1:L,Phase(index),P);
    
    dA_dc(1,k) = polyval(polyder(polyMag),k);
    dphi_dc(1,k) = polyval(polyder(polyPhase),k);
    
    D = exp(j*Phase(1,k))*(dA_dc(1,k) + j*(Mag(1,k)*dphi_dc(1,k)));
    
    Dvec(k,1) = abs(D); 
    Dvec(k,2) = phase(D); 
end

% ------------------------------------------------------------
% CALCULATE DISPLACEMENT 
% ------------------------------------------------------------    
Nloop = []; x = [0 C2]; codon = 0; InstPhase = [];
choices = []; codonsUsed = [];
n = 1; %Nucleotide number
k = 1; %Counter
while 1
    if (k > numcodons-1); break; end;
    this = seq(n:n+2); next = seq(n+1:n+3);
    choice = probe(x(k), this, next); 
    switch choice
        case {0,2}; codon = this; n = n+3;
        case 1; codon = next; n = n+1;
    end
    Nloop(k) = nloopcalcify(codon);
    codonsUsed = [codonsUsed ':' codon];
    choices = [choices choice];
   
    % Teflon code.
    phi_signal(1,k) = Dvec(k,2); x_temp = x(k);
    for wt=1:Nloop(k)
        phi_dx = ((pi/3)*x_temp) - phi_sp;
        dx = -C1*Dvec(k,1)*sin(phi_signal(k) + phi_dx);
        x_temp = dx + x_temp;
    end
	InstPhase(k) = (180/pi)*(phi_signal(k) + phi_dx);
        if(InstPhase(k)<0); InstPhase(k) = InstPhase(k) + 360; end
    x(k+1) = x_temp;

	% Check for stop codon.
    if (strmatch(codon, ['uag'; 'uaa'; 'uga'])); break; end;
    k = k+1;
end
%disp(choices);
%disp(codonsUsed);