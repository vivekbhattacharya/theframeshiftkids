% FCALCMPX: "Function to calculate magnitude, phase, displacement"
% Calculates the cummulative magnitude, phase and the
% displacement vector based on the frameshift model. 
% Aggregates all useful and desirable features into one function. 
% 
% 1. Allows for spacing between the tail and the A-site. 
% 2. Starts calculating wait-time from codon #2 onwards (since the start codon
%    is locked into the P-site at translation initiation)
% 3. Returns the differential vectors (Dvec) and the number of wait-cycles for
%    each codon (Nloop). 
% 
% USAGE: [Mag,Phase,InstPhase,x,Dvec,Nloop] = fcalcmpx(seq,signal,phi_sp,Names,TAV,C1,C2,Nstop,spc)
% seq = string containing the RNA sequence in lower case
% signal = row vector containing the free energy signal
% phi_sp = species-specific phase angle, typically the mean of all the
%          individual signal phase angles
% Names,TAV = pre-calculated tRNA availability
% C1 = step-size of the incremental displacement
% C2 = estimate of initial displacement
% Nstop = number of loops at the stop codon
% spc = number of codons between the 5'-end of the 16S rRNA tail and
%       the A-site
% 
% Note: 
% 1. This function does not check if the signal length is a multiple of 3
% 2. Code for calculating displacement is taken from xmodel7.m
%
% Return values:
% Mag, Phase: Arrays of data for the polar plot
% x: Displacement array

function [Mag,Phase,InstPhase,x,Dvec,Nloop] = fcalcmpx(seq,signal,phi_sp,Names,TAV,C1,C2,Nstop,spc)

% Make sure sequence and signal are of the same length!
if length(seq)~=length(signal)
    error('Sequence and signal are NOT of the same length!');
end

% ------------------------------------------------------------
% CALCULATE CUMMULATIVE MAGNITUDE AND PHASE
% ------------------------------------------------------------    
[Mag, Phase, Err] = cumm_mag_phase(signal);
I1 = find(Err(1,:));
if ~isempty(I1)
    fprintf(1,'\nMagnitude negative at %d indices:',length(I1));
end
I2 = find(Err(2,:)); 
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
    index = [x:x+2];
    polyMag = polyfit(1:L,Mag(index),P);
    polyPhase = polyfit(1:L,Phase(index),P);
    
    dA_dc(1,k) = polyval(polyder(polyMag),k);
    dphi_dc(1,k) = polyval(polyder(polyPhase),k);
    
    D = exp(j*Phase(1,k))*(dA_dc(1,k) + j*(Mag(1,k)*dphi_dc(1,k)));
    
    Dvec(k,1) = abs(D); 
    % Dvec(k,1) = sqrt(dA_dc(1,k)^2+(Mag(1,k)*dphi_dc(1,k))^2);
    Dvec(k,2) = phase(D); 
    % Dvec(k,2) = Phase(1,k) + atan2(Mag(1,k)*dphi_dc(1,k),dA_dc(1,k));   
end

% ------------------------------------------------------------
% CALCULATE DISPLACEMENT 
% ------------------------------------------------------------    
Nloop = []; x = [0 C2]; codon = []; InstPhase = [];
figure;
shift = 0;
codons = 0;
for k=2:numcodons-1
    % Choose appropriate codon, depending on the specified spacing, and
    % calculate nloop accordingly
    initial = 3*(k-1) + 3*spc;
    try
        codon=seq(initial+1+shift:initial+3+shift);
        other_codon=seq(initial+4+shift:initial+6+shift);
    catch
        break;
    end;
    codons = [codons codon];
    disp(codons);
%   if abs(x(1,k))<1
%       codon=seq(initial+1:initial+3);
%       other_codon=seq(initial+4:initial+6);
%   elseif x(1,k)<-1
%       codon=seq(initial+0:initial+2);
%       other_codon=seq(initial+3:initial+5);
%   elseif x(1,k)>1
%       codon=seq(initial+2:initial+4);
%       other_codon=seq(initial+5:initial+7);
%   end
    Nloop(k)=nloopcalc(codon,0,1,Names,TAV,Nstop);
    other_nloop=nloopcalc(other_codon,0,1,Names,TAV,Nstop);
    
    Nloop_floor = ceil(Nloop(k));
    other_Nloop_floor = ceil(other_nloop);
    
    real_loops = 2^(1/Nloop_floor);
    real_loops = real_loops / (real_loops - 1);
    
    other_real_loops = 2^(1/other_Nloop_floor);
    other_real_loops = other_real_loops / (other_real_loops - 1);
    
    phi_signal(1,k) = Dvec(k,2); x_temp = x(1,k); 
    
    phi_dx = ((pi/3)*x_temp)-phi_sp;
    P_fail_abc = 1;
    P_fail_bcd = 1;
    for wt=1:Nstop
        my_x_temp=x_temp-2*shift;
        weight_abc = cos(my_x_temp*pi/4)^4;             % Window Function
        weight_bcd = sin(my_x_temp*pi/4)^4;
        P_temp_abc = 1/real_loops*weight_abc;
        P_temp_bcd = 1/other_real_loops*weight_bcd;
        P_temp_fail_abc = 1-P_temp_abc;
        P_temp_fail_bcd = 1-P_temp_bcd;
        P_fail_abc = P_fail_abc * P_temp_fail_abc;
        P_fail_bcd = P_fail_bcd * P_temp_fail_bcd;
        P_abc = 1-P_fail_abc;
        P_bcd = 1-P_fail_bcd;
%       P_abc = 1-(((1-1/real_loops)*cos(x_temp*pi/4)^2))^wt;
%       P_bcd = 1-(((1-1/other_real_loops)*sin(x_temp*pi/4)^2))^wt;
        P_reloop =1 - P_abc - P_bcd;
        P = [P_abc P_bcd P_reloop];
        
        
        %disp(P);
        if ((P_reloop < P_abc) || (P_reloop < P_bcd))
            %disp(['life is complete' codon other_codon num2str(k)]);
            if max(P) == P_bcd
                shift = shift + 1;
            end;
            break;
        end;
        phi_dx = ((pi/3)*x_temp)-phi_sp;
        dx = -C1*Dvec(k,1)*sin(phi_signal(1,k) + phi_dx); % Correct
        x_temp = dx + x_temp;
    end

    InstPhase(k) = (180/pi)*(phi_signal(1,k) + phi_dx);
    if InstPhase(k)<0, InstPhase(k) = InstPhase(k) + 360; end
    x(1,k+1) = x_temp;  
    disp(x_temp);
end
disp(shift);
hold off
xlabel('Codon number, k');
ylabel('Total angle');