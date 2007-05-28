function [Phase,x,diffx] = displacement(seq,Nstop,spc,Phase,numcodons,Dvec)

phi_sp=-30*(pi/180); initialx = 0.1;
C1 = 0.005; C2 = initialx;
% ------------------------------------------------------------
% CALCULATE DISPLACEMENT 
% ------------------------------------------------------------    
Nloop = []; x = [0 C2]; codon = []; InstPhase = [];
ants = [];

shift = 0; codons = 0;
for k=2:numcodons-1
    % Choose appropriate codon, depending on the specified spacing, and
    % calculate nloop accordingly
    initial = 3*(k-1) + 3*spc;
    try
        codon=seq(initial+1+shift:initial+3+shift);
        other_codon=seq(initial+2+shift:initial+4+shift);
    catch
        break;
    end;
    codons = [codons codon];
    Nloop(k) = nloopcalcify(codon);
    other_nloop = nloopcalcify(other_codon);
    
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
        my_x_temp = x_temp-2*shift;
        weight_abc = cos(my_x_temp*pi/4)^10;             % Window Function
        weight_bcd = sin(my_x_temp*pi/4)^10;
        P_temp_fail_abc = 1- (1/real_loops*weight_abc);
        P_temp_fail_bcd = 1- (1/other_real_loops*weight_bcd);
        P_fail_abc = P_fail_abc * P_temp_fail_abc;
        P_fail_bcd = P_fail_bcd * P_temp_fail_bcd;
        P_abc = 1-P_fail_abc;
        P_bcd = 1-P_fail_bcd;
        P_reloop =1 - P_abc - P_bcd;
        P = [P_abc P_bcd P_reloop];
        
        r = rand; % Mersenne Twister
        if (P_reloop < P_abc) || (P_reloop < P_bcd)
            if(r < P_abc)
                break;
            elseif (r < P_abc + P_bcd)
                shift = shift + 1;
                ants = [ants ' ' codon ',' num2str(k) ';'];
                break;
            end;
        end;

        phi_dx = ((pi/3)*x_temp)-phi_sp;
        dx = -C1*Dvec(k,1)*sin(phi_signal(1,k) + phi_dx);
        x_temp = dx + x_temp;
    end
    
    InstPhase(k) = (180/pi)*(phi_signal(1,k) + phi_dx);
    if InstPhase(k)<0, InstPhase(k) = InstPhase(k) + 360; end
    x(1,k+1) = x_temp;
end

for k=1:length(Phase)
    if Phase(k)<0; Phase(k)=Phase(k)+(2*pi); end
end
for k=1:length(x)-1
    diffx(k)=x(k+1)-x(k);
end

%% Counters for megaunity
global shoals sands;
sands = sands + 1; disp(ants);
if(strcmp(ants, ' uga,25;')), shoals = shoals + 1; end;