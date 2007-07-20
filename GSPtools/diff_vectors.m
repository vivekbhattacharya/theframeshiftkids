function [Dvec, Phase] = diff_vectors(Mag, Phase)

% ------------------------------------------------------------
% CALCULATE DIFFERENTIAL VECTORS
% ------------------------------------------------------------    
L = 3; upper = length(Mag);
for k=1:upper
    x = min(max(1,k-1),upper-2);
    index = [x:x+2];
    dA_dc(1,k) = fakeslope(1:L,Mag(index));
    dphi_dc(1,k) = fakeslope(1:L,Phase(index));
    
    D = exp(j*Phase(1,k))*(dA_dc(1,k) + j*(Mag(1,k)*dphi_dc(1,k)));
    
    Dvec(k,1) = abs(D); 
    Dvec(k,2) = angle(D);
end

for k=1:length(Phase)
    if Phase(k)<0; Phase(k)=Phase(k)+(2*pi); end
end

clear global sands shoals;