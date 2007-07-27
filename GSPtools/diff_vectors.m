function [Dvec, Phase] = diff_vectors(Mag, Phase)

% ------------------------------------------------------------
% CALCULATE DIFFERENTIAL VECTORS
% ------------------------------------------------------------    
L = 3; upper = length(Mag);
for k=1:upper
    x = min(max(1,k-1),upper-2);
    
    index = [x:x+2];
    a = fakeslope(Mag(index));
    b = fakeslope(Phase(index));
    
    D = exp(j*Phase(k))*(a + j*Mag(k)*b);
    
    Dvec(k,1) = abs(D); 
    Dvec(k,2) = angle(D);
end

for k=1:length(Phase)
    if Phase(k)<0; Phase(k)=Phase(k)+(2*pi); end
end

clear global sands shoals;