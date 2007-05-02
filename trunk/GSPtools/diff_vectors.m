function [Dvec] = diff_vectors(Mag, Phase, numcodons)

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