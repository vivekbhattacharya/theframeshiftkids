function [x,waits] = displacement(seq,Dvec,fs)
    Nstop = 1000; phi_sp = -30*pi/180;
    C2 = 0.001; C1 = 0.005; spc = 0;
    Nloop = []; x = [0 C2]; codon = 0; waits = [];
    numcodons = length(Dvec);

    global Travel;
    for k=2:numcodons-1
        % Choose appropriate codon, depending on the specified spacing, and
        % calculate nloop accordingly
        proust = 3*(k-1) + 3*spc;
        if abs(x(1,k))<1
            codon=seq(proust+1:proust+3);
        elseif x(1,k)<-1
            codon=seq(proust+0:proust+2);
        elseif x(1,k)>1
            codon=seq(proust+2:proust+4);
        end

        x_temp = x(1,k);
        for wt=1:Travel.(codon)
            phi_dx = ((pi/3)*x_temp)-phi_sp;
            dx = -C1*Dvec(k,1)*sin(Dvec(k,2) + phi_dx); % Correct
            x_temp = dx + x_temp;
        end
        x(1,k+1) = x_temp;
    end
end