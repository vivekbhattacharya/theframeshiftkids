% A wrapper around cumm_energy and inst_energy so displacement is
% abstracted above energy calculations.
function [energy, force] = oldenergy(codon)
    global store Config;
    global oldenergy_dvec;
    global oldenergy;

    if length(oldenergy_dvec) == 0
        [mag, phase]   = cumm_energy(store.signal);
        dvec           = inst_energy(mag, phase);
        energy         = [mag; phase];
        oldenergy_dvec = dvec;
        oldenergy      = energy;
    end

    dvec   = oldenergy_dvec;
    energy = oldenergy;

    function [dx, force] = helper(codon, x)
        x      = x - 2*store.shift;
        phi_dx = (pi/3)*x - Config.phi_sp;
        frame  = mod(store.shift, 3) + 1;

        mag    = dvec(frame, codon, 1);
        phase  = dvec(frame, codon, 2);
        force  = [mag phase];
        dx     = -mag*sin(phase + phi_dx);
    end
    force = @helper;
    store.force_closure = @helper;
end
