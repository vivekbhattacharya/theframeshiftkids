% Computes the energy value at displacement x.
function y = energize(self, energy, codon, x)
    global Config;
    shift = get(self, 'shift');

    x      = x - 2*shift;
    phi_dx = (pi/3)*x - Config.phi_sp;
    frame  = mod(shift, 3) + 1;

    mag    = self.dvec(frame, codon, 1);
    phase  = self.dvec(frame, codon, 2);

    y = -mag*cos((pi/3)*x + phase - Config.phi_sp);
end
