% Computes the energy value at displacement x.
function y = energize(self, energy, x)
    global Config;
    shift = get(self, 'shift');

    x      = x - 2*shift;
    mag    = energy(1);
    phase  = energy(2);

    y = mag*cos((pi/3)*x + phase - Config.phi_sp);
end
