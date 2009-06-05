% Computes the energy value at displacement x.
function y = energize(self, energy, x)
    global Config;
    shift = get(self, 'shift');
    x = x - 3/pi*Config.phi_sp;
    y = polyval(energy, x - 2*shift);
end
