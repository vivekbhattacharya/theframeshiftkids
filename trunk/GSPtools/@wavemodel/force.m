% A method for the class model, part of the loop() protocol. Returns
% the force at the given position, given the codon number and energy
% vector.
function dx = force(self, energy, codon, x)
    global Config;
    shift = get(self.model, 'shift');

    % Energy parameters from energy() method.
    M     = energy(1);
    theta = energy(2);

    x  = x - 2 * shift;
    dx = M * sin((pi/3) * x + theta - Config.phi_sp);
end
