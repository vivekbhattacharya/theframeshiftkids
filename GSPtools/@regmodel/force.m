% A method for the @regmodel class, part of the loop() protocol.
% Returns the force at the given position, given the codon number and
% energy vector according to Dr. Ponnala's register model.
function dx = force(self, shift, codon, x)
    global Config;

    x      = x - 2*shift;
    phi_dx = (pi/3)*x - Config.phi_sp;
    frame  = mod(shift, 3) + 1;

    mag    = self.dvec(frame, codon, 1);
    phase  = self.dvec(frame, codon, 2);
    dx     = -mag*sin(phase + phi_dx);
end
