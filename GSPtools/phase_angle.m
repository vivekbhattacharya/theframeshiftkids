% Graphically represents the first 30 codons of the phasor's angle and
% magnitude as documented in the mechanics paper by Ponnala, et al.
% Takes the vector from diff_vector as the only argument. Returns 1000
% words.
function phase_angle(file)
    clear global Config;
    global Config;
    config();

    signal = get_signal(file);
    [a, b] = cumm_energy(signal);
    dvec = inst_energy(a, b);

    max = 31;

    figure(10);
    mag = dvec(1, 1:max);
    phase = dvec(2, 1:max);
    phase = phase - Config.phi_sp;

    % Polar vector derivative of the energy to produce force.
    index = 1:(max-1);
    dx = mag(index+1) .* cos(phase(index+1)) - mag(index) .* cos(phase(index));
    dy = mag(index+1) .* sin(phase(index+1)) - mag(index) .* sin(phase(index));
    dmag = sqrt(dx .^ 2 + dy .^ 2);
    dphase = atan2(dy, dx) * 180/pi;

    subplot(2, 1, 1);
    plot(1:length(dmag), dmag);
    xlabel('Codon');
    ylabel('Magnitude');

    subplot(2, 1, 2);
    plot(1:length(dphase), dphase);
    xlabel('Codon');
    ylabel('Phase');
end
