% Graphically represents the first 30 codons of the phasor's angle and
% magnitude as documented in the mechanics paper by Ponnala, et al.
% Takes the vector from diff_vector as the only argument. Returns 1000
% words.
function phase_angle(file)
    clear global Config;
    config();

    signal = get_signal(file);
    [mag, theta] = cumm_energy(signal);
    dvec = inst_energy(mag, theta);

    figure(10);
    mag = dvec(1, :);
    phase = dvec(2, :);

    mag = mag(1:30);
    phase = phase(1:30) * 180/pi;

    subplot(2, 1, 1);
    plot(1:length(mag), mag);
    xlabel('Codon');
    ylabel('Magnitude');

    subplot(2, 1, 2);
    plot(1:length(phase), phase);
    xlabel('Codon');
    ylabel('Phase');
end