% Graphically represents the first 30 codons of the phasor's angle and
% magnitude as documented in the mechanics paper by Ponnala, et al.
% Takes the vector from diff_vector as the only argument. Returns 1000
% words.
function phase_angle(dvec)
    figure(10);
    mag = dvec(:, 1);
    mag = mag(1:30);
    phase = dvec(:, 2);
    phase = phase(1:30);

    subplot(2, 1, 1);
    plot(1:length(mag), mag);
    xlabel('Codon');
    ylabel('Magnitude');
    
    subplot(2, 1, 2);
    plot(1:length(phase), phase);
    xlabel('Codon');
    ylabel('Phase');
end