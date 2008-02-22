function quasiforce(file)
    config();
    global Travel Names;
    if isempty(Travel), load Travel2.mat; end

    [signal, seq] = get_signal(file);
    [mag, theta] = cumm_mag_phase(signal);
    [dvec, theta] = diff_vectors(mag, theta);
    [x, wts] = displacement(seq(13:end), dvec, [25]);
    
    theta = theta * 180/pi;
    figure(1);
    % displacement uses theta in a way to make an entire codon's worth of
    % thetas disappear.
    % scatter(theta(1:25), wts(1:25), 49, 'filled');
    plot(theta(1:25), wts(1:25), 'LineWidth', 2);
    grid; xlabel('Angles (degrees)'); ylabel('Weights');
end