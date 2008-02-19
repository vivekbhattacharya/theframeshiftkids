% What the hell? This plots force?
function quasiforce(file)
    [signal, seq] = get_signal(file);
    [mag, theta] = cumm_mag_phase(signal);
    [dvec, theta] = diff_vectors(mag, theta);
    [x, wts] = displacement(seq(13:end), dvec, []);
    
    theta = theta * 180/pi;
    figure(1); plot(theta(1:length(wts)), wts, 'LineWidth', 2);
    grid; xlabel('Angles (degrees)'); ylabel('Weights');
end