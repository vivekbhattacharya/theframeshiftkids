function quasiforce(folder)
    config();
    global Travel;
    Travel = load_travel();

    classify(folder, 'sourcefource!', @helper, 'preparation');
    function helper(path, filename, image)
        disp(filename);

        h = graph(filename);
        [folder, file, ext] = fileparts(image);
        saveas(h, image, 'png');
        saveas(h, fullfile(folder, [file '.fig']), 'fig');
    end
end

function [h] = graph(file)
    [signal, seq] = get_signal(file);
    [mag, theta] = cumm_mag_phase(signal);
    [dvec, theta] = diff_vectors(mag, theta);
    [x, wts] = displacement(seq(13:end), dvec, [25]);

    theta = theta * 180/pi;
    h = figure(1);
    % displacement uses theta in a way to make an entire codon's worth of
    % thetas disappear.
    % scatter(theta(1:25), wts(1:25), 49, 'filled');
    plot3(theta(1:25), wts(1:25), 1:25, '-.r*', 'LineWidth', 2);
    campos([-320 -910 73]); grid;
    xlabel('Angles (degrees)'); ylabel('Weights');
end
