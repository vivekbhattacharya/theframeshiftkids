% Given a file or a folder of files, polardeviation
% outputs the deviation of the first few codons phasors'
% angles from the species angle, thus attempting to measure
% protein yield.
function polardeviation(folder)
    clear Config;
    config();
    global Config;

    % `mod` ensures 0 < phi_sp < 2*pi
    phi_sp = mod(Config.phi_sp, 2*pi);

    classify(folder, '', @helper, 'preparation');
    function helper(path, filename, image)
        theta = sowseal_surprise(path);
        sigma = sqrt(sum((theta - phi_sp) .^ 2));
        fprintf('%s: %g\n', filename, sigma);
    end
    fprintf('\n');

    function [theta] = sowseal_surprise(file)
        signal = get_signal(file);
        energy = cumm_energy(signal);
        theta = energy(2, 20:100);
    end
end
