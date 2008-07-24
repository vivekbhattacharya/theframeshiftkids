% Plots the p-values from applying moving window of free energy signal
% values to `periodogram.m` and saves the image to the
% moving_periodogram subdirectory within the folder you passed as a
% parameter.
%
% Usage: moving_periodogram(folder of genes, size of moving window)
function moving_periodogram(folder, interval)
    config;
    classify(folder, 'moving_periodogram', @helper, 'preparation');

    function helper(fullpath, file, image)
        disp(file);
        signal = get_signal(fullpath);

        max = length(signal) - interval + 1;
        pvals = zeros(1, max);

        for i = 1:max
            pvals(i) = plot_period(signal(i:i + interval - 1), 0);
        end

        h = figure(1);
        set(h, 'Visible', 'Off');
        grid; title(file);
        plot(-log10(pvals));

        myylim = ylim;
        axis([xlim 0 myylim(2)]);
        xlabel('Initial Nucleotide Number');
        ylabel('-log10(p-value)');
        saveas(h, image, 'png');
    end
end


