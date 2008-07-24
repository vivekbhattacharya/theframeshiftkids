function plot_moving_period(folder, interval)
    config;
    classify(folder, 'plot_moving_period', @helper, 'preparation');

    function helper(fullpath, file, image)
        disp(file);
        signal = get_signal(fullpath);

        max = length(signal) - interval + 1;
        pvals = zeros(1, max);

        for i = 1:max,
            pvals(i) = plot_period(signal(i:i + interval - 1), 0);
        end;

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


