function plot_moving_period(folder, interval)
    classify(folder, 'plot_moving_period', @helper);
    
    function helper(plot_period, n, file, image)
        disp(file);

        global Config;
        config;
        signal = get_signal(file);

        pvals = [];

        max_iter = length(signal) - interval + 1;

        for i = 1:max_iter,
            [H, pval] = plot_period(file, 0.05, 0, signal(i:i + interval - 1));
            pvals = [pvals pval];
        end;

        h = figure(1); set(h, 'Renderer', 'OpenGL');
            set(h, 'Visible', 'Off');
            plot(-log10(pvals));
            myylim = ylim;
            axis([xlim 0 myylim(2)]); grid; title(file);
            xlabel('Initial Nucleotide Number'); ylabel('Negative Log of p-Value');
        saveas(h, image, 'png');
    end
end
            
        