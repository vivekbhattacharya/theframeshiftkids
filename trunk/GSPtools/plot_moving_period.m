function plot_moving_period(file, interval)

    global Config;
    config;
    signal = get_signal(file);

    pvals = [];
    
    max_iter = length(signal) - interval + 1;
    
    for i = 1:max_iter,
        [H, pval] = plot_period(file, 0.05, 0, signal(i:i + interval - 1));
        pvals = [pvals pval];
    end;
    
    plot(-log10(pvals));
        myylim = ylim;
        axis([xlim 0 myylim(2)]);
        xlabel('Initial Nucleotide Number');
        ylabel('Negative Log of p-Value');
        title(file);
        