% Draws an animation that shows (1) the sine-wave fit done by Ponnala
% et al.'s research at each codon; (2) the codon window; and (3) the
% probabilities of frameshifting. Takes the filename of a sequence,
% like unity!
function draw_sine(file)
    global Travel store;
    Travel = load_travel();
    [signal, seq] = get_signal(file);
    [mag, phase] = cumm_energy(signal);
    dvec = inst_energy(mag, phase);

    mag = dvec(1, :);
    phase = dvec(2, :);
    figure(1000);
   
    old_codon_n = 0;
    wc = 1;

    function helper(x0, probs, codon_n)
        num_shift = length(store.anthill) - length(store.termites);
        
        clf
        subplot(3, 1, 1);
        draw_one_sine(mag(codon_n), phase(codon_n), file);
        ybig = max(3, ceil(mag(codon_n)));
        line([x0 x0], [mag(codon_n)*sin((2/6)*pi*(x0 - phase(codon_n))), mag(codon_n)*sin((2/6)*pi*(x0 - phase(codon_n))) + 0.4*ybig], 'LineWidth', 2);
        
        subplot(3, 1, 3);
        axis([-6 6 0 1]);
        codon_text = ['Codon: ' num2str(codon_n)];
        wc_text = ['Wait: ' num2str(wc)];
        text(4.1, 0.85, codon_text);
        text(4.1, 0.65, wc_text);
        line([-2 + 2 * num_shift, -2 + 2 * num_shift], [0, probs(1)], 'LineWidth', 2);
        line([0 + 2 * num_shift, 0 + 2 * num_shift], [0, probs(2)], 'LineWidth', 2);
        line([2 + 2 * num_shift, 2 + 2 * num_shift], [0, probs(3)], 'LineWidth', 2);
        
        subplot(3, 1, 2);
        draw_sequence((codon_n-1)*3 + 2, seq(13:end));
        fuzzy_bunny1 = xlim; fuzzy_bunny2 = ylim;
        if codon_n ~= old_codon_n
            rectangle('Position', [x0 - 3, 0.5, 6, 1.0], 'LineWidth', 2, 'EdgeColor', 'g');
            pause;
            wc = 0;
        else
            rectangle('Position', [x0 - 3, 0.5, 6, 1.0], 'LineWidth', 2, 'EdgeColor', 'r');
        end
              
        wc = wc + 1;
        pause(0.01);
        
        old_codon_n = codon_n;
    end    

    displacement(seq(13:end), dvec, [], @helper);
end

% Draws one sine wave like it says.
function draw_one_sine(mag, phase, file)
    x = -6:0.1:6;
    y = mag*sin((2/6)*pi*(x - phase));
    
    ybig = max(3, ceil(mag)); 
    change_axes = (ybig ~= 3);

    plot(x, y);
        axis([xlim -ybig ybig]);
        xlabel('Position');
        ylabel('Magnitude');
        title(file);
        if change_axes
            text(2.2, 0.76*ybig, 'NOTE: LARGER AXES', 'Color', 'r');
        end
end
