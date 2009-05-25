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

    function helper(x0, probs, codon_n)
        probs
        num_shift = length(store.anthill) - length(store.termites);
        
        clf;
        subplot(3, 1, 1);
        draw_one_sine(mag(codon_n), phase(codon_n), file);

        subplot(3, 1, 3);
        axis([-6 6 0 1]);
        rectangle('Position', [1 2 3 4]);
        line([-2 + 2 * num_shift, -2 + 2 * num_shift], [0, probs(1)], 'LineWidth', 2);
        line([0 + 2 * num_shift, 0 + 2 * num_shift], [0, probs(2)], 'LineWidth', 2);
        line([2 + 2 * num_shift, 2 + 2 * num_shift], [0, probs(3)], 'LineWidth', 2);
        
        subplot(3, 1, 2);
        draw_sequence((codon_n-1)*3 + 2, seq(13:end));
        fuzzy_bunny1 = xlim; fuzzy_bunny2 = ylim;
        if codon_n ~= old_codon_n
            rectangle('Position', [x0 - 3, 0.5, 6, 1.0], 'LineWidth', 2, 'EdgeColor', 'g');
            pause;
        else
            rectangle('Position', [x0 - 3, 0.5, 6, 1.0], 'LineWidth', 2, 'EdgeColor', 'r');
        end
                
        pause(0.01);
        
        old_codon_n = codon_n;
    end    

    displacement(seq(13:end), dvec, [], @helper);
end

% Draws one sine wave like it says.
function draw_one_sine(mag, phase, file)
    x = -6:0.1:6;
    y = mag*sin((2/6)*pi*(x - phase));
    plot(x, y);
        axis([xlim -3 3]);
        xlabel('Position');
        ylabel('Magnitude');
        title(file);
end
