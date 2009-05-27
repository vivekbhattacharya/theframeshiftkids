% Draws an animation that shows (1) the sine-wave fit done by Ponnala
% et al.'s research at each codon; (2) the codon window; and (3) the
% probabilities of frameshifting. Takes the filename of a sequence,
% like unity!
function draw_sine(file)
    global Travel store Config;
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
        mag_i     = mag(codon_n);
        phase_i   = phase(codon_n);
        clf;

        % Draw the signal.
        subplot(3, 1, 1);
        hold on;
        draw_one_sine(mag_i, phase_i, file);
        ymax = max(3, ceil(mag_i));

        % Position marker.
        x = -mag_i*cos((1/3)*pi*x0 + phase_i - Config.phi_sp);
        h = plot(x0, x, 'ko', ...
                 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        hold off;

        % Draw the probabilities.
        subplot(3, 1, 3);
        hold on;
        bar([-2 0 2] + 2*num_shift, probs);
        axis([-6 6 0 1]);
        xlabel('Position');
        ylabel('Probability');
        text(4.1, 0.85, sprintf('Codon: %g', codon_n));
        text(4.1, 0.65, sprintf('Wait: %g', wc));
        hold off;

        % Draw the mRNA sequence window.
        subplot(3, 1, 2);
        hold on;
        draw_sequence((codon_n - 1)*3 + 2, seq(13:end));
        h = rectangle('Position', [x0 - 3, 0.5, 6, 1.0], ...
                      'LineWidth', 2, 'EdgeColor', 'r');

        hold off;

        % Check to see if we're at a new codon.
        if codon_n ~= old_codon_n
            set(h, 'EdgeColor', 'g');
            wc = 0;
            pause;
        end

        wc = wc + 1;
        old_codon_n = codon_n;
        pause(0.005);
    end

    displacement(seq(13:end), dvec, [], @helper);
end

% Draws one sine wave like it says.
function draw_one_sine(mag, phase, file)
    global Config;
    x = -6:0.1:6;
    y = -mag*cos((1/3)*pi*x + phase - Config.phi_sp);
    ymax = max(3, ceil(mag));

    plot(x, y);
    axis([xlim -ymax ymax]);
    ylabel('Magnitude');
    title(file);
    if ymax ~= 3
        text(2.2, 0.76*ymax, 'NOTE: LARGER AXES', 'Color', 'r');
    end
end

% Draws the codon window given the nucleotide offset, the mRNA string,
% and chunky: the displacement vector from displacement().
function draw_sequence(index, seq)
    % Number of nucleotides to show minus one.
    len = 6;
    % Where to start the nucleotides on the axis.
    start = -2;
    % Left and right boundaries for text.
    min = index + start;
    max = min + len;
    % y-value to plot the text.
    y = 1;

    axis([-6 6 (y - 1) (y + 1)]);
    stop = length(seq);
    for i = min:max
        if (i < 1) || (i > stop)
            continue
        end

        % (i - min) shifts the interval to 0:len.
        % ((i - min) + start) shifts the interval to -2:len-2.
        % ((i - min) + start)*2 takes into account the idea that 2
        % spaces = 1 nucleotide.
        h = text((i - min + start)*2, y, upper(seq(i)));
        set(h, 'FontWeight', 'bold');
        set(h, 'FontSize', 14);
        set(h, 'FontName', 'Consolas');
        set(h, 'HorizontalAlignment', 'center');
    end
end
