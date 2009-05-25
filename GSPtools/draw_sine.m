% Draws an animation that shows (1) the sine-wave fit done by Ponnala
% et al.'s research at each codon; (2) the codon window; and (3) the
% probabilities of frameshifting. Takes the filename of a sequence,
% like unity!
function draw_sine(file)
    global Travel;
    Travel = load_travel();
    [signal, seq] = get_signal(file);
    [mag, phase] = cumm_energy(signal);
    dvec = inst_energy(mag, phase);

    mag = dvec(1, :);
    phase = dvec(2, :);
    figure(1000);
   
    old_codon_n = 0;

    function helper(x0, probs, codon_n)
        if codon_n ~= old_codon_n, pause; end;
        
        clf;
        axis([xlim -3 3]);
        xlabel('Position');
        ylabel('Magnitude');
        title(file);
        subplot(2, 1, 1);
        draw_one_sine(mag(codon_n), phase(codon_n));

        subplot(2, 1, 2);
        draw_sequence((codon_n-1)*3 + 2, seq(13:end));
        rectangle('Position', [x0 - 3, 0.5, 6, 1.0]);
        pause(0.01);
        
        old_codon_n = codon_n;
    end    

    displacement(seq(13:end), dvec, [], @helper);
end

% Draws one sine wave like it says.
function draw_one_sine(mag, phase)
    x = -6:0.1:6;
    y = mag*sin((2/3)*pi*(x - phase));
    plot(x, y);
end
