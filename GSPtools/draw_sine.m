% Draws an animation that shows (1) the sine-wave fit done by Ponnala
% et al.'s research at each codon; (2) the codon window; and (3) the
% probabilities of frameshifting. Takes the sequence string without
% the leader sequence, the vector from inst_energy(), and chunky: the
% displacement vector from displacement().
function draw_sine(seq, dvec, chunky)
    mag = dvec(1, :);
    phase = dvec(2, :);
    figure(1000);
    axis([xlim -3 3]);
    xlabel('Position');
    ylabel('Magnitude');
    title('O Spacious Skies');

    for i = 1:length(mag)
        clf;
        subplot(2, 1, 1);
        draw_one_sine(mag(i), phase(i));

        subplot(2, 1, 2);
        draw_sequence((i-1)*3 + 1, seq, chunky);
        pause
    end
end

% Draws one sine wave like it says.
function draw_one_sine(mag, phase)
    x = -4:0.1:4;
    y = mag*sin((2/3)*pi*(x - phase));
    plot(x, y);
end
