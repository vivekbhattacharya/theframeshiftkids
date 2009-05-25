% Draws the codon window given the nucleotide offset, the mRNA string,
% and chunky: the displacement vector from displacement().
function draw_sequence(index, seq, chunky)
    % Take me on. Take ... on ... me!
    hold on;

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
        % ((i - min) + start)*2 takes into account the idea that 2 spaces = 1
        % nucleotide.
        text((i - min + start)*2, y, seq(i));
    end
    hold off;
end
