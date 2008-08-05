% Plots displacment error bars for gene.txt, saving the plot to
% erroneous\gene.png in the work folder using an N sample size for
% each gene (txt or FASTA) in the work folder.
%
% Usage: erroneous('C:\work folder', N)
function erroneous(folder, limit)
    classify(folder, 'erroneous', @helper);

    function helper(displacement, n, file, image)
        disp(file);

        upper = ceil(n);
        a = zeros(limit, upper);
        for i=1:limit
            x = displacement([]);
            maxiderm = length(x);
            a(i, :) = x;
            a(find(x > maxiderm)) = -3.14;
        end
        [avg, stddev] = exmean(a);

        max_dom = 50*(ceil(length(x)/50));

        h = figure(1);
        set(h, 'Visible', 'off');
        errorbar(1:upper, avg, stddev);
        grid;
        xlabel('Codon number'); ylabel('Displacement');
        axis([0, max_dom, min(-3, floor(min(x))), max(3, ceil(max(x)))]);
        title(file);
        saveas(h, image, 'png');
    end
end
