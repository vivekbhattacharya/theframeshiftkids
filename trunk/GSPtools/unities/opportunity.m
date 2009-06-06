% Superimposes N iterations of the displacement plot for a directory
% of genes, saving that file to directory/opportunity/gene.png.
%
% > opportunity('/path/to/genes', N)
function opportunity(folder, limit)
    classify(folder, 'opportunity', [], @helper);

    function helper(model, file, png)
        plots = cell(1, limit*2);
        disp(file);

        for i=1:limit
            [bye, x] = displacement(model);
            plots(i*2 - 1) = {1:length(x)};
            plots(i*2) = {x};
        end

        top    = max(3, ceil(max(x)));
        bottom = min(-3, floor(min(x)));
        upper  = 50*(ceil(length(x)/50));

        fprintf('Saving to %s\n', png);

        h = figure(3000);
        set(h, 'Visible', 'off');
        plot(0,0);
        plot(plots{:});

        grid;
        title(file);
        xlabel('Codon Number');
        ylabel('Displacement');
        axis([0 upper bottom top]);
        saveas(h, png, 'png');
    end
end
