% Superimposes N iterations of the displacement
% plot for gene.txt, saving that file to
% opportunity\gene.png in the work folder for
% each gene (txt or FASTA) in the work folder.
%
% Usage: opportunity('C:\work folder', N)
function opportunity(folder, limit)
    classify(folder, 'opportunity', @helper);

    function helper(displacement, n, filename, image)
        plots = cell(1, limit*2);
        disp(filename);

        for i=1:limit
            x = displacement([]);
            plots(i*2 - 1) = {1:length(x)};
            plots(i*2) = {x};
        end

        top = max(3, ceil(max(x))); bottom = min(-3, floor(min(x)));
        max_dom = 50*(ceil(length(x)/50));

        h = figure(1); set(h, 'Renderer', 'OpenGL');
            set(h, 'Visible', 'off');
            plot(0,0); plot(plots{:}); title(filename);
            grid; xlabel('Codon Number'); ylabel('x(k)');
            axis([0, max_dom, bottom, top]);
        saveas(h, image, 'png');
    end
end
