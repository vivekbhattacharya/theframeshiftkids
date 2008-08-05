% Draws a polar plot for each gene in the work folder
% (txt or FASTA), saving that plot to polarplot\gene.png
% in the work folder.
%
% polarplot('c:\work folder')
function polarplot(folder)
    classify(folder, 'polarplot', @helper, 'preparation');
    function helper(path, filename, image)
        disp(filename);

        h = figure(2);
        set(h, 'Visible', 'off');
        [d, p] = walrus_surprise(path);
        p([]);
        title(filename);
        set(h, 'Visible', 'off');
        saveas(h, image, 'png');
    end
    disp(' ');
end
