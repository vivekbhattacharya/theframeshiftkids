% Draws a polar plot for each gene in the work folder
% (txt or FASTA), saving that plot to alopeciaunity\gene.png
% in the work folder.
%
% alopeciaunity('c:\work folder')
function alopeciaunity(folder)
    classify(folder, 'alopeciaunity', @helper, 'preparation');
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