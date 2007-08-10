% ------------------------------------------------
% Superimposes 5 iterations and saves that plot
% to a file for all files in a given folder where
% each file contains a gene sequence
% 
% Usage: opportunity('C:\folder', 10)
% ------------------------------------------------
function opportunity(folder, limit)
    classify(folder, 'thepictureshow', @helper);
    
    function helper(displacement, n, filename, image)
        plots = cell(1, limit*2);
        disp(filename);
        
        for i=1:limit
            x = displacement({}, {});
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