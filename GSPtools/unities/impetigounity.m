function impetigounity(folder, times, varargin)
    l = length(varargin);
    x_bound = -1; y_bound = -1;
    if l >= 1, x_bound = varargin{1}; end;
    if l >= 2, y_bound = varargin{2}; end;

    classify(folder, 'jail', @plot_average);
    function plot_average(displacement, n, filename, image)
        disp(filename);
        waitress = zeros(times, 0);
        for i = 1:times
            [x, waits] = displacement({}, {});
            waitress(i, 1:length(waits)) = waits;
        end
        tip = mean(waitress);
        
        viveklikes = x_bound;
        tospin = y_bound;
        
        if x_bound < 0, viveklikes = length(tip); end;
        if y_bound < 0, tospin = ylim;
        else tospin = [0 y_bound];
        end
        
        h = figure(1); set(h, 'Renderer', 'OpenGL');
            plot(0,0); bar(1:viveklikes, tip(1:viveklikes)); title(filename);
            axis([0 viveklikes+1 tospin]);
            grid; xlabel('Codon Number'); ylabel('Wait Time');
        saveas(h, image, 'png');
    end
end